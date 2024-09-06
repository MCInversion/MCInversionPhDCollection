#include "EvolverUtilsCommon.h"

#include <numeric>

#include "pmp/SurfaceMesh.h"
#include "pmp/ManifoldCurve2D.h"
#include "pmp/algorithms/Remeshing.h"
#include "pmp/algorithms/DifferentialGeometry.h"
#include "geometry/IcoSphereBuilder.h"

//#include <Spectra/SymEigsSolver.h>
//#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Eigen/SparseCore>

#include "geometry/MeshAnalysis.h"


CoVolumeStats AnalyzeMeshCoVolumes(pmp::SurfaceMesh& mesh, const AreaFunction& areaFunction)
{
	// vertex property for co-volume measures.
	if (!mesh.has_vertex_property(coVolMeasureVertexPropertyName))
		mesh.add_vertex_property<pmp::Scalar>(coVolMeasureVertexPropertyName);
	auto vCoVols = mesh.get_vertex_property<pmp::Scalar>(coVolMeasureVertexPropertyName);

	CoVolumeStats stats{};
	for (const auto v : mesh.vertices())
	{
		const double measure = areaFunction(mesh, v);
		vCoVols[v] = static_cast<pmp::Scalar>(measure);

		if (stats.Max < measure) stats.Max = measure;
		if (stats.Min > measure) stats.Min = measure;
		stats.Mean += measure;
	}
	stats.Mean /= static_cast<double>(mesh.n_vertices());
	return stats;
}

/// \brief The power of the stabilizing scale factor.
constexpr float SCALE_FACTOR_POWER = 1.0f / 2.0f;
/// \brief the reciprocal value of how many times the surface area element shrinks during evolution.
constexpr float INV_SHRINK_FACTOR = 5.0f;

float GetStabilizationScalingFactor(const double& timeStep, const float& icoRadius, const unsigned int& icoSubdiv, const float& stabilizationFactor)
{
	const unsigned int expectedVertexCount = (N_ICO_EDGES_0 * static_cast<unsigned int>(pow(4, icoSubdiv) - 1) + 3 * N_ICO_VERTS_0) / 3;
	const float weighedIcoRadius = icoRadius;
	const float expectedMeanCoVolArea = stabilizationFactor * (4.0f * static_cast<float>(M_PI) * weighedIcoRadius * weighedIcoRadius / static_cast<float>(expectedVertexCount));
	return pow(static_cast<float>(timeStep) / expectedMeanCoVolArea * INV_SHRINK_FACTOR, SCALE_FACTOR_POWER);
}

std::string InterpretSolverErrorCode(const Eigen::ComputationInfo& cInfo)
{
	if (cInfo == Eigen::Success)
		return "Eigen::Success";

	if (cInfo == Eigen::NumericalIssue)
		return "Eigen::NumericalIssue";

	if (cInfo == Eigen::NoConvergence)
		return "Eigen::NoConvergence";

	return "Eigen::InvalidInput";
}

pmp::vec3 ComputeTangentialUpdateVelocityAtVertex(const pmp::SurfaceMesh& mesh, const pmp::Vertex& v, const pmp::vec3& vNormal, const float& weight)
{
	pmp::vec3 result{};
	if (!mesh.has_vertex_property("v:normal"))
		return result;

	auto vCirculator = mesh.vertices(v);
	for (const auto w : vCirculator)
	{
		const auto e0 = mesh.position(w) - mesh.position(v);
		const auto e1 = mesh.position(*(++vCirculator)) - mesh.position(v);
		--vCirculator;
		const auto e0Normalized = pmp::normalize(e0);
		const auto e1Normalized = pmp::normalize(e1);
		const auto eDot = pmp::dot(e0Normalized, e1Normalized);
		result += (1.0f + eDot) * (e0 + e1);
	}
	result *= weight / static_cast<float>(mesh.valence(v));
	const auto resultDotNormal = pmp::dot(result, vNormal);
	return (result - resultDotNormal * vNormal);
}

pmp::vec2 ComputeTangentialUpdateVelocityAtVertex(const pmp::ManifoldCurve2D& curve, const pmp::Vertex& v, const pmp::vec2& vNormal, const float& weight)
{
	pmp::vec2 result{};
	if (!curve.has_vertex_property("v:normal") || curve.is_isolated(v) || curve.is_boundary(v))
		return result;

	const auto [vPrev, vNext] = curve.vertices(v);
	const auto e0Vec = curve.position(vPrev) - curve.position(v);
	const auto e1Vec = curve.position(vNext) - curve.position(v);
	result += 0.5f * (e0Vec + e1Vec);
	const auto resultDotNormal = pmp::dot(result, vNormal);
	return (result - resultDotNormal * vNormal);
}

/// \brief a unit speed of a shrink-wrapping sphere sufficiently far away from target.
constexpr double BASE_DISTANCE_MULTIPLIER = 1.0;

/// \brief a factor by which distance field variance is multiplied.
constexpr double BASE_DISTANCE_VARIANCE = 1.0;

AdvectionDiffusionParameters PreComputeAdvectionDiffusionParams(const double& distanceMax, const double& targetMinDimension)
{
	// the radius within which target object starts slowing down shrink wrapping.
	const double distanceVariance = BASE_DISTANCE_VARIANCE; //* 0.125 * targetMinDimension * targetMinDimension;
	const double distanceMultiplier = BASE_DISTANCE_MULTIPLIER / (1.0 - exp(-distanceMax * distanceMax / distanceVariance));
	return { distanceMultiplier, distanceVariance, distanceMultiplier, 1.0 };
}

bool IsRemeshingNecessary(const CoVolumeStats& stats, const double& tStep)
{
	return stats.Max > 1.2 * tStep;
}

constexpr float JACOBIAN_COND_MIN = 1.0f;
constexpr float JACOBIAN_COND_MAX = 1.5f;

bool IsRemeshingNecessary(const std::vector<float>& equilateralJacobianConditionNumbers)
{
	float minVal = FLT_MAX;
	float maxVal = -FLT_MAX;
	for (const auto& val : equilateralJacobianConditionNumbers)
	{
		if (val < minVal) minVal = val;
		if (val > maxVal) maxVal = val;
	}
	return !((minVal > JACOBIAN_COND_MIN && minVal < JACOBIAN_COND_MAX) &&
		     (maxVal > JACOBIAN_COND_MIN && maxVal < JACOBIAN_COND_MAX));
}

bool IsNonFeatureRemeshingNecessary(const pmp::SurfaceMesh& mesh)
{
	float minVal = FLT_MAX;
	float maxVal = -FLT_MAX;
	const auto vQualityProp = mesh.get_vertex_property<float>("v:equilateralJacobianCondition");
	const auto vIsFeature = mesh.get_vertex_property<bool>("v:feature");

	for (const auto v : mesh.vertices())
	{
		if (vIsFeature[v])
			continue;

		const auto& val = vQualityProp[v];
		if (val < minVal) minVal = val;
		if (val > maxVal) maxVal = val;
	}
	return !((minVal > JACOBIAN_COND_MIN && minVal < JACOBIAN_COND_MAX) &&
		(maxVal > JACOBIAN_COND_MIN && maxVal < JACOBIAN_COND_MAX));
}

// ------------------------------------------------------------------------------------

namespace
{
	[[nodiscard]] double ComputeSpectralRadius(const SparseMatrix& mat)
	{
		//// Define the operation for Spectra (for symmetric matrices)
		//Spectra::SparseSymMatProd<double> op(mat);

		//// Construct the eigen solver object, requesting the largest magnitude eigenvalue
		//Spectra::SymEigsSolver solver(op, 1, 6);

		//// Initialize and compute using the correct sorting rule
		//solver.init();
		//int nconv = solver.compute(Spectra::SortRule::LargestAlge);

		//// Check if the computation was successful
		//if (solver.info() != Spectra::CompInfo::Successful)
		//{
		//	throw std::runtime_error("Spectra solver did not converge!");
		//}

		//// Retrieve the eigenvalues
		//Eigen::VectorXd eigenvalues = solver.eigenvalues();

		//// Return the largest absolute eigenvalue (spectral radius)
		//return std::abs(eigenvalues[0]);
		return 1.0;
	}
	
} // anonymous namespace

bool IsRemeshingNecessary(const SparseMatrix& lswMatrix)
{
	return ComputeSpectralRadius(lswMatrix) >= 1.0;
}

// ------------------------------------------------------------------------------------

namespace
{
	[[nodiscard]] std::map<pmp::Vertex, float> ComputeLocalDensities(const pmp::ManifoldCurve2D& curve, float edgeLength)
	{
		// TODO: debug this! Perhaps we should use the actual co-volumes
		std::map<pmp::Vertex, float> densities;

		for (auto v : curve.vertices())
		{
			float localDensity = 0.0f;
			int count = 0;

			// Get the outgoing and incoming edges
			auto [e_out, e_in] = curve.edges(v);

			// Sum the lengths of adjacent edges
			if (e_out.is_valid())
			{
				localDensity += curve.edge_length(e_out);
				count++;
			}

			if (e_in.is_valid())
			{
				localDensity += curve.edge_length(e_in);
				count++;
			}

			// If there are adjacent edges, compute the local density
			if (count > 0)
			{
				localDensity /= count;
				densities[v] = 1.0f / localDensity; // Local density is the inverse of average edge length
			}
			else
			{
				densities[v] = 1.0f / edgeLength; // Default to edgeLength for isolated vertices
			}
		}

		return densities;
	}

	[[nodiscard]] float ComputeMeanDensity(const std::map<pmp::Vertex, float>& densities)
	{
		const float totalDensity = std::accumulate(densities.begin(), densities.end(), 0.0f,
			[](float sum, const std::pair<pmp::Vertex, float>& p) { return sum + p.second; });

		return totalDensity / densities.size();
	}

	[[nodiscard]] float ComputeDensityVariance(const std::map<pmp::Vertex, float>& densities, float meanDensity)
	{
		float variance = 0.0f;
		for (const auto& [v, density] : densities)
		{
			variance += (density - meanDensity) * (density - meanDensity);
		}

		return variance / densities.size();
	}
}

bool IsRemeshingNecessary(const pmp::ManifoldCurve2D& curve, const pmp::AdaptiveRemeshingSettings& remeshingSettings)
{
	// Step 1: Compute local densities
	const auto densities = ComputeLocalDensities(curve, (remeshingSettings.MinEdgeLength + remeshingSettings.MaxEdgeLength) / 2.0f);

	// Step 2: Compute mean density and its variance
	const float meanDensity = ComputeMeanDensity(densities);
	const float densityVariance = ComputeDensityVariance(densities, meanDensity);

	// Step 3: Check if density variance is within the acceptable range
	if (densityVariance > remeshingSettings.ApproxError)
	{
		return true;
	}

	// Step 4: Check if any edge lengths are outside the allowed range
	for (const auto& [vertex, density] : densities)
	{
		const float edgeLength = 1.0f / density;
		if (edgeLength < remeshingSettings.MinEdgeLength || edgeLength > remeshingSettings.MaxEdgeLength)
		{
			return true;
		}
	}

	// If the variance and edge lengths are within the acceptable range, remeshing is not necessary
	return false;
}


bool IsRemeshingNecessary(const pmp::SurfaceMesh& mesh, const pmp::AdaptiveRemeshingSettings& remeshingSettings)
{
	// TODO: fill in
}

// ------------------------------------------------------------------------------------

bool ShouldDetectFeatures(const std::vector<float>& distancePerVertexValues)
{
	const float minDist = *std::min(distancePerVertexValues.begin(), distancePerVertexValues.end());
	const float maxDist = *std::max(distancePerVertexValues.begin(), distancePerVertexValues.end());
	return minDist < 0.9f * maxDist;
	//return true;
}

/// \brief A set of time percentages for edge length adjustment for remeshing.
//static std::unordered_set<unsigned int> ADJUSTMENT_TIME_PERCENTAGES{
//	5, 10, 20, 50//, 60, 80
//};

/// \brief A set of time indices for edge length adjustment for remeshing.
static std::unordered_set<unsigned int> ADJUSTMENT_TIME_INDICES{
	3, 10, 20, 50
};

void SetRemeshingAdjustmentTimeIndices(const std::unordered_set<unsigned int>& valuesSet)
{
	ADJUSTMENT_TIME_INDICES = valuesSet;
}

bool ShouldAdjustRemeshingLengths(const unsigned int& ti /*, const unsigned int& NSteps*/)
{
	//const auto timePercentage = static_cast<unsigned int>(static_cast<float>(ti) / static_cast<float>(NSteps) * 100);
	//return ADJUSTMENT_TIME_PERCENTAGES.contains(timePercentage);
	return ADJUSTMENT_TIME_INDICES.contains(ti);
}

void AdjustRemeshingLengths(const float& decayFactor, float& minEdgeLength, float& maxEdgeLength, float& approxError)
{
	minEdgeLength *= decayFactor;
	maxEdgeLength *= decayFactor;
	approxError = 0.1f * (minEdgeLength + maxEdgeLength);
}

pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromIcoSphere(const std::shared_ptr<pmp::SurfaceMesh>& icosphere, float radius, const pmp::Point& center)
{
	if (!icosphere)
	{
		throw std::invalid_argument("CollectRemeshingSettingsFromIcoSphere: icosphere == nullptr!\n");
	}

	if (radius < FLT_EPSILON)
	{
		throw std::invalid_argument("CollectRemeshingSettingsFromIcoSphere: radius < FLT_EPSILON!\n");
	}

	pmp::AdaptiveRemeshingSettings settings;

	float minEdgeLength = std::numeric_limits<float>::max();
	float maxEdgeLength = std::numeric_limits<float>::lowest();
	float totalEdgeLength = 0.0f;
	float maxDeviation = 0.0f;

	// Calculate edge lengths
	for (const auto e : icosphere->edges())
	{
		float edgeLength = icosphere->edge_length(e);
		minEdgeLength = std::min(minEdgeLength, edgeLength);
		maxEdgeLength = std::max(maxEdgeLength, edgeLength);
		totalEdgeLength += edgeLength;
	}

	// Calculate the maximum deviation using the centroids of the faces
	for (const auto f : icosphere->faces())
	{
		const pmp::Point faceCentroid = centroid(*icosphere, f);
		const float distanceToCenter = norm(faceCentroid - center);

		// Calculate the deviation from the ideal spherical radius
		float deviation = std::abs(distanceToCenter - radius);
		maxDeviation = std::max(maxDeviation, deviation);
	}

	settings.MinEdgeLength = minEdgeLength;
	settings.MaxEdgeLength = maxEdgeLength;
	settings.ApproxError = maxDeviation;  // Use the maximum deviation as the approximation error
	settings.NRemeshingIterations = 10;
	settings.NTangentialSmoothingIters = 6;
	settings.UseProjection = true;

	return settings;
}

pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromMesh(const std::shared_ptr<pmp::SurfaceMesh>& mesh)
{
	if (!mesh)
	{
		throw std::invalid_argument("CollectRemeshingSettingsFromMesh: mesh == nullptr!\n");
	}

	pmp::AdaptiveRemeshingSettings settings;

	float minEdgeLength = std::numeric_limits<float>::max();
	float maxEdgeLength = std::numeric_limits<float>::lowest();

	// Calculate edge lengths
	for (const auto e : mesh->edges())
	{
		float edgeLength = mesh->edge_length(e);
		minEdgeLength = std::min(minEdgeLength, edgeLength);
		maxEdgeLength = std::max(maxEdgeLength, edgeLength);
	}

	// Calculate the maximum quadric approximation error across all vertices
	float maxDeviation = -FLT_MAX;
	for (const auto v : mesh->vertices())
	{
		float deviation = Geometry::CalculateQuadricApproximationErrorAtVertex(*mesh, v);
		maxDeviation = std::max(maxDeviation, deviation);
	}

	settings.MinEdgeLength = minEdgeLength;
	settings.MaxEdgeLength = maxEdgeLength;
	settings.ApproxError = maxDeviation;  // Use the maximum quadric approximation error as the approximation error
	settings.NRemeshingIterations = 10;
	settings.NTangentialSmoothingIters = 6;
	settings.UseProjection = true;

	return settings;
}

pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromCircleCurve(const std::shared_ptr<pmp::ManifoldCurve2D>& circlePolyline, float radius, const pmp::Point2& center)
{
	if (!circlePolyline)
	{
		throw std::invalid_argument("CollectRemeshingSettingsFromCircleCurve: circlePolyline == nullptr!\n");
	}

	if (radius < FLT_EPSILON)
	{
		throw std::invalid_argument("CollectRemeshingSettingsFromCircleCurve: radius < FLT_EPSILON!\n");
	}

	pmp::AdaptiveRemeshingSettings settings;

	float minEdgeLength = std::numeric_limits<float>::max();
	float maxEdgeLength = std::numeric_limits<float>::lowest();
	float totalEdgeLength = 0.0f;
	float maxDeviation = 0.0f;

	// Calculate edge lengths
	for (const auto e : circlePolyline->edges())
	{
		float edgeLength = circlePolyline->edge_length(e);
		minEdgeLength = std::min(minEdgeLength, edgeLength);
		maxEdgeLength = std::max(maxEdgeLength, edgeLength);
		totalEdgeLength += edgeLength;
	}

	// Calculate the maximum deviation using the centroids of the edges
	for (const auto e : circlePolyline->edges())
	{
		const pmp::Point2 faceCentroid = centroid(*circlePolyline, e);
		const float distanceToCenter = norm(faceCentroid - center);

		// Calculate the deviation from the ideal spherical radius
		float deviation = std::abs(distanceToCenter - radius);
		maxDeviation = std::max(maxDeviation, deviation);
	}

	settings.MinEdgeLength = minEdgeLength;
	settings.MaxEdgeLength = maxEdgeLength;
	settings.ApproxError = maxDeviation;  // Use the maximum deviation as the approximation error
	settings.NRemeshingIterations = 5;
	settings.NTangentialSmoothingIters = 4;
	settings.UseProjection = true;

	return settings;
}

pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromCurve(const std::shared_ptr<pmp::ManifoldCurve2D>& curve)
{
	if (!curve)
	{
		throw std::invalid_argument("CollectRemeshingSettingsFromCurve: mesh == nullptr!\n");
	}

	pmp::AdaptiveRemeshingSettings settings;

	float minEdgeLength = std::numeric_limits<float>::max();
	float maxEdgeLength = std::numeric_limits<float>::lowest();

	// Calculate edge lengths
	for (const auto e : curve->edges())
	{
		float edgeLength = curve->edge_length(e);
		minEdgeLength = std::min(minEdgeLength, edgeLength);
		maxEdgeLength = std::max(maxEdgeLength, edgeLength);
	}

	// Calculate the maximum quadric approximation error across all vertices
	float maxDeviation = -FLT_MAX;
	for (const auto v : curve->vertices())
	{
		float deviation = Geometry::CalculateCircularApproximationErrorAtVertex(*curve, v);
		maxDeviation = std::max(maxDeviation, deviation);
	}

	settings.MinEdgeLength = minEdgeLength;
	settings.MaxEdgeLength = maxEdgeLength;
	settings.ApproxError = maxDeviation;  // Use the maximum circular approximation error as the approximation error
	settings.NRemeshingIterations = 10;
	settings.NTangentialSmoothingIters = 6;
	settings.UseProjection = true;

	return settings;
}
