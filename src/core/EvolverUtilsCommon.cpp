#include "EvolverUtilsCommon.h"

#include <numeric>
#include <ranges>

#include "pmp/SurfaceMesh.h"
#include "pmp/ManifoldCurve2D.h"
#include "pmp/algorithms/DifferentialGeometry.h"
#include "pmp/algorithms/Features.h"

#include "geometry/IcoSphereBuilder.h"

//#include <Spectra/SymEigsSolver.h>
//#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Eigen/SparseCore>



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
	stats.Mean /= static_cast<pmp::Scalar>(mesh.n_vertices());
	return stats;
}

/// \brief The power of the stabilizing scale factor.
constexpr pmp::Scalar SCALE_FACTOR_POWER = 1.0 / 2.0;
/// \brief the reciprocal value of how many times the surface area element shrinks during evolution.
constexpr pmp::Scalar INV_SHRINK_FACTOR = 5.0;

pmp::Scalar GetStabilizationScalingFactor(const double& timeStep, const pmp::Scalar& icoRadius, const unsigned int& icoSubdiv, const pmp::Scalar& stabilizationFactor)
{
	const unsigned int expectedVertexCount = (N_ICO_EDGES_0 * static_cast<unsigned int>(pow(4, icoSubdiv) - 1) + 3 * N_ICO_VERTS_0) / 3;
	const pmp::Scalar weighedIcoRadius = icoRadius;
	const pmp::Scalar expectedMeanCoVolArea = stabilizationFactor * (4.0 * static_cast<pmp::Scalar>(M_PI) * weighedIcoRadius * weighedIcoRadius / static_cast<pmp::Scalar>(expectedVertexCount));
	return pow(static_cast<pmp::Scalar>(timeStep) / expectedMeanCoVolArea * INV_SHRINK_FACTOR, SCALE_FACTOR_POWER);
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

pmp::vec3 ComputeTangentialUpdateVelocityAtVertex(const pmp::SurfaceMesh& mesh, const pmp::Vertex& v, const pmp::vec3& vNormal, const pmp::Scalar& weight)
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
		result += (1.0 + eDot) * (e0 + e1);
	}
	result *= weight / static_cast<pmp::Scalar>(mesh.valence(v));
	const auto resultDotNormal = pmp::dot(result, vNormal);
	return (result - resultDotNormal * vNormal);
}

pmp::vec2 ComputeTangentialUpdateVelocityAtVertex(const pmp::ManifoldCurve2D& curve, const pmp::Vertex& v, const pmp::vec2& vNormal, const pmp::Scalar& weight)
{
	pmp::vec2 result{};
	if (!curve.has_vertex_property("v:normal") || curve.is_isolated(v) || curve.is_boundary(v))
		return result;

	const auto [vPrev, vNext] = curve.vertices(v);
	const auto e0Vec = curve.position(vPrev) - curve.position(v);
	const auto e1Vec = curve.position(vNext) - curve.position(v);
	result += 0.5 * (e0Vec + e1Vec);
	const auto resultDotNormal = pmp::dot(result, vNormal);
	return (result - resultDotNormal * vNormal);
}

void EvaluateCurvatureBasedFeatures(pmp::SurfaceMesh& mesh, const pmp::Scalar& curvatureAngle, const pmp::Scalar& curvatureFactor, const bool& excludeEdgesWithoutBothFeaturePts)
{
	pmp::Features featuresDetector(mesh);
	featuresDetector.detect_vertices_with_high_curvature(curvatureAngle, curvatureFactor, excludeEdgesWithoutBothFeaturePts);
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

bool IsRemeshingNecessary(const std::vector<pmp::Scalar>& equilateralJacobianConditionNumbers)
{
	pmp::Scalar minVal = FLT_MAX;
	pmp::Scalar maxVal = -FLT_MAX;
	for (const auto& val : equilateralJacobianConditionNumbers)
	{
		if (val < minVal) minVal = val;
		if (val > maxVal) maxVal = val;
	}
	return !((minVal > Geometry::JACOBIAN_COND_MIN && minVal < Geometry::JACOBIAN_COND_MAX) &&
		     (maxVal > Geometry::JACOBIAN_COND_MIN && maxVal < Geometry::JACOBIAN_COND_MAX));
}

bool IsNonFeatureRemeshingNecessary(const pmp::SurfaceMesh& mesh)
{
	pmp::Scalar minVal = FLT_MAX;
	pmp::Scalar maxVal = -FLT_MAX;
	const auto vQualityProp = mesh.get_vertex_property<pmp::Scalar>("v:equilateralJacobianCondition");
	const auto vIsFeature = mesh.get_vertex_property<pmp::Scalar>("v:feature");

	for (const auto v : mesh.vertices())
	{
		if (vIsFeature[v])
			continue;

		const auto& val = vQualityProp[v];
		if (val < minVal) minVal = val;
		if (val > maxVal) maxVal = val;
	}
	return !((minVal > Geometry::JACOBIAN_COND_MIN && minVal < Geometry::JACOBIAN_COND_MAX) &&
		(maxVal > Geometry::JACOBIAN_COND_MIN && maxVal < Geometry::JACOBIAN_COND_MAX));
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

bool IsRemeshingNecessary(const pmp::ManifoldCurve2D& curve, const pmp::AdaptiveRemeshingSettings& remeshingSettings)
{	
	for (const auto e : curve.edges())
	{
		const auto edgeLength = curve.edge_length(e);
		if (edgeLength < remeshingSettings.MinEdgeLength || edgeLength > remeshingSettings.MaxEdgeLength)
			return true;
	}
	return false;
}

bool IsRemeshingNecessary(const pmp::SurfaceMesh& mesh, const pmp::AdaptiveRemeshingSettings& remeshingSettings, const AreaFunction& areaFunc)
{
	for (const auto e : mesh.edges())
	{
		const auto edgeLength = mesh.edge_length(e);
		if (edgeLength < remeshingSettings.MinEdgeLength || edgeLength > remeshingSettings.MaxEdgeLength)
			return true;
	}
	return false;
}

// ------------------------------------------------------------------------------------

bool IsRemeshingNecessary(const pmp::SurfaceMesh& mesh, const Geometry::FaceQualityFunction& qualityFunc, const Geometry::FaceQualityRange& qualityRange)
{	
	return std::ranges::any_of(mesh.faces(), [&](const auto& f)
		{
			const auto qualityVal = qualityFunc(mesh, f);
			return !qualityRange(qualityVal);
		});
}

// ------------------------------------------------------------------------------------

bool ShouldDetectFeatures(const std::vector<pmp::Scalar>& distancePerVertexValues)
{
	const pmp::Scalar minDist = *std::min(distancePerVertexValues.begin(), distancePerVertexValues.end());
	const pmp::Scalar maxDist = *std::max(distancePerVertexValues.begin(), distancePerVertexValues.end());
	return minDist < 0.9 * maxDist;
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

std::unordered_set<unsigned int>& GetRemeshingAdjustmentTimeIndices()
{
	return ADJUSTMENT_TIME_INDICES;
}

bool ShouldAdjustRemeshingLengths(const unsigned int& ti /*, const unsigned int& NSteps*/)
{
	//const auto timePercentage = static_cast<unsigned int>(static_cast<pmp::Scalar>(ti) / static_cast<pmp::Scalar>(NSteps) * 100);
	//return ADJUSTMENT_TIME_PERCENTAGES.contains(timePercentage);
	return ADJUSTMENT_TIME_INDICES.contains(ti);
}

void AdjustRemeshingLengths(const pmp::Scalar& decayFactor, pmp::Scalar& minEdgeLength, pmp::Scalar& maxEdgeLength, pmp::Scalar& approxError)
{
	minEdgeLength *= decayFactor;
	maxEdgeLength *= decayFactor;
	approxError = 0.1 * (minEdgeLength + maxEdgeLength);
}

pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromIcoSphere_OLD(unsigned int subdiv, pmp::Scalar radius, const ManifoldAdaptiveRemeshingParams& remeshingParams)
{
	if (radius < FLT_EPSILON)
	{
		throw std::invalid_argument("CollectRemeshingSettingsFromIcoSphere_OLD: radius < FLT_EPSILON!\n");
	}

	pmp::AdaptiveRemeshingSettings settings;

	constexpr pmp::Scalar baseIcoHalfAngle = 2.0 * static_cast<pmp::Scalar>(M_PI) / 10.0;
	settings.MinEdgeLength = remeshingParams.MinEdgeMultiplier * 2.0 * radius * 
		sin(baseIcoHalfAngle * pow(2.0, -1.0 * static_cast<pmp::Scalar>(subdiv))); // from icosahedron edge length
	settings.MaxEdgeLength = 4.0 * settings.MinEdgeLength;
	settings.ApproxError = 0.25 * (settings.MinEdgeLength + settings.MaxEdgeLength);

	settings.NRemeshingIterations = remeshingParams.NRemeshingIters;
	settings.NTangentialSmoothingIters = remeshingParams.NTanSmoothingIters;
	settings.UseProjection = remeshingParams.UseBackProjection;

	return settings;
}

pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromIcoSphere(const std::shared_ptr<pmp::SurfaceMesh>& icosphere, pmp::Scalar radius, const pmp::Point& center)
{
	if (!icosphere)
	{
		throw std::invalid_argument("CollectRemeshingSettingsFromIcoSphere: icosphere == nullptr!\n");
	}

	if (radius < FLT_EPSILON)
	{
		throw std::invalid_argument("CollectRemeshingSettingsFromIcoSphere: radius < FLT_EPSILON!\n");
	}

	// TODO: debug this experimental implementation

	pmp::AdaptiveRemeshingSettings settings;

	pmp::Scalar minEdgeLength = std::numeric_limits<pmp::Scalar>::max();
	pmp::Scalar maxEdgeLength = std::numeric_limits<pmp::Scalar>::lowest();
	pmp::Scalar totalEdgeLength = 0.0;
	pmp::Scalar maxDeviation = 0.0;

	// Calculate edge lengths
	for (const auto e : icosphere->edges())
	{
		pmp::Scalar edgeLength = icosphere->edge_length(e);
		minEdgeLength = std::min(minEdgeLength, edgeLength);
		maxEdgeLength = std::max(maxEdgeLength, edgeLength);
		totalEdgeLength += edgeLength;
	}

	// Calculate the maximum deviation using the centroids of the faces
	for (const auto f : icosphere->faces())
	{
		const pmp::Point faceCentroid = centroid(*icosphere, f);
		const pmp::Scalar distanceToCenter = norm(faceCentroid - center);

		// Calculate the deviation from the ideal spherical radius
		pmp::Scalar deviation = std::abs(distanceToCenter - radius);
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

	pmp::Scalar minEdgeLength = std::numeric_limits<pmp::Scalar>::max();
	pmp::Scalar maxEdgeLength = std::numeric_limits<pmp::Scalar>::lowest();

	// Calculate edge lengths
	for (const auto e : mesh->edges())
	{
		pmp::Scalar edgeLength = mesh->edge_length(e);
		minEdgeLength = std::min(minEdgeLength, edgeLength);
		maxEdgeLength = std::max(maxEdgeLength, edgeLength);
	}

	// Calculate the maximum quadric approximation error across all vertices
	pmp::Scalar maxDeviation = -FLT_MAX;
	for (const auto v : mesh->vertices())
	{
		pmp::Scalar deviation = Geometry::CalculateQuadricApproximationErrorAtVertex(*mesh, v);
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

pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromCircleCurve(const std::shared_ptr<pmp::ManifoldCurve2D>& circlePolyline, pmp::Scalar radius, const pmp::Point2& center)
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

	pmp::Scalar minEdgeLength = std::numeric_limits<pmp::Scalar>::max();
	pmp::Scalar maxEdgeLength = std::numeric_limits<pmp::Scalar>::lowest();
	pmp::Scalar totalEdgeLength = 0.0;
	pmp::Scalar maxDeviation = 0.0;

	// Calculate edge lengths
	for (const auto e : circlePolyline->edges())
	{
		pmp::Scalar edgeLength = circlePolyline->edge_length(e);
		minEdgeLength = std::min(minEdgeLength, edgeLength);
		maxEdgeLength = std::max(maxEdgeLength, edgeLength);
		totalEdgeLength += edgeLength;
	}

	// Calculate the maximum deviation using the centroids of the edges
	for (const auto e : circlePolyline->edges())
	{
		const pmp::Point2 faceCentroid = centroid(*circlePolyline, e);
		const pmp::Scalar distanceToCenter = norm(faceCentroid - center);

		// Calculate the deviation from the ideal spherical radius
		pmp::Scalar deviation = std::abs(distanceToCenter - radius);
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
		throw std::invalid_argument("CollectRemeshingSettingsFromCurve: curve == nullptr!\n");
	}

	pmp::AdaptiveRemeshingSettings settings;

	pmp::Scalar minEdgeLength = std::numeric_limits<pmp::Scalar>::max();
	pmp::Scalar maxEdgeLength = std::numeric_limits<pmp::Scalar>::lowest();

	// Calculate edge lengths
	for (const auto e : curve->edges())
	{
		pmp::Scalar edgeLength = curve->edge_length(e);
		minEdgeLength = std::min(minEdgeLength, edgeLength);
		maxEdgeLength = std::max(maxEdgeLength, edgeLength);
	}

	// Calculate the maximum quadric approximation error across all vertices
	pmp::Scalar maxDeviation = -FLT_MAX;
	for (const auto v : curve->vertices())
	{
		pmp::Scalar deviation = Geometry::CalculateCircularApproximationErrorAtVertex(*curve, v);
		maxDeviation = std::max(maxDeviation, deviation);
	}

	settings.MinEdgeLength = minEdgeLength;
	settings.MaxEdgeLength = maxEdgeLength;
	settings.ApproxError = maxDeviation;  // Use the maximum circular approximation error as the approximation error
	settings.NRemeshingIterations = 2;
	settings.NTangentialSmoothingIters = 5;
	settings.UseProjection = true;

	return settings;
}

//
// ----------------------------------------------------------------------------------------
//

std::pair<
	std::optional<pmp::VertexProperty<pmp::Vertex>>, 
	std::optional<pmp::VertexProperty<pmp::Vertex>>> 
	GetNearestGapBoundaryVertices(
		pmp::ManifoldCurve2D& curve, 
		const std::shared_ptr<Geometry::ScalarGrid2D>& targetDistanceField,
		const std::vector<std::shared_ptr<Geometry::ScalarGrid2D>>& manifoldDistanceFields, 
		const ScalarGridInterpolationFunction2D& interpFunc, const NormalActivationSettings& settings)
{
	std::optional<pmp::VertexProperty<pmp::Vertex>> vForwardGapBoundary{ std::nullopt };
	std::optional<pmp::VertexProperty<pmp::Vertex>> vBackwardGapBoundary{ std::nullopt };

	if (curve.is_empty() || !settings.On)
		return { vForwardGapBoundary, vBackwardGapBoundary };

	// check for boundaries
	if (!curve.is_closed())
		return { vForwardGapBoundary, vBackwardGapBoundary };

	auto vGap0 = Geometry::GetVerticesWithinMinDistance(curve, { targetDistanceField },
		settings.ManifoldCriticalRadius, "v:target_activated", interpFunc);
	auto vGap1 = Geometry::GetVerticesWithinMinDistance(curve, manifoldDistanceFields,
		settings.ManifoldCriticalRadius, "v:manifold_activated", interpFunc);

	auto vGap = !curve.has_vertex_property("v:gap_activated") ?
		curve.vertex_property<bool>("v:gap_activated", false) : curve.get_vertex_property<bool>("v:gap_activated");
	for (const auto v : curve.vertices())
		vGap[v] = !vGap0[v] && vGap1[v];

	if (std::ranges::all_of(vGap.vector(), [](const auto& item) { return !item; }) ||
		std::ranges::all_of(vGap.vector(), [](const auto& item) { return item; }))
	{
		// no transition points
		return { vForwardGapBoundary, vBackwardGapBoundary };
	}

	return GetNearestGapBoundaryVertices(curve, vGap, settings);
}

std::pair<
	std::optional<pmp::VertexProperty<pmp::Vertex>>,
	std::optional<pmp::VertexProperty<pmp::Vertex>>>
	GetNearestGapBoundaryVertices(pmp::ManifoldCurve2D& curve,
		const pmp::VertexProperty<bool>& vGap, const NormalActivationSettings& settings)
{
	std::optional<pmp::VertexProperty<pmp::Vertex>> vNextBound{ std::nullopt };
	std::optional<pmp::VertexProperty<pmp::Vertex>> vPrevBound{ std::nullopt };

	if (curve.is_empty() || !settings.On)
		return { vNextBound, vPrevBound };

	// check for boundaries
	if (!curve.is_closed())
		return { vNextBound, vPrevBound };

	// collect vertices into a flat vector V[0..N-1]
	std::vector<pmp::Vertex> V;
	V.reserve(curve.n_vertices());
	{
		pmp::Vertex v = pmp::Vertex(0);
		do {
			V.push_back(v);
			auto e = curve.edge_from(v);
			v = curve.to_vertex(e);
		} while (v != V.front());
	}
	int N = int(V.size());
	unsigned int n = settings.NPointsFromCriticalBound;

	// extend the gap flags by the stride n
	std::vector<bool> inGap(N, 0);
	for (int i = 0; i < N; ++i) {
		if (vGap[V[i]]) {
			// mark [i-n .. i+n] mod N
			for (int d = -int(n); d <= int(n); ++d) {
				int j = (i + d + N) % N;
				inGap[j] = true;
			}
		}
	}

	// make or get our properties
	vPrevBound = !curve.has_vertex_property("v:prev_boundary") ?
		curve.vertex_property<pmp::Vertex>("v:prev_boundary", pmp::Vertex{}) : curve.get_vertex_property<pmp::Vertex>("v:prev_boundary");
	vNextBound = !curve.has_vertex_property("v:next_boundary") ?
		curve.vertex_property<pmp::Vertex>("v:next_boundary", pmp::Vertex{}) : curve.get_vertex_property<pmp::Vertex>("v:next_boundary");

	// 4) Scan inGap to build a list of segments [s,e]
	std::vector<std::pair<int, int>> segs;
	{
		int i = 0;
		while (i < N) {
			while (i < N && !inGap[i]) ++i;
			if (i >= N) break;
			int s = i;
			while (i < N && inGap[i]) ++i;
			int e = i - 1;
			segs.emplace_back(s, e);
		}
		if (segs.size() >= 2) {
			auto& first = segs.front();
			auto& last = segs.back();
			if (first.first == 0 && last.second == N - 1) {
				first.second = first.second + N;
				first.first = last.first;
				segs.pop_back();
			}
		}
	}

	// 5) Allocate result buffers and fill by segment
	std::vector<pmp::Vertex> prevB(N), nextB(N);
	for (auto se : segs) {
		int s = se.first;
		int e = se.second;
		int sMod = (s % N + N) % N;
		int eMod = (e % N + N) % N;
		int prevIdx = (sMod - int(n) + N) % N;
		int nextIdx = (eMod + int(n)) % N;
		for (int j = s; j <= e; ++j) {
			int idx = (j % N + N) % N;
			prevB[idx] = V[prevIdx];
			nextB[idx] = V[nextIdx];
		}
	}

	// 6) Write back into the curve properties
	for (int i = 0; i < N; ++i) {
		(*vPrevBound)[V[i]] = prevB[i];
		(*vNextBound)[V[i]] = nextB[i];
	}

	return { vNextBound, vPrevBound };
}

ImplicitBezierCurveVertexInfo CalculateBezierVertexInfo(
	const pmp::Point2& vNormalPrev, const pmp::Point2& vNormalNext, 
	const pmp::Scalar& vCurrentArcLength, const pmp::Scalar& vPrevArcLength, const pmp::Scalar& vNextArcLength,
	const NormalActivationSettings& settings)
{
	ImplicitBezierCurveVertexInfo result;

	// get tangents by rotating the normals by +/- 90 degrees
	const auto perpendicularRotPositive = pmp::rotation_matrix_2d(pmp::Point2{}, M_PI_2);
	const auto perpendicularRotNegative = pmp::inverse(perpendicularRotPositive);

	const auto T0 = pmp::affine_transform(perpendicularRotPositive, vNormalPrev);
	const auto T1 = pmp::affine_transform(perpendicularRotNegative, vNormalNext);

	// compute normalized parameter s in [0,1]
	double totalLen = static_cast<double>(vNextArcLength - vPrevArcLength);
	double s = 0.0;
	if (totalLen > 0.0)
	{
		s = (static_cast<double>(vCurrentArcLength - vPrevArcLength) / totalLen);
		s = std::clamp(s, 0.0, 1.0);
	}

	// cubic Bezier basis functions
	double o = 1.0 - s;
	double alpha = o * o * o + 3.0 * s * o * o;   // (1-s)^3 + 3 s (1-s)^2
	double beta = 3.0 * s * s * o + s * s * s;   // 3 s^2 (1-s) + s^3

	// build the two interior control point offsets 
	double h = settings.BezierWeightCoefficient * totalLen; // how far apart the tangents reach, per settings
	double gamma = 3.0 * s * o * o * h;  // weight for T0
	double delta = 3.0 * s * s * o * h;  // weight for T1

	const auto C = gamma * T0 - delta * T1;

	// soften to restore strict diagonal dominance
	double lambda = settings.BezierSofteningCoefficient;

	result.PrevGapBoundaryWeight = lambda * alpha;
	result.NextGapBoundaryWeight = lambda * beta;
	result.BezierRhs = lambda * C;

	return result;
}

