#include "EvolverUtilsCommon.h"

#include "pmp/SurfaceMesh.h"
#include "geometry/IcoSphereBuilder.h"

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