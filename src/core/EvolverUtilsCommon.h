#pragma once

#include "pmp/Types.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <unordered_set>

#include "geometry/MeshAnalysis.h"
#include "pmp/ManifoldCurve2D.h"


namespace pmp
{
	class SurfaceMesh;
	class ManifoldCurve2D;
	class Vertex;
	struct AdaptiveRemeshingSettings;
}

/// \brief a stats wrapper for co-volume measures affecting the stability of the finite volume method.
struct CoVolumeStats
{
	double Mean{ 0.0 };
	double Max{ -DBL_MAX };
	double Min{ DBL_MAX };
};

const std::string coVolMeasureVertexPropertyName{ "v:coVolumeMeasure" };

/// \brief a co-volume area evaluator.
using AreaFunction = std::function<double(const pmp::SurfaceMesh&, pmp::Vertex)>;

/**
 * \brief Analyzes the stats of co-volumes around each mesh vertex, and creates a vertex property for the measure values.
 * \param mesh             input mesh.
 * \param areaFunction     function to evaluate vertex co-volume area.
 * \return co-volume stats.
 */
[[nodiscard]] CoVolumeStats AnalyzeMeshCoVolumes(pmp::SurfaceMesh& mesh, const AreaFunction& areaFunction);

/*
 * \brief Computes scaling factor for stabilizing the finite volume method on assumed spherical surface meshes based on time step.
 * \param timeStep               time step size.
 * \param icoRadius              radius of an evolving geodesic icosahedron.
 * \param icoSubdiv              subdivision level of an evolving geodesic icosahedron.
 * \param stabilizationFactor    a multiplier for stabilizing mean co-volume area.
 * \return scaling factor for mesh and scalar grid.
 */
[[nodiscard]] float GetStabilizationScalingFactor(const double& timeStep, const float& icoRadius, const unsigned int& icoSubdiv, const float& stabilizationFactor = 1.0f);

/// \brief identifier for sparse matrix.
using SparseMatrix = Eigen::SparseMatrix<double>;

/// \brief A utility for converting Eigen::ComputationInfo to a string message.
[[nodiscard]] std::string InterpretSolverErrorCode(const Eigen::ComputationInfo& cInfo);

/**
 * \brief Computes angle-based tangential update velocity for a mesh vertex with a given weight.
 * \param mesh       a surface mesh to whom the vertex belongs.
 * \param v          vertex handle of a vertex where the tangential velocity is to be computed.
 * \param vNormal    normal to vertex with handle v.
 * \param weight     weight of the velocity vector.
 * \return tangential velocity vector.
 */
[[nodiscard]] pmp::vec3 ComputeTangentialUpdateVelocityAtVertex(const pmp::SurfaceMesh& mesh, const pmp::Vertex& v, const pmp::vec3& vNormal, const float& weight = 1.0f);

/**
 * \brief Computes tangential update velocity for a mesh vertex with a given weight.
 * \param mesh       a manifold curve to whom the vertex belongs.
 * \param v          vertex handle of a vertex where the tangential velocity is to be computed.
 * \param vNormal    normal to vertex with handle v.
 * \param weight     weight of the velocity vector.
 * \return tangential velocity vector.
 */
[[nodiscard]] pmp::vec2 ComputeTangentialUpdateVelocityAtVertex(const pmp::ManifoldCurve2D& curve, const pmp::Vertex& v, const pmp::vec2& vNormal, const float& weight = 1.0f);

// ======================================================================================================================

/// \brief An enumerator for the choice of mesh Laplacian scheme [Meyer, Desbrun, Schroder, Barr, 2003].
enum class [[nodiscard]] MeshLaplacian
{
	Voronoi = 0, //>! the finite volume for mesh Laplacian is generated from a true Voronoi neighborhood of each vertex.
	Barycentric = 1 //>! the finite volume for mesh Laplacian is generated from face barycenters.
};

/// \brief a list of triangle metrics to be computed.
using TriangleMetrics = std::vector<std::string>;

/**
 * \brief A wrapper for surface evolution advection-diffusion parameters.
 * \struct AdvectionDiffusionParameters
 */
struct AdvectionDiffusionParameters
{
	// constants for weight function:
	// C1 * (1 - exp(d^2 / C2)), where d is distance
	double MCFMultiplier{ 1.0 }; // C1
	double MCFVariance{ 1.0 }; // C2

	// constants for weight function:
	// D1 * d * ((-grad(d) . N) - D2 * sqrt(1 - (grad(d) . N)^2)), where d is distance, and N unit normal.
	double AdvectionMultiplier{ 1.0 }; // D1
	double AdvectionSineMultiplier{ 1.0 }; // D2

	bool MCFSupportPositive{ true }; //>! if true, the diffusion weight function has non-zero support for positive values only.
	bool AdvectionSupportPositive{ true }; //>! if true, the advection weight function has non-zero support for positive values only.
};

/**
 * \brief An enumerator for the type of feature detection function used during evolution.
 * \enum FeatureDetectionType
 */
enum class [[nodiscard]] FeatureDetectionType
{
	Angle = 0, //>! edges with dihedral angle larger than a given threshold value are marked as features.
	AngleWithinBounds = 1, //>! edges with dihedral angle within a given range are marked as features.
	PrincipalCurvatures = 2, //>! vertices with too much imbalance in principal curvatures are marked as features.
	MeanCurvature = 3 //>! features are preferred for vertices with high positive mean curvature.
};

/**
 * \brief A wrapper for parameters related to mesh topology adjustments (remeshing etc.).
 * \struct MeshTopologySettings
 */
struct MeshTopologySettings
{
	bool FixSelfIntersections{ true }; //>! if true, self-intersecting faces within the evolving surface will be removed, and the holes will be patched by pmp::HoleFilling.
	float MinEdgeMultiplier{ 0.14f }; //>! multiplier for minimum edge length in adaptive remeshing.
	double RemeshingStartTimeFactor{ 0.1 }; //>! the fraction of total time steps after which remeshing should take place.
	float EdgeLengthDecayFactor{ 0.98f }; //>! decay factor for minimum (and consequently maximum) edge length.
	double RemeshingSizeDecayStartTimeFactor{ 0.2 }; //>! decay of edge length bounds should take place after (this value) * NSteps of evolution.
	unsigned int StepStrideForEdgeDecay{ 5 }; //>! the number of steps after which edge length bound decay takes place.
	double FeatureDetectionStartTimeFactor{ 0.4 }; //>! feature detection becomes relevant after (this value) * NSteps.
	unsigned int NRemeshingIters{ 2 }; //>! the number of iterations for pmp::Remeshing.
	unsigned int NTanSmoothingIters{ 5 }; //>! the number of tangential smoothing iterations for pmp::Remeshing.
	bool UseBackProjection{ true }; //>! if true surface kd-tree back-projection will be used for pmp::Remeshing.

	FeatureDetectionType FeatureType{ FeatureDetectionType::MeanCurvature }; //>! type of feature detection function.
	double MinDihedralAngle{ 1.0 * M_PI_2 * 180.0 }; //>! critical dihedral angle for feature detection
	double MaxDihedralAngle{ 2.0 * M_PI_2 * 180.0 }; //>! critical dihedral angle for feature detection
	float PrincipalCurvatureFactor{ 2.0f }; //>! vertices with |Kmax| > \p principalCurvatureFactor * |Kmin| are marked as feature.
	float CriticalMeanCurvatureAngle{ 1.0f * static_cast<float>(M_PI_2) }; //>! vertices with curvature angles smaller than this value are feature vertices. 
	bool ExcludeEdgesWithoutBothFeaturePts{ false }; //>! if true, edges with only one vertex detected as feature will not be marked as feature.
};

/**
 * \brief Precomputes parameters for advection-diffusion model within (Iso)SurfaceEvolver.
 * \param distanceMax             maximum effective distance from target (affects diffusion term weight).
 * \param targetMinDimension      minimum target size (affects diffusion term variance).
 * \return desired advection-diffusion parameters.
 */
[[nodiscard]] AdvectionDiffusionParameters PreComputeAdvectionDiffusionParams(const double& distanceMax, const double& targetMinDimension);

/// \brief Evaluates whether remeshing is necessary from the co-volume stats and time step.
[[nodiscard]] bool IsRemeshingNecessary(const CoVolumeStats& stats, const double& tStep);

/// \brief Evaluates whether remeshing is necessary from the condition number metric for equilateral triangles.
[[nodiscard]] bool IsRemeshingNecessary(const std::vector<float>& equilateralJacobianConditionNumbers);

/// \brief Evaluates whether remeshing is necessary from the condition number metric for equilateral triangles that do not have a feature vertex.
[[nodiscard]] bool IsNonFeatureRemeshingNecessary(const pmp::SurfaceMesh& mesh);

/// \brief Evaluates whether remeshing is necessary from the system matrix spectral radius
[[nodiscard]] bool IsRemeshingNecessary(const SparseMatrix& lswMatrix);

/// \brief Evaluates whether remeshing is necessary for a manifold curve.
///	Uses midpoint co-volumes around each vertex.
[[nodiscard]] bool IsRemeshingNecessary(const pmp::ManifoldCurve2D& curve, const pmp::AdaptiveRemeshingSettings& remeshingSettings);

/// \brief Evaluates whether remeshing is necessary for a surface mesh based on vertex density.
///	\param[in] mesh                  input mesh.
///	\param[in] remeshingSettings     settings to be used.
///	\param[in] areaFunc              control volume measure evaluation.
///	\return true if remeshing is necessary.
[[nodiscard]] bool IsRemeshingNecessary(const pmp::SurfaceMesh& mesh, const pmp::AdaptiveRemeshingSettings& remeshingSettings, const AreaFunction& areaFunc);

/// \brief Evaluates whether remeshing is necessary for a surface mesh.
[[nodiscard]] bool IsRemeshingNecessary(const pmp::SurfaceMesh& mesh, const Geometry::FaceQualityFunction& qualityFunc, const Geometry::FaceQualityRange& qualityRange);

/// \brief A (one-time) evaluation whether the distance to target reaches a lower bound.
///	\param distancePerVertexValues    a vector of distance values on the evolving surface.
///	\return true if the conditions for feature detection are satisfied.
[[nodiscard]] bool ShouldDetectFeatures(const std::vector<float>& distancePerVertexValues);

/// \brief Sets a static container for time indices for a particular evolver setup.
void SetRemeshingAdjustmentTimeIndices(const std::unordered_set<unsigned int>& valuesSet);

/// 
/// \brief Evaluates whether the target edge lengths for adaptive remeshing should be decreased.
/// \param ti               time index from 1 to NSteps.
/// \return true if the adjustment should take place.
///
[[nodiscard]] bool ShouldAdjustRemeshingLengths(const unsigned int& ti /*, const unsigned int& NSteps*/);

///
/// \brief Adjusts edge lengths for adaptive remeshing.
///	\param decayFactor      decay factor from [0, 1] for edge length.
///	\param minEdgeLength    the minimum edge length to be adjusted.
///	\param maxEdgeLength    the maximum edge length to be adjusted.
///	\param approxError      approximation error to be adjusted.
///
void AdjustRemeshingLengths(const float& decayFactor, float& minEdgeLength, float& maxEdgeLength, float& approxError);

/**
 * \brief Logs manifolds which need remeshing.
 * \tparam ManifoldType  either pmp::ManifoldCurve2D or pmp::SurfaceMesh.
 */
template<typename ManifoldType>
class ManifoldsToRemeshTracker
{
public:
	/// \brief Add a manifold that requires remeshing
	void AddManifold(ManifoldType* manifoldPtr)
	{
		if (!manifoldPtr)
		{
			std::cerr << "ManifoldsToRemeshTracker::AddManifold: attempting to log a null manifold!\n";
			return;
		}
		m_Manifolds.push_back(manifoldPtr);
	}

	/// \brief Get all manifolds that need remeshing
	[[nodiscard]] std::vector<ManifoldType*> GetManifoldsToRemesh() const
	{
		return m_Manifolds;
	}

	/// \brief Clear all stored manifolds (reset for a new time step)
	void Reset()
	{
		m_Manifolds.clear();
	}

private:
	// TODO: use something like: bool m_IsReset{ true };
	std::vector<ManifoldType*> m_Manifolds;
};


[[nodiscard]] pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromIcoSphere_OLD(unsigned int subdiv, float radius, float minEdgeMultiplier = 0.14f);

/**
 * \brief A utility for computing the edge sizing and error limits from an icosphere mesh.
 * \param icosphere      the input icosphere surface.
 * \param radius         the radius of the smooth sphere which the icosphere's supposed to approximate.
 * \param center         the center of the smooth sphere which the icosphere's supposed to approximate.
 * \return the result AdaptiveRemeshingSettings.
 */
[[nodiscard]] pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromIcoSphere(
	const std::shared_ptr<pmp::SurfaceMesh>& icosphere,	float radius, const pmp::Point& center);

/// \brief  A utility for computing the edge sizing and error limits from an arbitrary mesh.
[[nodiscard]] pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromMesh(const std::shared_ptr<pmp::SurfaceMesh>& mesh);

/**
 * \brief A utility for computing the edge sizing and error limits from an icosphere mesh.
 * \param circlePolyline    the input circle polyline curve.
 * \param radius            the radius of the smooth circle which the circlePolyline's supposed to approximate.
 * \param center            the center of the smooth circle which the circlePolyline's supposed to approximate.
 * \return the result AdaptiveRemeshingSettings.
 */
[[nodiscard]] pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromCircleCurve(
	const std::shared_ptr<pmp::ManifoldCurve2D>& circlePolyline, float radius, const pmp::Point2& center);

/// \brief  A utility for computing the edge sizing and error limits from an arbitrary manifold curve.
[[nodiscard]] pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromCurve(const std::shared_ptr<pmp::ManifoldCurve2D>& curve);


/**
 * \brief Logs settings for manifold remeshing.
 * \tparam ManifoldType  either pmp::ManifoldCurve2D or pmp::SurfaceMesh.
 */
template<typename ManifoldType>
class ManifoldRemeshingSettingsWrapper
{
public:
	/// \brief Overload the [] operator for assigning and retrieving map values
	pmp::AdaptiveRemeshingSettings& operator[](ManifoldType* manifold)
	{
		return m_ManifoldSettings[manifold];
	}

	/// \brief Overload the [] operator for constant access (read-only)
	const pmp::AdaptiveRemeshingSettings& operator[](const ManifoldType* manifold) const
	{
		return m_ManifoldSettings.at(manifold);
	}

private:
	std::map<ManifoldType*, pmp::AdaptiveRemeshingSettings> m_ManifoldSettings{};
};

/**
 * \brief Logs initial sphere manifold settings.
 * \tparam ManifoldType  either pmp::ManifoldCurve2D or pmp::SurfaceMesh.
 * \tparam SphereType    either Sphere2D or Sphere3D.
 */
template<typename ManifoldType, typename SphereType>
class InitialSphereSettingsWrapper
{
public:
	/// \brief Overload the [] operator for assigning and retrieving map values
	SphereType& operator[](ManifoldType* manifold)
	{
		return m_ManifoldSettings[manifold];
	}

	/// \brief Overload the [] operator for constant access (read-only)
	SphereType& operator[](const ManifoldType* manifold) const
	{
		return m_ManifoldSettings.at(manifold);
	}

	/// \brief Find the maximum radius among all stored SphereType objects
	[[nodiscard]] pmp::Scalar MaxRadius() const
	{
		if (m_ManifoldSettings.empty())
		{
			return -1.0f; // Return a default value if the map is empty
		}

		// Use std::max_element to find the SphereType with the maximum radius
		auto maxIt = std::max_element(m_ManifoldSettings.begin(), m_ManifoldSettings.end(),
			[](const auto& lhs, const auto& rhs) {
				return lhs.second.Radius < rhs.second.Radius;
			});

		return maxIt->second.Radius;
	}

	/// \brief Find the minimum radius among all stored SphereType objects
	[[nodiscard]] pmp::Scalar MinRadius() const
	{
		if (m_ManifoldSettings.empty())
		{
			return -1.0f; // Return a default value if the map is empty
		}

		// Use std::max_element to find the SphereType with the maximum radius
		auto maxIt = std::min_element(m_ManifoldSettings.begin(), m_ManifoldSettings.end(),
			[](const auto& lhs, const auto& rhs) {
				return lhs.second.Radius < rhs.second.Radius;
			});

		return maxIt->second.Radius;
	}

private:
	std::map<ManifoldType*, SphereType> m_ManifoldSettings{};
};