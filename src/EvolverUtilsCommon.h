#pragma once

#include "pmp/MatVec.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>

namespace pmp
{
	class SurfaceMesh;
	class Vertex;
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

/// \brief A (one-time) evaluation whether the distance to target reaches a lower bound.
///	\param distancePerVertexValues    a vector of distance values on the evolving surface.
///	\return true if the conditions for feature detection are satisfied.
[[nodiscard]] bool ShouldDetectFeatures(const std::vector<float>& distancePerVertexValues);

/// 
/// \brief Evaluates whether the target edge lengths for adaptive remeshing should be decreased.
/// \param ti               time index from 1 to NSteps.
/// \param NSteps           the maximum number of time steps. 
/// \return true if the adjustment should take place.
///
[[nodiscard]] bool ShouldAdjustRemeshingLengths(const unsigned int& ti, const unsigned int& NSteps);

///
/// \brief Adjusts edge lengths for adaptive remeshing.
///	\param decayFactor      decay factor from [0, 1] for edge length.
///	\param minEdgeLength    the minimum edge length to be adjusted.
///	\param maxEdgeLength    the maximum edge length to be adjusted.
///	\param approxError      approximation error to be adjusted.
///
void AdjustRemeshingLengths(const float& decayFactor, float& minEdgeLength, float& maxEdgeLength, float& approxError);