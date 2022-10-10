#pragma once

#include "pmp/SurfaceMesh.h"
#include "pmp/algorithms/DifferentialGeometry.h"

#include "geometry/Grid.h"

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
	float MinEdgeMultiplier{ 0.17f }; //>! multiplier for minimum edge length in adaptive remeshing.
	double RemeshingStartTimeFactor{ 0.1 }; //>! the fraction of total time steps after which remeshing should take place.
	float EdgeLengthDecayFactor{ 0.97f }; //>! decay factor for minimum (and consequently maximum) edge length.
	double RemeshingSizeDecayStartTimeFactor{ 0.2 }; //>! decay of edge length bounds should take place after (this value) * NSteps of evolution.
	unsigned int StepStrideForEdgeDecay{ 5 }; //>! the number of steps after which edge length bound decay takes place.
	double FeatureDetectionStartTimeFactor{ 0.4 }; //>! feature detection becomes relevant after (this value) * NSteps.
	unsigned int NRemeshingIters{ 2 }; //>! the number of iterations for pmp::Remeshing.
	unsigned int NTanSmoothingIters{ 5 }; //>! the number of tangential smoothing iterations for pmp::Remeshing.
	bool UseBackProjection{ true }; //>! if true surface kd-tree back-projection will be used for pmp::Remeshing.

	FeatureDetectionType FeatureType{ FeatureDetectionType::MeanCurvature }; //>! type of feature detection function.
	double MinDihedralAngle{ 1.0 * M_PI_2 * 180.0 }; //>! critical dihedral angle for feature detection
	double MaxDihedralAngle{ 2.0 * M_PI_2 * 180.0 }; //>! critical dihedral angle for feature detection
	float PrincipalCurvatureFactor{ 4.0f }; //>! vertices with |Kmax| > \p principalCurvatureFactor * |Kmin| are marked as feature.
	float CriticalMeanCurvatureAngle{ 1.0f * static_cast<float>(M_PI_2) }; //>! vertices with curvature angles smaller than this value are feature vertices. 
	bool ExcludeEdgesWithoutBothFeaturePts{ false }; //>! if true, edges with only one vertex detected as feature will not be marked as feature.
};

/// \brief An enumerator for the choice of mesh Laplacian scheme [Meyer, Desbrun, Schroder, Barr, 2003].
enum class [[nodiscard]] MeshLaplacian
{
	Voronoi = 0, //>! the finite volume for mesh Laplacian is generated from a true Voronoi neighborhood of each vertex.
	Barycentric = 1 //>! the finite volume for mesh Laplacian is generated from face barycenters.
};

/// \brief a list of triangle metrics to be computed.
using TriangleMetrics = std::vector<std::string>;

/**
 * \brief A wrapper for surface evolution settings.
 * \struct SurfaceEvolutionSettings
 */
struct SurfaceEvolutionSettings
{
	std::string ProcedureName{}; //>! name for the evolution procedure.

	unsigned int NSteps{ 20 };   //>! number of time steps for surface evolution.
	double TimeStep{ 0.01 };     //>! time step size.
	double FieldIsoLevel{ 0.0 }; //>! target level of the scalar field (e.g. zero distance to target mesh).
	unsigned int IcoSphereSubdivisionLevel{ 3 }; //>! subdivision level of an evolving ico sphere surface.

	AdvectionDiffusionParameters ADParams{}; //>! parameters for the advection-diffusion model.
	MeshTopologySettings TopoParams{}; //>! parameters for mesh topology adjustments.

	float MinTargetSize{ 1.0f }; //>! minimum size of the target mesh bounding box.
	float MaxTargetSize{ 1.0f }; //>! maximum size of the target mesh bounding box.
	pmp::vec3 TargetOrigin{}; //>! origin of the evolution's target.

	bool ExportSurfacePerTimeStep{ false }; //>! whether to export evolving surface for each time step.
	bool ExportResultSurface{ true }; //>! whether to export resulting evolving surface.
	std::string OutputPath{}; //>! path where output surfaces are to be exported.

	MeshLaplacian LaplacianType{}; //>! type of mesh Laplacian.
	TriangleMetrics TriMetrics{}; //>! list of triangle metrics to be computed.

	float TangentialVelocityWeight{ 0.0f }; //>! the weight of tangential velocity update vector for each time step.
	bool DoRemeshing{ true }; //>! if true, adaptive remeshing will be performed after the first 10-th of time steps.
	bool DoFeatureDetection{ true }; //>! if true, feature detection will take place prior to remeshing.
	bool IdentityForBoundaryVertices{ true }; //>! if true, boundary vertices give rise to: updated vertex = previous vertex.
	bool IdentityForFeatureVertices{ false }; //>! if true, feature vertices give rise to: updated vertex = previous vertex.

	double MaxFractionOfVerticesOutOfBounds{ 0.02 }; //>! fraction of vertices allowed to be out of bounds (because it will be decimated).
};

/**
 * \brief A utility for evolving surfaces within a scalar field.
 * \class SurfaceEvolver
 */
class SurfaceEvolver
{
public:
	/**
	 * \brief Constructor. Initialize with a given scalar field environment.
	 * \param field                    pre-computed or loaded scalar field environment.
	 * \param fieldExpansionFactor     the factor by which target bounds are expanded (multiplying original bounds min dimension).
	 * \param settings                 surface evolution settings.
	 */
	SurfaceEvolver(const Geometry::ScalarGrid& field, const float& fieldExpansionFactor, const SurfaceEvolutionSettings& settings);

	/**
	 * \brief Main functionality.
	 */
	void Evolve();

private:
	/**
	 * \brief Preprocess for evolution, i.e.: generate m_EvolvingSurface, and transform both m_Field and m_EvolvingSurface for stabilization.
	 */
	void Preprocess();

	// ----------------------------------------------------------------

	/** // TODO: it's a design smell. Try to use a function object with a relevant number of parameters, or a variable-parameter function.
	 * \brief Performs a feature detection function based on given type.
	 * \param type    feature detection function type.
	 * \return number of edges detected as features.
	 */
	[[nodiscard]] size_t DetectFeatures(const FeatureDetectionType& type) const;

	/**
	 * \brief Weight function for Laplacian flow term, inspired by [Huska, Medla, Mikula, Morigi 2021].
	 * \param distanceAtVertex          the value of distance from evolving mesh vertex to target mesh.
	 * \return weight function value.
	 */
	[[nodiscard]] double LaplacianDistanceWeightFunction(const double& distanceAtVertex) const;

	/**
	 * \brief Weight function for advection flow term, inspired by [Huska, Medla, Mikula, Morigi 2021].
	 * \param distanceAtVertex          the value of distance from evolving mesh vertex to target mesh.
	 * \param negDistanceGradient       negative gradient vector of distance field at vertex position.
	 * \param vertexNormal              unit normal to vertex.
	 * \return weight function value.
	 */
	[[nodiscard]] double AdvectionDistanceWeightFunction(const double& distanceAtVertex,
		const pmp::dvec3& negDistanceGradient, const pmp::Point& vertexNormal) const;

	// ----------------------------------------------------------------

	/**
	 * \brief Writes m_EvolvingSurface using m_OutputMeshExtension.
	 * \param tId                     index of the current time step.
	 * \paran isResult                if true, a different "connecting name" is chosen for resulting surface after all time steps are completed.
	 * \param transformToOriginal     if true, m_TransformToOriginal matrix is used to transform stabilized geometry back to original.
	 */
	void ExportSurface(const unsigned int& tId, const bool& isResult = false, const bool& transformToOriginal = true) const;

	/**
	 * \brief Computes triangle metrics interpolated to vertices according to list m_EvolSettings.TriMetrics.
	 */
	void ComputeTriangleMetrics() const;

	// ----------------------------------------------------------------

	SurfaceEvolutionSettings m_EvolSettings{}; //>! settings.

	std::shared_ptr<Geometry::ScalarGrid> m_Field{ nullptr }; //>! scalar field environment.
	std::shared_ptr<pmp::SurfaceMesh> m_EvolvingSurface{ nullptr }; //>! (stabilized) evolving surface.

	float m_ExpansionFactor{ 0.0f }; //>! the factor by which target bounds are expanded (multiplying original bounds min dimension).
	pmp::Scalar m_StartingSurfaceRadius{ 1.0f }; //>! radius of the starting surface.
	pmp::Scalar m_ScalingFactor{ 1.0f }; //>! stabilization scaling factor value.

	std::function<pmp::ImplicitLaplaceInfo(const pmp::SurfaceMesh&, pmp::Vertex)> m_ImplicitLaplacianFunction{}; //>! a Laplacian function chosen from parameter MeshLaplacian.
	std::function<double(const pmp::SurfaceMesh&, pmp::Vertex)> m_LaplacianAreaFunction{}; //>! a Laplacian area function chosen from parameter MeshLaplacian.
	std::function<size_t(const pmp::Scalar&, const bool&)> m_FeatureFunction{}; //>! a function for detecting mesh features.

	// export
	std::string m_OutputMeshExtension = ".vtk"; //>! extension of the exported mesh geometry.
	pmp::mat4 m_TransformToOriginal{}; //>! transformation matrix from stabilized surface to original size (for export).

};

/**
 * \brief Reports SurfaceEvolver's input to a given stream.
 * \param evolSettings    settings for SurfaceEvolver.
 * \param os              output stream.
 */
void ReportInput(const SurfaceEvolutionSettings& evolSettings, std::ostream& os);

/**
 * \brief Precomputes parameters for advection-diffusion model within SurfaceEvolver.
 * \param distanceMax             maximum effective distance from target (affects diffusion term weight).
 * \param targetMinDimension      minimum target size (affects diffusion term variance).
 * \return desired advection-diffusion parameters.
 */
[[nodiscard]] AdvectionDiffusionParameters PreComputeAdvectionDiffusionParams(const double& distanceMax, const double& targetMinDimension);
