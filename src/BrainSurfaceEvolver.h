#pragma once

#include "geometry/Grid.h"
#include "pmp/SurfaceMesh.h"
#include "pmp/algorithms/DifferentialGeometry.h"

// [Smith, 2002]: Smith, S. M., Fast Robust Automated Brain Extraction, Human Brain Mapping, 2002

/**
 * \brief A wrapper for the constants (with some notes of their use) from the model of [Smith, 2002].
 * \struct BE_ThresholdSettings
 */
struct BE_ThresholdSettings
{
	unsigned int MinIntensitySearchDepth{ 7 }; //>! (d1) [mm] "determines how far into the brain the minimum intensity is searched for."
	unsigned int MaxIntensitySearchDepth{ 3 }; //>! (d2) [mm] "determines how far into the brain the maximum intensity is searched for."
	double Threshold2ndPercentile; // (t2) second percentile of image intensities.
	double Threshold98thPercentile; // (t98) 98th percentile of image intensities.
	// '[t] ...attempts to distinguish between brain matter and background (because bone appears dark in most MR images, "background" is taken to include bone).'
	double ThresholdEffective; // (t) 10th of the distance from t2 to t98:  t = t2 + 0.1 * (t98 - t2),
	double ThresholdEffectiveMedian; // (tm) median intensity value evaluated from the count of voxels within t2 and t98

	double BetMainParam{ 0.826450318154212 }; // (bt) magic constant: 'bet_main_parameter' = pow(0.5, 0.275)

	// ===================================================================================================================================================================================
	// NOTES:
	// "t2, tm, and t are used to limit the effect of very dark or very bright voxels, and t is included in the maximum intensity search to limit the effect of very bright voxels."
	// Imin = max(t2, min(tm, I[0],..., I[d1])), where list { I[0],..., I[d1] } is nearest-neighbor-interpolated from evolving surface normal direction N.
	// Imax = min(tm, max(t, I[0],..., I[d2])), where list { I[0],..., I[d2] } is nearest-neighbor-interpolated from evolving surface normal direction N.
	// t1 = (Imax - t2) * bt + t2
	// f3 = 2 * (Imin - t1) / (Imax - t2)
	// ===================================================================================================================================================================================
};

/**
 * \brief A wrapper for (pre-computed) ico sphere parameters
 * \struct BE_IcoSphereSettings
 */
struct BE_IcoSphereSettings
{
	pmp::vec3 Center{};
	float Radius{ 1.0f };
};

/**
 * \brief Curvature bounds for the evolving surface.
 * \struct BE_CurvatureSettings
 */
struct BE_CurvatureSettings
{
	double MinCurvature{ -1.0 };
	double MaxCurvature{ 1.0 };
};

/**
 * \brief A wrapper for parameters related to mesh topology adjustments (remeshing etc.).
 * \struct BE_MeshTopologySettings
 */
struct BE_MeshTopologySettings
{
	float MinEdgeMultiplier{ 0.07f }; //>! multiplier for minimum edge length in adaptive remeshing.
	double RemeshingStartTimeFactor{ 0.1 }; //>! the fraction of total time steps after which remeshing should take place.
	float EdgeLengthDecayFactor{ 0.97f }; //>! decay factor for minimum (and consequently maximum) edge length.
	double RemeshingSizeDecayStartTimeFactor{ 0.2 }; //>! decay of edge length bounds should take place after (this value) * NSteps of evolution.
	unsigned int StepStrideForEdgeDecay{ 5 }; //>! the number of steps after which edge length bound decay takes place.
	double FeatureDetectionStartTimeFactor{ 0.4 }; //>! feature detection becomes relevant after (this value) * NSteps.
	unsigned int NRemeshingIters{ 2 }; //>! the number of iterations for pmp::Remeshing.
	unsigned int NTanSmoothingIters{ 5 }; //>! the number of tangential smoothing iterations for pmp::Remeshing.
	bool UseBackProjection{ false }; //>! if true surface kd-tree back-projection will be used for pmp::Remeshing.

	double MinDihedralAngle{ 0.9 * M_PI_2 * 180.0 }; //>! critical dihedral angle for feature detection
	double MaxDihedralAngle{ 1.9 * M_PI_2 * 180.0 }; //>! critical dihedral angle for feature detection
};

/// \brief An enumerator for the choice of mesh Laplacian scheme [Meyer, Desbrun, Schroder, Barr, 2003].
enum class [[nodiscard]] BE_MeshLaplacian
{
	Voronoi = 0, //>! the finite volume for mesh Laplacian is generated from a true Voronoi neighborhood of each vertex.
	Barycentric = 1 //>! the finite volume for mesh Laplacian is generated from face barycenters.
};

/// \brief a list of triangle metrics to be computed.
using TriangleMetrics = std::vector<std::string>;

/**
 * \brief A wrapper for all settings for BrainSurfaceEvolver.
 */
struct BrainExtractionSettings
{
	std::string ProcedureName = ""; //>! name for the evolution procedure.

	unsigned int NSteps{ 20 };   //>! number of time steps for surface evolution.
	double TimeStep{ 0.01 };     //>! time step size.
	unsigned int IcoSphereSubdivisionLevel{ 3 }; //>! subdivision level of an evolving ico sphere surface.

	BE_CurvatureSettings CurvatureParams{}; //>! evolving surface curvature bounds.
	BE_ThresholdSettings ThresholdSettings{}; //>! settings for detecting relevant contours
	BE_IcoSphereSettings IcoSphereSettings{}; //>! settings for placing the initial ico-sphere surface in space
	BE_MeshTopologySettings TopoParams{}; //>! parameters for mesh topology adjustments.

	bool ExportSurfacePerTimeStep{ false }; //>! whether to export evolving surface for each time step.
	bool ExportResultSurface{ true }; //>! whether to export resulting evolving surface.
	std::string OutputPath{}; //>! path where output surfaces are to be exported.

	BE_MeshLaplacian LaplacianType{}; //>! type of mesh Laplacian.
	TriangleMetrics TriMetrics{}; //>! list of triangle metrics to be computed.

	bool DoRemeshing{ true }; //>! if true, adaptive remeshing will be performed after the first 10-th of time steps.
	bool DoFeatureDetection{ true }; //>! if true, feature detection will take place prior to remeshing.
	bool IdentityForBoundaryVertices{ true }; //>! if true, boundary vertices give rise to: updated vertex = previous vertex.
	bool IdentityForFeatureVertices{ false }; //>! if true, feature vertices give rise to: updated vertex = previous vertex.

	double MaxFractionOfVerticesOutOfBounds{ 0.02 }; //>! fraction of vertices allowed to be out of bounds (because it will be decimated).
};

/**
 * \brief A utility for extracting brain surfaces from CT scans.
 * \class BrainSurfaceEvolver
 */
class BrainSurfaceEvolver
{
public:
	/**
	 * \brief Constructor. Initialize with a given scalar field environment.
	 * \param field                    pre-computed or loaded scalar field environment.
	 * \param settings                 surface evolution settings.
	 */
	BrainSurfaceEvolver(const Geometry::ScalarGrid& field, const BrainExtractionSettings& settings);

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

	BrainExtractionSettings m_EvolSettings{}; //>! settings.

	std::shared_ptr<Geometry::ScalarGrid> m_Field{ nullptr }; //>! scalar field environment.
	std::shared_ptr<pmp::SurfaceMesh> m_EvolvingSurface{ nullptr }; //>! (stabilized) evolving surface.

	pmp::Scalar m_StartingSurfaceRadius{ 1.0f }; //>! radius of the starting surface.
	pmp::Scalar m_ScalingFactor{ 1.0f }; //>! stabilization scaling factor value.

	std::function<pmp::ImplicitLaplaceInfo(const pmp::SurfaceMesh&, pmp::Vertex)> m_ImplicitLaplacianFunction{}; //>! a Laplacian function chosen from parameter MeshLaplacian.
	std::function<double(const pmp::SurfaceMesh&, pmp::Vertex)> m_LaplacianAreaFunction{}; //>! a Laplacian area function chosen from parameter MeshLaplacian.
	
	// export
	pmp::mat4 m_TransformToOriginal{}; //>! transformation matrix from stabilized surface to original size (for export).
	std::string m_OutputMeshExtension = ".vtk"; //>! extension of the exported mesh geometry.
};