#pragma once

#include "EvolverUtilsCommon.h"
#include "geometry/Grid.h"

#include "pmp/SurfaceMesh.h"
#include "pmp/algorithms/DifferentialGeometry.h"

/**
 * \brief A wrapper for surface evolution settings.
 * \struct ConvexHullSurfaceEvolutionSettings
 */
struct ConvexHullSurfaceEvolutionSettings
{
	std::string ProcedureName{}; //>! name for the evolution procedure.

	unsigned int NSteps{ 20 };   //>! number of time steps for surface evolution.
	double TimeStep{ 0.01 };     //>! time step size.
	double FieldIsoLevel{ 0.0 }; //>! target level of the scalar field (e.g. zero distance to target mesh).

	unsigned int NVoxelsPerMinDimension{ 40 }; //>! fundamental discretization setting for distance field grid resolution.

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

	// TODO: clear the unnecessary settings

	float TangentialVelocityWeight{ 0.0f }; //>! the weight of tangential velocity update vector for each time step.
	bool DoRemeshing{ true }; //>! if true, adaptive remeshing will be performed after the first 10-th of time steps.
	bool DoFeatureDetection{ true }; //>! if true, feature detection will take place prior to remeshing.
	bool IdentityForBoundaryVertices{ true }; //>! if true, boundary vertices give rise to: updated vertex = previous vertex.
	bool IdentityForFeatureVertices{ false }; //>! if true, feature vertices give rise to: updated vertex = previous vertex.
	double MaxFractionOfVerticesOutOfBounds{ 0.02 }; //>! fraction of vertices allowed to be out of bounds (because it will be decimated).
};

class ConvexHullEvolver
{
public:
    // Constructor
    ConvexHullEvolver(const std::vector<pmp::Point>& pointCloud, const ConvexHullSurfaceEvolutionSettings& settings);

    // Main functionality
    void Evolve();

private:
	/**
	 * \brief Preprocess for evolution, i.e.: generate m_EvolvingSurface, and transform both m_Field and m_EvolvingSurface for stabilization.
	 */
    void Preprocess();

    /**
     * \brief Constructs the starting surface as a convex hull of the input point cloud using the quickhull algorithm.
     */
    void ConstructConvexHull();

    /**
     * \brief Computes the (unsigned) distance field to the input point cloud.
     */
    void ComputeDistanceField();

    // Members
    std::vector<pmp::Point> m_PointCloud;
	ConvexHullSurfaceEvolutionSettings m_EvolSettings; //>! settings.

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
 * \brief Reports ConvexHullEvolver's input to a given stream.
 * \param evolSettings    settings for ConvexHullEvolver.
 * \param os              output stream.
 */
void ReportInput(const ConvexHullSurfaceEvolutionSettings& evolSettings, std::ostream& os);