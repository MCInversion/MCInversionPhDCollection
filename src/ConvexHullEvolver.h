#pragma once

#include "EvolverUtilsCommon.h"
#include "geometry/Grid.h"

#include "pmp/SurfaceMesh.h"
#include "pmp/algorithms/Remeshing.h"
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

	// ================================================================

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

	// ================================================================

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

    // Members
    std::vector<pmp::Point> m_PointCloud;
	ConvexHullSurfaceEvolutionSettings m_EvolSettings; //>! settings.

	std::shared_ptr<Geometry::ScalarGrid> m_Field{ nullptr }; //>! scalar field environment.
    std::shared_ptr<pmp::SurfaceMesh> m_EvolvingSurface{ nullptr }; //>! (stabilized) evolving surface.
	std::shared_ptr<pmp::Remeshing> m_Remesher{ nullptr };  //>! a remesher which keeps the evolving surface with its vlocked_ info for initial convex hull vertices.

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
void ReportCHEvolverInput(const ConvexHullSurfaceEvolutionSettings& evolSettings, std::ostream& os);

/**
 * \brief Computes scaling factor for stabilizing the finite volume method on assumed convex hull surface meshes based on time step.
 * \param timeStep               time step size.
 * \param convexHullMesh         a remeshed convex hull mesh.
 * \param areaFunction           a Laplacian area function measuring the actual co-volume per vertex (Barycentric or Voronoi).
 * \param stabilizationFactor    a multiplier for stabilizing mean co-volume area.
 * \return scaling factor for mesh and scalar grid.
 */
[[nodiscard]] float GetConvexHullStabilizationScalingFactor(const double& timeStep, pmp::SurfaceMesh& convexHullMesh, const AreaFunction& areaFunction, const float& stabilizationFactor = 1.0f);