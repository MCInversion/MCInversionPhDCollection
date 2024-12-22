#pragma once

#include "pmp/SurfaceMesh.h"
#include "pmp/algorithms/DifferentialGeometry.h"

#include "geometry/Grid.h"

#include "EvolverUtilsCommon.h"

/**
 * \brief A wrapper for iso-surface evolution settings.
 * \struct IsoSurfaceEvolutionSettings
 */
struct IsoSurfaceEvolutionSettings
{
	std::string ProcedureName{}; //>! name for the evolution procedure.

	unsigned int NSteps{ 20 };   //>! number of time steps for surface evolution.
	double TimeStep{ 0.01 };     //>! time step size.
	double FieldIsoLevel{ 0.0 }; //>! target level of the scalar field (e.g. zero distance to target mesh).
	double FieldIsoLevelOffset{ 1.0 }; //! value by which the shrink-wrapped surface is offset from the target surface.

	pmp::Scalar ReSampledGridCellSize{ 1.0 }; //! cell size of the grid from which an isosurface will be generated.

	AdvectionDiffusionParameters ADParams{}; //>! parameters for the advection-diffusion model.
	MeshTopologySettings TopoParams{}; //>! parameters for mesh topology adjustments.

	bool ExportSurfacePerTimeStep{ false }; //>! whether to export evolving surface for each time step.
	bool ExportResultSurface{ true }; //>! whether to export resulting evolving surface.
	std::string OutputPath{}; //>! path where output surfaces are to be exported.

	MeshLaplacian LaplacianType{}; //>! type of mesh Laplacian.
	TriangleMetrics TriMetrics{}; //>! list of triangle metrics to be computed.

	pmp::Scalar TangentialVelocityWeight{ 0.0 }; //>! the weight of tangential velocity update vector for each time step.
	bool DoRemeshing{ true }; //>! if true, adaptive remeshing will be performed after the first 10-th of time steps.
	bool DoFeatureDetection{ true }; //>! if true, feature detection will take place prior to remeshing.
	bool IdentityForBoundaryVertices{ true }; //>! if true, boundary vertices give rise to: updated vertex = previous vertex.
	bool IdentityForFeatureVertices{ false }; //>! if true, feature vertices give rise to: updated vertex = previous vertex.

	double MaxFractionOfVerticesOutOfBounds{ 0.02 }; //>! fraction of vertices allowed to be out of bounds (because it will be decimated).
};

/**
 * \brief A utility for evolving iso-surfaces within a scalar field.
 * \class IsoSurfaceEvolver
 */
class IsoSurfaceEvolver
{
public:
	/**
	 * \brief Constructor. Initialize with a given scalar field environment.
	 * \param field                    pre-computed or loaded scalar field environment.
	 * \param fieldExpansionFactor     the factor by which target bounds are expanded (multiplying original bounds min dimension).
	 * \param settings                 surface evolution settings.
	 */
	IsoSurfaceEvolver(const Geometry::ScalarGrid& field, const pmp::Scalar& fieldExpansionFactor, const IsoSurfaceEvolutionSettings& settings);

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

	IsoSurfaceEvolutionSettings m_EvolSettings{}; //>! settings.

	std::shared_ptr<Geometry::ScalarGrid> m_Field{ nullptr }; //>! scalar field environment.
	std::shared_ptr<pmp::SurfaceMesh> m_EvolvingSurface{ nullptr }; //>! (stabilized) evolving surface.

	pmp::Scalar m_ExpansionFactor{ 0.0 }; //>! the factor by which target bounds are expanded (multiplying original bounds min dimension).
	pmp::Scalar m_ScalingFactor{ 1.0 }; //>! stabilization scaling factor value.

	std::function<pmp::ImplicitLaplaceInfo(const pmp::SurfaceMesh&, pmp::Vertex)> m_ImplicitLaplacianFunction{}; //>! a Laplacian function chosen from parameter MeshLaplacian.
	std::function<double(const pmp::SurfaceMesh&, pmp::Vertex)> m_LaplacianAreaFunction{}; //>! a Laplacian area function chosen from parameter MeshLaplacian.
	std::function<size_t(const pmp::Scalar&, const bool&)> m_FeatureFunction{}; //>! a function for detecting mesh features.

	// export
	std::string m_OutputMeshExtension = ".vtk"; //>! extension of the exported mesh geometry.
	pmp::mat4 m_TransformToOriginal{}; //>! transformation matrix from stabilized surface to original size (for export).

};

/**
 * \brief Reports IsoSurfaceEvolver's input to a given stream.
 * \param evolSettings    settings for IsoSurfaceEvolver.
 * \param os              output stream.
 */
void ReportInput(const IsoSurfaceEvolutionSettings& evolSettings, std::ostream& os);

/**
 * \brief Computes scaling factor for stabilizing the finite volume method on assumed isosurface meshes with face size close to voxel size.
 * \param timeStep               time step size.
 * \param cellSize               size of the voxel.
 * \param stabilizationFactor    a multiplier for stabilizing mean co-volume area.
 * \return scaling factor for mesh and scalar grid.
 */
[[nodiscard]] pmp::Scalar GetStabilizationScalingFactor(const double& timeStep, const pmp::Scalar& cellSize, const pmp::Scalar& stabilizationFactor = 1.0);
