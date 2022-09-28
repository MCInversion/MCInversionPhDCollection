#pragma once

#include "pmp/SurfaceMesh.h"
#include "geometry/Grid.h"

/**
 * \brief A wrapper for surface evolution advectoion-diffusion parameters.
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
 * \brief A wrapper for surface evolution settings.
 * \struct SurfaceEvolutionSettings
 */
struct SurfaceEvolutionSettings
{
	std::string ProcedureName = "";    //>! name for the evolution procedure.

	unsigned int NSteps{ 20 };   //>! number of time steps for surface evolution.
	double TimeStep{ 0.01 };     //>! time step size.
	double FieldIsoLevel{ 0.0 }; //>! target level of the scalar field (e.g. zero distance to target mesh).
	unsigned int IcoSphereSubdivisionLevel{ 3 }; //>! subdivision level of an evolving ico sphere surface.

	AdvectionDiffusionParameters ADParams{}; //>! parameters for the advection-diffusion model.

	float MinTargetSize{ 1.0f }; //>! minimum size of the target mesh bounding box.
	float MaxTargetSize{ 1.0f }; //>! maximum size of the target mesh bounding box.
	pmp::vec3 TargetOrigin{}; //>! origin of the evolution's target.

	bool ExportSurfacePerTimeStep{ false }; //>! whether to export evolving surface for each time step.
	bool ExportResultSurface{ true }; //>! whether to export resulting evolving surface.
	std::string OutputPath{}; //>! path where output surfaces are to be exported.
	bool DoRemeshing{ true };
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
	 * \param field           pre-computed or loaded scalar field environment.
	 * \param settings        surface evolution settings.
	 */
	SurfaceEvolver(const Geometry::ScalarGrid& field, const SurfaceEvolutionSettings& settings)
		: m_EvolSettings(settings), m_Field(std::make_shared<Geometry::ScalarGrid>(field))
	{ }

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

	// ----------------------------------------------------------------

	SurfaceEvolutionSettings m_EvolSettings{}; //>! settings.

	std::shared_ptr<Geometry::ScalarGrid> m_Field; //>! scalar field environment.
	std::shared_ptr<pmp::SurfaceMesh> m_EvolvingSurface{ nullptr }; //>! (stabilized) evolving surface.

	pmp::Scalar m_StartingSurfaceRadius{ 1.0f }; //>! radius of the starting surface.
	pmp::Scalar m_ScalingFactor{ 1.0f }; //>! stabilization scaling factor value.
	pmp::mat4 m_TransformToOriginal{}; //>! transformation matrix from stabilized surface to original size (for export).

	std::string m_OutputMeshExtension = ".vtk"; //>! extension of the exported mesh geometry.
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
