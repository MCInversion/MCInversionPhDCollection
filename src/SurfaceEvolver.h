#pragma once

#include "pmp/SurfaceMesh.h"
#include "geometry/Grid.h"


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

	bool ExportSurfacePerTimeStep{ false }; //>! whether to export evolving surface for each time step.
	std::string OutputPath{}; //>! path where output surfaces are to be exported.
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

	SurfaceEvolutionSettings m_EvolSettings{}; //>! settings.

	std::shared_ptr<Geometry::ScalarGrid> m_Field; //>! scalar field environment.
	std::shared_ptr<pmp::SurfaceMesh> m_EvolvingSurface{ nullptr }; //>! (stabilized) evolving surface.

	pmp::mat4 m_TransformToOriginal{}; //>! transformation matrix from stabilized surface to original size (for export).

	std::string m_OutputMeshExtension = ".vtk"; //>! extension of the exported mesh geometry.
};

/**
 * \brief Reports SurfaceEvolver's input to a given stream.
 * \param evolSettings    settings for SurfaceEvolver.
 * \param os              output stream.
 */
void ReportInput(const SurfaceEvolutionSettings& evolSettings, std::ostream& os);
