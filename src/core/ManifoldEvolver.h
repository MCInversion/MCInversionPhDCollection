#pragma once

#include "geometry/GeometryAdapters.h"

#include "EvolverUtilsCommon.h"
#include "ManifoldEvolver.h"

namespace Geometry
{
	// TODO: replace dummy fwd decl with actual class
	class ManifoldGeometryAdapter;

	class GeometryAdapter; // for point clouds \Gamma

	class ScalarFieldAdapter; // for distance fields
	class VectorFieldAdapter; // for normalized negative gradients of distance fields

}

struct AmbientFieldSettings
{
    float FieldExpansionFactor{ 1.0f }; //>! the factor by which target bounds are expanded (multiplying original bounds min dimension).
    unsigned int NVoxelsPerMinDimension{ 20 }; //>! the number of voxels per smallest dimension of the resulting scalar grid.
};

/**
 * \brief A wrapper for manifold evolution settings.
 * \struct ManifoldEvolutionSettings
 */
struct ManifoldEvolutionSettings
{
	std::string ProcedureName{}; //>! name for the evolution procedure.

	unsigned int Dimension{ 2 }; //>! the dimension of the evolving manifold(s). 1 for curves, 2 for surfaces.

	// TODO: standardize for circles and icospheres, e.g.: the base circle could be a regular pentagon
	unsigned int LevelOfDetail{ 3 }; //>! Level of detail. Reflected in the number of vertices used during discretization.

	unsigned int NSteps{ 20 };   //>! number of time steps for surface evolution.
	double TimeStep{ 0.01 };     //>! time step size.

	AdvectionDiffusionParameters ADParams{}; //>! parameters for the advection-diffusion model.

	bool ExportPerTimeStep{ false };  //>! whether to export evolving manifold for each time step.
	bool ExportResult{ true }; //>! whether to export resulting evolving manifold.

	std::string OutputPath{}; //>! path where output manifolds are to be exported.

	bool DoRemeshing{ true }; //>! if true, adaptive remeshing will be performed after the first 10-th of time steps.

	double MaxFractionOfVerticesOutOfBounds{ 0.02 }; //>! fraction of vertices allowed to be out of bounds (because it will be decimated).
};

/**
 * \brief A class for evolving general manifolds (curves and surfaces) in Euclidean space.
 */
class ManifoldEvolver
{
public:

	ManifoldEvolver(ManifoldEvolutionSettings settings, std::shared_ptr<Geometry::GeometryAdapter> targetGeom = nullptr)
        : m_Settings(std::move(settings)) //, m_TargetGeometryAdapter(std::make_shared<Geometry::GeometryAdapter>(targetGeom))
	{
	}

    /**
     * \brief Evolves the manifold over the specified number of time steps.
     */
    void Evolve();

private:

    void PreprocessFieldFromTargetGeometry();

    void InitializeManifolds();

    /**
     * \brief Performs a single evolution step.
     * \param step      The current time step index.
     */
    void PerformEvolutionStep(unsigned int step);

    /**
     * \brief Performs adaptive remeshing if required.
     */
    void AdaptiveRemeshing();

    /**
     * \brief Exports the current state of the evolving manifold.
     * \param step       Index of the current time step.
     */
    void ExportCurrentState(unsigned int step);

    /**
     * \brief Exports the final result of the evolving manifold.
     */
    void ExportFinalResult();

    //
    // ======== Members ===============
    //

    std::shared_ptr<Geometry::GeometryAdapter> m_TargetGeometryAdapter{ nullptr }; //>! target geometry (mesh/point cloud) representing the spatial data for the reconstructed manifold.

    std::shared_ptr<Geometry::ManifoldGeometryAdapter> m_OuterManifoldAdapter{ nullptr }; //>! evolving outer manifold which evolves inward, attempting to shrink-wrap m_TargetGeometryAdapter if it's defined.
    std::vector<std::shared_ptr<Geometry::ManifoldGeometryAdapter>> m_InnerManifoldAdapters{}; //>! evolving inner manifolds which evolves outward

    std::shared_ptr<Geometry::ScalarFieldAdapter> m_Field{ nullptr };
    std::shared_ptr<Geometry::VectorFieldAdapter> m_Gradient{ nullptr };

    ManifoldEvolutionSettings m_Settings;

    pmp::Scalar m_ScalingFactor{ 1.0f }; //>! stabilization scaling factor value.

    // export
    std::string m_OutputExtension = ".vtk"; //>! extension of the exported geometry.
    pmp::mat4 m_TransformToOriginal{}; //>! transformation matrix from stabilized surface to original size (for export).
};