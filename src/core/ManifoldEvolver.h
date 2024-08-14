#pragma once

#include "geometry/GeometryAdapters.h"
#include "geometry/Grid.h"

#include "EvolverUtilsCommon.h"

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
    bool UseSemiImplicit{ true }; //>! whether to use a more numerically stable, but computationally costly numerical scheme (has to solve a linear system for each coord for each time step).
    unsigned int LevelOfDetail{ 3 }; //>! Level of detail. Reflected in the number of vertices used during discretization. Standardized for circles and icospheres, e.g.: the base circle could be a regular pentagon

    unsigned int NSteps{ 20 };   //>! number of time steps for surface evolution.
    double TimeStep{ 0.01 };     //>! time step size.

    AdvectionDiffusionParameters ADParams{}; //>! parameters for the advection-diffusion model.

    bool ExportPerTimeStep{ false };  //>! whether to export evolving manifold for each time step.
    bool ExportResult{ true }; //>! whether to export resulting evolving manifold.
    bool ExportInnerManifolds{ false }; //>! if true also m_InnerManifoldAdapters for each time step. 

    std::string OutputPath{}; //>! path where output manifolds are to be exported.

    bool DoRemeshing{ true }; //>! if true, adaptive remeshing will be performed after the first 10-th of time steps.

    double MaxFractionOfVerticesOutOfBounds{ 0.02 }; //>! fraction of vertices allowed to be out of bounds (because it will be decimated).
};

class ManifoldEvolutionStrategy
{
public:
    virtual ~ManifoldEvolutionStrategy() = default;

    virtual void Preprocess() = 0;

    virtual void PerformEvolutionStep(unsigned int step) = 0;

    virtual [[nodiscard]] bool ShouldRemesh() = 0;

    virtual void Remesh() = 0;

    virtual void ExportCurrentState(unsigned int step) = 0;

    virtual void ExportFinalResult() = 0;
private:

    std::function<void(unsigned int /* step */)> m_Integrate; //>! numerical integration function (clearly, derived classes have different coefficients/matrices).
};

class ManifoldCurveEvolutionStrategy : public ManifoldEvolutionStrategy
{
public:
    ManifoldCurveEvolutionStrategy(std::shared_ptr<std::vector<pmp::Point2>> targetPointCloud = nullptr)
	    : m_TargetPointCloud(std::move(targetPointCloud))
    {
    }

private:

    std::shared_ptr<std::vector<pmp::Point2>> m_TargetPointCloud{ nullptr }; //>! target point cloud geometry representing the spatial data for the reconstructed manifold.

    std::shared_ptr<pmp::ManifoldCurve2D> m_OuterCurve{ nullptr }; //>! evolving outer manifold which evolves inward, attempting to shrink-wrap m_TargetGeometryAdapter if it's defined.
    std::vector<std::shared_ptr<pmp::ManifoldCurve2D>> m_InnerCurves{}; //>! evolving inner manifolds which evolves outward

    std::shared_ptr<Geometry::ScalarGrid2D> m_DistanceField{ nullptr }; //>! the computed distance field of m_TargetPointCloud on a 2D scalar grid.
    std::shared_ptr<Geometry::VectorGrid2D> m_DFNegNormalizedGradient{ nullptr }; //>! the normalized negative gradient of m_DistanceField.
};

class ManifoldSurfaceEvolutionStrategy : public ManifoldEvolutionStrategy
{
public:
    ManifoldSurfaceEvolutionStrategy(std::shared_ptr<std::vector<pmp::Point>> targetPointCloud = nullptr)
        : m_TargetPointCloud(std::move(targetPointCloud))
    {
    }

private:

    std::shared_ptr<std::vector<pmp::Point>> m_TargetPointCloud{ nullptr }; //>! target point cloud geometry representing the spatial data for the reconstructed manifold.

    std::shared_ptr<pmp::SurfaceMesh> m_OuterSurface{ nullptr }; //>! evolving outer manifold which evolves inward, attempting to shrink-wrap m_TargetGeometryAdapter if it's defined.
    std::vector<std::shared_ptr<pmp::SurfaceMesh>> m_InnerSurfaces{ nullptr }; //>! evolving inner manifolds which evolves outward

    std::shared_ptr<Geometry::ScalarGrid> m_DistanceField{ nullptr }; //>! the computed distance field of m_TargetPointCloud on a 3D scalar grid.
    std::shared_ptr<Geometry::VectorGrid> m_DFNegNormalizedGradient{ nullptr }; //>! the normalized negative gradient of m_DistanceField.
};

/**
 * \brief A class for evolving general manifolds (curves and surfaces) in Euclidean space.
 */
class ManifoldEvolver
{
public:

	ManifoldEvolver(ManifoldEvolutionSettings settings, std::unique_ptr<ManifoldEvolutionStrategy> strategy)
        : m_Settings(std::move(settings)),
		  m_Strategy(std::move(strategy))
	{
	}

    /**
     * \brief Evolves the manifold over the specified number of time steps.
     */
    void Evolve() const
	{
        m_Strategy->Preprocess();

        for (unsigned int step = 0; step < m_Settings.NSteps; ++step)
        {
            m_Strategy->PerformEvolutionStep(step);

            if (m_Settings.ExportPerTimeStep)
            {
                m_Strategy->ExportCurrentState(step);
            }

            if (m_Settings.DoRemeshing && m_Strategy->ShouldRemesh())
            {
                m_Strategy->Remesh();
            }
        }

        if (m_Settings.ExportResult) 
        {
            m_Strategy->ExportFinalResult();
        }
    }

private:

    //
    // ======== Members ===============
    //

    // Key settings
    ManifoldEvolutionSettings m_Settings;

    // Evolution strategy
    std::unique_ptr<ManifoldEvolutionStrategy> m_Strategy{ nullptr };
};