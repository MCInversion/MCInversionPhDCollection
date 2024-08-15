#pragma once

#include "geometry/GeometryAdapters.h"
#include "geometry/Grid.h"

#include "EvolverUtilsCommon.h"

//
// ===============================================================================================
//                                          Settings
// -----------------------------------------------------------------------------------------------
//

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
    unsigned int LevelOfDetail{ 3 }; //>! Level of detail. Reflected in the number of vertices used during discretization. Standardized for circles and icospheres, e.g.: the base circle is a regular pentagon, and the sphere an icosahedron.

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

//
// ===============================================================================================
//                                     Evolver strategies
// -----------------------------------------------------------------------------------------------
//

/**
 * \brief An interface for the internal functionality of the ManifoldEvolver handling specific geometry types.
 * \class ManifoldEvolutionStrategy
 */
class ManifoldEvolutionStrategy
{
public:
    virtual ~ManifoldEvolutionStrategy() = default;

    /**
     * \brief Preprocess for evolution, i.e.: construct the evolving manifolds, and transform the target data's distance field, and the DF's normalized neg gradient for stabilization.
     */
    virtual void Preprocess() = 0;

    /**
     * \brief Performs a single step of manifold evolution from the configuration at previous time step.
     */
    virtual void PerformEvolutionStep(unsigned int step) = 0;

    /**
     * \brief Verifies whether the tessellation quality of evolving manifold(s) deteriorated so much that remeshing is necessary.
     */
    virtual [[nodiscard]] bool ShouldRemesh() = 0;

    /**
     * \brief Perform remeshing on the evolving manifold(s).
     */
    virtual void Remesh() = 0;

    /**
     * \brief Exports the current state of evolving manifold(s).
     */
    virtual void ExportCurrentState(unsigned int step) = 0;

    /**
     * \brief Exports the final state of evolving manifold(s).
     */
    virtual void ExportFinalResult() = 0;

private:

    std::function<void(unsigned int /* step */)> m_Integrate; //>! numerical integration function (clearly, derived classes have different coefficients/matrices).
};

/**
 * \brief The internal functionality of the ManifoldEvolver handling pmp::ManifoldCurve2D and 2D target set and its distance represented by 2D scalar and vector field.
 * \class ManifoldCurveEvolutionStrategy
 */
class ManifoldCurveEvolutionStrategy : public ManifoldEvolutionStrategy
{
public:
    explicit ManifoldCurveEvolutionStrategy(std::shared_ptr<std::vector<pmp::Point2>> targetPointCloud = nullptr)
	    : m_TargetPointCloud(std::move(targetPointCloud))
    {
    }

    /**
     * \brief Preprocess for evolution, i.e.: construct the evolving manifolds, and transform the target data's distance field, and the DF's normalized neg gradient for stabilization.
     */
    void Preprocess() override;

    /**
     * \brief Performs a single step of manifold evolution from the configuration at previous time step.
     */
    void PerformEvolutionStep(unsigned int step) override;

    /**
     * \brief Verifies whether the tessellation quality of evolving manifold(s) deteriorated so much that remeshing is necessary.
     */
    [[nodiscard]] bool ShouldRemesh() override;

    /**
     * \brief Perform remeshing on the evolving manifold(s).
     */
    void Remesh() override;

    /**
     * \brief Exports the current state of evolving manifold(s).
     */
    void ExportCurrentState(unsigned int step) override;

    /**
     * \brief Exports the final state of evolving manifold(s).
     */
    void ExportFinalResult() override;

    /// \brief get outer curve as a pmp::ManifoldCurve2D
    [[nodiscard]] std::shared_ptr<pmp::ManifoldCurve2D> GetResultOuterCurve() const
    {
        return m_OuterCurve;
    }

    /// \brief get inner curves as pmp::ManifoldCurve2D instances
    [[nodiscard]] std::vector<std::shared_ptr<pmp::ManifoldCurve2D>> GetResultInnerCurves() const
    {
        return m_InnerCurves;
    }

protected:

    [[nodiscard]] std::shared_ptr<pmp::ManifoldCurve2D>& GetOuterCurve()
	{
        return m_OuterCurve;
    }

    [[nodiscard]] std::vector<std::shared_ptr<pmp::ManifoldCurve2D>>& GetInnerCurves()
    {
        return m_InnerCurves;
    }

private:

    std::shared_ptr<std::vector<pmp::Point2>> m_TargetPointCloud{ nullptr }; //>! target point cloud geometry representing the spatial data for the reconstructed manifold.

    std::shared_ptr<pmp::ManifoldCurve2D> m_OuterCurve{ nullptr }; //>! evolving outer manifold which evolves inward, attempting to shrink-wrap m_TargetGeometryAdapter if it's defined.
    std::vector<std::shared_ptr<pmp::ManifoldCurve2D>> m_InnerCurves{}; //>! evolving inner manifolds which evolves outward

    std::shared_ptr<Geometry::ScalarGrid2D> m_DistanceField{ nullptr }; //>! the computed distance field of m_TargetPointCloud on a 2D scalar grid.
    std::shared_ptr<Geometry::VectorGrid2D> m_DFNegNormalizedGradient{ nullptr }; //>! the normalized negative gradient of m_DistanceField.
};

/**
 * \brief An internal functionality of the ManifoldEvolver handling custom pmp::ManifoldCurve2Ds and 2D target set and its distance represented by 2D scalar and vector field.
 * \class CustomManifoldCurveEvolutionStrategy
 */
class CustomManifoldCurveEvolutionStrategy : public ManifoldCurveEvolutionStrategy
{
public:
    /**
     * \brief Constructor that accepts custom outer and inner curves.
     * \param outerCurve Custom outer curve.
     * \param innerCurves Vector of custom inner curves.
     * \param targetPointCloud target point cloud.
     */
    CustomManifoldCurveEvolutionStrategy(
        pmp::ManifoldCurve2D outerCurve,
        std::vector<pmp::ManifoldCurve2D> innerCurves = {},
        std::shared_ptr<std::vector<pmp::Point2>> targetPointCloud = nullptr)
        : ManifoldCurveEvolutionStrategy(targetPointCloud)
    {
        GetOuterCurve() = std::make_shared<pmp::ManifoldCurve2D>(std::move(outerCurve));
        for (auto& c : innerCurves)
            GetInnerCurves().push_back(std::make_shared<pmp::ManifoldCurve2D>(std::move(c)));
    }

    /**
     * \brief Special overridden preprocessing method which omits construction of inner/outer curves.
     */
    void Preprocess() override;
};

/**
 * \brief The internal functionality of the ManifoldEvolver handling pmp::SurfaceMesh and 3D target set and its distance represented by 3D scalar and vector field.
 * \class ManifoldSurfaceEvolutionStrategy
 */
class ManifoldSurfaceEvolutionStrategy : public ManifoldEvolutionStrategy
{
public:
    explicit ManifoldSurfaceEvolutionStrategy(std::shared_ptr<std::vector<pmp::Point>> targetPointCloud = nullptr)
        : m_TargetPointCloud(std::move(targetPointCloud))
    {
    }

    /**
     * \brief Preprocess for evolution, i.e.: construct the evolving manifolds, and transform the target data's distance field, and the DF's normalized neg gradient for stabilization.
     */
    void Preprocess() override;

    /**
     * \brief Performs a single step of manifold evolution from the configuration at previous time step.
     */
    void PerformEvolutionStep(unsigned int step) override;

    /**
     * \brief Verifies whether the tessellation quality of evolving manifold(s) deteriorated so much that remeshing is necessary.
     */
    [[nodiscard]] bool ShouldRemesh() override;

    /**
     * \brief Perform remeshing on the evolving manifold(s).
     */
    void Remesh() override;

    /**
     * \brief Exports the current state of evolving manifold(s).
     */
    void ExportCurrentState(unsigned int step) override;

    /**
     * \brief Exports the final state of evolving manifold(s).
     */
    void ExportFinalResult() override;

    /// \brief get outer surface as a pmp::SurfaceMesh
    [[nodiscard]] std::shared_ptr<pmp::SurfaceMesh> GetResultOuterSurface() const
    {
        return m_OuterSurface;
    }

    /// \brief get inner surfaces as pmp::SurfaceMesh instances
    [[nodiscard]] std::vector<std::shared_ptr<pmp::SurfaceMesh>> GetResultInnerSurfaces() const
    {
        return m_InnerSurfaces;
    }

private:

    std::shared_ptr<std::vector<pmp::Point>> m_TargetPointCloud{ nullptr }; //>! target point cloud geometry representing the spatial data for the reconstructed manifold.

    std::shared_ptr<pmp::SurfaceMesh> m_OuterSurface{ nullptr }; //>! evolving outer manifold which evolves inward, attempting to shrink-wrap m_TargetGeometryAdapter if it's defined.
    std::vector<std::shared_ptr<pmp::SurfaceMesh>> m_InnerSurfaces{ nullptr }; //>! evolving inner manifolds which evolves outward

    std::shared_ptr<Geometry::ScalarGrid> m_DistanceField{ nullptr }; //>! the computed distance field of m_TargetPointCloud on a 3D scalar grid.
    std::shared_ptr<Geometry::VectorGrid> m_DFNegNormalizedGradient{ nullptr }; //>! the normalized negative gradient of m_DistanceField.
};

//
// ===============================================================================================
//                                   Main evolver interface 
// -----------------------------------------------------------------------------------------------
//

/**
 * \brief A class for evolving general manifolds (curves and surfaces) in Euclidean space.
 */
class ManifoldEvolver
{
public:

	ManifoldEvolver(ManifoldEvolutionSettings settings, std::shared_ptr<ManifoldEvolutionStrategy> strategy)
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

    const std::shared_ptr<ManifoldEvolutionStrategy>& GetStrategy() const
	{
        return m_Strategy;
	}

private:

    //
    // ======== Members ===============
    //

    // Key settings
    ManifoldEvolutionSettings m_Settings;

    // Evolution strategy
    std::shared_ptr<ManifoldEvolutionStrategy> m_Strategy{ nullptr };
};