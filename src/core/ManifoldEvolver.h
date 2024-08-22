#pragma once

#include "pmp/algorithms/DifferentialGeometry.h"
#include "pmp/MatVec.h"

#include "geometry/GeometryAdapters.h"
#include "geometry/Grid.h"
#include "geometry/GridUtil.h"

#include "EvolverUtilsCommon.h"

//
// ===============================================================================================
//                                     Evolver functions
// -----------------------------------------------------------------------------------------------
//

/// \brief eps(d) = 1, eps(d) = C_1 * (1 - exp(-d^2/C_2)), or any other control function for the Laplacian term in the equation. 
using CurvatureCtrlFunction = std::function<double(double /* distance */)>;

const CurvatureCtrlFunction TRIVIAL_EPSILON = [](double /* distance */) { return 1.0; };
const CurvatureCtrlFunction STANDARD_EPSILON = [](double distance) { return 1.0 - exp(-distance * distance); };

/// \brief eta(d) = 0, eta(d) = D_1 * d * ((grad d . N) - D_2 * sqrt(1 - (grad d . N)^2), or any other control function for the advection term in the equation.
using AdvectionCtrlFunction = std::function<double(double /* distance */, double /* negGradDotNormal */)>;

const AdvectionCtrlFunction TRIVIAL_ETA = [](double /* distance */, double /* negGradDotNormal */) { return 0.0; };
const AdvectionCtrlFunction STANDARD_ETA = [](double distance, double negGradDotNormal)
{
    if (distance >= Geometry::DEFAULT_SCALAR_GRID_INIT_VAL)
        return 0.0;
    return distance * (negGradDotNormal - sqrt(1.0 - negGradDotNormal * negGradDotNormal));
};

/// \brief Numerical integration function specific to the scheme and dimension.
using NumericalStepIntegrateFunction = std::function<void(unsigned int /* step */)>;

/// \brief Functions for interpolating 2D scalar grids
using ScalarGridInterpolationFunction2D = std::function<double(const pmp::Point2& /* pos */, const Geometry::ScalarGrid2D& /* grid */)>;
/// \brief Functions for interpolating 3D scalar grids
using ScalarGridInterpolationFunction3D = std::function<double(const pmp::Point& /* pos */, const Geometry::ScalarGrid& /* grid */)>;
/// \brief Functions for interpolating 2D vector grids
using VectorGridInterpolationFunction2D = std::function<pmp::dvec2(const pmp::Point2& /* pos */, const Geometry::VectorGrid2D& /* vecGrid */)>;
/// \brief Functions for interpolating 3D vector grids
using VectorGridInterpolationFunction3D = std::function<pmp::dvec3(const pmp::Point& /* pos */, const Geometry::VectorGrid& /* vecGrid */)>;

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
    bool UseInnerManifolds{ true }; //>! whether to construct outward-evolving inner manifolds.
    bool UseSemiImplicit{ true }; //>! whether to use a more numerically stable, but computationally costly numerical scheme (has to solve a linear system for each coord for each time step).
    unsigned int LevelOfDetail{ 3 }; //>! Level of detail. Reflected in the number of vertices used during discretization. Standardized for circles and icospheres, e.g.: the base circle is a regular pentagon, and the sphere an icosahedron.

    double TimeStep{ 0.01 };     //>! time step size.
    bool UseStabilizationViaScaling{ true }; //>! whether to scale all evolving geometries to match with the time step.

    bool UseLinearGridInterpolation{ true }; //>! whether to use d-linear interpolation for scalar and vector fields where d is the grid dimension.
    CurvatureCtrlFunction Epsilon{ TRIVIAL_EPSILON }; //>! control function for the curvature term.
    AdvectionCtrlFunction Eta{ TRIVIAL_ETA }; //>! control function for the advection term.

    bool AdvectionInteractWithOtherManifolds{ false }; //>! whether to use the minimum from the distances to all manifolds including target data.

    AmbientFieldSettings FieldSettings{}; //>! the settings for the construction of the ambient fields.

    double TangentialVelocityWeight{ 0.05 }; //>! the weight of tangential velocity update vector for each time step.

    double MaxFractionOfVerticesOutOfBounds{ 0.02 }; //>! fraction of vertices allowed to be out of bounds (because it will be decimated).
};

/**
 * \brief A wrapper for the global manifold evolution settings.
 * \struct GlobalManifoldEvolutionSettings
 */
struct GlobalManifoldEvolutionSettings
{
    std::string ProcedureName{}; //>! name for the evolution procedure.

    unsigned int NSteps{ 20 };   //>! number of time steps for surface evolution.

    bool ExportPerTimeStep{ false };  //>! whether to export evolving manifold for each time step.
    bool ExportResult{ true }; //>! whether to export resulting evolving manifold.
    bool ExportInnerManifolds{ false }; //>! if true also m_InnerManifoldAdapters for each time step. 
    bool ExportTargetDistanceFieldAsImage{ false }; //>! if true the m_DistanceField in the active strategy will be exported as a png/vti image.

    std::string OutputPath{}; //>! path where output manifolds are to be exported.

    bool DoRemeshing{ true }; //>! if true, adaptive remeshing will be performed after the first 10-th of time steps.
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

    /// \brief Base constructor from given settings.
    explicit ManifoldEvolutionStrategy(ManifoldEvolutionSettings settings)
	    : m_Settings(settings)
    {
    }

    virtual ~ManifoldEvolutionStrategy() = default;

    /**
     * \brief Preprocess for evolution, i.e.: construct the evolving manifolds, and transform the target data's distance field, and the DF's normalized neg gradient for stabilization.
     */
    virtual void Preprocess() = 0;

    /**
     * \brief Performs a single step of manifold evolution from the configuration at previous time step.
     */
    virtual void PerformEvolutionStep(unsigned int stepId) = 0;

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
    virtual void ExportCurrentState(unsigned int step, const std::string& baseOutputFilename) = 0;

    /**
     * \brief Exports the final state of evolving manifold(s).
     */
    virtual void ExportFinalResult(const std::string& baseOutputFilename) = 0;

    /**
	 * \brief Exports the target point cloud distance field as image data (if defined).
	 */
    virtual void ExportTargetDistanceFieldAsImage(const std::string& baseOutputFilename) = 0;

    /// \brief Settings getter.
    ManifoldEvolutionSettings& GetSettings()
    {
        return m_Settings;
    }

    /// \brief Calculate the distance fields to evolving manifolds if there should be an interaction between these manifolds.
    virtual void ComputeVariableDistanceFields() = 0;

protected:

    /// \brief A getter for the curvature term control function (epsilon).
    CurvatureCtrlFunction& GetCurvatureTermControlFunction()
    {
        return m_Settings.Epsilon;
    }

    /// \brief A getter for the advection control function (eta).
    AdvectionCtrlFunction& GetAdvectionCtrlFunction()
    {
        return m_Settings.Eta;
    }

    /// \brief A getter for the stabilization scaling factor.
    pmp::Scalar& GetScalingFactor()
    {
        return m_ScalingFactor;
    }

    /// \brief Transform all of the geometries so that numerical stability is ensured.
    /// \param[in] outerRadius   radius of outer sphere to calculate the approx co-volume measure.
    /// \param[in] stabilizationFactor     a multiplier for stabilizing mean co-volume measure.
    virtual void StabilizeGeometries(float outerRadius, float stabilizationFactor = 1.0f) = 0;

    /// \brief Prepares the property arrays for all evolving manifolds.
    virtual void PrepareManifoldProperties() = 0;

    /// \brief A getter for the numerical integration step function.
    NumericalStepIntegrateFunction& GetIntegrate()
    {
        return m_Integrate;
    }

    /// \brief A placeholder for the semi-implicit integration method.
    virtual void SemiImplicitIntegrationStep(unsigned int step) = 0;

    /// \brief A placeholder for the explicit integration method.
    virtual void ExplicitIntegrationStep(unsigned int step) = 0;

    /// \brief A getter for the standardized cell size for the variable distance fields
    pmp::Scalar& GetFieldCellSize()
    {
        return m_FieldCellSize;
    }

private:

    ManifoldEvolutionSettings m_Settings{};       //>! settings for the evolution strategy.
    NumericalStepIntegrateFunction m_Integrate{}; //>! numerical integration function (clearly, derived classes have different coefficients/matrices).

    pmp::Scalar m_ScalingFactor{ 1.0f }; //>! the computed scaling factor for numerical stabilization.
    pmp::Scalar m_FieldCellSize{ 0.1f }; //>! standardized cell size for the variable distance fields.
};

/**
 * \brief The internal functionality of the ManifoldEvolver handling pmp::ManifoldCurve2D and 2D target set and its distance represented by 2D scalar and vector field.
 * \class ManifoldCurveEvolutionStrategy
 */
class ManifoldCurveEvolutionStrategy : public ManifoldEvolutionStrategy
{
public:
	/**
	 * \brief A base constructor for curve evolution strategy.
	 * \param settings           evolution strategy settings.
	 * \param targetPointCloud   target point cloud data.
	 */
	explicit ManifoldCurveEvolutionStrategy(ManifoldEvolutionSettings settings, 
	                                        std::shared_ptr<std::vector<pmp::Point2>> targetPointCloud = nullptr)
	    : ManifoldEvolutionStrategy(settings),
        m_TargetPointCloud(std::move(targetPointCloud))
    {
        if (GetSettings().UseSemiImplicit)
        {
            GetIntegrate() = [this](unsigned int step)
            {
	            SemiImplicitIntegrationStep(step);
            };
        }
        else
        {
            GetIntegrate() = [this](unsigned int step)
            {
	            ExplicitIntegrationStep(step);
            };
        }

        if (GetSettings().UseLinearGridInterpolation)
        {
            m_ScalarInterpolate = Geometry::BilinearInterpolateScalarValue;
            m_VectorInterpolate = Geometry::BilinearInterpolateVectorValue;
        }
        else
        {
            m_ScalarInterpolate = Geometry::GetNearestNeighborScalarValue2D;
            m_VectorInterpolate = Geometry::GetNearestNeighborVectorValue2D;
        }
    }

    /**
     * \brief Preprocess for evolution, i.e.: construct the evolving manifolds, and transform the target data's distance field, and the DF's normalized neg gradient for stabilization.
     */
    void Preprocess() override;

    /**
     * \brief Performs a single step of manifold evolution from the configuration at previous time step.
     */
    void PerformEvolutionStep(unsigned int stepId) override;

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
    void ExportCurrentState(unsigned int step, const std::string& baseOutputFilename) override;

    /**
     * \brief Exports the final state of evolving manifold(s).
     */
    void ExportFinalResult(const std::string& baseOutputFilename) override;

    /**
     * \brief Exports the target point cloud distance field as png image data (if defined).
     */
    void ExportTargetDistanceFieldAsImage(const std::string& baseOutputFilename) override;

    /// \brief A specialized external getter which also transforms the stabilized outer curve to its original scale.
    [[nodiscard]] std::shared_ptr<pmp::ManifoldCurve2D> GetOuterCurveInOrigScale() const;

    /// \brief A specialized external getter which also transforms the stabilized inner curves to its original scale.
    [[nodiscard]] std::vector<std::shared_ptr<pmp::ManifoldCurve2D>> GetInnerCurvesInOrigScale() const;

    /// \brief Calculate m_OuterCurveDistanceField and m_InnerCurvesDistanceFields if there should be an interaction between these manifolds.
    void ComputeVariableDistanceFields() override;

protected:

    /// \brief Semi-implicit integration method step.
    void SemiImplicitIntegrationStep(unsigned int step) override;

    /// \brief Explicit integration method step.
    void ExplicitIntegrationStep(unsigned int step) override;

    /// \brief Prepares the property arrays for all evolving manifolds.
    void PrepareManifoldProperties() override;

    /// \brief A getter for the outer curve
    std::shared_ptr<pmp::ManifoldCurve2D>& GetOuterCurve()
	{
        return m_OuterCurve;
    }

    /// \brief A getter for the vector of inner curves
    std::vector<std::shared_ptr<pmp::ManifoldCurve2D>>& GetInnerCurves()
    {
        return m_InnerCurves;
    }

    /// \brief A const getter for the outer curve
    [[nodiscard]] const std::shared_ptr<pmp::ManifoldCurve2D>& GetOuterCurve() const
    {
        return m_OuterCurve;
    }

    /// \brief A const getter for the vector of inner curves
    [[nodiscard]] const std::vector<std::shared_ptr<pmp::ManifoldCurve2D>>& GetInnerCurves() const
    {
        return m_InnerCurves;
    }

    /// \brief A getter for the inverse stabilization transformation matrix.
    pmp::mat3& GetTransformToOriginal()
    {
        return m_TransformToOriginal;
    }

    /// \brief A getter for the distance field scalar grid.
    std::shared_ptr<Geometry::ScalarGrid2D>& GetDistanceField()
    {
        return m_DistanceField;
    }

    /// \brief A getter for the distance field normalized negative gradient vector grid.
    std::shared_ptr<Geometry::VectorGrid2D>& GetDFNegNormalizedGradient()
    {
        return m_DFNegNormalizedGradient;
    }

    /// \brief Calculate the m_DistanceField and m_DFNegNormalizedGradient.
    /// \return triple { minTargetSize, maxTargetSize, targetBoundsCenter }.
    [[nodiscard]] std::tuple<float, float, pmp::Point2> ComputeAmbientFields();

    /// \brief Construct m_OuterCurve and m_InnerCurves from settings.
    /// \param[in] minTargetSize        minimal size of the target data bounding box. Used for computing the radius of the outer manifold.
    /// \param[in] maxTargetSize        maximal size of the target data bounding box. Used for computing the radius of the outer manifold.
    /// \param[in] targetBoundsCenter   the center of the target data bounding box. Used for proper centering the initial outer manifold.
    /// \return radius of outer sphere.
    [[nodiscard]] float ConstructInitialManifolds(float minTargetSize, float maxTargetSize, const pmp::Point2& targetBoundsCenter);

    /// \brief Transform all of the geometries so that numerical stability is ensured.
    /// \param[in] outerRadius             radius of outer sphere to calculate the approx co-volume measure.
    /// \param[in] stabilizationFactor     a multiplier for stabilizing mean co-volume measure.
    void StabilizeGeometries(float outerRadius, float stabilizationFactor = 1.0f) override;

    /// \brief A getter for the scalar grid interpolator function.
    ScalarGridInterpolationFunction2D& GetScalarInterpolate()
    {
        return m_ScalarInterpolate;
    }

    /// \brief A getter for the vector grid interpolator function.
    VectorGridInterpolationFunction2D GetVectorInterpolate()
    {
        return m_VectorInterpolate;
    }

    /// \brief A getter for the test box for numerical validity.
    pmp::BoundingBox2& GetEvolBox()
    {
        return m_EvolBox;
    }

    /// \brief A getter for the distance field to outer curve
    std::shared_ptr<Geometry::ScalarGrid2D>& GetOuterCurveDistanceField()
    {
        return m_OuterCurveDistanceField;
    }

    /// \brief A getter for the distance fields of inner curves
    std::vector<std::shared_ptr<Geometry::ScalarGrid2D>>& GetInnerCurvesDistanceFields()
    {
        return m_InnerCurvesDistanceFields;
    }

private:

    std::shared_ptr<std::vector<pmp::Point2>> m_TargetPointCloud{ nullptr }; //>! target point cloud geometry representing the spatial data for the reconstructed manifold.

    std::shared_ptr<pmp::ManifoldCurve2D> m_OuterCurve{ nullptr }; //>! evolving outer manifold which evolves inward, attempting to shrink-wrap m_TargetGeometryAdapter if it's defined.
    std::vector<std::shared_ptr<pmp::ManifoldCurve2D>> m_InnerCurves{}; //>! evolving inner manifolds which evolves outward

    std::shared_ptr<Geometry::ScalarGrid2D> m_DistanceField{ nullptr }; //>! the computed distance field of m_TargetPointCloud on a 2D scalar grid.
    std::shared_ptr<Geometry::VectorGrid2D> m_DFNegNormalizedGradient{ nullptr }; //>! the normalized negative gradient of m_DistanceField.

    // TODO: investigate performance
    std::shared_ptr<Geometry::ScalarGrid2D> m_OuterCurveDistanceField{ nullptr }; //>! the updated distance field of the evolving outer surface.
    std::vector<std::shared_ptr<Geometry::ScalarGrid2D>> m_InnerCurvesDistanceFields{}; //>! the updated distance fields of the evolving inner surfaces.

    pmp::mat3 m_TransformToOriginal = pmp::mat3::identity(); //>! a transformation matrix to transform the stabilized geometry back to its original scale.

    ScalarGridInterpolationFunction2D m_ScalarInterpolate{}; //>! a parametrizeable function for interpolating values within Geometry::ScalarGrid2D.
    VectorGridInterpolationFunction2D m_VectorInterpolate{};  //>! a parametrizeable function for interpolating vector values within Geometry::VectorGrid2D.

    pmp::BoundingBox2 m_EvolBox{}; //>! the test box for numerical validity of the evolution.
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
     * \param settings            evolution strategy settings.
     * \param outerCurve          Custom outer curve.
     * \param innerCurves         Vector of custom inner curves.
     * \param targetPointCloud    target point cloud.
     */
    explicit CustomManifoldCurveEvolutionStrategy(
        ManifoldEvolutionSettings settings,
        pmp::ManifoldCurve2D& outerCurve,
        std::vector<pmp::ManifoldCurve2D>& innerCurves,
        std::shared_ptr<std::vector<pmp::Point2>> targetPointCloud = nullptr)
        : ManifoldCurveEvolutionStrategy(settings, std::move(targetPointCloud))
    {
        GetOuterCurve() = std::make_shared<pmp::ManifoldCurve2D>(outerCurve);
        for (auto& c : innerCurves)
            GetInnerCurves().emplace_back(std::make_shared<pmp::ManifoldCurve2D>(c));
    }

    /**
     * \brief Special overridden preprocessing method which omits construction of inner/outer curves.
     */
    void Preprocess() override;

private:

    /// \brief Verifies whether all custom inner curves are contained within the custom outer curve if not, throw std::invalid_argument
    [[nodiscard]] bool HasValidInnerOuterManifolds() const;

    /// \brief Computes the full range of co-volume sizes (lengths) to help compute the stabilization scaling factor.
    [[nodiscard]] std::pair<float, float> CalculateCoVolumeRange() const;

    /// \brief Transform all of the geometries so that numerical stability is ensured.
	/// \param[in] minLength               the minimum length of a 1D co-volume within the custom curves.
	/// \param[in] maxLength               the maximum length of a 1D co-volume within the custom curves.
	/// \param[in] stabilizationFactor     a multiplier for stabilizing mean co-volume measure.
    void StabilizeCustomGeometries(float minLength, float maxLength, float stabilizationFactor = 1.0f);
};

/**
 * \brief The internal functionality of the ManifoldEvolver handling pmp::SurfaceMesh and 3D target set and its distance represented by 3D scalar and vector field.
 * \class ManifoldSurfaceEvolutionStrategy
 */
class ManifoldSurfaceEvolutionStrategy : public ManifoldEvolutionStrategy
{
public:
    explicit ManifoldSurfaceEvolutionStrategy(
        ManifoldEvolutionSettings settings, MeshLaplacian laplacianType,
        std::shared_ptr<std::vector<pmp::Point>> targetPointCloud = nullptr)
        : ManifoldEvolutionStrategy(std::move(settings)),
		m_TargetPointCloud(std::move(targetPointCloud))
    {
        m_LaplacianAreaFunction = (laplacianType == MeshLaplacian::Barycentric ?
                pmp::barycentric_area : pmp::voronoi_area);
        m_ExplicitLaplacianFunction = (laplacianType == MeshLaplacian::Barycentric ?
            pmp::laplace_barycentric : pmp::laplace_voronoi);
        m_ImplicitLaplacianFunction = (laplacianType == MeshLaplacian::Barycentric ?
                pmp::laplace_implicit_barycentric : pmp::laplace_implicit_voronoi);

        if (GetSettings().UseSemiImplicit)
        {
            GetIntegrate() = [this](unsigned int step)
            {
	            SemiImplicitIntegrationStep(step);
            };
        }
        else
        {
            GetIntegrate() = [this](unsigned int step)
            {
	            ExplicitIntegrationStep(step);
            };
        }

        if (GetSettings().UseLinearGridInterpolation)
        {
            m_ScalarInterpolate = Geometry::TrilinearInterpolateScalarValue;
            m_VectorInterpolate = Geometry::TrilinearInterpolateVectorValue;
        }
        else
        {
            m_ScalarInterpolate = Geometry::GetNearestNeighborScalarValue;
            m_VectorInterpolate = Geometry::GetNearestNeighborVectorValue;
        }
    }

    /**
     * \brief Preprocess for evolution, i.e.: construct the evolving manifolds, and transform the target data's distance field, and the DF's normalized neg gradient for stabilization.
     */
    void Preprocess() override;

    /**
     * \brief Performs a single step of manifold evolution from the configuration at previous time step.
     */
    void PerformEvolutionStep(unsigned int stepId) override;

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
    void ExportCurrentState(unsigned int step, const std::string& baseOutputFilename) override;

    /**
     * \brief Exports the final state of evolving manifold(s).
     */
    void ExportFinalResult(const std::string& baseOutputFilename) override;

    /**
     * \brief Exports the target point cloud distance field as vti image data (if defined).
     */
    void ExportTargetDistanceFieldAsImage(const std::string& baseOutputFilename) override;

    /// \brief A specialized external getter which also transforms the stabilized outer surface to its original scale.
    [[nodiscard]] std::shared_ptr<pmp::SurfaceMesh> GetOuterSurfaceInOrigScale() const;

    /// \brief A specialized external getter which also transforms the stabilized inner surfaces to its original scale.
    [[nodiscard]] std::vector<std::shared_ptr<pmp::SurfaceMesh>> GetInnerSurfacesInOrigScale() const;

    /// \brief Calculate m_OuterSurfaceDistanceField and m_InnerSurfacesDistanceFields if there should be an interaction between these manifolds.
    void ComputeVariableDistanceFields() override;

protected:

    /// \brief Semi-implicit integration method step.
    void SemiImplicitIntegrationStep(unsigned int step) override;

    /// \brief Explicit integration method step.
    void ExplicitIntegrationStep(unsigned int step) override;

    /// \brief Prepares the property arrays for all evolving manifolds.
    void PrepareManifoldProperties() override;

    /// \brief A getter for the outer surface
    std::shared_ptr<pmp::SurfaceMesh>& GetOuterSurface()
    {
        return m_OuterSurface;
    }

    /// \brief A getter for the vector of inner surfaces
    std::vector<std::shared_ptr<pmp::SurfaceMesh>>& GetInnerSurfaces()
    {
        return m_InnerSurfaces;
    }

    /// \brief A const getter for the outer surface
    [[nodiscard]] const std::shared_ptr<pmp::SurfaceMesh>& GetOuterSurface() const
    {
        return m_OuterSurface;
    }

    /// \brief A const getter for the vector of inner surfaces
    [[nodiscard]] const std::vector<std::shared_ptr<pmp::SurfaceMesh>>& GetInnerSurfaces() const
    {
        return m_InnerSurfaces;
    }

    /// \brief Compute fields m_DistanceField and m_DFNegNormalizedGradient.
    /// \return triple { minTargetSize, maxTargetSize, targetBoundsCenter }.
    [[nodiscard]] std::tuple<float, float, pmp::Point> ComputeAmbientFields();

    /// \brief Construct m_OuterSurface and m_InnerSurfaces from settings.
    /// \param[in] minTargetSize        minimal size of the target data bounding box. Used for computing the radius of the outer manifold.
    /// \param[in] maxTargetSize        maximal size of the target data bounding box. Used for computing the radius of the outer manifold.
    /// \param[in] targetBoundsCenter   the center of the target data bounding box. Used for proper centering the initial outer manifold.
    /// \return radius of outer sphere.
    [[nodiscard]] float ConstructInitialManifolds(float minTargetSize, float maxTargetSize, const pmp::Point& targetBoundsCenter);

    /// \brief Transform all of the geometries so that numerical stability is ensured.
    /// \param[in] outerRadius             radius of outer sphere to calculate the approx co-volume measure.
    /// \param[in] stabilizationFactor     a multiplier for stabilizing mean co-volume measure.
    void StabilizeGeometries(float outerRadius, float stabilizationFactor = 1.0f) override;

    /// \brief A getter for the inverse stabilization transformation matrix.
    pmp::mat4& GetTransformToOriginal()
    {
        return m_TransformToOriginal;
    }

    /// \brief A getter for the distance field scalar grid.
    std::shared_ptr<Geometry::ScalarGrid>& GetDistanceField()
    {
        return m_DistanceField;
    }

    /// \brief A getter for the distance field normalized negative gradient vector grid.
    std::shared_ptr<Geometry::VectorGrid>& GetDFNegNormalizedGradient()
    {
        return m_DFNegNormalizedGradient;
    }

    /// \brief A getter for the area function of the chosen Laplacian scheme.
    AreaFunction& GetLaplacianAreaFunction()
    {
        return m_LaplacianAreaFunction;
    }

    /// \brief A const getter for the area function of the chosen Laplacian scheme.
    [[nodiscard]] const AreaFunction& GetLaplacianAreaFunction() const
    {
        return m_LaplacianAreaFunction;
    }

    /// \brief A getter for the scalar interpolator function.
    ScalarGridInterpolationFunction3D& GetScalarInterpolate()
    {
        return m_ScalarInterpolate;
    }

    /// \brief A getter for the vector interpolator function.
    VectorGridInterpolationFunction3D& GetVectorInterpolate()
    {
        return m_VectorInterpolate;
    }

    /// \brief A getter for the test box for numerical validity.
    pmp::BoundingBox& GetEvolBox()
    {
        return m_EvolBox;
    }

    /// \brief A getter for the distance field to outer surface
    std::shared_ptr<Geometry::ScalarGrid>& GetOuterSurfraceDistanceField()
    {
        return m_OuterSurfaceDistanceField;
    }

    /// \brief A getter for the distance fields of inner surfaces
    std::vector<std::shared_ptr<Geometry::ScalarGrid>>& GetInnerSurfacesDistanceFields()
    {
        return m_InnerSurfacesDistanceFields;
    }

private:

    std::shared_ptr<std::vector<pmp::Point>> m_TargetPointCloud{ nullptr }; //>! target point cloud geometry representing the spatial data for the reconstructed manifold.

    std::shared_ptr<pmp::SurfaceMesh> m_OuterSurface{ nullptr }; //>! evolving outer manifold which evolves inward, attempting to shrink-wrap m_TargetGeometryAdapter if it's defined.
    std::vector<std::shared_ptr<pmp::SurfaceMesh>> m_InnerSurfaces{}; //>! evolving inner manifolds which evolves outward

    std::shared_ptr<Geometry::ScalarGrid> m_DistanceField{ nullptr }; //>! the computed distance field of m_TargetPointCloud on a 3D scalar grid.
    std::shared_ptr<Geometry::VectorGrid> m_DFNegNormalizedGradient{ nullptr }; //>! the normalized negative gradient of m_DistanceField.

    // TODO: investigate performance
    std::shared_ptr<Geometry::ScalarGrid> m_OuterSurfaceDistanceField{ nullptr }; //>! the updated distance field of the evolving outer surface.
    std::vector<std::shared_ptr<Geometry::ScalarGrid>> m_InnerSurfacesDistanceFields{}; //>! the updated distance fields of the evolving inner surfaces.

    pmp::mat4 m_TransformToOriginal = pmp::mat4::identity(); //>! a transformation matrix to transform the stabilized geometry back to its original scale.

    AreaFunction m_LaplacianAreaFunction{}; //>! a function for calculating co-volume areas (see [Meyer, Desbrun, Schroder, Barr, 2003])

    ScalarGridInterpolationFunction3D m_ScalarInterpolate{}; //>! a parametrizeable function for interpolating values within Geometry::ScalarGrid2D.
    VectorGridInterpolationFunction3D m_VectorInterpolate{};  //>! a parametrizeable function for interpolating vector values within Geometry::VectorGrid2D.

    std::function<pmp::ImplicitLaplaceInfo(const pmp::SurfaceMesh& /* mesh */, pmp::Vertex /* v */)> m_ImplicitLaplacianFunction{}; //>! a Laplacian function chosen from parameter laplacianType.
    std::function<pmp::Point(const pmp::SurfaceMesh& /* mesh */, pmp::Vertex /* v */)> m_ExplicitLaplacianFunction{}; //>! a Laplacian function chosen from parameter laplacianType.

    pmp::BoundingBox m_EvolBox{}; //>! the test box for numerical validity of the evolution.
};

/**
 * \brief An internal functionality of the ManifoldEvolver handling custom pmp::ManifoldCurve2Ds and 2D target set and its distance represented by 2D scalar and vector field.
 * \class CustomManifoldSurfaceEvolutionStrategy
 */
class CustomManifoldSurfaceEvolutionStrategy : public ManifoldSurfaceEvolutionStrategy
{
public:
    /**
     * \brief Constructor that accepts custom outer and inner surfaces.
     * \param settings         evolution strategy settings.
     * \param laplacianType    the type of Laplacian area element used during calculations.
     * \param outerSurface     Custom outer surface.
     * \param innerSurfaces    Vector of custom inner surfaces.
     * \param targetPointCloud target point cloud.
     */
    CustomManifoldSurfaceEvolutionStrategy(
        ManifoldEvolutionSettings settings, 
        MeshLaplacian laplacianType,
        pmp::SurfaceMesh& outerSurface,
        std::vector<pmp::SurfaceMesh>& innerSurfaces,
        std::shared_ptr<std::vector<pmp::Point>> targetPointCloud = nullptr)
        : ManifoldSurfaceEvolutionStrategy(settings, laplacianType, std::move(targetPointCloud))
    {
        GetOuterSurface() = std::make_shared<pmp::SurfaceMesh>(outerSurface);
        for (auto& s : innerSurfaces)
            GetInnerSurfaces().emplace_back(std::make_shared<pmp::SurfaceMesh>(s));
    }

    /**
     * \brief Special overridden preprocessing method which omits construction of inner/outer curves.
     */
    void Preprocess() override;

private:

    /// \brief Verifies whether all custom inner surfaces are contained within the custom outer surface if not, throw std::invalid_argument
    [[nodiscard]] bool HasValidInnerOuterManifolds() const;

    /// \brief Computes the full range of co-volume sizes (lengths) to help compute the stabilization scaling factor.
    [[nodiscard]] std::pair<float, float> CalculateCoVolumeRange() const;

    /// \brief Transform all of the geometries so that numerical stability is ensured.
    /// \param[in] minArea                 the minimum area of a 2D co-volume within the custom surfaces.
    /// \param[in] maxArea                 the maximum area of a 2D co-volume within the custom surfaces.
    /// \param[in] stabilizationFactor     a multiplier for stabilizing mean co-volume measure.
    void StabilizeCustomGeometries(float minArea, float maxArea, float stabilizationFactor = 1.0f);
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

	ManifoldEvolver(GlobalManifoldEvolutionSettings settings, std::shared_ptr<ManifoldEvolutionStrategy> strategy)
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
        if (m_Settings.ExportTargetDistanceFieldAsImage)
        {
            m_Strategy->ExportTargetDistanceFieldAsImage(m_Settings.OutputPath + m_Settings.ProcedureName);
        }

        if (m_Settings.ExportPerTimeStep)
        {
            m_Strategy->ExportCurrentState(0, m_Settings.OutputPath + m_Settings.ProcedureName);
        }

        for (unsigned int step = 1; step <= m_Settings.NSteps; ++step)
        {
            m_Strategy->PerformEvolutionStep(step);

            if (m_Settings.DoRemeshing && m_Strategy->ShouldRemesh())
            {
                m_Strategy->Remesh();
            }

            m_Strategy->ComputeVariableDistanceFields();

            if (m_Settings.ExportPerTimeStep)
            {
                m_Strategy->ExportCurrentState(step, m_Settings.OutputPath + m_Settings.ProcedureName);
            }
        }

        if (m_Settings.ExportResult) 
        {
            m_Strategy->ExportFinalResult(m_Settings.OutputPath + m_Settings.ProcedureName);
        }
    }

    /// \brief A getter for this evolver's strategy.
    [[nodiscard]] std::shared_ptr<ManifoldEvolutionStrategy>& GetStrategy()
	{
        return m_Strategy;
	}

private:

    //
    // ======== Members ===============
    //

    // Key settings
    GlobalManifoldEvolutionSettings m_Settings;

    // Evolution strategy
    std::shared_ptr<ManifoldEvolutionStrategy> m_Strategy{ nullptr };
};