#pragma once

#include "pmp/algorithms/DifferentialGeometry.h"
#include "pmp/algorithms/Remeshing.h"
#include "pmp/algorithms/ArcLengthCalculator.h"

#include "geometry/Grid.h"
#include "geometry/GridUtil.h"

#include "InscribedManifold.h"
#include "EvolverUtilsCommon.h"

//
// ===============================================================================================
//                                          Settings
// -----------------------------------------------------------------------------------------------
//

/// \brief the smallest allowed number of vertices in a manifold curve.
constexpr unsigned int N_CIRCLE_VERTS_0{ 5 };

/**
 * \brief Settings for printing out debug info into a specified log file.
 * \struct DiagnosticSettings
 */
struct DiagnosticSettings
{
    bool LogOuterManifoldEpsilon{ false }; //>! whether to log values of epsilon (curvature control function) for each vertex and each time step for the outer manifold.
    bool LogInnerManifoldsEpsilon{ false }; //>! whether to log values of epsilon (curvature control function) for each vertex and each time step for each inner manifold.
    bool LogOuterManifoldEta{ false };  //>! whether to log values of eta (advection control function) for each vertex and each time step for the outer manifold.
    bool LogInnerManifoldsEta{ false }; //>! whether to log values of eta (advection control function) for each vertex and each time step for each inner manifold.

    //bool LogOuterManifoldDisplacements{ false };
    //bool LogInnerManifoldsDisplacements{ false };

    // numerical diagnostics
    bool LogOuterManifoldErrors{ false };
    bool LogInnerManifoldsErrors{ false };

    bool LogOuterManifoldXErrors{ false };
    bool LogOuterManifoldYErrors{ false };
    bool LogInnerManifoldsXErrors{ false };
    bool LogInnerManifoldsYErrors{ false };
};

/**
 * \brief A wrapper for manifold evolution settings.
 * \struct ManifoldEvolutionSettings
 */
struct ManifoldEvolutionSettings
{
    bool UseInnerManifolds{ true }; //>! whether to construct outward-evolving inner manifolds.
    bool UseOuterManifolds{ true }; //>! whether to construct inward-evolving outer manifolds.

	bool UseSemiImplicit{ true }; //>! whether to use a more numerically stable, but computationally costly numerical scheme (has to solve a linear system for each coord for each time step).
    unsigned int LevelOfDetail{ 3 }; //>! Level of detail. Reflected in the number of vertices used during discretization. Standardized for circles and icospheres, e.g.: the base circle is a regular pentagon, and the sphere an icosahedron.

    double TimeStep{ 0.01 };     //>! time step size.
    bool UseStabilizationViaScaling{ true }; //>! whether to scale all evolving geometries to match with the time step.

    bool UseLinearGridInterpolation{ true }; //>! whether to use d-linear interpolation for scalar and vector fields where d is the grid dimension.

    DistanceSelectionType DistanceSelection{ DistanceSelectionType::PlainMinimum }; //>! The type of selection function when evaluating multiple interaction distances.
    double DistanceBlendingRadius{ 0.0 }; //>! the radius for blending interaction distances.
    double InteractionDistanceRatio{ 1.0 }; //>! the ratio between actual interaction distances (other manifold : target) to make them more-or-less equal.

    NormalActivationSettings NormalActivation; //>! settings for the normal activation procedure.

	CurvatureCtrlFunction OuterManifoldEpsilon{ TRIVIAL_EPSILON }; //>! control function for the curvature term of the outer manifold.
    AdvectionCtrlFunction OuterManifoldEta{ TRIVIAL_ETA }; //>! control function for the advection term of the outer manifold.
    RepulsionFunction     OuterManifoldRepulsion{ TRIVIAL_REPULSION }; //>! repulsion control function for the outer manifold.

    CurvatureCtrlFunction InnerManifoldEpsilon{ TRIVIAL_EPSILON }; //>! control function for the curvature term of the inner manifolds.
    AdvectionCtrlFunction InnerManifoldEta{ TRIVIAL_ETA }; //>! control function for the advection term of the inner manifolds.
    RepulsionFunction     InnerManifoldRepulsion{ TRIVIAL_REPULSION }; //>! repulsion control function for the inner manifolds.

    bool AdvectionInteractWithOtherManifolds{ false }; //>! whether to use the minimum from the distances to all manifolds including target data.

    AmbientFieldSettings FieldSettings{}; //>! the settings for the construction of the ambient fields.

    double TangentialVelocityWeight{ 0.05 }; //>! the weight of tangential velocity update vector for each time step.

    double MaxFractionOfVerticesOutOfBounds{ 0.02 }; //>! fraction of vertices allowed to be out of bounds (because it will be decimated).

    FaceQualitySettings QualitySettings{}; //>! settings for mesh quality evaluation
    ManifoldAdaptiveRemeshingParams RemeshingSettings{}; //>! settings for pmp::Remeshing.
    FeatureDetectionSettings FeatureSettings{}; //>! settings for feature detection.

    bool ExportVariableScalarFieldsDimInfo{ false }; //>! whether to export the dimensions of the variable scalar fields.
    bool ExportVariableVectorFieldsDimInfo{ false }; //>! whether to export the dimensions of the variable vector fields.

    DiagnosticSettings DiagSettings{}; //>! evolution diagnostics.
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
    std::unordered_set<unsigned int> RemeshingResizeTimeIds{}; //>! a list of time indices during which remeshing lengths are resized by a given factor.
    pmp::Scalar RemeshingResizeFactor{ 0.98 }; //>! a factor by which remeshing edge lengths are downsized when a particular time step (logged in RemeshingResizeTimeIds) is reached. DISCLAIMER: This value is not used if RemeshingResizeTimeIds is empty!

    bool DetectFeatures{ true }; //>! if true, vertices with critical curvature will be marked as feature and frozen with respect to remeshing to avoid evolution past the target.
};

/// \brief Stabilization weight param from [0, 1]
constexpr pmp::Scalar STABILIZATION_FACTOR{ 1.0 };

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
     * \brief Separate logging method for initializing value buffers for logging.
     */
    virtual void InitLogger(const std::string& baseOutputFileName) = 0;

    /**
     * \brief Separate method for passing the resulting logged data into the log file.
     * \param[in] omitLastTimeStep        a flag useful when an exception has been thrown during some time step.
     */
    virtual void SaveLog(bool omitLastTimeStep = false) = 0;

    /**
     * \brief Preprocess for evolution, i.e.: construct the evolving manifolds, and transform the target data's distance field, and the DF's normalized neg gradient for stabilization.
     */
    virtual void Preprocess() = 0;

    /**
     * \brief Cleanup procedures after evolution is finished.
     */
    virtual void Postprocess() = 0;

    /**
     * \brief Performs a single step of manifold evolution from the configuration at previous time step.
     */
    virtual void PerformEvolutionStep(unsigned int stepId) = 0;

    /**
     * \brief Perform remeshing on the evolving manifold(s).
     */
    virtual void Remesh() = 0;

    /**
     * \brief Resizes remeshing settings the evolving manifold(s) by a given factor.
     */
    virtual void ResizeRemeshingSettings(pmp::Scalar resizeFactor) = 0;

    /**
     * \brief Marks essential vertices that should not be displaced by remeshing as feature.
     */
    virtual void DetectFeatures() = 0;

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

    /**
     * \brief Initializing value buffers for logging a new time step.
     */
    virtual void InitNewTimeStepLog(unsigned int stepId) = 0;

    /// \brief A getter for the stabilization scaling factor.
    pmp::Scalar& GetScalingFactor()
    {
        return m_ScalingFactor;
    }

    /// \brief Transform all of the geometries so that numerical stability is ensured.
    /// \param[in] stabilizationFactor     a multiplier for stabilizing mean co-volume measure.
    virtual void StabilizeGeometries(pmp::Scalar stabilizationFactor = STABILIZATION_FACTOR) = 0;

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

    /// \brief Calculates and assigns remeshing settings to the strategy wrapper for remeshing settings.
    virtual void AssignRemeshingSettingsToEvolvingManifolds() = 0;

    /// \brief Verifies whether variable fields need to be defined.
    virtual [[nodiscard]] bool NeedsVariableFieldsCalculation() = 0;

    /// \brief Verifies whether target or variable fields need to be defined.
    virtual [[nodiscard]] bool NeedsFieldsCalculation() = 0;

    /// \brief Checks the diagnostic flags in m_Settings.
    [[nodiscard]] bool LogOuterManifoldValues() const;

    /// \brief Checks the diagnostic flags in m_Settings.
    [[nodiscard]] bool LogInnerManifoldValues() const;

    /// \brief Checks the diagnostic flags in m_Settings.
    [[nodiscard]] bool LogManifoldValues() const
    {
        return LogOuterManifoldValues() || LogInnerManifoldValues();
    }

    /// \brief Computes lower bounds for control functions from maximum available distance range after scaling by m_ScalingFactor.
    virtual void ComputeControlFunctionsLowerBounds() = 0;

private:

    ManifoldEvolutionSettings m_Settings{};       //>! settings for the evolution strategy.
    NumericalStepIntegrateFunction m_Integrate{}; //>! numerical integration function (clearly, derived classes have different coefficients/matrices).

    pmp::Scalar m_ScalingFactor{ 1.0 }; //>! the computed scaling factor for numerical stabilization.
    pmp::Scalar m_FieldCellSize{ 0.1 }; //>! standardized cell size for the variable distance fields.
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
	 * \param settings               evolution strategy settings.
	 * \param targetPointCloud       target point cloud data (optional).
     * \param targetDistanceField    precomputed target distance field (optional).
	 */
	explicit ManifoldCurveEvolutionStrategy(ManifoldEvolutionSettings settings, 
	                                        std::shared_ptr<std::vector<pmp::Point2>> targetPointCloud = nullptr,
                                            std::shared_ptr<Geometry::ScalarGrid2D> targetDistanceField = nullptr)
	    : ManifoldEvolutionStrategy(settings),
        m_TargetPointCloud(std::move(targetPointCloud)),
        m_DistanceField(std::move(targetDistanceField))
    {
        if (GetSettings().UseSemiImplicit)
        {
            GetIntegrate() = std::bind(&ManifoldCurveEvolutionStrategy::SemiImplicitIntegrationStep, this, std::placeholders::_1);
        }
        else
        {
            GetIntegrate() = std::bind(&ManifoldCurveEvolutionStrategy::ExplicitIntegrationStep, this, std::placeholders::_1);
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
     * \brief Separate logging method for initializing value buffers for logging.
     */
    void InitLogger(const std::string& baseOutputFileName) override;

    /**
     * \brief Separate method for passing the resulting logged data into the log file.
     * \param[in] omitLastTimeStep        a flag useful when an exception has been thrown during some time step.
     */
    void SaveLog(bool omitLastTimeStep = false) override;

    /**
     * \brief Preprocess for evolution, i.e.: construct the evolving manifolds, and transform the target data's distance field, and the DF's normalized neg gradient for stabilization.
     */
    void Preprocess() override;

    /**
     * \brief Cleanup procedures after evolution is finished.
     */
    void Postprocess() override;

    /**
     * \brief Performs a single step of manifold evolution from the configuration at previous time step.
     */
    void PerformEvolutionStep(unsigned int stepId) override;

    /**
     * \brief Perform remeshing on the evolving manifold(s) logged in the m_RemeshTracker.
     */
    void Remesh() override;

    /**
	 * \brief Resizes remeshing settings the evolving manifold(s) by a given factor.
	 */
    void ResizeRemeshingSettings(pmp::Scalar resizeFactor) override;

    /**
	 * \brief Marks essential vertices that should not be displaced by remeshing as feature.
	 */
    void DetectFeatures() override;

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

    /**
     * \brief Initializing value buffers for logging a new time step.
     */
    void InitNewTimeStepLog(unsigned int stepId) override;

    /// \brief Semi-implicit integration method step.
    void SemiImplicitIntegrationStep(unsigned int step) override;

    /// \brief Explicit integration method step.
    void ExplicitIntegrationStep(unsigned int step) override;

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

    /// \brief  A getter for the inverse stabilization scaling matrix.
    pmp::mat3& GetScaleFieldToOriginal()
	{
		return m_ScaleFieldToOriginal;
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
    [[nodiscard]] std::tuple<pmp::Scalar, pmp::Scalar, pmp::Point2> ComputeAmbientFields();

    /// \brief Construct m_OuterCurve and m_InnerCurves from settings.
    /// \param[in] minTargetSize        minimal size of the target data bounding box. Used for computing the radius of the outer manifold.
    /// \param[in] maxTargetSize        maximal size of the target data bounding box. Used for computing the radius of the outer manifold.
    /// \param[in] targetBoundsCenter   the center of the target data bounding box. Used for proper centering the initial outer manifold.
    void ConstructInitialManifolds(pmp::Scalar minTargetSize, pmp::Scalar maxTargetSize, const pmp::Point2& targetBoundsCenter);

    /// \brief Transform all of the geometries so that numerical stability is ensured.
    /// \param[in] stabilizationFactor     a multiplier for stabilizing mean co-volume measure.
    void StabilizeGeometries(pmp::Scalar stabilizationFactor = STABILIZATION_FACTOR) override;

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

    /// \brief A getter for the negative normalized gradient of the distance field to outer curve
    std::shared_ptr<Geometry::VectorGrid2D>& GetOuterCurveDFNegNormalizedGradient()
    {
        return m_OuterCurveDFNegNormalizedGradient;
    }

    /// \brief A getter for the distance fields of inner curves
    std::vector<std::shared_ptr<Geometry::ScalarGrid2D>>& GetInnerCurvesDistanceFields()
    {
        return m_InnerCurvesDistanceFields;
    }

    /// \brief A getter for the negative normalized gradients of the distance fields to inner curves
    std::vector<std::shared_ptr<Geometry::VectorGrid2D>>& GetInnerCurvesDFNegNormalizedGradients()
    {
        return m_InnerCurvesDFNegNormalizedGradients;
    }

    /// \brief A getter for the distance blending strategy.
    std::shared_ptr<BaseDistanceBlendingStrategy<pmp::dvec2>>& GetDistBlendStrategy()
    {
        return m_DistBlendStrategy;
    }

    /// \brief Calculates and assigns remeshing settings to the strategy wrapper for remeshing settings.
    void AssignRemeshingSettingsToEvolvingManifolds() override;

    /// \brief Initializes the arc length calculation objects for all available evolving curves.
    void InitializeArcLengthCalculation();

    /// \brief A getter for the remeshing settings wrapper
    ManifoldRemeshingSettingsWrapper<pmp::ManifoldCurve2D>& GetRemeshingSettings()
    {
        return m_RemeshingSettings;
    }

    /// \brief Verifies whether variable fields need to be defined.
    [[nodiscard]] bool NeedsVariableFieldsCalculation() override;

    /// \brief Verifies whether target or variable fields need to be defined.
    [[nodiscard]] bool NeedsFieldsCalculation() override;

    /// \brief A getter for the diagnostic logger
    VertexValueLogger<pmp::ManifoldCurve2D>& GetLogger()
    {
        return m_Logger;
    }

    /// \brief Computes lower bounds for control functions from maximum available distance range after scaling by m_ScalingFactor.
    void ComputeControlFunctionsLowerBounds() override;

private:

    std::shared_ptr<std::vector<pmp::Point2>> m_TargetPointCloud{ nullptr }; //>! target point cloud geometry representing the spatial data for the reconstructed manifold.

    std::shared_ptr<pmp::ManifoldCurve2D> m_OuterCurve{ nullptr }; //>! evolving outer manifold which evolves inward, attempting to shrink-wrap m_TargetGeometryAdapter if it's defined.
    std::vector<std::shared_ptr<pmp::ManifoldCurve2D>> m_InnerCurves{}; //>! evolving inner manifolds which evolves outward

    std::shared_ptr<Geometry::ScalarGrid2D> m_DistanceField{ nullptr }; //>! the computed distance field of m_TargetPointCloud on a 2D scalar grid.
    std::shared_ptr<Geometry::VectorGrid2D> m_DFNegNormalizedGradient{ nullptr }; //>! the normalized negative gradient of m_DistanceField.

    // TODO: investigate performance
    std::shared_ptr<Geometry::ScalarGrid2D> m_OuterCurveDistanceField{ nullptr }; //>! the updated distance field of the evolving outer manifold.
    std::shared_ptr<Geometry::VectorGrid2D> m_OuterCurveDFNegNormalizedGradient{ nullptr }; //>! the updated gradient of distance field to evolving outer manifold.
    std::vector<std::shared_ptr<Geometry::ScalarGrid2D>> m_InnerCurvesDistanceFields{}; //>! the updated distance fields of the evolving inner manifolds.
    std::vector<std::shared_ptr<Geometry::VectorGrid2D>> m_InnerCurvesDFNegNormalizedGradients{}; //>! the updated gradients of distance fields to evolving inner manifolds.

    pmp::mat3 m_TransformToOriginal = pmp::mat3::identity(); //>! a transformation matrix to transform the stabilized geometry back to its original scale and translate to the target center.
    pmp::mat3 m_ScaleFieldToOriginal = pmp::mat3::identity(); //>! a transformation matrix to transform the stabilized geometry back to its original scale.

    ScalarGridInterpolationFunction2D m_ScalarInterpolate{}; //>! a parametrizeable function for interpolating values within Geometry::ScalarGrid2D.
    VectorGridInterpolationFunction2D m_VectorInterpolate{};  //>! a parametrizeable function for interpolating vector values within Geometry::VectorGrid2D.

    std::shared_ptr<BaseDistanceBlendingStrategy<pmp::dvec2>> m_DistBlendStrategy = std::make_shared<PlainMinimumStrategy<pmp::dvec2>>(); //>! blending strategy for interaction distance.

    pmp::BoundingBox2 m_EvolBox{}; //>! the test box for numerical validity of the evolution.

    ManifoldsToRemeshTracker<pmp::ManifoldCurve2D> m_RemeshTracker{}; //>! a utility which logs curves that need remeshing.
    ManifoldRemeshingSettingsWrapper<pmp::ManifoldCurve2D> m_RemeshingSettings{}; //>! a wrapper with remeshing settings assigned to evolving manifolds
    InitialSphereSettingsWrapper<pmp::ManifoldCurve2D, Geometry::Circle2D> m_InitialSphereSettings{}; //>! a wrapper for the initial sphere settings for each manifold.

    VertexValueLogger<pmp::ManifoldCurve2D> m_Logger{}; //>! a utility for exporting the chosen vertex values to a file for debugging purposes.

    std::unordered_map<pmp::ManifoldCurve2D*, std::shared_ptr<pmp::EvolvingArcLengthCalculator>> m_ArcLengthCalculators{}; // Calculators for the arc length of vertex positions of the evolving curves.
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
     * \param settings               evolution strategy settings.
     * \param outerCurve             Custom outer curve (optional).
     * \param innerCurves            Vector of custom inner curves.
     * \param targetPointCloud       target point cloud (optional).
     * \param targetDistanceField    precomputed target distance field (optional).
     */
    explicit CustomManifoldCurveEvolutionStrategy(
        ManifoldEvolutionSettings settings,
        std::optional<pmp::ManifoldCurve2D> outerCurve,
        std::vector<pmp::ManifoldCurve2D>& innerCurves,
        std::shared_ptr<std::vector<pmp::Point2>> targetPointCloud = nullptr,
        std::shared_ptr<Geometry::ScalarGrid2D> targetDistanceField = nullptr)
        : ManifoldCurveEvolutionStrategy(settings, std::move(targetPointCloud), std::move(targetDistanceField))
    {
        GetOuterCurve() = outerCurve.has_value() ? std::make_shared<pmp::ManifoldCurve2D>(*outerCurve) : nullptr;
        for (auto& c : innerCurves)
            GetInnerCurves().emplace_back(std::make_shared<pmp::ManifoldCurve2D>(c));
    }

    /**
     * \brief Special overridden preprocessing method which omits construction of inner/outer curves.
     */
    void Preprocess() override;

private:

    /// \brief Calculates and assigns remeshing settings to the strategy wrapper for remeshing settings.
    void AssignRemeshingSettingsToEvolvingManifolds() override;

    /// \brief Verifies whether all custom inner curves are contained within the custom outer curve if not, throw std::invalid_argument
    [[nodiscard]] bool HasValidInnerOuterManifolds() const;

    /// \brief Computes the full range of co-volume sizes (lengths) to help compute the stabilization scaling factor.
    [[nodiscard]] std::pair<pmp::Scalar, pmp::Scalar> CalculateCoVolumeRange() const;

    /// \brief Transform all of the geometries so that numerical stability is ensured.
	/// \param[in] minLength               the minimum length of a 1D co-volume within the custom curves.
	/// \param[in] maxLength               the maximum length of a 1D co-volume within the custom curves.
	/// \param[in] stabilizationFactor     a multiplier for stabilizing mean co-volume measure.
    void StabilizeCustomGeometries(pmp::Scalar minLength, pmp::Scalar maxLength, pmp::Scalar stabilizationFactor = 1.0);
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
            GetIntegrate() = std::bind(&ManifoldSurfaceEvolutionStrategy::SemiImplicitIntegrationStep, this, std::placeholders::_1);
        }
        else
        {
            GetIntegrate() = std::bind(&ManifoldSurfaceEvolutionStrategy::ExplicitIntegrationStep, this, std::placeholders::_1);
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
     * \brief Separate logging method for initializing value buffers for logging.
     */
    void InitLogger(const std::string& baseOutputFileName) override;

    /**
     * \brief Separate method for passing the resulting logged data into the log file.
     * \param[in] omitLastTimeStep        a flag useful when an exception has been thrown during some time step.
     */
    void SaveLog(bool omitLastTimeStep = false) override;

    /**
     * \brief Preprocess for evolution, i.e.: construct the evolving manifolds, and transform the target data's distance field, and the DF's normalized neg gradient for stabilization.
     */
    void Preprocess() override;

    /**
     * \brief Cleanup procedures after evolution is finished.
     */
    void Postprocess() override;

    /**
     * \brief Performs a single step of manifold evolution from the configuration at previous time step.
     */
    void PerformEvolutionStep(unsigned int stepId) override;

    /**
     * \brief Perform remeshing on the evolving manifold(s) logged in the m_RemeshTracker.
     */
    void Remesh() override;

    /**
	 * \brief Resizes remeshing settings the evolving manifold(s) by a given factor.
	 */
    void ResizeRemeshingSettings(pmp::Scalar resizeFactor) override;

    /**
	 * \brief Marks essential vertices that should not be displaced by remeshing as feature.
	 */
    void DetectFeatures() override;

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

    /**
     * \brief Initializing value buffers for logging a new time step.
     */
    void InitNewTimeStepLog(unsigned int stepId) override;

    /// \brief Semi-implicit integration method step.
    void SemiImplicitIntegrationStep(unsigned int step) override;

    /// \brief Explicit integration method step.
    void ExplicitIntegrationStep(unsigned int step) override;

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
    [[nodiscard]] std::tuple<pmp::Scalar, pmp::Scalar, pmp::Point> ComputeAmbientFields();

    /// \brief Construct m_OuterSurface and m_InnerSurfaces from settings.
    /// \param[in] minTargetSize        minimal size of the target data bounding box. Used for computing the radius of the outer manifold.
    /// \param[in] maxTargetSize        maximal size of the target data bounding box. Used for computing the radius of the outer manifold.
    /// \param[in] targetBoundsCenter   the center of the target data bounding box. Used for proper centering the initial outer manifold.
    void ConstructInitialManifolds(pmp::Scalar minTargetSize, pmp::Scalar maxTargetSize, const pmp::Point& targetBoundsCenter);

    /// \brief Transform all of the geometries so that numerical stability is ensured.
    /// \param[in] stabilizationFactor     a multiplier for stabilizing mean co-volume measure.
    void StabilizeGeometries(pmp::Scalar stabilizationFactor = STABILIZATION_FACTOR) override;

    /// \brief A getter for the inverse stabilization transformation matrix.
    pmp::mat4& GetTransformToOriginal()
    {
        return m_TransformToOriginal;
    }

    /// \brief  A getter for the inverse stabilization scaling matrix.
    pmp::mat4& GetScaleFieldToOriginal()
	{
		return m_ScaleFieldToOriginal;
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

    /// \brief A getter for the normalized negative gradient of the distance field to outer surface
    std::shared_ptr<Geometry::VectorGrid>& GetOuterSurfaceDFNegNormalizedGradient()
    {
        return m_OuterSurfaceDFNegNormalizedGradient;
    }

    /// \brief A getter for the distance fields of inner surfaces
    std::vector<std::shared_ptr<Geometry::ScalarGrid>>& GetInnerSurfacesDistanceFields()
    {
        return m_InnerSurfacesDistanceFields;
    }

    /// \brief A getter for the normalized negative gradients of the distance fields to inner surfaces
    std::vector<std::shared_ptr<Geometry::VectorGrid>>& GetInnerSurfacesDFNegNormalizedGradients()
    {
        return m_InnerSurfacesDFNegNormalizedGradients;
    }

    /// \brief A getter for the distance blending strategy.
    std::shared_ptr<BaseDistanceBlendingStrategy<pmp::dvec3>>& GetDistBlendStrategy()
    {
        return m_DistBlendStrategy;
    }

    /// \brief Calculates and assigns remeshing settings to the strategy wrapper for remeshing settings.
    void AssignRemeshingSettingsToEvolvingManifolds() override;

    /// \brief A getter for the remeshing settings wrapper
    ManifoldRemeshingSettingsWrapper<pmp::SurfaceMesh>& GetRemeshingSettings()
    {
        return m_RemeshingSettings;
    }

    /// \brief Verifies whether variable fields need to be defined.
    [[nodiscard]] bool NeedsVariableFieldsCalculation() override;

    /// \brief Verifies whether target or variable fields need to be defined.
    [[nodiscard]] bool NeedsFieldsCalculation() override;

    /// \brief A getter for the diagnostic logger
    VertexValueLogger<pmp::SurfaceMesh>& GetLogger()
    {
        return m_Logger;
    }

    /// \brief Computes lower bounds for control functions from maximum available distance range after scaling by m_ScalingFactor.
    void ComputeControlFunctionsLowerBounds() override;

private:

    std::shared_ptr<std::vector<pmp::Point>> m_TargetPointCloud{ nullptr }; //>! target point cloud geometry representing the spatial data for the reconstructed manifold.

    std::shared_ptr<pmp::SurfaceMesh> m_OuterSurface{ nullptr }; //>! evolving outer manifold which evolves inward, attempting to shrink-wrap m_TargetGeometryAdapter if it's defined.
    std::vector<std::shared_ptr<pmp::SurfaceMesh>> m_InnerSurfaces{}; //>! evolving inner manifolds which evolves outward

    std::shared_ptr<Geometry::ScalarGrid> m_DistanceField{ nullptr }; //>! the computed distance field of m_TargetPointCloud on a 3D scalar grid.
    std::shared_ptr<Geometry::VectorGrid> m_DFNegNormalizedGradient{ nullptr }; //>! the normalized negative gradient of m_DistanceField.

    // TODO: investigate performance
    std::shared_ptr<Geometry::ScalarGrid> m_OuterSurfaceDistanceField{ nullptr }; //>! the updated distance field of the evolving outer manifold.
    std::shared_ptr<Geometry::VectorGrid> m_OuterSurfaceDFNegNormalizedGradient{ nullptr }; //>! the updated gradient of distance field to evolving outer manifold.
    std::vector<std::shared_ptr<Geometry::ScalarGrid>> m_InnerSurfacesDistanceFields{}; //>! the updated distance fields of the evolving inner manifolds.
    std::vector<std::shared_ptr<Geometry::VectorGrid>> m_InnerSurfacesDFNegNormalizedGradients{}; //>! the updated gradients of distance fields to evolving inner manifolds.

    pmp::mat4 m_TransformToOriginal = pmp::mat4::identity(); //>! a transformation matrix to transform the stabilized geometry back to its original scale and translate to the target center.
    pmp::mat4 m_ScaleFieldToOriginal = pmp::mat4::identity(); //>! a transformation matrix to transform the stabilized geometry back to its original scale.

    AreaFunction m_LaplacianAreaFunction{}; //>! a function for calculating co-volume areas (see [Meyer, Desbrun, Schroder, Barr, 2003])

    ScalarGridInterpolationFunction3D m_ScalarInterpolate{}; //>! a parametrizeable function for interpolating values within Geometry::ScalarGrid2D.
    VectorGridInterpolationFunction3D m_VectorInterpolate{};  //>! a parametrizeable function for interpolating vector values within Geometry::VectorGrid2D.

    std::shared_ptr<BaseDistanceBlendingStrategy<pmp::dvec3>> m_DistBlendStrategy = std::make_shared<PlainMinimumStrategy<pmp::dvec3>>(); //>! blending strategy for interaction distance.

    std::function<pmp::ImplicitLaplaceInfo(const pmp::SurfaceMesh& /* mesh */, pmp::Vertex /* v */)> m_ImplicitLaplacianFunction{}; //>! a Laplacian function chosen from parameter laplacianType.
    std::function<pmp::Point(const pmp::SurfaceMesh& /* mesh */, pmp::Vertex /* v */)> m_ExplicitLaplacianFunction{}; //>! a Laplacian function chosen from parameter laplacianType.

    pmp::BoundingBox m_EvolBox{}; //>! the test box for numerical validity of the evolution.

    ManifoldsToRemeshTracker<pmp::SurfaceMesh> m_RemeshTracker{}; //>! a utility which logs surfaces that need remeshing.
    ManifoldRemeshingSettingsWrapper<pmp::SurfaceMesh> m_RemeshingSettings{}; //>! a wrapper with remeshing settings assigned to evolving manifolds
    InitialSphereSettingsWrapper<pmp::SurfaceMesh, Geometry::Sphere3D> m_InitialSphereSettings{}; //>! a wrapper for the initial sphere settings for each manifold.

    VertexValueLogger<pmp::SurfaceMesh> m_Logger{}; //>! a utility for exporting the chosen vertex values to a file for debugging purposes.
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
     * \param outerSurface     Custom outer surface (optional).
     * \param innerSurfaces    Vector of custom inner surfaces.
     * \param targetPointCloud target point cloud.
     */
    CustomManifoldSurfaceEvolutionStrategy(
        ManifoldEvolutionSettings settings, 
        MeshLaplacian laplacianType,
        std::optional<pmp::SurfaceMesh> outerSurface,
        std::vector<pmp::SurfaceMesh>& innerSurfaces,
        std::shared_ptr<std::vector<pmp::Point>> targetPointCloud = nullptr)
        : ManifoldSurfaceEvolutionStrategy(settings, laplacianType, std::move(targetPointCloud))
    {
        GetOuterSurface() = outerSurface.has_value() ? std::make_shared<pmp::SurfaceMesh>(*outerSurface) : nullptr;
        for (auto& s : innerSurfaces)
            GetInnerSurfaces().emplace_back(std::make_shared<pmp::SurfaceMesh>(s));
    }

    /**
     * \brief Special overridden preprocessing method which omits construction of inner/outer curves.
     */
    void Preprocess() override;

private:

    /// \brief Calculates and assigns remeshing settings to the strategy wrapper for remeshing settings.
    void AssignRemeshingSettingsToEvolvingManifolds() override;

    /// \brief Verifies whether all custom inner surfaces are contained within the custom outer surface if not, throw std::invalid_argument
    [[nodiscard]] bool HasValidInnerOuterManifolds() const;

    /// \brief Computes the full range of co-volume sizes (lengths) to help compute the stabilization scaling factor.
    [[nodiscard]] std::pair<pmp::Scalar, pmp::Scalar> CalculateCoVolumeRange() const;

    /// \brief Transform all of the geometries so that numerical stability is ensured.
    /// \param[in] minArea                 the minimum area of a 2D co-volume within the custom surfaces.
    /// \param[in] maxArea                 the maximum area of a 2D co-volume within the custom surfaces.
    /// \param[in] stabilizationFactor     a multiplier for stabilizing mean co-volume measure.
    void StabilizeCustomGeometries(pmp::Scalar minArea, pmp::Scalar maxArea, pmp::Scalar stabilizationFactor = STABILIZATION_FACTOR);
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
        m_Strategy->InitLogger(m_Settings.OutputPath + m_Settings.ProcedureName);
        m_Strategy->Preprocess();

        if (m_Settings.ExportTargetDistanceFieldAsImage)
        {
            m_Strategy->ExportTargetDistanceFieldAsImage(m_Settings.OutputPath + m_Settings.ProcedureName);
        }

        if (m_Settings.ExportPerTimeStep)
        {
            m_Strategy->ExportCurrentState(0, m_Settings.OutputPath + m_Settings.ProcedureName);
        }

        try
        {
            for (unsigned int step = 1; step <= m_Settings.NSteps; ++step)
            {
                // Print progress to the console
                std::cout << "\rManifoldEvolver::Evolve ... Step: " << step << "/" << m_Settings.NSteps << std::flush;

                m_Strategy->PerformEvolutionStep(step);

                if (m_Settings.DetectFeatures)
                {
                    m_Strategy->DetectFeatures();
                }

                if (m_Settings.DoRemeshing)
                {
                    if (m_Settings.RemeshingResizeTimeIds.contains(step))
                    {
	                    m_Strategy->ResizeRemeshingSettings(m_Settings.RemeshingResizeFactor);
                    }
                    m_Strategy->Remesh();
                }

                m_Strategy->ComputeVariableDistanceFields();

                if (m_Settings.ExportPerTimeStep)
                {
                    m_Strategy->ExportCurrentState(step, m_Settings.OutputPath + m_Settings.ProcedureName);
                }
            }
        }
        catch (std::invalid_argument& ex)
        {
            m_Strategy->SaveLog(true);
            std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Exiting... < < < < < \n";
            return;
        }
        catch (std::runtime_error& ex)
        {
            m_Strategy->SaveLog(true);
            std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Exiting... < < < < < \n";
            return;
        }
        catch (...)
        {
            m_Strategy->SaveLog(true);
            std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an unknown exception! Exiting... < < < < < \n";
            return;
        }

        m_Strategy->Postprocess();
        m_Strategy->SaveLog();

        if (m_Settings.ExportResult) 
        {
            m_Strategy->ExportFinalResult(m_Settings.OutputPath + m_Settings.ProcedureName);
        }

        // Clear the progress message after the loop ends
        std::cout << "\rManifoldEvolver::Evolve ... Step: " << m_Settings.NSteps << "/" << m_Settings.NSteps << ". Done.\n";
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