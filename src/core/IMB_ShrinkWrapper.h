#pragma once

#include "pmp/SurfaceMesh.h"

#include "pmp/algorithms/DifferentialGeometry.h"
#include "pmp/algorithms/Remeshing.h"

#include "geometry/GeometryUtil.h"
#include "geometry/Grid.h"

#include "EvolverUtilsCommon.h"

// --------------------------------------------------------------------------------------------
/// \brief Settings for the IMB_ShrinkWrapper.
/// \struct IMB_ShrinkWrapperSettings
// --------------------------------------------------------------------------------------------
struct IMB_ShrinkWrapperSettings
{
	std::string ProcedureName{}; //>! name for the evolution procedure.

	unsigned int LevelOfDetail{ 3 }; //>! Level of detail. Reflected in the number of vertices used during discretization. Standardized for icosphere.
	double TimeStep{ 0.01 };     //>! time step size.

	MeshLaplacian LaplacianType{ MeshLaplacian::Barycentric }; //>! the type of triangle mesh implicit Laplacian scheme.

	bool UseLinearGridInterpolation{ true }; //>! whether to use d-linear interpolation for scalar and vector fields where d is the grid dimension.
	AmbientFieldSettings FieldSettings{}; //>! the settings for the construction of the ambient fields.
	ManifoldAdaptiveRemeshingParams RemeshingSettings{}; //>! settings for pmp::Remeshing.
	FaceQualitySettings QualitySettings{}; //>! settings for mesh quality evaluation

	unsigned int MaxSteps{ 200 }; //>! the maximum number of allowed evolution steps.

	pmp::Scalar PointActivationRadius{ 0.5 }; //>! radius within which the point becomes activated for normal estimation.
	double ActivatedPointPercentageThreshold{ 0.8 }; //>! the fraction of target points to be within PointActivationRadius for the evolution to terminate.

	CurvatureCtrlFunction Epsilon{ TRIVIAL_EPSILON }; //>! control function for the curvature term of the shrink-wrapping surface.
	AdvectionCtrlFunction Eta{ TRIVIAL_ETA }; //>! control function for the advection term of the shrink-wrapping surface.

	double TangentialVelocityWeight{ 0.05 }; //>! the weight of tangential velocity update vector for each time step.
	double MaxFractionOfVerticesOutOfBounds{ 0.02 }; //>! fraction of vertices allowed to be out of bounds (because it will be decimated).

	PreStepFunction PreStep; //>! optional pre-step callback.
};

// --------------------------------------------------------------------------------------------
/// \brief An object for performing shrink-wrapping evolution towards a given point cloud, until a given percentage of this target data is covered.
/// \class IMB_ShrinkWrapper
// --------------------------------------------------------------------------------------------
class IMB_ShrinkWrapper
{
public:

	// -----------------------------------------------------------------------------
	/// \brief Constructor.
	// -----------------------------------------------------------------------------
	explicit IMB_ShrinkWrapper(const IMB_ShrinkWrapperSettings& settings, const std::vector<pmp::Point>& points);

	// -----------------------------------------------------------------------------
	/// \brief Main functionality.
	// -----------------------------------------------------------------------------
	[[nodiscard]] std::optional<pmp::SurfaceMesh> Perform();

private:

	// -----------------------------------------------------------------------------------------------------------------
	/// \brief Compute fields m_DistanceField and m_DFNegNormalizedGradient.
	/// \param[in/out] minTargetSize           Minimum bbox size of the target point cloud (m_Points).
	/// \param[in/out] maxTargetSize           Maximum bbox size of the target point cloud (m_Points).
	/// \param[in/out] targetBoundsCenter      Center of the bbox of the target point cloud (m_Points). 
	/// \return true if calculation successful.
	// -----------------------------------------------------------------------------------------------------------------
	[[nodiscard]] bool ComputeAmbientFields(pmp::Scalar& minTargetSize, pmp::Scalar& maxTargetSize, pmp::Point& targetBoundsCenter);

	// -----------------------------------------------------------------------------------------------------------------
	/// \brief Constructs the initial ico-sphere surface from the given parameters.
	/// \param[in] minTargetSize           Minimum bbox size of the target point cloud (m_Points).
	/// \param[in] maxTargetSize           Maximum bbox size of the target point cloud (m_Points).
	/// \param[in] targetBoundsCenter      Center of the bbox of the target point cloud (m_Points). 
	// -----------------------------------------------------------------------------------------------------------------
	void ConstructInitialSurface(const pmp::Scalar& minTargetSize, const pmp::Scalar& maxTargetSize, const pmp::Point& targetBoundsCenter);

	// -----------------------------------------------------------------------------------------------
	/// \brief Evaluates and assigns remeshing settings for the particular initial surface.
	// -----------------------------------------------------------------------------------------------
	void AssignRemeshingSettings();

	//
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//

	/**
	 * \brief Preprocess for evolution, i.e.: construct the evolving manifolds, and transform the target data's distance field, and the DF's normalized neg gradient for stabilization.
	 * \return false if the process fails.
	 */
	[[nodiscard]] bool Preprocess();

	// -----------------------------------------------------------------------------
	/// \brief Transform all of the geometries so that numerical stability is ensured.
	/// \param[in] minTargetSize           Minimum bbox size of the target point cloud (m_Points).
	/// \param[in] stabilizationFactor     A multiplier for stabilizing mean co-volume measure.
	// -----------------------------------------------------------------------------
	void StabilizeGeometries(const pmp::Scalar& minTargetSize, const pmp::Scalar& stabilizationFactor = 1.0);

	/**
	 * \brief Perform remeshing on the evolving m_Surface.
	 */
	void Remesh();

	/**
	 * \brief Performs a single step of manifold evolution from the configuration at previous time step.
	 * \return false if an error is encountered.
	 */
	[[nodiscard]] bool PerformEvolutionStep(unsigned int stepId);

	//
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//

	// -----------------------------------------------------------------------------
	/// \brief Calculates the percentage of how many m_Points are wrapped by m_Surface within m_Settings.PointActivationRadius.
	// -----------------------------------------------------------------------------
	[[nodiscard]] double GetVertexCoveragePercentage() const;

	//
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//

	// -----------------------------------------------------------------------------
	/// An internal optional getter.
	// -----------------------------------------------------------------------------
	[[nodiscard]] std::optional<pmp::SurfaceMesh> GetSurface() const;

	//
	// ========================================================
	//

	IMB_ShrinkWrapperSettings m_Settings; //>! settings for this shrink-wrapper instance.

	const std::vector<pmp::Point>& m_Points; //>! Input point cloud.

	std::shared_ptr<Geometry::ScalarGrid> m_DistanceField{ nullptr }; //>! Distance field to m_Points.
	std::shared_ptr<Geometry::VectorGrid> m_DFNegNormalizedGradient{ nullptr }; //>! Normalized negative gradient of m_DistanceField.

	std::shared_ptr<pmp::SurfaceMesh> m_Surface{ nullptr }; //>! the evolving surface shrink-wrapping m_Points.

	pmp::mat4 m_TransformToOriginal = pmp::mat4::identity(); //>! a transformation matrix to transform the stabilized geometry back to its original scale and translate to the target center.

	AreaFunction m_LaplacianAreaFunction{}; //>! a function for calculating co-volume areas (see [Meyer, Desbrun, Schroder, Barr, 2003])
	std::function<pmp::ImplicitLaplaceInfo(const pmp::SurfaceMesh& /* mesh */, pmp::Vertex /* v */)> m_ImplicitLaplacianFunction{}; //>! a Laplacian function chosen from parameter laplacianType.

	pmp::Scalar m_ScalingFactor{ 1.0 }; //>! the computed scaling factor for numerical stabilization.

	ScalarGridInterpolationFunction3D m_ScalarInterpolate{}; //>! a parametrizeable function for interpolating values within Geometry::ScalarGrid2D.
	VectorGridInterpolationFunction3D m_VectorInterpolate{};  //>! a parametrizeable function for interpolating vector values within Geometry::VectorGrid2D.

	Geometry::Sphere3D m_InitialSphere; //>! parameters for the initial sphere.

	pmp::AdaptiveRemeshingSettings m_RemeshingSettingsCustom; //>! remeshing settings customized for the initial surface.
	bool m_ShouldRemesh{ false }; //>! a flag turned on whenever remeshing is necessary after a time step.

	pmp::BoundingBox m_EvolBox{}; //>! the test box for numerical validity of the evolution.
};