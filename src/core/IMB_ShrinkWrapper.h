#pragma once

#include "pmp/SurfaceMesh.h"

#include "pmp/algorithms/DifferentialGeometry.h"
#include "pmp/algorithms/Remeshing.h"

#include "geometry/Grid.h"

#include "EvolverUtilsCommon.h"


struct IMB_ShrinkWrapperSettings
{
	unsigned int LevelOfDetail{ 3 }; //>! Level of detail. Reflected in the number of vertices used during discretization. Standardized for icosphere.
	double TimeStep{ 0.01 };     //>! time step size.

	bool UseLinearGridInterpolation{ true }; //>! whether to use d-linear interpolation for scalar and vector fields where d is the grid dimension.
	AmbientFieldSettings FieldSettings{}; //>! the settings for the construction of the ambient fields.
	ManifoldAdaptiveRemeshingParams RemeshingSettings{}; //>! settings for pmp::Remeshing.
	FaceQualitySettings QualitySettings{}; //>! settings for mesh quality evaluation

	unsigned int MaxSteps{ 200 }; //>! the maximum number of evolution steps to prevent freezing

	pmp::Scalar PointActivationRadius{ 0.5 }; //>! radius within which the point becomes activated for normal estimation.
	double ActivatedPointPercentageThreshold{ 0.8 }; //>! the fraction of target points to be within PointActivationRadius for the evolution to terminate.
};

class IMB_ShrinkWrapper
{
public:

	// -----------------------------------------------------------------------------
	/// \brief Constructor.
	// -----------------------------------------------------------------------------
	explicit IMB_ShrinkWrapper(const IMB_ShrinkWrapperSettings& settings, const std::vector<pmp::Point>& points)
		: m_Settings(settings), m_Points(points) {}

	// -----------------------------------------------------------------------------
	/// \brief Main functionality.
	// -----------------------------------------------------------------------------
	[[nodiscard]] std::optional<pmp::SurfaceMesh> Perform();

private:

	/**
	 * \brief Preprocess for evolution, i.e.: construct the evolving manifolds, and transform the target data's distance field, and the DF's normalized neg gradient for stabilization.
	 */
	[[nodiscard]] bool Preprocess();

	// -----------------------------------------------------------------------------
	/// \brief Transform all of the geometries so that numerical stability is ensured.
	/// \param[in] stabilizationFactor     a multiplier for stabilizing mean co-volume measure.
	// -----------------------------------------------------------------------------
	void StabilizeGeometries(pmp::Scalar stabilizationFactor = 1.0);

	/**
	 * \brief Perform remeshing on the evolving manifold(s) logged in the m_RemeshTracker.
	 */
	[[nodiscard]] bool Remesh();

	/**
	 * \brief Performs a single step of manifold evolution from the configuration at previous time step.
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

	const IMB_ShrinkWrapperSettings& m_Settings;
	const std::vector<pmp::Point>& m_Points;

	std::shared_ptr<Geometry::ScalarGrid> m_DistanceField{ nullptr };
	std::shared_ptr<Geometry::VectorGrid> m_DFNegNormalizedGradient{ nullptr };

	std::shared_ptr<pmp::SurfaceMesh> m_Surface{ nullptr };

	pmp::mat4 m_TransformToOriginal = pmp::mat4::identity(); //>! a transformation matrix to transform the stabilized geometry back to its original scale and translate to the target center.

	std::function<pmp::ImplicitLaplaceInfo(const pmp::SurfaceMesh& /* mesh */, pmp::Vertex /* v */)> m_ImplicitLaplacianFunction{}; //>! a Laplacian function chosen from parameter laplacianType.

	pmp::Scalar m_ScalingFactor{ 1.0 }; //>! the computed scaling factor for numerical stabilization.

	ScalarGridInterpolationFunction3D m_ScalarInterpolate{}; //>! a parametrizeable function for interpolating values within Geometry::ScalarGrid2D.
	VectorGridInterpolationFunction3D m_VectorInterpolate{};  //>! a parametrizeable function for interpolating vector values within Geometry::VectorGrid2D.

	pmp::AdaptiveRemeshingSettings m_RemeshingSettings;
	bool m_ShouldRemesh{ false }; //>! a flag turned on whenever remeshing is necessary after a time step.

	pmp::BoundingBox2 m_EvolBox{}; //>! the test box for numerical validity of the evolution.
};