#pragma once

#include "pmp/SurfaceMesh.h"

#include "geometry/Grid.h"

#include "EvolverUtilsCommon.h"


struct IMB_ShrinkWrapperSettings
{
	unsigned int LevelOfDetail{ 3 }; //>! Level of detail. Reflected in the number of vertices used during discretization. Standardized for icosphere.
	double TimeStep{ 0.01 };     //>! time step size.

	bool UseLinearGridInterpolation{ true }; //>! whether to use d-linear interpolation for scalar and vector fields where d is the grid dimension.
	AmbientFieldSettings FieldSettings{}; //>! the settings for the construction of the ambient fields.
	ManifoldAdaptiveRemeshingParams RemeshingSettings{}; //>! settings for pmp::Remeshing.

	unsigned int MaxSteps{ 200 }; //>! the maximum number of evolution steps to prevent freezing

	pmp::Scalar PointActivationRadius{ 0.5 }; //>! radius within which the point becomes activated for normal estimation.
	double ActivatedPointPercentageThreshold{ 0.8 }; //>! the fraction of target points to be within PointActivationRadius for the evolution to terminate.
};

class IMB_ShrinkWrapper
{
public:

	explicit IMB_ShrinkWrapper(const IMB_ShrinkWrapperSettings& settings, const std::vector<pmp::Point>& points)
		: m_Settings(settings), m_Points(points) {}

	[[nodiscard]] std::optional<pmp::SurfaceMesh> Perform();

private:

	[[nodiscard]] bool Preprocess();

	const IMB_ShrinkWrapperSettings& m_Settings;
	const std::vector<pmp::Point>& m_Points;

	std::shared_ptr<Geometry::ScalarGrid> m_DistanceField{ nullptr };
	std::shared_ptr<Geometry::VectorGrid> m_DFNegNormalizedGradient{ nullptr };

	std::shared_ptr<pmp::SurfaceMesh> m_Surface{ nullptr };
};