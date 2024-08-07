#pragma once

#include "pmp/Types.h"
#include "pmp/ManifoldCurve2D.h"
#include "pmp/algorithms/CurveFactory.h"

#include "geometry/Grid.h"

/// \brief A wrapper for input data for the InscribedCircleBuilder
struct InscribedCircleInputData
{
	std::vector<pmp::Point2> Points{}; //>! the evaluated point cloud.
	std::shared_ptr<Geometry::ScalarGrid2D> DistanceField{ nullptr }; //>! a pointer to a pre-computed distance field to Points. 
};

/// \brief Parameters of a 2D circle.
struct Circle2D
{
	pmp::Point2 Center{};
	pmp::Scalar Radius{ 1.0f };
};

/// \brief A base utility to calculate the centers and radii of circles inscribed to a point cloud.
class InscribedCircleCalculator
{
public:
	virtual ~InscribedCircleCalculator() = default;

	/**
	 * \brief Estimates inscribed circles to a point cloud.
	 * \param data        input point cloud data.
	 * \param method      a calculation method to be used.
	 * \return a vector of circles.
	 */
	virtual [[nodiscard]] std::vector<Circle2D> Calculate(const InscribedCircleInputData& data) = 0;
};

/// \brief Calculates the centers and radii of circles inscribed to a point cloud using the naive approach:
/// Use the center of the bounding box of the input point cloud, and the distance to the closest point as radius.
class NaiveInscribedCircleCalculator : public InscribedCircleCalculator
{
public:
	/// \copydoc InscribedCircleCalculator::Calculate
	[[nodiscard]] std::vector<Circle2D> Calculate(const InscribedCircleInputData& data) override;
};

/// \brief Calculates the centers and radii of circles inscribed to a point cloud using the distance-field-based approach:
/// Find the approximate locations of local maxima of the point cloud distance field which will serve as centers, and the distance to the closest point as radii.
class DistanceFieldInscribedCircleCalculator : public InscribedCircleCalculator
{
public:
	/// \copydoc InscribedCircleCalculator::Calculate
	[[nodiscard]] std::vector<Circle2D> Calculate(const InscribedCircleInputData& data) override;
};

/**
 * \brief Tessellates the resulting curve from the computed radius and center.
 * \param circle        A parametric circle to be reconstructed.
 * \param nSegments     the number of segments to be used in the reconstruction.
 * \return The tessellation of the resulting circle (must be called after CalculateInnerCircleRadiusAndCenter).
 */
inline [[nodiscard]] pmp::ManifoldCurve2D ConstructCircle(const Circle2D& circle, size_t nSegments = 32)
{
	return pmp::CurveFactory::circle(circle.Center, circle.Radius, nSegments);
}