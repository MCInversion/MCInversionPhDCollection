#pragma once

#include "pmp/Types.h"
#include "pmp/ManifoldCurve2D.h"

namespace Geometry
{
	/// \brief A builder object for an inner circle curve to an arbitrary 2D geometry.
	class InnerCircleBuilder
	{
	public:
		/**
		 * \brief Solves an optimization problem to fit a sphere in the input point cloud.
		 * \param points        input point cloud.
		 * \return pair { radius, center } of the resulting circle.
		 */
		static [[nodiscard]] std::pair<pmp::Scalar, pmp::Point2> CalculateInnerCircleRadiusAndCenter(const std::vector<pmp::Point2>& points);

		/**
		 * \brief Tessellates the resulting curve from the computed radius and center.
		 * \return The tessellation of the resulting circle (must be called after CalculateInnerCircleRadiusAndCenter).
		 */
		static [[nodiscard]] pmp::ManifoldCurve2D ConstructCircle();

	private:

		pmp::Scalar m_Radius;
		pmp::Point2 m_Center;
	};
	
} // namespace Geometry