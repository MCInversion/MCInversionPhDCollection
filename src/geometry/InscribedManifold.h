#pragma once

#include "pmp/Types.h"
#include "pmp/ManifoldCurve2D.h"

namespace Geometry
{
	/// \brief Enumerator for choosing the calculation method for inner circle/sphere.
	enum class [[nodiscard]] InscribedSphereMethod
	{
		MinDistToBBoxCenter = 0 //>! this is a naive method to calculate the inner sphere.
	};


	/// \brief A builder object for an inscribed circle curve to an arbitrary 2D geometry.
	class InscribedCircleBuilder
	{
	public:
		/**
		 * \brief Solves an optimization problem to fit a sphere in the input point cloud.
		 * \param points        input point cloud.
		 * \return pair { radius, center } of the resulting circle.
		 */
		static [[nodiscard]] std::pair<pmp::Scalar, pmp::Point2> CalculateInscribedCircleRadiusAndCenter(const std::vector<pmp::Point2>& points);

		/**
		 * \brief Tessellates the resulting curve from the computed radius and center.
		 * \return The tessellation of the resulting circle (must be called after CalculateInnerCircleRadiusAndCenter).
		 */
		static [[nodiscard]] pmp::ManifoldCurve2D ConstructCircle();

	private:

		pmp::Scalar m_Radius; //>! the radius of the inner circle
		pmp::Point2 m_Center; //>! the center point of the inner circle
	};
	
} // namespace Geometry