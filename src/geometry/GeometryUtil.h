#pragma once

#include "pmp/MatVec.h"

namespace Geometry
{
	/**
	 * \brief Computes the squared distance from point to a triangle.
	 * \param vertices     list of (three) vertices of a triangle.
	 * \param point        point from which the distance is to be computed.
	 * \return squared distance from point to triangle.
	 */
	[[nodiscard]] double GetDistanceToTriangleSq(const std::vector<pmp::vec3>& vertices, const pmp::vec3& point);

	/**
	 * \brief An intersection test between a triangle and a box.
	 * \param vertices     list of (three) vertices of a triangle.
	 * \param boxCenter    center of the box to be queried.
	 * \param boxHalfSize  half-size of the box to be queried.
	 * \return true if the box intersects the triangle.
	 */
	[[nodiscard]] bool TriangleIntersectsBox(const std::vector<pmp::vec3>& vertices, const pmp::vec3& boxCenter, const pmp::vec3& boxHalfSize);
	
} // namespace Geometry