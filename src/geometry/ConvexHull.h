#pragma once

#include "pmp/SurfaceMesh.h"

namespace Geometry
{
	/**
	 * \brief Internal utility for filling a convex hull mesh with points.
	 */
	[[nodiscard]] pmp::SurfaceMesh ComputeConvexHullFromPoints(const std::vector<pmp::Point>& points, const float& distTolerance);
	
} // namespace Geometry