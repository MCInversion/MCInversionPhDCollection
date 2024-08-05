// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#pragma once

#include "pmp/Types.h"

namespace pmp {

//! \addtogroup algorithms
//! @{

//! Compute the distance of a point p to a 3D line segment given by points (v0,v1).
Scalar dist_point_line_segment(const Point& p, const Point& v0, const Point& v1,
                               Point& nearest_point);

//! Compute the distance of a point p to a 2D line segment given by 2D points (v0,v1).
Scalar dist_point_line_segment(const Point2& p, const Point2& v0, const Point2& v1,
    Point2& nearest_point);

//! Compute the distance of a point p to the triangle given by points (v0, v1, v2).
Scalar dist_point_triangle(const Point& p, const Point& v0, const Point& v1,
                           const Point& v2, Point& nearest_point);

//! @}

} // namespace pmp
