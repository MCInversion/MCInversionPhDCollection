// Copyright 2011-2021 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#pragma once

#include "pmp/ManifoldCurve2D.h"

namespace pmp {

    //! \addtogroup algorithms
    //! @{

    //! Factory class to generate different types of basic 1D manifold shapes.
    class CurveFactory
    {
    public:
        //! Generate a circular arc centered at \p center with radius \p radius, and with \p nSegments segments, start angle \p startAngle, and end angle \p endAngle
        static ManifoldCurve2D circle(
            const Point2& center = Point2{0.0f, 0.0f},
            Scalar radius = 1.0f,
            size_t nSegments = 6,
            Scalar startAngle = 0.0f, 
            Scalar endAngle = 2.0f * M_PI);
    };
} // namespace pmp