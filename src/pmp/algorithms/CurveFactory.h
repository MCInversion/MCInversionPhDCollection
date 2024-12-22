// Copyright 2011-2021 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#pragma once

#include "pmp/ManifoldCurve2D.h"

namespace pmp {

    /// \brief Vertices of an equilateral triangle with side 1 with barycenter at (0, 0).
    const std::vector BASE_EQUILATERAL_TRIANGLE_VERTS{ 
        Point2{-0.5, -sqrtf(3.0) / 6.0}, 
        Point2{0.5, -sqrtf(3.0) / 6.0}, 
        Point2{0.0, sqrtf(3.0) / 3.0} };

    //! \addtogroup algorithms
    //! @{

    //! Factory class to generate different types of basic 1D manifold shapes.
    class CurveFactory
    {
    public:
        //! Generate a circular arc centered at \p center with radius \p radius, and with \p nSegments segments, start angle \p startAngle, and end angle \p endAngle
        static ManifoldCurve2D circle(
            const Point2& center = Point2{0.0, 0.0},
            Scalar radius = 1.0,
            size_t nSegments = 6,
            Scalar startAngle = 0.0, 
            Scalar endAngle = 2.0 * M_PI);
        //! Generate a circular arc centered at \p center with radius \p radius, and with \p nSegments segments, start angle \p startAngle, and end angle \p endAngle deformed with a sine wave with \p amplitude and \p freq.
        static ManifoldCurve2D sine_deformed_circle(
            const Point2& center = Point2{ 0.0, 0.0 },
            Scalar radius = 1.0,
            size_t nSegments = 6,
            Scalar amplitude = 1.0,
            Scalar freq = 4,
            Scalar startAngle = 0.0,
            Scalar endAngle = 2.0 * M_PI
        );
        //! Generate a rectangular polyline centered at \p with sides \p sideX and \p sideY, and with \p nSegments, \p nChamfer, start angle \p startAngle, and end angle \p endAngle
        static ManifoldCurve2D rectangle(
            const Point2& center = Point2{ 0.0, 0.0 },
            Scalar sideX = 1.0,
            Scalar sideY = 1.0,
            size_t nSegments = 6,
            bool chamferCorners = true);
        //! Generate a resampled polyline from a given base polygon.
        static ManifoldCurve2D sampled_polygon(
            const std::vector<Point2>& polyVertices = BASE_EQUILATERAL_TRIANGLE_VERTS,
            size_t nSegments = 10,
            bool chamferCorners = true,
            bool closeLoop = true,
            bool checkOrientation = true);
        //! Generate a hyperellipse centered at \p center with radii \p radiusX and \p radiusY of degree 
        //! \p degree, and with \p nSegments segments, start angle \p startAngle, and end angle \p endAngle
        static ManifoldCurve2D hyper_ellipse(
			const Point2& center = Point2{ 0.0, 0.0 },
			Scalar radiusX = 1.0,
			Scalar radiusY = 1.0,
            size_t degree = 4,
			size_t nSegments = 6,
			Scalar startAngle = 0.0,
			Scalar endAngle = 2.0 * M_PI
		);
    };
} // namespace pmp