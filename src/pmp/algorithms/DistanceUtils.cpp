// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "pmp/algorithms/DistanceUtils.h"

#include <cmath>

#include <limits>

namespace pmp {

Scalar dist_point_line_segment(const Point& p, const Point& v0, const Point& v1,
                               Point& nearest_point)
{
    Point d1(p - v0);
    Point d2(v1 - v0);
    Point min_v(v0);
    Scalar t = dot(d2, d2);

    if (t > std::numeric_limits<Scalar>::min())
    {
        t = dot(d1, d2) / t;
        if (t > 1.0)
            d1 = p - (min_v = v1);
        else if (t > 0.0)
            d1 = p - (min_v = v0 + d2 * t);
    }

    nearest_point = min_v;
    return norm(d1);
}

Scalar dist_point_line_segment(const Point2& p, const Point2& v0, const Point2& v1, Point2& nearest_point)
{
    // Vector from v0 to v1
    pmp::vec2 v(v1[0] - v0[0], v1[1] - v0[1]);

    // Vector from v0 to the point p
    pmp::vec2 w(p[0] - v0[0], p[1] - v0[1]);

    // Calculate the projection of w onto v to find the closest point on the line segment
    Scalar c1 = w[0] * v[0] + w[1] * v[1];
    Scalar c2 = v[0] * v[0] + v[1] * v[1];

    // If the line segment is degenerate (i.e., v0 == v1), return the distance to v0
    if (std::abs(c2) < 1e-6)
    {
        nearest_point = v0;
        return sqrt(w[0] * w[0] + w[1] * w[1]);
    }

    Scalar t = c1 / c2;

    // Clamp t to the range [0, 1] to ensure the closest point is on the segment
    t = std::max(Scalar(0.0), std::min(Scalar(1.0), t));

    // Calculate the closest point on the line segment
    nearest_point[0] = v0[0] + t * v[0];
    nearest_point[1] = v0[1] + t * v[1];

    // Vector from the closest point to the point p
    pmp::vec2 distVec(p[0] - nearest_point[0], p[1] - nearest_point[1]);

    return sqrt(distVec[0] * distVec[0] + distVec[1] * distVec[1]);
}

Scalar dist_point_triangle(const Point& p, const Point& v0, const Point& v1,
                           const Point& v2, Point& nearest_point)
{
    Point v0v1 = v1 - v0;
    Point v0v2 = v2 - v0;
    Point n = cross(v0v1, v0v2); // not normalized !
    Scalar d = sqrnorm(n);

    // Check if the triangle is degenerated -> measure dist to line segments
    if (fabs(d) < std::numeric_limits<Scalar>::min())
    {
        Point query;
        Point nearest;

        auto distance = dist_point_line_segment(p, v0, v1, nearest);

        auto other = dist_point_line_segment(p, v1, v2, query);
        if (other < distance)
        {
            distance = other;
            nearest = query;
        }

        other = dist_point_line_segment(p, v2, v0, query);
        if (other < distance)
        {
            distance = other;
            nearest = query;
        }

        nearest_point = nearest;
        return distance;
    }

    Scalar inv_d = 1.0 / d;
    Point v1v2 = v2;
    v1v2 -= v1;
    Point v0p = p;
    v0p -= v0;
    Point t = cross(v0p, n);
    Scalar a = dot(t, v0v2) * -inv_d;
    Scalar b = dot(t, v0v1) * inv_d;
    Scalar s01, s02, s12;

    // Calculate the distance to an edge or a corner vertex
    if (a < 0)
    {
        s02 = dot(v0v2, v0p) / sqrnorm(v0v2);
        if (s02 < 0.0)
        {
            s01 = dot(v0v1, v0p) / sqrnorm(v0v1);
            if (s01 <= 0.0)
            {
                v0p = v0;
            }
            else if (s01 >= 1.0)
            {
                v0p = v1;
            }
            else
            {
                (v0p = v0) += (v0v1 *= s01);
            }
        }
        else if (s02 > 1.0)
        {
            s12 = dot(v1v2, (p - v1)) / sqrnorm(v1v2);
            if (s12 >= 1.0)
            {
                v0p = v2;
            }
            else if (s12 <= 0.0)
            {
                v0p = v1;
            }
            else
            {
                (v0p = v1) += (v1v2 *= s12);
            }
        }
        else
        {
            (v0p = v0) += (v0v2 *= s02);
        }
    }

    // Calculate the distance to an edge or a corner vertex
    else if (b < 0.0)
    {
        s01 = dot(v0v1, v0p) / sqrnorm(v0v1);
        if (s01 < 0.0)
        {
            s02 = dot(v0v2, v0p) / sqrnorm(v0v2);
            if (s02 <= 0.0)
            {
                v0p = v0;
            }
            else if (s02 >= 1.0)
            {
                v0p = v2;
            }
            else
            {
                (v0p = v0) += (v0v2 *= s02);
            }
        }
        else if (s01 > 1.0)
        {
            s12 = dot(v1v2, (p - v1)) / sqrnorm(v1v2);
            if (s12 >= 1.0)
            {
                v0p = v2;
            }
            else if (s12 <= 0.0)
            {
                v0p = v1;
            }
            else
            {
                (v0p = v1) += (v1v2 *= s12);
            }
        }
        else
        {
            (v0p = v0) += (v0v1 *= s01);
        }
    }

    // Calculate the distance to an edge or a corner vertex
    else if (a + b > 1.0)
    {
        s12 = dot(v1v2, (p - v1)) / sqrnorm(v1v2);
        if (s12 >= 1.0)
        {
            s02 = dot(v0v2, v0p) / sqrnorm(v0v2);
            if (s02 <= 0.0)
            {
                v0p = v0;
            }
            else if (s02 >= 1.0)
            {
                v0p = v2;
            }
            else
            {
                (v0p = v0) += (v0v2 *= s02);
            }
        }
        else if (s12 <= 0.0)
        {
            s01 = dot(v0v1, v0p) / sqrnorm(v0v1);
            if (s01 <= 0.0)
            {
                v0p = v0;
            }
            else if (s01 >= 1.0)
            {
                v0p = v1;
            }
            else
            {
                (v0p = v0) += (v0v1 *= s01);
            }
        }
        else
        {
            (v0p = v1) += (v1v2 *= s12);
        }
    }

    // Calculate the distance to an interior point of the triangle
    else
    {
        n *= (dot(n, v0p) * inv_d);
        (v0p = p) -= n;
    }

    nearest_point = v0p;
    v0p -= p;
    return norm(v0p);
}

} // namespace pmp
