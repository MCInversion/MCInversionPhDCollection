// Copyright 2013-2021 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#pragma once

#include <cassert>

#include "pmp/Types.h"

namespace pmp {

//! Simple class for representing a bounding box.
//! \ingroup core
class BoundingBox
{
public:
    //! Construct infinite/invalid bounding box.
    BoundingBox()
        : min_(std::numeric_limits<Scalar>::max()),
          max_(-std::numeric_limits<Scalar>::max())
    {
    }

    //! Construct from min and max points.
    BoundingBox(const Point& min, const Point& max) : min_(min), max_(max) {}

    //! Add point to the bounding box.
    BoundingBox& operator+=(const Point& p)
    {
        for (int i = 0; i < 3; ++i)
        {
            if (p[i] < min_[i])
                min_[i] = p[i];
            if (p[i] > max_[i])
                max_[i] = p[i];
        }
        return *this;
    }

    //! Add two bounding boxes.
    BoundingBox& operator+=(const BoundingBox& bb)
    {
        for (int i = 0; i < 3; ++i)
        {
            if (bb.min_[i] < min_[i])
                min_[i] = bb.min_[i];
            if (bb.max_[i] > max_[i])
                max_[i] = bb.max_[i];
        }
        return *this;
    }

    //! Get min point.
    Point& min() { return min_; }

    //! Get max point.
    Point& max() { return max_; }

    //! Get min point.
    [[nodiscard]] const Point& min() const { return min_; }

    //! Get max point.
    [[nodiscard]] const Point& max() const { return max_; }

    //! Get center point.
    Point center() const { return 0.5f * (min_ + max_); }

    //! Indicate if the bounding box is empty.
    bool is_empty() const
    {
        return (max_[0] < min_[0] || max_[1] < min_[1] || max_[2] < min_[2]);
    }

    //! Get the size of the bounding box.
    Scalar size() const
    {
        return is_empty() ? Scalar(0.0) : distance(max_, min_);
    }

    //! Get the intersection value of this box with another.
    [[nodiscard]] bool Intersects(const BoundingBox& other) const
    {        
        if (other.max_[0] < min_[0] || other.min_[0] > max_[0])
            return false;

        if (other.max_[1] < min_[1] || other.min_[1] > max_[1])
            return false;

        return (other.max_[2] >= min_[2] && other.min_[2] <= max_[2]);
    }

    //! Expand the size of the bounding box.
    void expand(const float& x, const float& y, const float& z)
    {
        assert(x >= 0.0f && y >= 0.0f && z >= 0.0f);
        min_ -= Point(x, y, z);
        max_ += Point(x, y, z);
    }

private:
    Point min_, max_;
};

} // namespace pmp
