// Copyright 2013-2021 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#pragma once

#include <cassert>

#include "pmp/Types.h"

namespace pmp {

//! Simple class for representing a 3D bounding box.
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

    //! Construct from a vector of points.
    explicit BoundingBox(const std::vector<Point>& pts)
        : min_(std::numeric_limits<Scalar>::max()),
	      max_(-std::numeric_limits<Scalar>::max())
    {
        for (int pId = 0; const auto& p : pts)
		{
            for (int i = 0; i < 3; ++i)
            {
                if (p[i] < min_[i])
                    min_[i] = p[i];
                if (p[i] > max_[i])
                    max_[i] = p[i];
            }

            ++pId;
        }
    }

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

    BoundingBox& operator+=(const std::vector<Point>& pts)
    {
        for (const auto& pt : pts)
            operator+=(pt);
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
    [[nodiscard]] Point center() const { return 0.5f * (min_ + max_); }

    //! Indicate if the bounding box is empty.
    [[nodiscard]] bool is_empty() const
    {
        return (max_[0] < min_[0] || max_[1] < min_[1] || max_[2] < min_[2]);
    }

    //! Get the size of the bounding box.
    [[nodiscard]] Scalar size() const
    {
        return is_empty() ? 0.0f : distance(max_, min_);
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

    //! Get the intersection box of this box with another.
    [[nodiscard]] BoundingBox Intersect(const BoundingBox& other) const
    {
        if (!Intersects(other))
            return {}; // empty intersection

        BoundingBox result;
        for (int i = 0; i < 3; i++)
        {
	        result.min_[i] = std::fmaxf(min_[i], other.min_[i]);
	        result.max_[i] = std::fminf(max_[i], other.max_[i]);
        }
        return result;
    }

    //! Verifies whether a point is contained within this bounding box.
    [[nodiscard]] bool Contains(const Point& pt) const
    {
        if (pt[0] < min_[0])
            return false;

        if (pt[0] > max_[0])
            return false;

        if (pt[1] < min_[1])
            return false;

        if (pt[1] > max_[1])
            return false;

        if (pt[2] < min_[2])
            return false;

        return (pt[2] <= max_[2]);
    }

    //! Expand the size of the bounding box.
    void expand(const float& x, const float& y, const float& z)
    {
        assert(x >= 0.0f && y >= 0.0f && z >= 0.0f);
        min_ -= Point(x, y, z);
        max_ += Point(x, y, z);
    }

    BoundingBox& operator*= (const mat4& mat)
    {
        min_ = affine_transform(mat, min_);
        max_ = affine_transform(mat, max_);
        return *this;
    }

    [[nodiscard]] Scalar volume() const
    {
        const auto dims = max_ - min_;
        return dims[0] * dims[1] * dims[2];
    }

    [[nodiscard]] Scalar width() const
    {
        return max_[0] - min_[0];
    }

    [[nodiscard]] Scalar depth() const
    {
        return max_[1] - min_[1];
    }

    [[nodiscard]] Scalar height() const
    {
        return max_[2] - min_[2];
    }

private:
    Point min_, max_;
};

//! Simple class for representing a 2D bounding box.
//! \ingroup core
class BoundingBox2
{
public:
    //! Construct infinite/invalid bounding box.
    BoundingBox2()
        : min_(std::numeric_limits<Scalar>::max()),
        max_(-std::numeric_limits<Scalar>::max())
    {
    }

    //! Construct from min and max points.
    BoundingBox2(const Point2& min, const Point2& max) : min_(min), max_(max) {}

    //! Construct from a vector of points.
    explicit BoundingBox2(const std::vector<Point>& pts)
        : min_(std::numeric_limits<Scalar>::max()),
        max_(-std::numeric_limits<Scalar>::max())
    {
        for (int pId = 0; const auto & p : pts)
        {
            for (int i = 0; i < 2; ++i)
            {
                if (p[i] < min_[i])
                    min_[i] = p[i];
                if (p[i] > max_[i])
                    max_[i] = p[i];
            }

            ++pId;
        }
    }

    //! Add point to the bounding box.
    BoundingBox2& operator+=(const Point2& p)
    {
        for (int i = 0; i < 2; ++i)
        {
            if (p[i] < min_[i])
                min_[i] = p[i];
            if (p[i] > max_[i])
                max_[i] = p[i];
        }
        return *this;
    }

    //! Add two bounding boxes.
    BoundingBox2& operator+=(const BoundingBox2& bb)
    {
        for (int i = 0; i < 2; ++i)
        {
            if (bb.min_[i] < min_[i])
                min_[i] = bb.min_[i];
            if (bb.max_[i] > max_[i])
                max_[i] = bb.max_[i];
        }
        return *this;
    }

    BoundingBox2& operator+=(const std::vector<Point2>& pts)
    {
        for (const auto& pt : pts)
            operator+=(pt);
        return *this;
    }

    //! Get min point.
    Point2& min() { return min_; }

    //! Get max point.
    Point2& max() { return max_; }

    //! Get min point.
    [[nodiscard]] const Point2& min() const { return min_; }

    //! Get max point.
    [[nodiscard]] const Point2& max() const { return max_; }

    //! Get center point.
    [[nodiscard]] Point2 center() const { return 0.5f * (min_ + max_); }

    //! Indicate if the bounding box is empty.
    [[nodiscard]] bool is_empty() const
    {
        return (max_[0] < min_[0] || max_[1] < min_[1]);
    }

    //! Get the size of the bounding box.
    [[nodiscard]] Scalar size() const
    {
        return is_empty() ? 0.0f : distance(max_, min_);
    }

    //! Get the intersection value of this box with another.
    [[nodiscard]] bool Intersects(const BoundingBox2& other) const
    {
        if (other.max_[0] < min_[0] || other.min_[0] > max_[0])
            return false;

        return (other.max_[1] >= min_[1] && other.min_[1] <= max_[1]);
    }

    //! Get the intersection box of this box with another.
    [[nodiscard]] BoundingBox2 Intersect(const BoundingBox2& other) const
    {
        if (!Intersects(other))
            return {}; // empty intersection

        BoundingBox2 result;
        for (int i = 0; i < 2; i++)
        {
            result.min_[i] = std::fmaxf(min_[i], other.min_[i]);
            result.max_[i] = std::fminf(max_[i], other.max_[i]);
        }
        return result;
    }

    //! Verifies whether a point is contained within this bounding box.
    [[nodiscard]] bool Contains(const Point2& pt) const
    {
        if (pt[0] < min_[0])
            return false;

        if (pt[0] > max_[0])
            return false;

        if (pt[1] < min_[1])
            return false;

        return (pt[1] <= max_[1]);
    }

    //! Expand the size of the bounding box.
    void expand(const float& x, const float& y)
    {
        assert(x >= 0.0f && y >= 0.0f);
        min_ -= Point2(x, y);
        max_ += Point2(x, y);
    }

    BoundingBox2& operator*= (const mat3& mat)
    {
        min_ = affine_transform(mat, min_);
        max_ = affine_transform(mat, max_);
        return *this;
    }

    [[nodiscard]] Scalar volume() const
    {
        const auto dims = max_ - min_;
        return dims[0] * dims[1];
    }

    [[nodiscard]] Scalar width() const
    {
        return max_[0] - min_[0];
    }

    [[nodiscard]] Scalar height() const
    {
        return max_[1] - min_[1];
    }

private:
    Point2 min_, max_;
};


inline [[nodiscard]] mat4 CalculateTransformMatrixBetweenBoxes(const BoundingBox& fromBox, const BoundingBox& toBox, const bool forceUniform = false)
{
    const auto translationVector = toBox.min() - fromBox.min();
    const mat4 translationMatrix = translation_matrix(translationVector);

    // Calculate scale factors for each dimension
    const float scaleX = toBox.width() / fromBox.width();
    const float scaleY = toBox.height() / fromBox.height();
    const float scaleZ = toBox.depth() / fromBox.depth();
    mat4 scalingMatrix;
    if (forceUniform)
    {
		const float uniformScale = (scaleX + scaleY + scaleZ) / 3.0f;
        scalingMatrix = scaling_matrix(uniformScale);	    
    }
    else
    {
        const vec3 scaleVec{ scaleX, scaleY, scaleZ };
        scalingMatrix = scaling_matrix(scaleVec);
    }

    // Combine translation and scaling
    return translationMatrix * scalingMatrix;
}

inline [[nodiscard]] mat3 CalculateTransformationMatrixBetweenBoxes2D(const BoundingBox2& fromBox, const BoundingBox2& toBox, const bool forceUniform = false)
{
    const auto translationVector = toBox.min() - fromBox.min();
    const mat3 translationMatrix = translation_matrix(translationVector);

    // Calculate scale factors for each dimension
    const float scaleX = toBox.width() / fromBox.width();
    const float scaleY = toBox.height() / fromBox.height();
    mat3 scalingMatrix;
    if (forceUniform)
    {
        const float uniformScale = (scaleX + scaleY) / 3.0f;
        scalingMatrix = scaling_matrix_2d(uniformScale);
    }
    else
    {
        const vec2 scaleVec{ scaleX, scaleY};
        scalingMatrix = scaling_matrix_2d(scaleVec);
    }

    // Combine translation and scaling
    return translationMatrix * scalingMatrix;
}

} // namespace pmp
