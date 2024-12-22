#include "gtest/gtest.h"

#include "geometry/GeometryUtil.h"

#include <pmp/Types.h>
#include <optional>

using namespace Geometry;
using namespace pmp;

// Testing macros
#define EXPECT_VEC3_NEAR(vec1, vec2, tol) \
    do { \
        EXPECT_NEAR((vec1)[0], (vec2)[0], tol); \
        EXPECT_NEAR((vec1)[1], (vec2)[1], tol); \
        EXPECT_NEAR((vec1)[2], (vec2)[2], tol); \
    } while(0)

#define EXPECT_VEC2_NEAR(vec1, vec2, tol) \
    do { \
        EXPECT_NEAR((vec1)[0], (vec2)[0], tol); \
        EXPECT_NEAR((vec1)[1], (vec2)[1], tol); \
    } while(0)

// ========= Distance utils ===========
TEST(DistanceUtilTests, GetDistanceToTriangleSq)
{
    // Arrange
    const std::vector triangle = { vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0) };
    const vec3 point(0.5, 0.5, 1);

    // Act
    const double distance = GetDistanceToTriangleSq(triangle, point);

    // Assert
    EXPECT_NEAR(distance, 1.0, 1e-6);
}

TEST(DistanceUtilTests, GetDistanceToTriangleSq_NonIntersectingPoint)
{
    // Arrange
    const std::vector triangle = { vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0) };
    const vec3 point(2, 2, 2);

    // Act
    const double distance = GetDistanceToTriangleSq(triangle, point);

    // Assert
    EXPECT_NEAR(distance, 8.5, 1e-6);
}

TEST(DistanceUtilTests, GetDistanceToLine2DSq)
{
    // Arrange
    const std::vector line = { Point2(0, 0), Point2(1, 0) };
    const vec2 point(0.5, 1);

    // Act
    const double distance = GetDistanceToLine2DSq(line, point);

    // Assert
    EXPECT_NEAR(distance, 1.0, 1e-6);
}

TEST(DistanceUtilTests, GetDistanceToLine2DSq_NonIntersectingPoint)
{
    // Arrange
    const std::vector line = { Point2(0, 0), Point2(1, 0) };
    const vec2 point(0.5, 2);

    // Act
    const double distance = GetDistanceToLine2DSq(line, point);

    // Assert
    EXPECT_NEAR(distance, 4.0, 1e-6); // distance should be 2^2
}


// ======== Intersection utils =============

TEST(IntersectionUtilTests, TriangleIntersectsBox)
{
    // Arrange
    const std::vector triangle = { vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0) };
    const vec3 boxCenter(0.5, 0.5, 0);
    const vec3 boxHalfSize(0.5, 0.5, 0.5);

    // Act/Assert
    EXPECT_TRUE(TriangleIntersectsBox(triangle, boxCenter, boxHalfSize));
}

TEST(IntersectionUtilTests, TriangleIntersectsBox_NoIntersection)
{
    // Arrange
    const std::vector triangle = { vec3(2, 2, 2), vec3(3, 2, 2), vec3(2, 3, 2) };
    const vec3 boxCenter(0.5, 0.5, 0);
    const vec3 boxHalfSize(0.5, 0.5, 0.5);

    // Act/Assert
    EXPECT_FALSE(TriangleIntersectsBox(triangle, boxCenter, boxHalfSize));
}

TEST(IntersectionUtilTests, TriangleIntersectsBox_AlmostTouching)
{
    // Arrange
    const std::vector triangle = { vec3(1.01, 0, 0), vec3(2.01, 0, 0), vec3(1.01, 1, 0) };
    const vec3 boxCenter(0.5, 0.5, 0);
    const vec3 boxHalfSize(0.5, 0.5, 0.5);

    // Act/Assert
    EXPECT_FALSE(TriangleIntersectsBox(triangle, boxCenter, boxHalfSize));
}

TEST(IntersectionUtilTests, Line2DIntersectsBox)
{
    // Arrange
    const std::vector line = { Point2(0, 0), Point2(1, 1) };
    const Point2 boxCenter(0.5, 0.5);
    const vec2 boxHalfSize(0.5, 0.5);

    // Act/Assert
	EXPECT_TRUE(Line2DIntersectsBox(line, boxCenter, boxHalfSize));
}

TEST(IntersectionUtilTests, Line2DIntersectsBox_NoIntersection)
{
    // Arrange
    const std::vector line = { Point2(2, 2), Point2(3, 3) };
    const Point2 boxCenter(0.5, 0.5);
    const vec2 boxHalfSize(0.5, 0.5);

    // Act/Assert
    EXPECT_FALSE(Line2DIntersectsBox(line, boxCenter, boxHalfSize));
}

TEST(IntersectionUtilTests, Line2DIntersectsBox_AlmostTouching)
{
    // Arrange
    const std::vector line = { Point2(1.01, 1.01), Point2(2, 2) };
    const Point2 boxCenter(0.5, 0.5);
    const vec2 boxHalfSize(0.5, 0.5);

    // Act/Assert
    EXPECT_FALSE(Line2DIntersectsBox(line, boxCenter, boxHalfSize));
}


TEST(IntersectionUtilTests, TriangleIntersectsTriangle)
{
    // Arrange
    const std::vector triangle1 = { vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0) };
    const std::vector triangle2 = { vec3(0.5, 0.5, 0), vec3(1.5, 0.5, 0), vec3(0.5, 1.5, 0) };

    // Act/Assert
	EXPECT_TRUE(TriangleIntersectsTriangle(triangle1, triangle2));
}

TEST(IntersectionUtilTests, TriangleIntersectsTriangle_NoIntersection)
{
    // Arrange
    const std::vector triangle1 = { vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0) };
    const std::vector triangle2 = { vec3(2, 2, 0), vec3(3, 2, 0), vec3(2, 3, 0) };

    // Act/Assert
    EXPECT_FALSE(TriangleIntersectsTriangle(triangle1, triangle2));
}

TEST(IntersectionUtilTests, TriangleIntersectsTriangle_AlmostTouching)
{
    // Arrange
    const std::vector triangle1 = { vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0) };
    const std::vector triangle2 = { vec3(1.01, 0.01, 0), vec3(2, 0.01, 0), vec3(1.01, 1, 0) };

    // Act/Assert
    EXPECT_FALSE(TriangleIntersectsTriangle(triangle1, triangle2));
}

// TODO: Investigate why this one returns large coordinates despite not returning nullopt
//TEST(IntersectionUtilTests, ComputeTriangleTriangleIntersectionLine_CoplanarTouching)
//{
//    // Arrange
//    const std::vector triangle1 = { vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0) };
//    const std::vector triangle2 = { vec3(0.4, 0.4, -0.5), vec3(1.5, 0.5, 0.5), vec3(0.5, 1.5, 0.5) };
//
//    // Act
//    const auto intersection = ComputeTriangleTriangleIntersectionLine(triangle1, triangle2);
//
//    // Assert
//    EXPECT_TRUE(intersection.has_value());
//    if (intersection) {
//        EXPECT_EQ(intersection->first, vec3(0.5, 0.5, 0));
//        EXPECT_EQ(intersection->second, vec3(1, 0, 0));
//    }
//}

TEST(IntersectionUtilTests, ComputeTriangleTriangleIntersectionLine)
{
    // Arrange
    const std::vector triangle1 = { vec3(0, 0, 0), vec3(2, 0, 0), vec3(0, 2, 0) };
    const std::vector triangle2 = { vec3(0.4, 0.4, -0.2), vec3(1.5, 0.5, 0.2), vec3(0.5, 1.5, 0.2) };

    // Act
	const auto intersection = ComputeTriangleTriangleIntersectionLine(triangle1, triangle2);

    // Assert
	EXPECT_TRUE(intersection.has_value());
    if (intersection) {
        EXPECT_VEC3_NEAR(intersection->first, vec3(0.95, 0.45, 0), 1e-6);
        EXPECT_VEC3_NEAR(intersection->second, vec3(0.45, 0.95, 0), 1e-6);
    }
}

TEST(IntersectionUtilTests, ComputeTriangleTriangleIntersectionLine_NoIntersection)
{
    // Arrange
    const std::vector triangle1 = { vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0) };
    const std::vector triangle2 = { vec3(2, 2, 0), vec3(3, 2, 0), vec3(2, 3, 0) };

    // Act
    const auto intersection = ComputeTriangleTriangleIntersectionLine(triangle1, triangle2);

    // Assert
    EXPECT_FALSE(intersection.has_value());
}

TEST(IntersectionUtilTests, ComputeTriangleTriangleIntersectionLine_AlmostTouching)
{
    // Arrange
    const std::vector triangle1 = { vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0) };
    const std::vector triangle2 = { vec3(1.01, 0.01, 0), vec3(2, 0.01, 0), vec3(1.01, 1, 0) };

    // Act
    const auto intersection = ComputeTriangleTriangleIntersectionLine(triangle1, triangle2);

    // Assert
    EXPECT_FALSE(intersection.has_value());
}
