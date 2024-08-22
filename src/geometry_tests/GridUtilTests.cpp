#include "gtest/gtest.h"
#include "geometry/Grid.h"
#include "geometry/GridUtil.h"

using namespace Geometry;
using namespace pmp;

constexpr float BOX_EPSILON = 1e-4f;

TEST(GridInterpolationTest, BilinearInterpolateScalarValue2D_case)
{
    // Setup for 2D ScalarGrid
    ScalarGrid2D grid2D(1.0f, pmp::BoundingBox2(
        pmp::vec2(0.0f + BOX_EPSILON, 0.0f + BOX_EPSILON),
        pmp::vec2(2.0f + BOX_EPSILON, 2.0f + BOX_EPSILON)));
    grid2D.Values() = {
        0, 1, 2,
        1, 2, 3,
        2, 3, 4 };

    pmp::vec2 samplePt1(0.5f, 0.5f);
    pmp::vec2 samplePt2(1.5f, 0.5f);
    pmp::vec2 samplePt3(0.5f, 1.5f);
    pmp::vec2 samplePt4(1.5f, 1.5f);

    double result1 = BilinearInterpolateScalarValue(samplePt1, grid2D);
    double result2 = BilinearInterpolateScalarValue(samplePt2, grid2D);
    double result3 = BilinearInterpolateScalarValue(samplePt3, grid2D);
    double result4 = BilinearInterpolateScalarValue(samplePt4, grid2D);

    EXPECT_NEAR(result1, 1.0, 1e-6);
    EXPECT_NEAR(result2, 2.0, 1e-6);
    EXPECT_NEAR(result3, 2.0, 1e-6);
    EXPECT_NEAR(result4, 3.0, 1e-6);
}

TEST(GridInterpolationTest, BilinearInterpolateVectorValue2D_case)
{
    // Setup for 2D VectorGrid
    VectorGrid2D vectorGrid2D(1.0f, pmp::BoundingBox2(
        pmp::vec2(0.0f + BOX_EPSILON, 0.0f + BOX_EPSILON),
        pmp::vec2(2.0f + BOX_EPSILON, 2.0f + BOX_EPSILON)));
    vectorGrid2D.ValuesX() = {
        0, 1, 2,
        1, 2, 3,
        2, 3, 4 };
    vectorGrid2D.ValuesY() = {
        6, 5, 4,
        5, 4, 3,
        4, 3, 2 };

    pmp::vec2 samplePt1(0.5f, 0.5f);
    pmp::vec2 samplePt2(1.5f, 0.5f);
    pmp::vec2 samplePt3(0.5f, 1.5f);
    pmp::vec2 samplePt4(1.5f, 1.5f);

    pmp::dvec2 result1 = BilinearInterpolateVectorValue(samplePt1, vectorGrid2D);
    pmp::dvec2 result2 = BilinearInterpolateVectorValue(samplePt2, vectorGrid2D);
    pmp::dvec2 result3 = BilinearInterpolateVectorValue(samplePt3, vectorGrid2D);
    pmp::dvec2 result4 = BilinearInterpolateVectorValue(samplePt4, vectorGrid2D);

    EXPECT_NEAR(result1[0], 1.0, 1e-6);
    EXPECT_NEAR(result1[1], 5.0, 1e-6);
    EXPECT_NEAR(result2[0], 2.0, 1e-6);
    EXPECT_NEAR(result2[1], 4.0, 1e-6);
    EXPECT_NEAR(result3[0], 2.0, 1e-6);
    EXPECT_NEAR(result3[1], 4.0, 1e-6);
    EXPECT_NEAR(result4[0], 3.0, 1e-6);
    EXPECT_NEAR(result4[1], 3.0, 1e-6);
}

TEST(GridInterpolationTest, TrilinearInterpolateScalarValue_case)
{
    // Setup for 3D ScalarGrid
    ScalarGrid grid3D(1.0f, pmp::BoundingBox(
        pmp::vec3(0.0f + BOX_EPSILON, 0.0f + BOX_EPSILON, 0.0f + BOX_EPSILON),
        pmp::vec3(2.0f + BOX_EPSILON, 2.0f + BOX_EPSILON, 2.0f + BOX_EPSILON)));
    grid3D.Values() = {
        0, 1, 2,
        1, 2, 3,
        2, 3, 4,

        1, 2, 3,
        2, 3, 4,
        3, 4, 5,

        2, 3, 4,
        3, 4, 5,
        4, 5, 6 };

    pmp::vec3 samplePt1(0.5f, 0.5f, 0.5f);
    pmp::vec3 samplePt2(1.5f, 0.5f, 0.5f);
    pmp::vec3 samplePt3(0.5f, 1.5f, 0.5f);
    pmp::vec3 samplePt4(0.5f, 0.5f, 1.5f);

    double result1 = TrilinearInterpolateScalarValue(samplePt1, grid3D);
    double result2 = TrilinearInterpolateScalarValue(samplePt2, grid3D);
    double result3 = TrilinearInterpolateScalarValue(samplePt3, grid3D);
    double result4 = TrilinearInterpolateScalarValue(samplePt4, grid3D);

    EXPECT_NEAR(result1, 1.5, 1e-6);
    EXPECT_NEAR(result2, 2.5, 1e-6);
    EXPECT_NEAR(result3, 2.5, 1e-6);
    EXPECT_NEAR(result4, 2.5, 1e-6);
}

TEST(GridInterpolationTest, TrilinearInterpolateVectorValue_case)
{
    // Setup for 3D VectorGrid
    VectorGrid vectorGrid3D(1.0f, pmp::BoundingBox(
        pmp::vec3(0.0f + BOX_EPSILON, 0.0f + BOX_EPSILON, 0.0f + BOX_EPSILON),
        pmp::vec3(2.0f + BOX_EPSILON, 2.0f + BOX_EPSILON, 2.0f + BOX_EPSILON)));
    vectorGrid3D.ValuesX() = { 0, 1, 2, 1, 2, 3, 2, 3, 4, 1, 2, 3, 2, 3, 4, 3, 4, 5, 2, 3, 4, 3, 4, 5, 4, 5, 6 };
    vectorGrid3D.ValuesY() = { 6, 5, 4, 5, 4, 3, 4, 3, 2, 5, 4, 3, 4, 3, 2, 3, 2, 1, 4, 3, 2, 3, 2, 1, 2, 1, 0 };
    vectorGrid3D.ValuesZ() = { 2, 3, 4, 3, 4, 5, 4, 5, 6, 3, 4, 5, 4, 5, 6, 5, 6, 7, 4, 5, 6, 5, 6, 7, 6, 7, 8 };

    pmp::vec3 samplePt1(0.5f, 0.5f, 0.5f);
    pmp::vec3 samplePt2(1.5f, 0.5f, 0.5f);
    pmp::vec3 samplePt3(0.5f, 1.5f, 0.5f);
    pmp::vec3 samplePt4(0.5f, 0.5f, 1.5f);

    pmp::dvec3 result1 = TrilinearInterpolateVectorValue(samplePt1, vectorGrid3D);
    pmp::dvec3 result2 = TrilinearInterpolateVectorValue(samplePt2, vectorGrid3D);
    pmp::dvec3 result3 = TrilinearInterpolateVectorValue(samplePt3, vectorGrid3D);
    pmp::dvec3 result4 = TrilinearInterpolateVectorValue(samplePt4, vectorGrid3D);

    EXPECT_NEAR(result1[0], 1.5, 1e-6);
    EXPECT_NEAR(result1[1], 4.5, 1e-6);
    EXPECT_NEAR(result1[2], 3.5, 1e-6);
    EXPECT_NEAR(result2[0], 2.5, 1e-6);
    EXPECT_NEAR(result2[1], 3.5, 1e-6);
    EXPECT_NEAR(result2[2], 4.5, 1e-6);
    EXPECT_NEAR(result3[0], 2.5, 1e-6);
    EXPECT_NEAR(result3[1], 3.5, 1e-6);
    EXPECT_NEAR(result3[2], 4.5, 1e-6);
    EXPECT_NEAR(result4[0], 2.5, 1e-6);
    EXPECT_NEAR(result4[1], 3.5, 1e-6);
    EXPECT_NEAR(result4[2], 4.5, 1e-6);
}

TEST(GridInterpolationTest, GetNearestNeighborScalarValue_case)
{
    // Setup for 3D ScalarGrid
    ScalarGrid grid3D(1.0f, pmp::BoundingBox(
        pmp::vec3(0.0f + BOX_EPSILON, 0.0f + BOX_EPSILON, 0.0f + BOX_EPSILON),
        pmp::vec3(2.0f + BOX_EPSILON, 2.0f + BOX_EPSILON, 2.0f + BOX_EPSILON)));
    grid3D.Values() = { 0, 1, 2, 1, 2, 3, 2, 3, 4, 1, 2, 3, 2, 3, 4, 3, 4, 5, 2, 3, 4, 3, 4, 5, 4, 5, 6 };

    pmp::vec3 samplePt(1.4f, 0.6f, 0.6f);
    double result = GetNearestNeighborScalarValue(samplePt, grid3D);
    EXPECT_NEAR(result, 3.0, 1e-6);
}

TEST(GridInterpolationTest, GetNearestNeighborVectorValue_case)
{
    // Setup for 3D VectorGrid
    VectorGrid vectorGrid3D(1.0f, pmp::BoundingBox(
        pmp::vec3(0.0f + BOX_EPSILON, 0.0f + BOX_EPSILON, 0.0f + BOX_EPSILON),
        pmp::vec3(2.0f + BOX_EPSILON, 2.0f + BOX_EPSILON, 2.0f + BOX_EPSILON)));
    vectorGrid3D.ValuesX() = { 0, 1, 2, 1, 2, 3, 2, 3, 4, 1, 2, 3, 2, 3, 4, 3, 4, 5, 2, 3, 4, 3, 4, 5, 4, 5, 6 };
    vectorGrid3D.ValuesY() = { 6, 5, 4, 5, 4, 3, 4, 3, 2, 5, 4, 3, 4, 3, 2, 3, 2, 1, 4, 3, 2, 3, 2, 1, 2, 1, 0 };
    vectorGrid3D.ValuesZ() = { 2, 3, 4, 3, 4, 5, 4, 5, 6, 3, 4, 5, 4, 5, 6, 5, 6, 7, 4, 5, 6, 5, 6, 7, 6, 7, 8 };

    pmp::vec3 samplePt(1.4f, 0.6f, 0.4f);
    pmp::dvec3 result = GetNearestNeighborVectorValue(samplePt, vectorGrid3D);
    EXPECT_NEAR(result[0], 2.0, 1e-6);
    EXPECT_NEAR(result[1], 4.0, 1e-6);
    EXPECT_NEAR(result[2], 4.0, 1e-6);
}

TEST(GridInterpolationTest, GetNearestNeighborScalarValue2D_case)
{
    // Setup for 2D ScalarGrid
    ScalarGrid2D grid2D(1.0f, pmp::BoundingBox2(
        pmp::vec2(0.0f + BOX_EPSILON, 0.0f + BOX_EPSILON),
        pmp::vec2(2.0f + BOX_EPSILON, 2.0f + BOX_EPSILON)));
    grid2D.Values() = {
        0, 1, 2,
        1, 2, 3,
        2, 3, 4 };

    pmp::vec2 samplePt(1.6f, 1.7f);
    double result = GetNearestNeighborScalarValue2D(samplePt, grid2D);
    EXPECT_NEAR(result, 4.0, 1e-6);
}

TEST(GridInterpolationTest, GetNearestNeighborVectorValue2D_case)
{
    // Setup for 2D VectorGrid
    VectorGrid2D vectorGrid2D(1.0f, pmp::BoundingBox2(
        pmp::vec2(0.0f + BOX_EPSILON, 0.0f + BOX_EPSILON),
        pmp::vec2(2.0f + BOX_EPSILON, 2.0f + BOX_EPSILON)));
    vectorGrid2D.ValuesX() = { 0, 1, 2, 1, 2, 3, 2, 3, 4 };
    vectorGrid2D.ValuesY() = { 6, 5, 4, 5, 4, 3, 4, 3, 2 };

    pmp::vec2 samplePt(1.4f, 0.6f);
    pmp::dvec2 result = GetNearestNeighborVectorValue2D(samplePt, vectorGrid2D);
    EXPECT_NEAR(result[0], 2.0, 1e-6);
    EXPECT_NEAR(result[1], 4.0, 1e-6);
}

// Local maxima

TEST(ScalarGridLocalMaximumTests3D, MaximumWithinCellPointsWithUnitRadius)
{
    // Arrange
    ScalarGrid grid(1.0f, BoundingBox(Point(-1, -1, -1), Point(1, 1, 1)));
    grid.Values() = {
        1, 2, 1,
        2, 3, 2,
        1, 2, 1,

        2, 3, 2,
        3, 5, 3,  // Maximum at center (0, 0, 0)
        2, 3, 2,

        1, 2, 1,
        2, 3, 2,
        1, 2, 1
    };
    const unsigned int ix = 1, iy = 1, iz = 1;

    // Act
    const auto localMax = FindLocalMaximumNearScalarGridCell(grid, ix, iy, iz, 1);

    // Assert
    ASSERT_TRUE(localMax.has_value());
    EXPECT_NEAR((*localMax)[0], 0.0f, 1e-5f);
    EXPECT_NEAR((*localMax)[1], 0.0f, 1e-5f);
    EXPECT_NEAR((*localMax)[2], 0.0f, 1e-5f);
}



TEST(ScalarGridLocalMaximumTests3D, MaximumOutsideOfCellPointsWithUnitRadius)
{
    // Arrange
    ScalarGrid grid(1.0f, BoundingBox(pmp::Point(-1, -1, -1), Point(1, 1, 1)));
    grid.Values() = {
        1, 2, 1,
        2, 3, 2,
        1, 2, 1,

        2, 3, 2,
        3, 4, 3,
        2, 3, 2,

        2, 3, 2,
        3, 5, 3,
        2, 3, 2
    };
    const unsigned int ix = 1, iy = 1, iz = 1;

    // Act
    const auto localMax = FindLocalMaximumNearScalarGridCell(grid, ix, iy, iz, 1);

    // Assert
    ASSERT_FALSE(localMax.has_value());
}

// TODO: tests with larger data
//TEST(ScalarGridLocalMaximumTests3D, MaximumWithinCellPointsWithRadius3)
//{
//
//}
//
//
//TEST(ScalarGridLocalMaximumTests3D, MaximumOutsideOfCellPointsWithRadius3)
//{
//
//}