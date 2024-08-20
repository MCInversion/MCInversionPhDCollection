#include "gtest/gtest.h"

#include "geometry/Grid.h"
#include "geometry/GridUtil.h"

using namespace Geometry;
using namespace pmp;

namespace
{
    class GridInterpolationTestFixture : public ::testing::Test 
    {
    protected:
        std::shared_ptr<ScalarGrid> grid3D{ nullptr };
        std::shared_ptr<ScalarGrid2D> grid2D{ nullptr };
        std::shared_ptr<VectorGrid> vectorGrid3D{ nullptr };
        std::shared_ptr<VectorGrid2D> vectorGrid2D{ nullptr };

        void SetUp() override 
        {
            // Setup for 3D ScalarGrid
            grid3D = std::make_shared<ScalarGrid>(1.0f, pmp::BoundingBox(pmp::vec3(0.0f, 0.0f, 0.0f), pmp::vec3(2.0f, 2.0f, 2.0f)));
            grid3D->Values() = { 0, 1, 2, 1, 2, 3, 2, 3, 4, 1, 2, 3, 2, 3, 4, 3, 4, 5, 2, 3, 4, 3, 4, 5, 4, 5, 6 };

            // Setup for 2D ScalarGrid
            grid2D = std::make_shared<ScalarGrid2D>(1.0f, pmp::BoundingBox2(pmp::vec2(0.0f, 0.0f), pmp::vec2(2.0f, 2.0f)));
            grid2D->Values() = { 0, 1, 2, 1, 2, 3, 2, 3, 4 };

            // Setup for 3D VectorGrid
            vectorGrid3D = std::make_shared<VectorGrid>(1.0f, pmp::BoundingBox(pmp::vec3(0.0f, 0.0f, 0.0f), pmp::vec3(2.0f, 2.0f, 2.0f)));
            vectorGrid3D->ValuesX() = { 0, 1, 2, 1, 2, 3, 2, 3, 4, 1, 2, 3, 2, 3, 4, 3, 4, 5, 2, 3, 4, 3, 4, 5, 4, 5, 6 };
            vectorGrid3D->ValuesY() = { 6, 5, 4, 5, 4, 3, 4, 3, 2, 5, 4, 3, 4, 3, 2, 3, 2, 1, 4, 3, 2, 3, 2, 1, 2, 1, 0 };
            vectorGrid3D->ValuesZ() = { 2, 3, 4, 3, 4, 5, 4, 5, 6, 3, 4, 5, 4, 5, 6, 5, 6, 7, 4, 5, 6, 5, 6, 7, 6, 7, 8 };

            // Setup for 2D VectorGrid
            vectorGrid2D = std::make_shared<VectorGrid2D>(1.0f, pmp::BoundingBox2(pmp::vec2(0.0f, 0.0f), pmp::vec2(2.0f, 2.0f)));
            vectorGrid2D->ValuesX() = { 0, 1, 2, 1, 2, 3, 2, 3, 4 };
            vectorGrid2D->ValuesY() = { 6, 5, 4, 5, 4, 3, 4, 3, 2 };
        }
    };

} // anonymous namespace

TEST_F(GridInterpolationTestFixture, BilinearInterpolateScalarValue2D_case)
{
    pmp::vec2 samplePt1(0.5f, 0.5f);
    pmp::vec2 samplePt2(1.5f, 0.5f);
    pmp::vec2 samplePt3(0.5f, 1.5f);
    pmp::vec2 samplePt4(1.5f, 1.5f);

    double result1 = BilinearInterpolateScalarValue(samplePt1, *grid2D);
    double result2 = BilinearInterpolateScalarValue(samplePt2, *grid2D);
    double result3 = BilinearInterpolateScalarValue(samplePt3, *grid2D);
    double result4 = BilinearInterpolateScalarValue(samplePt4, *grid2D);

    EXPECT_NEAR(result1, 1.0, 1e-6);
    EXPECT_NEAR(result2, 2.0, 1e-6);
    EXPECT_NEAR(result3, 2.0, 1e-6);
    EXPECT_NEAR(result4, 3.0, 1e-6);
}

TEST_F(GridInterpolationTestFixture, BilinearInterpolateVectorValue2D_case)
{
    pmp::vec2 samplePt1(0.5f, 0.5f);
    pmp::vec2 samplePt2(1.5f, 0.5f);
    pmp::vec2 samplePt3(0.5f, 1.5f);
    pmp::vec2 samplePt4(1.5f, 1.5f);

    pmp::dvec2 result1 = BilinearInterpolateVectorValue(samplePt1, *vectorGrid2D);
    pmp::dvec2 result2 = BilinearInterpolateVectorValue(samplePt2, *vectorGrid2D);
    pmp::dvec2 result3 = BilinearInterpolateVectorValue(samplePt3, *vectorGrid2D);
    pmp::dvec2 result4 = BilinearInterpolateVectorValue(samplePt4, *vectorGrid2D);

    EXPECT_NEAR(result1[0], 1.0, 1e-6);
    EXPECT_NEAR(result1[1], 5.0, 1e-6);
    EXPECT_NEAR(result2[0], 2.0, 1e-6);
    EXPECT_NEAR(result2[1], 4.0, 1e-6);
    EXPECT_NEAR(result3[0], 2.0, 1e-6);
    EXPECT_NEAR(result3[1], 4.0, 1e-6);
    EXPECT_NEAR(result4[0], 3.0, 1e-6);
    EXPECT_NEAR(result4[1], 3.0, 1e-6);
}

TEST_F(GridInterpolationTestFixture, TrilinearInterpolateScalarValue_case)
{
    pmp::vec3 samplePt1(0.5f, 0.5f, 0.5f);
    pmp::vec3 samplePt2(1.5f, 0.5f, 0.5f);
    pmp::vec3 samplePt3(0.5f, 1.5f, 0.5f);
    pmp::vec3 samplePt4(0.5f, 0.5f, 1.5f);

    double result1 = TrilinearInterpolateScalarValue(samplePt1, *grid3D);
    double result2 = TrilinearInterpolateScalarValue(samplePt2, *grid3D);
    double result3 = TrilinearInterpolateScalarValue(samplePt3, *grid3D);
    double result4 = TrilinearInterpolateScalarValue(samplePt4, *grid3D);

    EXPECT_NEAR(result1, 2.0, 1e-6);
    EXPECT_NEAR(result2, 3.0, 1e-6);
    EXPECT_NEAR(result3, 2.5, 1e-6);
    EXPECT_NEAR(result4, 3.0, 1e-6);
}

TEST_F(GridInterpolationTestFixture, TrilinearInterpolateVectorValue_case)
{
    pmp::vec3 samplePt1(0.5f, 0.5f, 0.5f);
    pmp::vec3 samplePt2(1.5f, 0.5f, 0.5f);
    pmp::vec3 samplePt3(0.5f, 1.5f, 0.5f);
    pmp::vec3 samplePt4(0.5f, 0.5f, 1.5f);

    pmp::dvec3 result1 = TrilinearInterpolateVectorValue(samplePt1, *vectorGrid3D);
    pmp::dvec3 result2 = TrilinearInterpolateVectorValue(samplePt2, *vectorGrid3D);
    pmp::dvec3 result3 = TrilinearInterpolateVectorValue(samplePt3, *vectorGrid3D);
    pmp::dvec3 result4 = TrilinearInterpolateVectorValue(samplePt4, *vectorGrid3D);

    EXPECT_NEAR(result1[0], 2.0, 1e-6);
    EXPECT_NEAR(result1[1], 4.0, 1e-6);
    EXPECT_NEAR(result1[2], 4.0, 1e-6);
    EXPECT_NEAR(result2[0], 3.0, 1e-6);
    EXPECT_NEAR(result2[1], 3.0, 1e-6);
    EXPECT_NEAR(result2[2], 4.5, 1e-6);
    EXPECT_NEAR(result3[0], 2.5, 1e-6);
    EXPECT_NEAR(result3[1], 4.5, 1e-6);
    EXPECT_NEAR(result3[2], 4.0, 1e-6);
    EXPECT_NEAR(result4[0], 3.0, 1e-6);
    EXPECT_NEAR(result4[1], 3.0, 1e-6);
    EXPECT_NEAR(result4[2], 5.0, 1e-6);
}

TEST_F(GridInterpolationTestFixture, GetNearestNeighborScalarValue_case)
{
    pmp::vec3 samplePt(1.4f, 0.6f, 0.6f);
    double result = GetNearestNeighborScalarValue(samplePt, *grid3D);
    EXPECT_NEAR(result, 3.0, 1e-6);
}

TEST_F(GridInterpolationTestFixture, GetNearestNeighborVectorValue_case)
{
    pmp::vec3 samplePt(1.4f, 0.6f, 0.6f);
    pmp::dvec3 result = GetNearestNeighborVectorValue(samplePt, *vectorGrid3D);
    EXPECT_NEAR(result[0], 3.0, 1e-6);
    EXPECT_NEAR(result[1], 3.0, 1e-6);
    EXPECT_NEAR(result[2], 4.5, 1e-6);
}

TEST_F(GridInterpolationTestFixture, GetNearestNeighborScalarValue2D_case)
{
    pmp::vec2 samplePt(1.4f, 0.6f);
    double result = GetNearestNeighborScalarValue2D(samplePt, *grid2D);
    EXPECT_NEAR(result, 2.0, 1e-6);
}

TEST_F(GridInterpolationTestFixture, GetNearestNeighborVectorValue2D_case)
{
    pmp::vec2 samplePt(1.4f, 0.6f);
    pmp::dvec2 result = GetNearestNeighborVectorValue2D(samplePt, *vectorGrid2D);
    EXPECT_NEAR(result[0], 2.0, 1e-6);
    EXPECT_NEAR(result[1], 4.0, 1e-6);
}
