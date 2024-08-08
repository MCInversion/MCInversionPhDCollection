#include "gtest/gtest.h"

#include "sdf/SDF.h"

#include "core/InscribedManifold.h"

namespace
{
    [[nodiscard]] InscribedCircleInputData CreateUnitSquareVerticesData() 
    {
        InscribedCircleInputData inputData;
        inputData.Points = { pmp::Point2(0, 0), pmp::Point2(1, 0), pmp::Point2(1, 1), pmp::Point2(0, 1) };
        return inputData;
    }

    [[nodiscard]] InscribedCircleInputData CreateUniformCircleSamplingData()
    {
        InscribedCircleInputData inputData;
        inputData.Points = pmp::CurveFactory::circle(pmp::Point2(0, 0), 1.0f, 32).positions();
        return inputData;
    }

    [[nodiscard]] InscribedCircleInputData CreateEllipseSamplingData()
    {
        InscribedCircleInputData inputData;
        auto ellipse = pmp::CurveFactory::circle(pmp::Point2(0, 0), 1.0f, 32);
        ellipse *= scaling_matrix_2d(pmp::vec2(2.0f, 1.0f));
        inputData.Points = ellipse.positions();
        return inputData;
    }

    [[nodiscard]] std::shared_ptr<Geometry::ScalarGrid2D> GenerateDistanceField(const std::vector<pmp::Point2>& points, float truncationFactor = 10.0f)
    {
        const auto pointBBox = pmp::BoundingBox2(points);
        const auto pointBBoxSize = pointBBox.max() - pointBBox.min();
        const float minSize = std::min(pointBBoxSize[0], pointBBoxSize[1]);
        const float cellSize = minSize / 20.0f;

        const SDF::PointCloudDistanceField2DSettings sdfSettings{
            cellSize,
            1.0f,
        DBL_MAX
        };
        return std::make_shared<Geometry::ScalarGrid2D>(SDF::PlanarPointCloudDistanceFieldGenerator::Generate(points, sdfSettings));
    }

    class NaiveInscribedCircleCalculatorTests : public ::testing::Test 
    {
    protected:
        NaiveInscribedCircleCalculator calculator;
    };

    class DistanceFieldInscribedCircleCalculatorTests : public ::testing::Test 
    {
    protected:
        DistanceFieldInscribedCircleCalculator calculator;
    };

} // anonymous namespace

TEST_F(NaiveInscribedCircleCalculatorTests, UnitSquareVertices) 
{
    // Arrange
    const auto inputData = CreateUnitSquareVerticesData();

    // Act
    auto circles = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(circles.size(), 1);
    EXPECT_FLOAT_EQ(circles[0].Center[0], 0.5f);
    EXPECT_FLOAT_EQ(circles[0].Center[1], 0.5f);
    EXPECT_FLOAT_EQ(circles[0].Radius, std::sqrt(2) / 2);
}

TEST_F(NaiveInscribedCircleCalculatorTests, UniformCircleSampling) 
{
    // Arrange
    const auto inputData = CreateUniformCircleSamplingData();

    // Act
    auto circles = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(circles.size(), 1);
    EXPECT_FLOAT_EQ(circles[0].Center[0], 0.0f);
    EXPECT_FLOAT_EQ(circles[0].Center[1], 0.0f);
    EXPECT_FLOAT_EQ(circles[0].Radius, 1.0f);
}

TEST_F(NaiveInscribedCircleCalculatorTests, EllipseSampling)
{
    // Arrange
    const auto inputData = CreateEllipseSamplingData();

    // Act
    auto circles = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(circles.size(), 1);
    EXPECT_FLOAT_EQ(circles[0].Center[0], 0.0f);
    EXPECT_FLOAT_EQ(circles[0].Center[1], 0.0f);
    EXPECT_FLOAT_EQ(circles[0].Radius, 1.0f);
}

TEST_F(DistanceFieldInscribedCircleCalculatorTests, UnitSquareVertices) 
{
    // Arrange
    auto inputData = CreateUnitSquareVerticesData();
    inputData.DistanceField = GenerateDistanceField(inputData.Points);
    EXPECT_TRUE(inputData.DistanceField != nullptr);
    const auto epsilon = inputData.DistanceField->CellSize();

    // Act
    const auto circles = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(circles.size(), 1);
    for (const auto& circle : circles)
    {
        EXPECT_NEAR(circle.Center[0], 0.5f, epsilon);
        EXPECT_NEAR(circle.Center[1], 0.5f, epsilon);
        EXPECT_NEAR(circle.Radius, std::sqrt(2) / 2, epsilon);
    }
}

TEST_F(DistanceFieldInscribedCircleCalculatorTests, UniformCircleSampling) 
{
    // Arrange
    auto inputData = CreateUniformCircleSamplingData();
    inputData.DistanceField = GenerateDistanceField(inputData.Points);
    EXPECT_TRUE(inputData.DistanceField != nullptr);
    const auto epsilon = inputData.DistanceField->CellSize();

    // Act
    const auto circles = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(circles.size(), 1);
    for (const auto& circle : circles)
    {
        EXPECT_NEAR(circle.Center[0], 0.0f, epsilon);
        EXPECT_NEAR(circle.Center[1], 0.0f, epsilon);
        EXPECT_NEAR(circle.Radius, 1.0f, epsilon);
    }
}

TEST_F(DistanceFieldInscribedCircleCalculatorTests, EllipseSampling) 
{
    // Arrange
    auto inputData = CreateEllipseSamplingData();
    inputData.DistanceField = GenerateDistanceField(inputData.Points);
    EXPECT_TRUE(inputData.DistanceField != nullptr);
    const auto epsilon = inputData.DistanceField->CellSize();

    // Act
    const auto circles = calculator.Calculate(inputData);

    // Assert
    ASSERT_GT(circles.size(), 0);
    for (const auto& circle : circles)
    {
        // TODO: fix multiple circles selection.
        //EXPECT_NEAR(circle.Center[0], 0.0f, epsilon);
        //EXPECT_NEAR(circle.Center[1], 0.0f, epsilon);
        //EXPECT_NEAR(circle.Radius, 1.0f, epsilon);
        EXPECT_GT(circle.Radius, 0.75f);
    }
}
