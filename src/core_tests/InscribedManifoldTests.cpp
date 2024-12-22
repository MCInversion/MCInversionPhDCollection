#include "gtest/gtest.h"

#include "sdf/SDF.h"

#include "geometry/GridUtil.h"

#include "core/InscribedManifold.h"
#include "core/ConversionUtils.h"

#include <filesystem>

// set up root directory
const std::filesystem::path fsRootPath = DROOT_DIR;
const auto fsDataDirPath = fsRootPath / "data\\";
const auto fsDataOutPath = fsRootPath / "output\\";
const std::string dataDirPath = fsDataDirPath.string();
const std::string dataOutPath = fsDataOutPath.string();

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
        inputData.Points = pmp::CurveFactory::circle(pmp::Point2(0, 0), 1.0, 32).positions();
        return inputData;
    }

    [[nodiscard]] InscribedCircleInputData CreateEllipseSamplingData()
    {
        InscribedCircleInputData inputData;
        auto ellipse = pmp::CurveFactory::circle(pmp::Point2(0, 0), 1.0, 32);
        ellipse *= scaling_matrix_2d(pmp::vec2(2.0, 1.0));
        inputData.Points = ellipse.positions();
        return inputData;
    }

    [[nodiscard]] std::shared_ptr<Geometry::ScalarGrid2D> GenerateDistanceField(const std::vector<pmp::Point2>& points, float truncationFactor = 10.0)
    {
        const auto pointBBox = pmp::BoundingBox2(points);
        const auto pointBBoxSize = pointBBox.max() - pointBBox.min();
        const pmp::Scalar minSize = std::min(pointBBoxSize[0], pointBBoxSize[1]);
        const pmp::Scalar cellSize = minSize / 20.0;

        const SDF::PointCloudDistanceField2DSettings sdfSettings{
            cellSize,
            1.0,
        DBL_MAX
        };
        return std::make_shared<Geometry::ScalarGrid2D>(SDF::PlanarPointCloudDistanceFieldGenerator::Generate(points, sdfSettings));
    }

    // --------------------------------------------------------------

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

    class HierarchicalDistanceFieldInscribedCircleCalculatorTests : public ::testing::Test
    {
    protected:
        HierarchicalDistanceFieldInscribedCircleCalculator calculator;
    };

    class ParticleSwarmDistanceFieldInscribedCircleCalculatorTests : public ::testing::Test
    {
    protected:
        ParticleSwarmDistanceFieldInscribedCircleCalculator calculator;
    };

    // --------------------------------------------------------------

    class NaiveInscribedSphereCalculatorTests : public ::testing::Test
    {
    protected:
        NaiveInscribedSphereCalculator calculator;
    };

    class DistanceFieldInscribedSphereCalculatorTests : public ::testing::Test
    {
    protected:
        DistanceFieldInscribedSphereCalculator calculator;
    };

    class HierarchicalDistanceFieldInscribedSphereCalculatorTests : public ::testing::Test
    {
    protected:
        HierarchicalDistanceFieldInscribedSphereCalculator calculator;
    };

    class ParticleSwarmDistanceFieldInscribedSphereCalculatorTests : public ::testing::Test
    {
    protected:
        ParticleSwarmDistanceFieldInscribedSphereCalculator calculator;
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
    EXPECT_FLOAT_EQ(circles[0].Center[0], 0.5);
    EXPECT_FLOAT_EQ(circles[0].Center[1], 0.5);
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
    EXPECT_FLOAT_EQ(circles[0].Center[0], 0.0);
    EXPECT_FLOAT_EQ(circles[0].Center[1], 0.0);
    EXPECT_FLOAT_EQ(circles[0].Radius, 1.0);
}

TEST_F(NaiveInscribedCircleCalculatorTests, EllipseSampling)
{
    // Arrange
    const auto inputData = CreateEllipseSamplingData();

    // Act
    auto circles = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(circles.size(), 1);
    EXPECT_FLOAT_EQ(circles[0].Center[0], 0.0);
    EXPECT_FLOAT_EQ(circles[0].Center[1], 0.0);
    EXPECT_FLOAT_EQ(circles[0].Radius, 1.0);
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
        EXPECT_NEAR(circle.Center[0], 0.5, epsilon);
        EXPECT_NEAR(circle.Center[1], 0.5, epsilon);
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
        EXPECT_NEAR(circle.Center[0], 0.0, epsilon);
        EXPECT_NEAR(circle.Center[1], 0.0, epsilon);
        EXPECT_NEAR(circle.Radius, 1.0, epsilon);
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
    ASSERT_EQ(circles.size(), 1);
    for (const auto& circle : circles)
    {
        EXPECT_NEAR(circle.Radius, 1.0, epsilon);
    }
}

TEST_F(HierarchicalDistanceFieldInscribedCircleCalculatorTests, UnitSquareVertices)
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
        EXPECT_NEAR(circle.Center[0], 0.5, epsilon);
        EXPECT_NEAR(circle.Center[1], 0.5, epsilon);
        EXPECT_NEAR(circle.Radius, std::sqrt(2) / 2, epsilon);
    }
}

TEST_F(HierarchicalDistanceFieldInscribedCircleCalculatorTests, UniformCircleSampling)
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
        EXPECT_NEAR(circle.Center[0], 0.0, epsilon);
        EXPECT_NEAR(circle.Center[1], 0.0, epsilon);
        EXPECT_NEAR(circle.Radius, 1.0, epsilon);
    }
}

TEST_F(HierarchicalDistanceFieldInscribedCircleCalculatorTests, EllipseSampling)
{
    // Arrange
    auto inputData = CreateEllipseSamplingData();
    inputData.DistanceField = GenerateDistanceField(inputData.Points);
    EXPECT_TRUE(inputData.DistanceField != nullptr);
    const auto epsilon = inputData.DistanceField->CellSize();

    // Act
    const auto circles = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(circles.size(), 1);
    for (const auto& circle : circles)
    {
        EXPECT_NEAR(circle.Radius, 1.0, epsilon);
    }
}

TEST_F(ParticleSwarmDistanceFieldInscribedCircleCalculatorTests, UnitSquareVertices)
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
        EXPECT_NEAR(circle.Center[0], 0.5, epsilon);
        EXPECT_NEAR(circle.Center[1], 0.5, epsilon);
        EXPECT_NEAR(circle.Radius, std::sqrt(2) / 2, epsilon);
    }
}

TEST_F(ParticleSwarmDistanceFieldInscribedCircleCalculatorTests, UniformCircleSampling)
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
        EXPECT_NEAR(circle.Center[0], 0.0, epsilon);
        EXPECT_NEAR(circle.Center[1], 0.0, epsilon);
        EXPECT_NEAR(circle.Radius, 1.0, epsilon);
    }
}

TEST_F(ParticleSwarmDistanceFieldInscribedCircleCalculatorTests, EllipseSampling)
{
    // Arrange
    auto inputData = CreateEllipseSamplingData();
    inputData.DistanceField = GenerateDistanceField(inputData.Points);
    EXPECT_TRUE(inputData.DistanceField != nullptr);
    const auto epsilon = 1.5 * inputData.DistanceField->CellSize();

    // Act
    const auto circles = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(circles.size(), 1);
    for (const auto& circle : circles)
    {
        EXPECT_NEAR(circle.Radius, 1.0, epsilon);
    }
}

TEST(ParticleSwarmDistanceFieldInscribedCircleCalculator_Tests, ProblematicIncompleteCircle)
{
    // Arrange
    pmp::ManifoldCurve2D targetCurve = pmp::CurveFactory::circle(pmp::Point2(0.0, 0.0), 0.75, 16);
    auto targetPts = targetCurve.positions();
    targetPts.erase(targetPts.begin());
    targetPts.erase(targetPts.begin());
    targetPts.erase(targetPts.begin());
    InscribedCircleInputData inputData;
    inputData.Points = targetPts;
    inputData.DistanceField = GenerateDistanceField(inputData.Points);

    // Act
    ParticleSwarmDistanceFieldInscribedCircleCalculator calculator;
    const auto circles = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(circles.size(), 1);
}

// =======================================================================================================================================

namespace
{
    [[nodiscard]] pmp::SurfaceMesh ConstructIcoSphere(const pmp::Point& center, const pmp::Scalar& radius, const unsigned int& subdiv)
    {
        Geometry::IcoSphereBuilder icoBuilder({ subdiv, radius });
        icoBuilder.BuildBaseData();
        icoBuilder.BuildPMPSurfaceMesh();
        if (center == pmp::Point(0, 0, 0))
            return icoBuilder.GetPMPSurfaceMeshResult();

        auto mesh = icoBuilder.GetPMPSurfaceMeshResult();
        const auto translationMatrix = translation_matrix(center);
        mesh *= translationMatrix;
        return mesh;
    }

    [[nodiscard]] InscribedSphereInputData CreateUnitCubeVerticesData()
    {
        InscribedSphereInputData data;
        data.Points = {
            pmp::Point(0, 0, 0),
            pmp::Point(1, 0, 0),
            pmp::Point(0, 1, 0),
            pmp::Point(0, 0, 1),
            pmp::Point(1, 1, 0),
            pmp::Point(1, 0, 1),
            pmp::Point(0, 1, 1),
            pmp::Point(1, 1, 1)
        };
        return data;
    }

    [[nodiscard]] InscribedSphereInputData CreateUniformSphereSamplingData()
    {
        pmp::SurfaceMesh sphereMesh = ConstructIcoSphere(pmp::Point(0, 0, 0), 1.0, 2);
        InscribedSphereInputData data;
        data.Points = sphereMesh.positions();
        return data;
    }

    [[nodiscard]] InscribedSphereInputData CreateEllipsoidSamplingData()
    {
        pmp::SurfaceMesh sphereMesh = ConstructIcoSphere(pmp::Point(0, 0, 0), 1.0, 2);
        constexpr pmp::Scalar a = 1.0;
        constexpr pmp::Scalar b = 1.5;
        constexpr pmp::Scalar c = 2.0;
        const auto scalingMat = scaling_matrix(pmp::vec3{a, b, c});
        sphereMesh *= scalingMat;
        InscribedSphereInputData data;
        data.Points = sphereMesh.positions();
        return data;
    }

    [[nodiscard]] std::shared_ptr<Geometry::ScalarGrid> GenerateDistanceField(const std::vector<pmp::Point>& points, pmp::Scalar truncationFactor = 10.0)
    {
        const auto pointBBox = pmp::BoundingBox(points);
        const auto pointBBoxSize = pointBBox.max() - pointBBox.min();
        const pmp::Scalar minSize = std::min({pointBBoxSize[0], pointBBoxSize[1], pointBBoxSize[2]});
        const pmp::Scalar cellSize = minSize / 20.0;

        const SDF::PointCloudDistanceFieldSettings sdfSettings{
            cellSize,
            1.0,
            DBL_MAX
        };
        return std::make_shared<Geometry::ScalarGrid>(SDF::PointCloudDistanceFieldGenerator::Generate(points, sdfSettings));
    }

} // anonymous namespace


TEST_F(NaiveInscribedSphereCalculatorTests, UnitCubeVertices)
{
    // Arrange
    const auto inputData = CreateUnitCubeVerticesData();

    // Act
    auto spheres = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(spheres.size(), 1);
    EXPECT_FLOAT_EQ(spheres[0].Center[0], 0.5);
    EXPECT_FLOAT_EQ(spheres[0].Center[1], 0.5);
    EXPECT_FLOAT_EQ(spheres[0].Center[2], 0.5);
    EXPECT_FLOAT_EQ(spheres[0].Radius, std::sqrt(3) / 2);  // Inscribed sphere in a unit cube
}

TEST_F(NaiveInscribedSphereCalculatorTests, UniformSphereSampling)
{
    // Arrange
    const auto inputData = CreateUniformSphereSamplingData();

    // Act
    auto spheres = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(spheres.size(), 1);
    EXPECT_FLOAT_EQ(spheres[0].Center[0], 0.0);
    EXPECT_FLOAT_EQ(spheres[0].Center[1], 0.0);
    EXPECT_FLOAT_EQ(spheres[0].Center[2], 0.0);
    EXPECT_FLOAT_EQ(spheres[0].Radius, 1.0);
}

TEST_F(NaiveInscribedSphereCalculatorTests, EllipsoidSampling)
{
    // Arrange
    const auto inputData = CreateEllipsoidSamplingData();

    // Act
    auto spheres = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(spheres.size(), 1);
    EXPECT_FLOAT_EQ(spheres[0].Center[0], 0.0);
    EXPECT_FLOAT_EQ(spheres[0].Center[1], 0.0);
    EXPECT_FLOAT_EQ(spheres[0].Center[2], 0.0);
    EXPECT_FLOAT_EQ(spheres[0].Radius, 1.0);
}

TEST_F(DistanceFieldInscribedSphereCalculatorTests, UnitCubeVertices)
{
    // Arrange
    auto inputData = CreateUnitCubeVerticesData();
    inputData.DistanceField = GenerateDistanceField(inputData.Points);
    EXPECT_TRUE(inputData.DistanceField != nullptr);
    const auto epsilon = inputData.DistanceField->CellSize();

    // Act
    const auto spheres = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(spheres.size(), 1);
    for (const auto& sphere : spheres)
    {
        EXPECT_NEAR(sphere.Center[0], 0.5, epsilon);
        EXPECT_NEAR(sphere.Center[1], 0.5, epsilon);
        EXPECT_NEAR(sphere.Center[2], 0.5, epsilon);
        EXPECT_NEAR(sphere.Radius, std::sqrt(3) / 2, epsilon);
    }
}

TEST_F(DistanceFieldInscribedSphereCalculatorTests, UniformSphereSampling)
{
    // Arrange
    auto inputData = CreateUniformSphereSamplingData();
    inputData.DistanceField = GenerateDistanceField(inputData.Points);
    EXPECT_TRUE(inputData.DistanceField != nullptr);
    const auto epsilon = inputData.DistanceField->CellSize();

    //const auto subDF = Geometry::ExtractSubGrid(*inputData.DistanceField, 28, 28, 28, 31, 31, 31);
    //const auto subDF = Geometry::ExtractSubGrid(*inputData.DistanceField, 27, 28, 28, 30, 31, 31);
    //ExportToVTI(dataOutPath + "\\core_tests\\UniformSphereSampling_subDF", subDF);

    // Act
    const auto spheres = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(spheres.size(), 1);
    for (const auto& sphere : spheres)
    {
        EXPECT_NEAR(sphere.Center[0], 0.0, epsilon);
        EXPECT_NEAR(sphere.Center[1], 0.0, epsilon);
        EXPECT_NEAR(sphere.Center[2], 0.0, epsilon);
        EXPECT_NEAR(sphere.Radius, 1.0, epsilon);
    }
}

TEST_F(DistanceFieldInscribedSphereCalculatorTests, EllipsoidSampling)
{
    // Arrange
    auto inputData = CreateEllipsoidSamplingData();
    inputData.DistanceField = GenerateDistanceField(inputData.Points);
    EXPECT_TRUE(inputData.DistanceField != nullptr);
    const auto epsilon = inputData.DistanceField->CellSize();

    //ExportToVTI(dataOutPath + "\\core_tests\\EllipsoidSampling_SDF", *inputData.DistanceField);

    // Act
    const auto spheres = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(spheres.size(), 1);
    for (const auto& sphere : spheres)
    {
        EXPECT_NEAR(sphere.Radius, 1.0, epsilon);
    }
}

TEST_F(HierarchicalDistanceFieldInscribedSphereCalculatorTests, UnitCubeVertices)
{
    // Arrange
    auto inputData = CreateUnitCubeVerticesData();
    inputData.DistanceField = GenerateDistanceField(inputData.Points);
    EXPECT_TRUE(inputData.DistanceField != nullptr);
    const auto epsilon = inputData.DistanceField->CellSize();

    // Act
    const auto spheres = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(spheres.size(), 1);
    for (const auto& sphere : spheres)
    {
        EXPECT_NEAR(sphere.Center[0], 0.5, epsilon);
        EXPECT_NEAR(sphere.Center[1], 0.5, epsilon);
        EXPECT_NEAR(sphere.Center[2], 0.5, epsilon);
        EXPECT_NEAR(sphere.Radius, std::sqrt(3) / 2, epsilon);
    }
}

TEST_F(HierarchicalDistanceFieldInscribedSphereCalculatorTests, UniformSphereSampling)
{
    // Arrange
    auto inputData = CreateUniformSphereSamplingData();
    inputData.DistanceField = GenerateDistanceField(inputData.Points);
    EXPECT_TRUE(inputData.DistanceField != nullptr);
    const auto epsilon = inputData.DistanceField->CellSize();

    //ExportToVTI(dataOutPath + "\\core_tests\\UniformSphereSampling_SDF", *inputData.DistanceField);

    // Act
    const auto spheres = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(spheres.size(), 1);
    for (const auto& sphere : spheres)
    {
        EXPECT_NEAR(sphere.Center[0], 0.0, epsilon);
        EXPECT_NEAR(sphere.Center[1], 0.0, epsilon);
        EXPECT_NEAR(sphere.Center[2], 0.0, epsilon);
        EXPECT_NEAR(sphere.Radius, 1.0, epsilon);
    }
}

TEST_F(HierarchicalDistanceFieldInscribedSphereCalculatorTests, EllipsoidSampling)
{
    // Arrange
    auto inputData = CreateEllipsoidSamplingData();
    inputData.DistanceField = GenerateDistanceField(inputData.Points);
    EXPECT_TRUE(inputData.DistanceField != nullptr);
    const auto epsilon = inputData.DistanceField->CellSize();

    // Act
    const auto spheres = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(spheres.size(), 1);
    for (const auto& sphere : spheres)
    {
        EXPECT_NEAR(sphere.Radius, 1.0, epsilon);
    }
}

TEST_F(ParticleSwarmDistanceFieldInscribedSphereCalculatorTests, UnitCubeVertices)
{
    // Arrange
    auto inputData = CreateUnitCubeVerticesData();
    inputData.DistanceField = GenerateDistanceField(inputData.Points);
    EXPECT_TRUE(inputData.DistanceField != nullptr);
    const auto epsilon = inputData.DistanceField->CellSize();

    // Act
    const auto spheres = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(spheres.size(), 1);
    for (const auto& sphere : spheres)
    {
        EXPECT_NEAR(sphere.Center[0], 0.5, epsilon);
        EXPECT_NEAR(sphere.Center[1], 0.5, epsilon);
        EXPECT_NEAR(sphere.Center[2], 0.5, epsilon);
        EXPECT_NEAR(sphere.Radius, std::sqrt(3) / 2, epsilon);
    }
}

TEST_F(ParticleSwarmDistanceFieldInscribedSphereCalculatorTests, UniformSphereSampling)
{
    // Arrange
    auto inputData = CreateUniformSphereSamplingData();
    inputData.DistanceField = GenerateDistanceField(inputData.Points);
    EXPECT_TRUE(inputData.DistanceField != nullptr);
    const auto epsilon = inputData.DistanceField->CellSize();

    // Act
    const auto spheres = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(spheres.size(), 1);
    for (const auto& sphere : spheres)
    {
        EXPECT_NEAR(sphere.Center[0], 0.0, epsilon);
        EXPECT_NEAR(sphere.Center[1], 0.0, epsilon);
        EXPECT_NEAR(sphere.Center[2], 0.0, epsilon);
        EXPECT_NEAR(sphere.Radius, 1.0, epsilon);
    }
}

TEST_F(ParticleSwarmDistanceFieldInscribedSphereCalculatorTests, EllipsoidSampling)
{
    // Arrange
    auto inputData = CreateEllipsoidSamplingData();
    inputData.DistanceField = GenerateDistanceField(inputData.Points);
    EXPECT_TRUE(inputData.DistanceField != nullptr);
    const auto epsilon = 1.5 * inputData.DistanceField->CellSize();

    // Act
    const auto spheres = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(spheres.size(), 1);
    for (const auto& sphere : spheres)
    {
        EXPECT_NEAR(sphere.Radius, 1.0, epsilon);
    }
}

TEST(ParticleSwarmDistanceFieldInscribedSphereCalculator_Tests, ProblematicIncompleteSphere)
{
    // Arrange
    pmp::SurfaceMesh targetSurface = ConstructIcoSphere(pmp::Point(0.0, 0.0, 0.0), 0.75, 1);
    auto targetPts = targetSurface.positions();
    // Remove a few points to simulate incompleteness
    targetPts.erase(targetPts.begin(), targetPts.begin() + 3);
    InscribedSphereInputData inputData;
    inputData.Points = targetPts;
    inputData.DistanceField = GenerateDistanceField(inputData.Points);

    // Act
    ParticleSwarmDistanceFieldInscribedSphereCalculator calculator;
    const auto spheres = calculator.Calculate(inputData);

    // Assert
    ASSERT_EQ(spheres.size(), 1);
}


