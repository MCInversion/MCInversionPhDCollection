#include "gtest/gtest.h"

#include "sdf/SDF.h"

#include "pmp/BoundingBox.h"

#include "geometry/Grid.h"
#include "geometry/GeometryUtil.h"
#include "geometry/GeometryAdapters.h"

using namespace SDF;
using namespace Geometry;
using namespace pmp;

TEST(DistanceField2DTests, PlanarDistanceFieldGenerator_SimpleBaseCurve)
{
    // Arrange
    const std::vector curveVertices = { Point2(0, 0), Point2(1, 0), Point2(1, 1), Point2(0, 1) };
    BaseCurveGeometryData curveData;
    curveData.Vertices = curveVertices;
    curveData.EdgeIndices = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
    const BaseCurveAdapter curveAdapter(std::make_shared<BaseCurveGeometryData>(curveData));

    const auto curveBBox = curveAdapter.GetBounds();
    const auto curveBBoxSize = curveBBox.max() - curveBBox.min();
    const float minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
    const float cellSize = minSize / 10.0f;

    const DistanceField2DSettings sdfSettings{
        cellSize,
        1.0f,
        0.1,
        KDTreeSplitType::Center,
        SignComputation2D::PixelFloodFill,
        PreprocessingType2D::Quadtree
    };

    // Act
    const auto sdf = PlanarDistanceFieldGenerator::Generate(curveAdapter, sdfSettings);

    // Assert
    ASSERT_TRUE(sdf.IsValid());
    EXPECT_EQ(sdf.Dimensions().Nx, 10);
    EXPECT_EQ(sdf.Dimensions().Ny, 10);
}

TEST(DistanceField2DTests, PlanarDistanceFieldGenerator_SimpleManifoldCurve)
{
    // Arrange
    const std::vector curveVertices = { Point2(0, 0), Point2(1, 0), Point2(1, 1), Point2(0, 1) };
    ManifoldCurve2D curve;
    std::vector<Vertex> vertices;
    std::ranges::transform(curveVertices, std::back_inserter(vertices), [&curve](const auto v) { return curve.add_vertex(v); });
    for (size_t i = 0; i < vertices.size() - 1; ++i) { curve.new_edge(vertices[i], vertices[i + 1]); }
	const ManifoldCurve2DAdapter curveAdapter(std::make_shared<ManifoldCurve2D>(curve));

    const auto curveBBox = curveAdapter.GetBounds();
    const auto curveBBoxSize = curveBBox.max() - curveBBox.min();
    const float minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
    const float cellSize = minSize / 10.0f;

    const DistanceField2DSettings sdfSettings{
        cellSize,
        1.0f,
        0.1,
        KDTreeSplitType::Center,
        SignComputation2D::PixelFloodFill,
        PreprocessingType2D::Quadtree
    };

    // Act
    const auto sdf = PlanarDistanceFieldGenerator::Generate(curveAdapter, sdfSettings);

    // Assert
    ASSERT_TRUE(sdf.IsValid());
    EXPECT_EQ(sdf.Dimensions().Nx, 10);
    EXPECT_EQ(sdf.Dimensions().Ny, 10);
}

TEST(DistanceField2DTests, PlanarPointCloudDistanceFieldGenerator_SimplePointCloud)
{
    // Arrange
    const std::vector points = { Point2(0, 0), Point2(1, 0), Point2(1, 1), Point2(0, 1) };
    const auto pointBBox = BoundingBox2(Point2(0, 0), Point2(1, 1));
    const auto pointBBoxSize = pointBBox.max() - pointBBox.min();
    const float minSize = std::min(pointBBoxSize[0], pointBBoxSize[1]);
    const float cellSize = minSize / 10.0f;

    const PointCloudDistanceField2DSettings sdfSettings{
        cellSize,
        1.0f,
        0.1
    };

    // Act
    const auto sdf = PlanarPointCloudDistanceFieldGenerator::Generate(points, sdfSettings);

    // Assert
    ASSERT_TRUE(sdf.IsValid());
    EXPECT_EQ(sdf.Dimensions().Nx, 10);
    EXPECT_EQ(sdf.Dimensions().Ny, 10);
}

// Add more tests to verify edge cases and different settings
