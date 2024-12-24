#include "gtest/gtest.h"

#include "sdf/SDF.h"

#include "pmp/BoundingBox.h"

#include "geometry/Grid.h"
#include "geometry/GridUtil.h"
#include "geometry/GeometryUtil.h"
#include "geometry/GeometryAdapters.h"
#include "geometry/GeometryIOUtils.h"
#include "pmp/algorithms/CurveFactory.h"

using namespace SDF;
using namespace Geometry;
using namespace pmp;

// Hashes of the tested distance fields
constexpr size_t HASH_square_SDF{ 11463448448184433784 };
constexpr size_t HASH_openSquare_SDF{ 9360073553033082223 };
constexpr size_t HASH_squarePts_SDF{ 12327925769871880303 };

#define EXPORT_FIELDS_AND_SOURCES true

#if EXPORT_FIELDS_AND_SOURCES

#include <filesystem>
// set up root directory
const std::filesystem::path fsRootPath = DROOT_DIR;
const auto fsDataDirPath = fsRootPath / "data\\";
const auto fsDataOutPath = fsRootPath / "output\\";
const std::string dataDirPath = fsDataDirPath.string();
const std::string dataOutPath = fsDataOutPath.string();

static void TestExportDFAndPointCloud(const Geometry::ScalarGrid2D& df, const std::vector<Point2>& pts, const std::string& absFileName)
{
    ExportScalarGridDimInfo2D(absFileName + ".gdim2d", df);
    constexpr double colorMapPlotScaleFactor = 1.0; // scale the distance field color map down to show more detail
    ExportScalarGrid2DToPNG(absFileName + ".png", df,
        Geometry::BilinearInterpolateScalarValue,
        //Geometry::GetNearestNeighborScalarValue2D,
        10, 10, RAINBOW_TO_WHITE_MAP * colorMapPlotScaleFactor);
    if (!Export2DPointCloudToPLY(pts, absFileName + ".ply"))
        std::cerr << "TestExportDFAndPointCloud: failed Export2DPointCloudToPLY!\n";
}

#define EXPORT_FIELD_AND_PT_CLOUD(df, pts, name) \
    do { \
        const auto fileName = dataOutPath + "/sdf_tests/" + (name); \
        TestExportDFAndPointCloud(df, pts, fileName); \
    } while (0)

static void TestExportDFAndCurve(const Geometry::ScalarGrid2D& df, const ManifoldCurve2D& curve, const std::string& absFileName)
{
    ExportScalarGridDimInfo2D(absFileName + ".gdim2d", df);
    constexpr double colorMapPlotScaleFactor = 1.0; // scale the distance field color map down to show more detail
    ExportScalarGrid2DToPNG(absFileName + ".png", df,
        Geometry::BilinearInterpolateScalarValue,
        //Geometry::GetNearestNeighborScalarValue2D,
        10, 10, RAINBOW_TO_WHITE_MAP * colorMapPlotScaleFactor);
    if (!ExportManifoldCurve2DToPLY(curve, absFileName + ".ply"))
        std::cerr << "TestExportDFAndCurve: failed ExportManifoldCurve2DToPLY!\n";
}

#define EXPORT_FIELD_AND_CURVE(df, curve, name) \
    do { \
        const auto fileName = dataOutPath + "/sdf_tests/" + (name); \
        TestExportDFAndCurve(df, curve, fileName); \
    } while (0)

static void TestExportDFAndBaseCurve(const Geometry::ScalarGrid2D& df, const Geometry::BaseCurveGeometryData& curve, const std::string& absFileName)
{
    ExportScalarGridDimInfo2D(absFileName + ".gdim2d", df);
    constexpr double colorMapPlotScaleFactor = 1.0; // scale the distance field color map down to show more detail
    ExportScalarGrid2DToPNG(absFileName + ".png", df,
        Geometry::BilinearInterpolateScalarValue,
        //Geometry::GetNearestNeighborScalarValue2D,
        10, 10, RAINBOW_TO_WHITE_MAP * colorMapPlotScaleFactor);
    if (!ExportBaseCurveGeometryDataToPLY(curve, absFileName + ".ply"))
        std::cerr << "TestExportDFAndBaseCurve: failed ExportBaseCurveGeometryDataToPLY!\n";
}

#define EXPORT_FIELD_AND_BASE_CURVE(df, curve, name) \
    do { \
        const auto fileName = dataOutPath + "/sdf_tests/" + (name); \
        TestExportDFAndBaseCurve(df, curve, fileName); \
    } while (0)

#else
#define EXPORT_FIELD_AND_PT_CLOUD(df, pts, name) do {} while (0)
#define EXPORT_FIELD_AND_CURVE(df, curve, name) do {} while (0)
#define EXPORT_FIELD_AND_BASE_CURVE(df, curve, name) do {} while (0)
#endif

TEST(DistanceField2DTests, PlanarDistanceFieldGenerator_SimpleClosedBaseCurve)
{
    // Arrange
    const std::vector curveVertices = { Point2(0, 0), Point2(1, 0), Point2(1, 1), Point2(0, 1) };
    BaseCurveGeometryData curveData;
    curveData.Vertices = curveVertices;
    curveData.EdgeIndices = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
    const BaseCurveAdapter curveAdapter(std::make_shared<BaseCurveGeometryData>(curveData));

    const auto curveBBox = curveAdapter.GetBounds();
    const auto curveBBoxSize = curveBBox.max() - curveBBox.min();
    const Scalar minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
    const Scalar cellSize = minSize / 10.0;

    const DistanceField2DSettings sdfSettings{
        cellSize,
        1.0,
        DBL_MAX,
        KDTreeSplitType::Center,
        SignComputation2D::PixelFloodFill,
        PreprocessingType2D::Quadtree
    };

    // Act
    const auto sdf = PlanarDistanceFieldGenerator::Generate(curveAdapter, sdfSettings);
    EXPORT_FIELD_AND_BASE_CURVE(sdf, curveAdapter.Get(), "SimpleClosedBaseCurve");

    // Assert
    ASSERT_TRUE(sdf.IsValid());
    EXPECT_EQ(HashScalarGrid(sdf), HASH_square_SDF);
    ASSERT_EQ(sdf.Dimensions().Nx, 31);
    ASSERT_EQ(sdf.Dimensions().Ny, 31);
    auto index = [&sdf](size_t i, size_t j) {
        return i + j * sdf.Dimensions().Nx;
    };
    EXPECT_TRUE(sdf.Values()[index(14, 14)] < 0.0);
    EXPECT_NEAR(sdf.Values()[index(1, 1)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(0, 29)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 29)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 0)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(0, 14)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(15, 0)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 14)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(14, 29)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(14, 14)], -0.5, cellSize);
}

TEST(DistanceField2DTests, PlanarDistanceFieldGenerator_SimpleClosedBaseCurveUnsigned)
{
    // Arrange
    const std::vector curveVertices = { Point2(0, 0), Point2(1, 0), Point2(1, 1), Point2(0, 1) };
    BaseCurveGeometryData curveData;
    curveData.Vertices = curveVertices;
    curveData.EdgeIndices = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
    const BaseCurveAdapter curveAdapter(std::make_shared<BaseCurveGeometryData>(curveData));

    const auto curveBBox = curveAdapter.GetBounds();
    const auto curveBBoxSize = curveBBox.max() - curveBBox.min();
    const Scalar minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
    const Scalar cellSize = minSize / 10.0;

    const DistanceField2DSettings sdfSettings{
        cellSize,
        1.0,
        DBL_MAX,
        KDTreeSplitType::Center,
        SignComputation2D::None,
        PreprocessingType2D::Quadtree
    };

    // Act
    const auto sdf = PlanarDistanceFieldGenerator::Generate(curveAdapter, sdfSettings);
    EXPORT_FIELD_AND_BASE_CURVE(sdf, curveAdapter.Get(), "SimpleClosedBaseCurveUnsigned");

    // Assert
    ASSERT_TRUE(sdf.IsValid());
    EXPECT_EQ(HashScalarGrid(sdf), HASH_square_SDF);
    ASSERT_EQ(sdf.Dimensions().Nx, 31);
    ASSERT_EQ(sdf.Dimensions().Ny, 31);
    auto index = [&sdf](size_t i, size_t j) {
        return i + j * sdf.Dimensions().Nx;
    };
    EXPECT_TRUE(sdf.Values()[index(14, 14)] < 0.0);
    EXPECT_NEAR(sdf.Values()[index(1, 1)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(0, 29)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 29)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 0)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(0, 14)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(15, 0)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 14)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(14, 29)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(14, 14)], -0.5, cellSize);
}

TEST(DistanceField2DTests, PlanarDistanceFieldGenerator_SimpleClosedBaseCurveUnsignedInFirstQuadrant)
{
    // Arrange
    const vec2 translateToFirstQuadrant{ 2.5, 3.7 };
    const std::vector curveVertices = {
    	Point2(0, 0) + translateToFirstQuadrant,
    	Point2(1, 0) + translateToFirstQuadrant,
    	Point2(1, 1) + translateToFirstQuadrant,
    	Point2(0, 1) + translateToFirstQuadrant
    };
    BaseCurveGeometryData curveData;
    curveData.Vertices = curveVertices;
    curveData.EdgeIndices = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
    const BaseCurveAdapter curveAdapter(std::make_shared<BaseCurveGeometryData>(curveData));

    const auto curveBBox = curveAdapter.GetBounds();
    const auto curveBBoxSize = curveBBox.max() - curveBBox.min();
    const Scalar minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
    const Scalar cellSize = minSize / 10.0;

    const DistanceField2DSettings sdfSettings{
        cellSize,
        1.0,
        DBL_MAX,
        KDTreeSplitType::Center,
        SignComputation2D::None,
        PreprocessingType2D::Quadtree
    };

    // Act
    const auto sdf = PlanarDistanceFieldGenerator::Generate(curveAdapter, sdfSettings);
    EXPORT_FIELD_AND_BASE_CURVE(sdf, curveAdapter.Get(), "SimpleClosedBaseCurveUnsignedInFirstQuadrant");

    // Assert
    ASSERT_TRUE(sdf.IsValid());
}

TEST(DistanceField2DTests, PlanarDistanceFieldGenerator_SimpleClosedManifoldCurveUnsignedInFirstQuadrant)
{
    // Arrange
    const vec2 translateToFirstQuadrant{ 2.5, 3.7 };
    const std::vector curveVertices = {
        Point2(0, 0) + translateToFirstQuadrant,
        Point2(1, 0) + translateToFirstQuadrant,
        Point2(1, 1) + translateToFirstQuadrant,
        Point2(0, 1) + translateToFirstQuadrant
    };
    ManifoldCurve2D curve;
    std::vector<Vertex> vertices;
    std::ranges::transform(curveVertices, std::back_inserter(vertices), [&curve](const auto v) { return curve.add_vertex(v); });
    for (size_t i = 0; i < vertices.size(); ++i) { curve.new_edge(vertices[i], vertices[(i + 1) % vertices.size()]); }
    const ManifoldCurve2DAdapter curveAdapter(std::make_shared<ManifoldCurve2D>(curve));

    const auto curveBBox = curveAdapter.GetBounds();
    const auto curveBBoxSize = curveBBox.max() - curveBBox.min();
    const Scalar minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
    const Scalar cellSize = minSize / 10.0;

    const DistanceField2DSettings sdfSettings{
        cellSize,
        1.0,
        DBL_MAX,
        KDTreeSplitType::Center,
        SignComputation2D::None,
        PreprocessingType2D::Quadtree
    };

    // Act
    const auto sdf = PlanarDistanceFieldGenerator::Generate(curveAdapter, sdfSettings);
    EXPORT_FIELD_AND_CURVE(sdf, curveAdapter.Get(), "SimpleClosedManifoldCurveUnsignedInFirstQuadrant");

    // Assert
    ASSERT_TRUE(sdf.IsValid());
}

TEST(DistanceField2DTests, PlanarDistanceFieldGenerator_SimpleOpenBaseCurve)
{
    // Arrange
    const std::vector curveVertices = { Point2(0, 0), Point2(1, 0), Point2(1, 1), Point2(0, 1) };
    BaseCurveGeometryData curveData;
    curveData.Vertices = curveVertices;
    curveData.EdgeIndices = { {0, 1}, {1, 2}, {2, 3} };
    const BaseCurveAdapter curveAdapter(std::make_shared<BaseCurveGeometryData>(curveData));

    const auto curveBBox = curveAdapter.GetBounds();
    const auto curveBBoxSize = curveBBox.max() - curveBBox.min();
    const Scalar minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
    const Scalar cellSize = minSize / 10.0;

    const DistanceField2DSettings sdfSettings{
        cellSize,
        1.0,
        DBL_MAX,
        KDTreeSplitType::Center,
        SignComputation2D::None,
        PreprocessingType2D::Quadtree
    };

    // Act
    const auto sdf = PlanarDistanceFieldGenerator::Generate(curveAdapter, sdfSettings);
    EXPORT_FIELD_AND_BASE_CURVE(sdf, curveAdapter.Get(), "SimpleOpenBaseCurve");

    // Assert
    ASSERT_TRUE(sdf.IsValid());
    EXPECT_EQ(HashScalarGrid(sdf), HASH_openSquare_SDF);
    ASSERT_EQ(sdf.Dimensions().Nx, 31);
    ASSERT_EQ(sdf.Dimensions().Ny, 31);
    auto index = [&sdf](size_t i, size_t j) {
        return i + j * sdf.Dimensions().Nx;
    };
    EXPECT_FALSE(sdf.Values()[index(14, 14)] < 0.0);
    EXPECT_NEAR(sdf.Values()[index(1, 1)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(0, 29)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 29)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 0)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(14, 0)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 14)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(14, 29)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(14, 14)], 0.5, cellSize);
}

TEST(DistanceField2DTests, PlanarDistanceFieldGenerator_SimpleClosedManifoldCurve)
{
    // Arrange
    const std::vector curveVertices = { Point2(0, 0), Point2(1, 0), Point2(1, 1), Point2(0, 1) };
    ManifoldCurve2D curve;
    std::vector<Vertex> vertices;
    std::ranges::transform(curveVertices, std::back_inserter(vertices), [&curve](const auto v) { return curve.add_vertex(v); });
    for (size_t i = 0; i < vertices.size(); ++i) { curve.new_edge(vertices[i], vertices[(i + 1) % vertices.size()]); }
	const ManifoldCurve2DAdapter curveAdapter(std::make_shared<ManifoldCurve2D>(curve));

    const auto curveBBox = curveAdapter.GetBounds();
    const auto curveBBoxSize = curveBBox.max() - curveBBox.min();
    const Scalar minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
    const Scalar cellSize = minSize / 10.0;

    const DistanceField2DSettings sdfSettings{
        cellSize,
        1.0,
        DBL_MAX,
        KDTreeSplitType::Center,
        SignComputation2D::PixelFloodFill,
        PreprocessingType2D::Quadtree
    };

    // Act
    const auto sdf = PlanarDistanceFieldGenerator::Generate(curveAdapter, sdfSettings);
    EXPORT_FIELD_AND_CURVE(sdf, curveAdapter.Get(), "SimpleClosedManifoldCurve");

    // Assert
    ASSERT_TRUE(sdf.IsValid());
    EXPECT_EQ(HashScalarGrid(sdf), HASH_square_SDF);
    ASSERT_EQ(sdf.Dimensions().Nx, 31);
    ASSERT_EQ(sdf.Dimensions().Ny, 31);
    auto index = [&sdf](size_t i, size_t j) {
        return i + j * sdf.Dimensions().Nx;
    };
    EXPECT_TRUE(sdf.Values()[index(14, 14)] < 0.0);
    EXPECT_NEAR(sdf.Values()[index(1, 1)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(0, 29)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 29)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 0)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(0, 14)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(15, 0)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 14)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(14, 29)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(14, 14)], -0.5, cellSize);
}

TEST(DistanceField2DTests, PlanarDistanceFieldGenerator_SimpleOpenManifoldCurve)
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
    const Scalar minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
    const Scalar cellSize = minSize / 10.0;

    const DistanceField2DSettings sdfSettings{
        cellSize,
        1.0,
        DBL_MAX,
        KDTreeSplitType::Center,
        SignComputation2D::None,
        PreprocessingType2D::Quadtree
    };

    // Act
    const auto sdf = PlanarDistanceFieldGenerator::Generate(curveAdapter, sdfSettings);
    EXPORT_FIELD_AND_CURVE(sdf, curveAdapter.Get(), "SimpleOpenManifoldCurve");

    // Assert
    ASSERT_TRUE(sdf.IsValid());
    EXPECT_EQ(HashScalarGrid(sdf), HASH_openSquare_SDF);
    ASSERT_EQ(sdf.Dimensions().Nx, 31);
    ASSERT_EQ(sdf.Dimensions().Ny, 31);
    auto index = [&sdf](size_t i, size_t j) {
        return i + j * sdf.Dimensions().Nx;
    };
    EXPECT_FALSE(sdf.Values()[index(14, 14)] < 0.0);
    EXPECT_NEAR(sdf.Values()[index(1, 1)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(0, 29)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 29)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 0)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(14, 0)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 14)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(14, 29)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(14, 14)], 0.5, cellSize);
}

TEST(DistanceField2DTests, PlanarPointCloudDistanceFieldGenerator_SimplePointCloud)
{
    // Arrange
    const std::vector points = { Point2(0, 0), Point2(1, 0), Point2(1, 1), Point2(0, 1) };
    const auto pointBBox = BoundingBox2(points);
    const auto pointBBoxSize = pointBBox.max() - pointBBox.min();
    const Scalar minSize = std::min(pointBBoxSize[0], pointBBoxSize[1]);
    const Scalar cellSize = minSize / 10.0;

    const PointCloudDistanceField2DSettings sdfSettings{
        cellSize,
        1.0,
        DBL_MAX
    };

    // Act
    const auto sdf = PlanarPointCloudDistanceFieldGenerator::Generate(points, sdfSettings);
    EXPORT_FIELD_AND_PT_CLOUD(sdf, points, "SimplePointCloud");

    // Assert
    ASSERT_TRUE(sdf.IsValid());
    EXPECT_EQ(HashScalarGrid(sdf), HASH_squarePts_SDF);
    ASSERT_EQ(sdf.Dimensions().Nx, 31);
    ASSERT_EQ(sdf.Dimensions().Ny, 31);
    auto index = [&sdf](size_t i, size_t j) {
        return i + j * sdf.Dimensions().Nx;
    };
    EXPECT_FALSE(sdf.Values()[index(14, 14)] < 0.0);
    EXPECT_NEAR(sdf.Values()[index(0, 0)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(0, 29)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 29)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 0)], sqrt(2.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(0, 14)], sqrt(5.0) / 2.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(14, 0)], sqrt(5.0) / 2.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 14)], sqrt(5.0) / 2.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(14, 29)], sqrt(5.0) / 2.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(29, 29)], sqrt(2.0), cellSize);
}

TEST(DistanceField2DTests, PlanarDistanceFieldGenerator_MoreComplexPointCloud)
{
    // Arrange
    std::vector points{
        Point2{39.507142, 14.544772},
        Point2{46.104261, 5.542495},
        Point2{61.36143, 4.906308},
        Point2{68.282948, 13.11281},
        Point2{66.153916, 31.426095},
        Point2{69.924933, 39.754365},
        Point2{111.082270, 39.723965},
        Point2{117.795930, 46.528945},
        Point2{117.765430, 66.419005},
        Point2{113.358230, 72.215125},
        Point2{89.030514, 71.743755},
        Point2{82.788235, 77.337235},
        Point2{87.565954, 122.613040},
        Point2{80.332640, 129.222760},
        Point2{68.372952, 128.366110},
        Point2{61.451434, 122.392570},
        Point2{57.089489, 85.394595},
        Point2{36.929786, 84.297475},
        Point2{10.835265, 83.544745},
        Point2{3.558908, 76.519305},
        Point2{3.558908, 57.450225},
        Point2{10.584357, 47.664785},
        Point2{32.413427, 48.166595},
        Point2{40.865615, 41.985195}
    };
    const auto sampledCurve = CurveFactory::sampled_polygon(points, 50, false);
    points = sampledCurve.positions();
    const auto pointBBox = BoundingBox2(points);
    const auto pointBBoxSize = pointBBox.max() - pointBBox.min();
    const Scalar minSize = std::min(pointBBoxSize[0], pointBBoxSize[1]);
    const Scalar cellSize = minSize / 10.0;

    const PointCloudDistanceField2DSettings sdfSettings{
        cellSize,
        1.0,
        DBL_MAX
    };

    // Act
    const auto sdf = PlanarPointCloudDistanceFieldGenerator::Generate(points, sdfSettings);
    EXPORT_FIELD_AND_PT_CLOUD(sdf, points, "MoreComplexPointCloud");

    // Assert
    ASSERT_TRUE(sdf.IsValid());
}

TEST(DistanceField2DTests, PlanarDistanceFieldGenerator_ManifoldCircleCurve)
{
    // Arrange
    //const ManifoldCurve2DAdapter curveAdapter(std::make_shared<ManifoldCurve2D>(CurveFactory::circle(Point2{ 53.669357, 34.419353 }, 13.217741, 10)));
    const ManifoldCurve2DAdapter curveAdapter(std::make_shared<ManifoldCurve2D>(CurveFactory::circle(Point2{ 0, 0 }, 13.217741, 10)));
	const auto curveBBox = curveAdapter.GetBounds();
    const auto curveBBoxSize = curveBBox.max() - curveBBox.min();
    const Scalar minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
    const Scalar cellSize = minSize / 10.0;

    const DistanceField2DSettings sdfSettings{
        cellSize,
        1.0,
        DBL_MAX,
        KDTreeSplitType::Center,
        SignComputation2D::None,
        PreprocessingType2D::Quadtree
    };

    // Act
    const auto sdf = PlanarDistanceFieldGenerator::Generate(curveAdapter, sdfSettings);
    EXPORT_FIELD_AND_CURVE(sdf, curveAdapter.Get(), "ManifoldCircleCurve");

    // Assert
    ASSERT_TRUE(sdf.IsValid());
}

TEST(DistanceField2DTests, PlanarDistanceFieldGenerator_CirclePointCloud)
{
    // Arrange
    const auto sampledCurve = CurveFactory::circle(Point2{ 0, 0 }, 13.217741, 10);
    const auto points = sampledCurve.positions();
    const auto pointBBox = BoundingBox2(points);
    const auto pointBBoxSize = pointBBox.max() - pointBBox.min();
    const Scalar minSize = std::min(pointBBoxSize[0], pointBBoxSize[1]);
    const Scalar cellSize = minSize / 10.0;

    const PointCloudDistanceField2DSettings sdfSettings{
        cellSize,
        1.0,
        DBL_MAX
    };

    // Act
    const auto sdf = PlanarPointCloudDistanceFieldGenerator::Generate(points, sdfSettings);
    EXPORT_FIELD_AND_PT_CLOUD(sdf, points, "CirclePointCloud");

    // Assert
    ASSERT_TRUE(sdf.IsValid());
}
TEST(DistanceField2DTests, PlanarDistanceFieldGenerator_MoreComplexPolygonalManifoldCurve)
{
    // Arrange
    const std::vector points{
        Point2{39.507142, 14.544772},
        Point2{46.104261, 5.542495},
        Point2{61.36143, 4.906308},
        Point2{68.282948, 13.11281},
        Point2{66.153916, 31.426095},
        Point2{69.924933, 39.754365},
        Point2{111.082270, 39.723965},
        Point2{117.795930, 46.528945},
        Point2{117.765430, 66.419005},
        Point2{113.358230, 72.215125},
        Point2{89.030514, 71.743755},
        Point2{82.788235, 77.337235},
        Point2{87.565954, 122.613040},
        Point2{80.332640, 129.222760},
        Point2{68.372952, 128.366110},
        Point2{61.451434, 122.392570},
        Point2{57.089489, 85.394595},
        Point2{36.929786, 84.297475},
        Point2{10.835265, 83.544745},
        Point2{3.558908, 76.519305},
        Point2{3.558908, 57.450225},
        Point2{10.584357, 47.664785},
        Point2{32.413427, 48.166595},
        Point2{40.865615, 41.985195}
    };
    const ManifoldCurve2DAdapter curveAdapter(std::make_shared<ManifoldCurve2D>(CurveFactory::sampled_polygon(points, 50, false)));
    const auto curveBBox = curveAdapter.GetBounds();
    const auto curveBBoxSize = curveBBox.max() - curveBBox.min();
    const Scalar minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
    const Scalar cellSize = minSize / 10.0;

    const DistanceField2DSettings sdfSettings{
        cellSize,
        1.0,
        DBL_MAX,
        KDTreeSplitType::Center,
        SignComputation2D::None,
        PreprocessingType2D::Quadtree
    };

    // Act
    const auto sdf = PlanarDistanceFieldGenerator::Generate(curveAdapter, sdfSettings);
    EXPORT_FIELD_AND_CURVE(sdf, curveAdapter.Get(), "MoreComplexPolygonalManifoldCurve");

    // Assert
    ASSERT_TRUE(sdf.IsValid());
}
