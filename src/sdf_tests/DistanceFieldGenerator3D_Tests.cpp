#include "gtest/gtest.h"

#include "sdf/SDF.h"

#include "pmp/BoundingBox.h"

#include "geometry/Grid.h"
#include "geometry/GridUtil.h"
#include "geometry/GeometryUtil.h"
#include "geometry/GeometryAdapters.h"
#include "geometry/GeometryConversionUtils.h"

using namespace SDF;
using namespace Geometry;
using namespace pmp;

// Hashes of the tested distance fields
constexpr size_t HASH_fullTetra_SDF{ 14803762804065774602 };
constexpr size_t HASH_tetraWithoutBottomFace_SDF{ 8980635575871867015 };
constexpr size_t HASH_tetraPoints_SDF{ 10659542913674947031 };

TEST(DistanceField3DTests, DistanceFieldGenerator_BaseWatertightTetrahedron)
{
    // Arrange
    const std::vector tetraVertices = { Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1) };
    BaseMeshGeometryData meshData;
    meshData.Vertices = tetraVertices;
    meshData.PolyIndices = { {0, 1, 3}, {1, 0, 2}, {0, 3, 2}, {1, 2, 3} };
    const BaseMeshAdapter meshAdapter(std::make_shared<BaseMeshGeometryData>(meshData));

    const auto meshBBox = meshAdapter.GetBounds();
    const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
    const float minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
    const float cellSize = minSize / 10.0f;

    const DistanceFieldSettings sdfSettings{
        cellSize,
        1.0f,
        DBL_MAX,
        KDTreeSplitType::Center,
        SignComputation::VoxelFloodFill,
        BlurPostprocessingType::None,
        PreprocessingType::Octree
    };

    // Act
    const auto sdf = DistanceFieldGenerator::Generate(meshAdapter, sdfSettings);

    // Assert
    ASSERT_TRUE(sdf.IsValid());
    EXPECT_EQ(HashScalarGrid(sdf), HASH_fullTetra_SDF);
    ASSERT_EQ(sdf.Dimensions().Nx, 30);
    ASSERT_EQ(sdf.Dimensions().Ny, 30);
    ASSERT_EQ(sdf.Dimensions().Nz, 30);
    auto index = [&sdf](size_t i, size_t j, size_t k) {
        return i + j * sdf.Dimensions().Nx + k * sdf.Dimensions().Nx * sdf.Dimensions().Ny;
    };
    EXPECT_TRUE(sdf.Values()[index(12, 12, 12)] < 0.0);
    EXPECT_NEAR(sdf.Values()[index(1, 1, 1)], sqrt(3.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(1, 12, 12)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 1, 12)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 12, 1)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 12, 0)], sqrt(GetDistanceToTriangleSq(std::vector{ tetraVertices[0], tetraVertices[1], tetraVertices[2] }, Point(0.25f, 0.25f, -1.0f))), cellSize);
    EXPECT_NEAR(sdf.Values()[index(0, 12, 12)], sqrt(GetDistanceToTriangleSq(std::vector{ tetraVertices[0], tetraVertices[3], tetraVertices[2] }, Point(-1.0f, 0.25f, 0.25f))), cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 0, 12)], sqrt(GetDistanceToTriangleSq(std::vector{ tetraVertices[0], tetraVertices[1], tetraVertices[3] }, Point(0.25f, -1.0f, 0.25f))), cellSize);
    EXPECT_NEAR(sdf.Values()[index(20, 20, 20)], sqrt(GetDistanceToTriangleSq(std::vector{ tetraVertices[1], tetraVertices[2], tetraVertices[3] }, Point(1.0f, 1.0f, 1.0f))), cellSize);
}

TEST(DistanceField3DTests, DistanceFieldGenerator_BaseOpenTetrahedron)
{
    // Arrange
    const std::vector tetraVertices = { Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1) };
    BaseMeshGeometryData meshData;
    meshData.Vertices = tetraVertices;
    meshData.PolyIndices = { {0, 1, 3}, {0, 3, 2}, {1, 2, 3} };
    const BaseMeshAdapter meshAdapter(std::make_shared<BaseMeshGeometryData>(meshData));

    const auto meshBBox = meshAdapter.GetBounds();
    const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
    const float minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
    const float cellSize = minSize / 10.0f;

    const DistanceFieldSettings sdfSettings{
        cellSize,
        1.0f,
        DBL_MAX,
        KDTreeSplitType::Center,
        SignComputation::None,
        BlurPostprocessingType::None,
        PreprocessingType::Octree
    };

    // Act
    const auto sdf = DistanceFieldGenerator::Generate(meshAdapter, sdfSettings);

    // Assert
    ASSERT_TRUE(sdf.IsValid());
    EXPECT_EQ(HashScalarGrid(sdf), HASH_tetraWithoutBottomFace_SDF);
    ASSERT_EQ(sdf.Dimensions().Nx, 30);
    ASSERT_EQ(sdf.Dimensions().Ny, 30);
    ASSERT_EQ(sdf.Dimensions().Nz, 30);
    auto index = [&sdf](size_t i, size_t j, size_t k) {
        return i + j * sdf.Dimensions().Nx + k * sdf.Dimensions().Nx * sdf.Dimensions().Ny;
    };
    EXPECT_FALSE(sdf.Values()[index(12, 12, 12)] < 0.0);
    EXPECT_NEAR(sdf.Values()[index(1, 1, 1)], sqrt(3.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(1, 12, 12)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 1, 12)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 12, 1)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(0, 12, 12)], sqrt(GetDistanceToTriangleSq(std::vector{ tetraVertices[0], tetraVertices[3], tetraVertices[2] }, Point(-1.0f, 0.25f, 0.25f))), cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 0, 12)], sqrt(GetDistanceToTriangleSq(std::vector{ tetraVertices[0], tetraVertices[1], tetraVertices[3] }, Point(0.25f, -1.0f, 0.25f))), cellSize);
    EXPECT_NEAR(sdf.Values()[index(20, 20, 20)], sqrt(GetDistanceToTriangleSq(std::vector{ tetraVertices[1], tetraVertices[2], tetraVertices[3] }, Point(1.0f, 1.0f, 1.0f))), cellSize);
}

TEST(DistanceField3DTests, DistanceFieldGenerator_PMPWatertightTetrahedron)
{
    // Arrange
    const std::vector tetraVertices = { Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1) };
    BaseMeshGeometryData meshData;
    meshData.Vertices = tetraVertices;
    meshData.PolyIndices = { {0, 1, 3}, {1, 0, 2}, {0, 3, 2}, {1, 2, 3} };
    const PMPSurfaceMeshAdapter meshAdapter(std::make_shared<SurfaceMesh>(ConvertBufferGeomToPMPSurfaceMesh(meshData)));

    const auto meshBBox = meshAdapter.GetBounds();
    const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
    const float minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
    const float cellSize = minSize / 10.0f;

    const DistanceFieldSettings sdfSettings{
        cellSize,
        1.0f,
        DBL_MAX,
        KDTreeSplitType::Center,
        SignComputation::VoxelFloodFill,
        BlurPostprocessingType::None,
        PreprocessingType::Octree
    };

    // Act
    const auto sdf = DistanceFieldGenerator::Generate(meshAdapter, sdfSettings);

    // Assert
    ASSERT_TRUE(sdf.IsValid());
    EXPECT_EQ(HashScalarGrid(sdf), HASH_fullTetra_SDF);
    ASSERT_EQ(sdf.Dimensions().Nx, 30);
    ASSERT_EQ(sdf.Dimensions().Ny, 30);
    ASSERT_EQ(sdf.Dimensions().Nz, 30);
    auto index = [&sdf](size_t i, size_t j, size_t k) {
        return i + j * sdf.Dimensions().Nx + k * sdf.Dimensions().Nx * sdf.Dimensions().Ny;
    };
    EXPECT_TRUE(sdf.Values()[index(12, 12, 12)] < 0.0);
    EXPECT_NEAR(sdf.Values()[index(1, 1, 1)], sqrt(3.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(1, 12, 12)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 1, 12)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 12, 1)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 12, 0)], sqrt(GetDistanceToTriangleSq(std::vector{ tetraVertices[0], tetraVertices[1], tetraVertices[2] }, Point(0.25f, 0.25f, -1.0f))), cellSize);
    EXPECT_NEAR(sdf.Values()[index(0, 12, 12)], sqrt(GetDistanceToTriangleSq(std::vector{ tetraVertices[0], tetraVertices[3], tetraVertices[2] }, Point(-1.0f, 0.25f, 0.25f))), cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 0, 12)], sqrt(GetDistanceToTriangleSq(std::vector{ tetraVertices[0], tetraVertices[1], tetraVertices[3] }, Point(0.25f, -1.0f, 0.25f))), cellSize);
    EXPECT_NEAR(sdf.Values()[index(20, 20, 20)], sqrt(GetDistanceToTriangleSq(std::vector{ tetraVertices[1], tetraVertices[2], tetraVertices[3] }, Point(1.0f, 1.0f, 1.0f))), cellSize);
}

TEST(DistanceField3DTests, DistanceFieldGenerator_PMPOpenTetrahedron)
{
    // Arrange
    const std::vector tetraVertices = { Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1) };
    BaseMeshGeometryData meshData;
    meshData.Vertices = tetraVertices;
    meshData.PolyIndices = { {0, 1, 3}, {0, 3, 2}, {1, 2, 3} };
    const PMPSurfaceMeshAdapter meshAdapter(std::make_shared<SurfaceMesh>(ConvertBufferGeomToPMPSurfaceMesh(meshData)));

    const auto meshBBox = meshAdapter.GetBounds();
    const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
    const float minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
    const float cellSize = minSize / 10.0f;

    const DistanceFieldSettings sdfSettings{
        cellSize,
        1.0f,
        DBL_MAX,
        KDTreeSplitType::Center,
        SignComputation::None,
        BlurPostprocessingType::None,
        PreprocessingType::Octree
    };

    // Act
    const auto sdf = DistanceFieldGenerator::Generate(meshAdapter, sdfSettings);

    // Assert
    ASSERT_TRUE(sdf.IsValid());
    EXPECT_EQ(HashScalarGrid(sdf), HASH_tetraWithoutBottomFace_SDF);
    ASSERT_EQ(sdf.Dimensions().Nx, 30);
    ASSERT_EQ(sdf.Dimensions().Ny, 30);
    ASSERT_EQ(sdf.Dimensions().Nz, 30);
    auto index = [&sdf](size_t i, size_t j, size_t k) {
        return i + j * sdf.Dimensions().Nx + k * sdf.Dimensions().Nx * sdf.Dimensions().Ny;
    };
    EXPECT_FALSE(sdf.Values()[index(12, 12, 12)] < 0.0);
    EXPECT_NEAR(sdf.Values()[index(1, 1, 1)], sqrt(3.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(1, 12, 12)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 1, 12)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 12, 1)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(0, 12, 12)], sqrt(GetDistanceToTriangleSq(std::vector{ tetraVertices[0], tetraVertices[3], tetraVertices[2] }, Point(-1.0f, 0.25f, 0.25f))), cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 0, 12)], sqrt(GetDistanceToTriangleSq(std::vector{ tetraVertices[0], tetraVertices[1], tetraVertices[3] }, Point(0.25f, -1.0f, 0.25f))), cellSize);
    EXPECT_NEAR(sdf.Values()[index(20, 20, 20)], sqrt(GetDistanceToTriangleSq(std::vector{ tetraVertices[1], tetraVertices[2], tetraVertices[3] }, Point(1.0f, 1.0f, 1.0f))), cellSize);
}

TEST(DistanceField3DTests, PointCloudDistanceFieldGenerator_TetrahedronPointCloud)
{
    // Arrange
    const std::vector points = { Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1) };
    const auto pointBBox = BoundingBox(points);
    const auto pointBBoxSize = pointBBox.max() - pointBBox.min();
    const float minSize = std::min({ pointBBoxSize[0], pointBBoxSize[1], pointBBoxSize[2] });
    const float cellSize = minSize / 10.0f;

    const PointCloudDistanceFieldSettings sdfSettings{
        cellSize,
        1.0f,
        DBL_MAX,
        BlurPostprocessingType::None
    };

    // Act
    const auto sdf = PointCloudDistanceFieldGenerator::Generate(points, sdfSettings);

    // Assert
    ASSERT_TRUE(sdf.IsValid());
    EXPECT_EQ(HashScalarGrid(sdf), HASH_tetraPoints_SDF);
    ASSERT_EQ(sdf.Dimensions().Nx, 30);
    ASSERT_EQ(sdf.Dimensions().Ny, 30);
    ASSERT_EQ(sdf.Dimensions().Nz, 30);
    auto index = [&sdf](size_t i, size_t j, size_t k) {
        return i + j * sdf.Dimensions().Nx + k * sdf.Dimensions().Nx * sdf.Dimensions().Ny;
    };
    EXPECT_FALSE(sdf.Values()[index(12, 12, 12)] < 0.0);
    EXPECT_NEAR(sdf.Values()[index(1, 1, 1)], sqrt(3.0), cellSize);
    EXPECT_NEAR(sdf.Values()[index(1, 12, 12)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 1, 12)], 1.0, cellSize);
    EXPECT_NEAR(sdf.Values()[index(12, 12, 1)], 1.0, cellSize);
}