
#include "gtest/gtest.h"

#include "geometry/Grid.h"
#include "geometry/GridUtil.h"

using namespace Geometry;
using namespace pmp;

constexpr pmp::Scalar BOX_EPSILON = 1e-4;

TEST(GridInterpolationTest, BilinearInterpolateScalarValue2D_case)
{
    // Setup for 2D ScalarGrid
    ScalarGrid2D grid2D(1.0, BoundingBox2(
        vec2(0.0 + BOX_EPSILON, 0.0 + BOX_EPSILON),
        vec2(2.0 - BOX_EPSILON, 2.0 - BOX_EPSILON)));
    grid2D.Values() = {
        0, 1, 2,
        1, 2, 3,
        2, 3, 4 };
    const auto [Nx, Ny] = grid2D.Dimensions();
    ASSERT_EQ(Nx * Ny, grid2D.Values().size());

    const vec2 samplePt1(0.5, 0.5);
    const vec2 samplePt2(1.5, 0.5);
    const vec2 samplePt3(0.5, 1.5);
    const vec2 samplePt4(1.5, 1.5);

    const double result1 = BilinearInterpolateScalarValue(samplePt1, grid2D);
    const double result2 = BilinearInterpolateScalarValue(samplePt2, grid2D);
    const double result3 = BilinearInterpolateScalarValue(samplePt3, grid2D);
    const double result4 = BilinearInterpolateScalarValue(samplePt4, grid2D);

    EXPECT_NEAR(result1, 1.0, 1e-6);
    EXPECT_NEAR(result2, 2.0, 1e-6);
    EXPECT_NEAR(result3, 2.0, 1e-6);
    EXPECT_NEAR(result4, 3.0, 1e-6);
}

TEST(GridInterpolationTest, BilinearInterpolateVectorValue2D_case)
{
    // Setup for 2D VectorGrid
    VectorGrid2D vectorGrid2D(1.0, BoundingBox2(
        vec2(0.0 + BOX_EPSILON, 0.0 + BOX_EPSILON),
        vec2(2.0 - BOX_EPSILON, 2.0 - BOX_EPSILON)));
    vectorGrid2D.ValuesX() = {
        0, 1, 2,
        1, 2, 3,
        2, 3, 4 };
    vectorGrid2D.ValuesY() = {
        6, 5, 4,
        5, 4, 3,
        4, 3, 2 };
    const auto [Nx, Ny] = vectorGrid2D.Dimensions();
    ASSERT_EQ(Nx * Ny, vectorGrid2D.ValuesX().size());
    ASSERT_EQ(Nx * Ny, vectorGrid2D.ValuesY().size());
    const vec2 samplePt1(0.5, 0.5);
    const vec2 samplePt2(1.5, 0.5);
    const vec2 samplePt3(0.5, 1.5);
    const vec2 samplePt4(1.5, 1.5);

    const dvec2 result1 = BilinearInterpolateVectorValue(samplePt1, vectorGrid2D);
    const dvec2 result2 = BilinearInterpolateVectorValue(samplePt2, vectorGrid2D);
    const dvec2 result3 = BilinearInterpolateVectorValue(samplePt3, vectorGrid2D);
    const dvec2 result4 = BilinearInterpolateVectorValue(samplePt4, vectorGrid2D);

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
    ScalarGrid grid3D(1.0, BoundingBox(
        vec3(0.0 + BOX_EPSILON, 0.0 + BOX_EPSILON, 0.0 + BOX_EPSILON),
        vec3(2.0 - BOX_EPSILON, 2.0 - BOX_EPSILON, 2.0 - BOX_EPSILON)));
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
    const auto [Nx, Ny, Nz] = grid3D.Dimensions();
    ASSERT_EQ(Nx * Ny * Nz, grid3D.Values().size());

    const vec3 samplePt1(0.5, 0.5, 0.5);
    const vec3 samplePt2(1.5, 0.5, 0.5);
    const vec3 samplePt3(0.5, 1.5, 0.5);
    const vec3 samplePt4(0.5, 0.5, 1.5);

    const double result1 = TrilinearInterpolateScalarValue(samplePt1, grid3D);
    const double result2 = TrilinearInterpolateScalarValue(samplePt2, grid3D);
    const double result3 = TrilinearInterpolateScalarValue(samplePt3, grid3D);
    const double result4 = TrilinearInterpolateScalarValue(samplePt4, grid3D);

    EXPECT_NEAR(result1, 1.5, 1e-6);
    EXPECT_NEAR(result2, 2.5, 1e-6);
    EXPECT_NEAR(result3, 2.5, 1e-6);
    EXPECT_NEAR(result4, 2.5, 1e-6);
}

TEST(GridInterpolationTest, TrilinearInterpolateVectorValue_case)
{
    // Setup for 3D VectorGrid
    VectorGrid vectorGrid3D(1.0, BoundingBox(
        vec3(0.0 + BOX_EPSILON, 0.0 + BOX_EPSILON, 0.0 + BOX_EPSILON),
        vec3(2.0 - BOX_EPSILON, 2.0 - BOX_EPSILON, 2.0 - BOX_EPSILON)));
    vectorGrid3D.ValuesX() = { 0, 1, 2, 1, 2, 3, 2, 3, 4, 1, 2, 3, 2, 3, 4, 3, 4, 5, 2, 3, 4, 3, 4, 5, 4, 5, 6 };
    vectorGrid3D.ValuesY() = { 6, 5, 4, 5, 4, 3, 4, 3, 2, 5, 4, 3, 4, 3, 2, 3, 2, 1, 4, 3, 2, 3, 2, 1, 2, 1, 0 };
    vectorGrid3D.ValuesZ() = { 2, 3, 4, 3, 4, 5, 4, 5, 6, 3, 4, 5, 4, 5, 6, 5, 6, 7, 4, 5, 6, 5, 6, 7, 6, 7, 8 };
    const auto [Nx, Ny, Nz] = vectorGrid3D.Dimensions();
    ASSERT_EQ(Nx * Ny * Nz, vectorGrid3D.ValuesX().size());
    ASSERT_EQ(Nx * Ny * Nz, vectorGrid3D.ValuesY().size());
    ASSERT_EQ(Nx * Ny * Nz, vectorGrid3D.ValuesZ().size());
    const vec3 samplePt1(0.5, 0.5, 0.5);
    const vec3 samplePt2(1.5, 0.5, 0.5);
    const vec3 samplePt3(0.5, 1.5, 0.5);
    const vec3 samplePt4(0.5, 0.5, 1.5);

    const dvec3 result1 = TrilinearInterpolateVectorValue(samplePt1, vectorGrid3D);
    const dvec3 result2 = TrilinearInterpolateVectorValue(samplePt2, vectorGrid3D);
    const dvec3 result3 = TrilinearInterpolateVectorValue(samplePt3, vectorGrid3D);
    const dvec3 result4 = TrilinearInterpolateVectorValue(samplePt4, vectorGrid3D);

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
    ScalarGrid grid3D(1.0, BoundingBox(
        vec3(0.0 + BOX_EPSILON, 0.0 + BOX_EPSILON, 0.0 + BOX_EPSILON),
        vec3(2.0 - BOX_EPSILON, 2.0 - BOX_EPSILON, 2.0 - BOX_EPSILON)));
    grid3D.Values() = { 0, 1, 2, 1, 2, 3, 2, 3, 4, 1, 2, 3, 2, 3, 4, 3, 4, 5, 2, 3, 4, 3, 4, 5, 4, 5, 6 };
    const auto [Nx, Ny, Nz] = grid3D.Dimensions();
    ASSERT_EQ(Nx * Ny * Nz, grid3D.Values().size());

    vec3 samplePt(1.4, 0.6, 0.6);
    double result = GetNearestNeighborScalarValue(samplePt, grid3D);
    EXPECT_NEAR(result, 3.0, 1e-6);
}

TEST(GridInterpolationTest, GetNearestNeighborVectorValue_case)
{
    // Setup for 3D VectorGrid
    VectorGrid vectorGrid3D(1.0, BoundingBox(
        vec3(0.0 + BOX_EPSILON, 0.0 + BOX_EPSILON, 0.0 + BOX_EPSILON),
        vec3(2.0 - BOX_EPSILON, 2.0 - BOX_EPSILON, 2.0 - BOX_EPSILON)));
    vectorGrid3D.ValuesX() = { 0, 1, 2, 1, 2, 3, 2, 3, 4, 1, 2, 3, 2, 3, 4, 3, 4, 5, 2, 3, 4, 3, 4, 5, 4, 5, 6 };
    vectorGrid3D.ValuesY() = { 6, 5, 4, 5, 4, 3, 4, 3, 2, 5, 4, 3, 4, 3, 2, 3, 2, 1, 4, 3, 2, 3, 2, 1, 2, 1, 0 };
    vectorGrid3D.ValuesZ() = { 2, 3, 4, 3, 4, 5, 4, 5, 6, 3, 4, 5, 4, 5, 6, 5, 6, 7, 4, 5, 6, 5, 6, 7, 6, 7, 8 };
    const auto [Nx, Ny, Nz] = vectorGrid3D.Dimensions();
    ASSERT_EQ(Nx * Ny * Nz, vectorGrid3D.ValuesX().size());
    ASSERT_EQ(Nx * Ny * Nz, vectorGrid3D.ValuesY().size());
    ASSERT_EQ(Nx * Ny * Nz, vectorGrid3D.ValuesZ().size());

    vec3 samplePt(1.4, 0.6, 0.4);
    dvec3 result = GetNearestNeighborVectorValue(samplePt, vectorGrid3D);
    EXPECT_NEAR(result[0], 2.0, 1e-6);
    EXPECT_NEAR(result[1], 4.0, 1e-6);
    EXPECT_NEAR(result[2], 4.0, 1e-6);
}

TEST(GridInterpolationTest, GetNearestNeighborScalarValue2D_case)
{
    // Setup for 2D ScalarGrid
    ScalarGrid2D grid2D(1.0, BoundingBox2(
        vec2(0.0 + BOX_EPSILON, 0.0 + BOX_EPSILON),
        vec2(2.0 - BOX_EPSILON, 2.0 - BOX_EPSILON)));
    grid2D.Values() = {
        0, 1, 2,
        1, 2, 3,
        2, 3, 4 };
    const auto [Nx, Ny] = grid2D.Dimensions();
    ASSERT_EQ(Nx * Ny, grid2D.Values().size());

    vec2 samplePt(1.6, 1.7);
    double result = GetNearestNeighborScalarValue2D(samplePt, grid2D);
    EXPECT_NEAR(result, 4.0, 1e-6);
}

TEST(GridInterpolationTest, GetNearestNeighborVectorValue2D_case)
{
    // Setup for 2D VectorGrid
    VectorGrid2D vectorGrid2D(1.0, BoundingBox2(
        vec2(0.0 + BOX_EPSILON, 0.0 + BOX_EPSILON),
        vec2(2.0 - BOX_EPSILON, 2.0 - BOX_EPSILON)));
    vectorGrid2D.ValuesX() = { 0, 1, 2, 1, 2, 3, 2, 3, 4 };
    vectorGrid2D.ValuesY() = { 6, 5, 4, 5, 4, 3, 4, 3, 2 };
    const auto [Nx, Ny] = vectorGrid2D.Dimensions();
    ASSERT_EQ(Nx * Ny, vectorGrid2D.ValuesX().size());
    ASSERT_EQ(Nx * Ny, vectorGrid2D.ValuesY().size());

    vec2 samplePt(1.4, 0.6);
    dvec2 result = GetNearestNeighborVectorValue2D(samplePt, vectorGrid2D);
    EXPECT_NEAR(result[0], 2.0, 1e-6);
    EXPECT_NEAR(result[1], 4.0, 1e-6);
}

// Local maxima

TEST(ScalarGridLocalMaximumTests2D, DFSlice3x3ContainingLocalMax)
{
	// Arrange
    ScalarGrid2D grid(1.0, BoundingBox2(
        Point2(-1 + BOX_EPSILON, -1 + BOX_EPSILON),
        Point2(1 - BOX_EPSILON, 1 - BOX_EPSILON)));
    grid.Values() = {
        0.852389336, 0.899750948, 0.885099947,
        0.899742365, 0.957528651 /* this one's a discrete max */, 0.952537298,
        0.900355637, 0.952537298, 0.955666423
    };
    const auto [Nx, Ny] = grid.Dimensions();
    ASSERT_EQ(Nx * Ny, grid.Values().size());
    const auto discreteMaxIter = std::ranges::max_element(grid.Values());
    const auto discreteMaxGridIndex = std::distance(grid.Values().begin(), discreteMaxIter);
    const unsigned int ix = discreteMaxGridIndex % Nx;
    const unsigned int iy = (discreteMaxGridIndex / Nx) % Ny;
    const auto epsilon = 0.7 * grid.CellSize();

    // Act
    const auto localMax = FindLocalMaximumNearScalarGridCell(grid, ix, iy, 1);

    // Assert
    ASSERT_TRUE(localMax.has_value());
    EXPECT_NEAR((*localMax)[0], 0.0, epsilon);
    EXPECT_NEAR((*localMax)[1], 0.0, epsilon);
}

TEST(ScalarGridLocalMaximumTests3D, MaximumWithinCellPointsWithUnitRadius)
{
    // Arrange
    ScalarGrid grid(1.0, BoundingBox(
        Point(-0.5, -0.5, -0.5),
        Point(0.5, 0.5, 0.5)));
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
    const auto [Nx, Ny, Nz] = grid.Dimensions();
    ASSERT_EQ(Nx * Ny * Nz, grid.Values().size());
    const unsigned int ix = 1, iy = 1, iz = 1;

    // Act
    const auto localMax = FindLocalMaximumNearScalarGridCell(grid, ix, iy, iz, 1);

    // Assert
    ASSERT_TRUE(localMax.has_value());
    EXPECT_NEAR((*localMax)[0], 0.0, 0.5);
    EXPECT_NEAR((*localMax)[1], 0.0, 0.5);
    EXPECT_NEAR((*localMax)[2], 0.0, 0.5);
}



TEST(ScalarGridLocalMaximumTests3D, MaximumOutsideOfCellPointsWithUnitRadius)
{
    // Arrange
    ScalarGrid grid(1.0, BoundingBox(
        Point(-0.5, -0.5, -0.5),
        Point(0.5, 0.5, 0.5)));
    grid.Values() = {
        1, 2, 1,
        2, 3, 2,
        1, 2, 1,

        2, 2, 2,
        3, 2, 3,
        2, 3, 2,

        2, 3, 3,
        3, 3, 3,
        2, 3, 2
    };
    const auto [Nx, Ny, Nz] = grid.Dimensions();
    ASSERT_EQ(Nx * Ny * Nz, grid.Values().size());
    const unsigned int ix = 1, iy = 1, iz = 1;

    // Act
    const auto localMax = FindLocalMaximumNearScalarGridCell(grid, ix, iy, iz, 1);

    // Assert
    ASSERT_FALSE(localMax.has_value());
}

TEST(ScalarGridLocalMaximumTests3D, DistanceFieldMaximumOnA4x4x4Grid)
{
    // Arrange
    ScalarGrid grid(1.0,
        BoundingBox(Point(-1 + BOX_EPSILON, -1 + BOX_EPSILON, -1 + BOX_EPSILON),
            Point(2 - BOX_EPSILON, 2 - BOX_EPSILON, 2 - BOX_EPSILON)));
    grid.Values() = {
		0.81023629286903087, 0.85238934794879040, 0.83914933112888346, 0.80354128818732840,
        0.85239291668881567, 0.89974235459188201, 0.88511633497951381, 0.84865449608501553,
        0.83913910333340569, 0.90035563805955754, 0.89886805137362868, 0.84998801142285851,
        0.80352713017607513, 0.85438082927200598, 0.85216718931308788, 0.81597085291517268,
        0.85241283968808357, 0.89975094099893516, 0.90036883229000197, 0.85439444458819092,
        0.89979855218883731, 0.95752867446787882 /* this one's a discrete max */, 0.95256121730760279, 0.90031943442697238,
        0.88509994314652107, 0.95253730050092644, 0.95567193715342036, 0.89611810084265564,
        0.84862831557833784, 0.90029482979158448, 0.89886538570375674, 0.86223548739688782,
        0.83913910333340569, 0.88509994358805388, 0.89886236419191556, 0.85215848859126753,
        0.90035561966136124, 0.95253728894960732, 0.95566664730415296, 0.89885937553849304,
        0.89886214466203296, 0.95566642738693652, 0.94619437949884622, 0.88845914342289634,
        0.84998556810887471, 0.89611624156026271, 0.88845911192517624, 0.85689880931903240,
        0.80352713017607513, 0.84862831731183797, 0.84998591379318611, 0.81596602250509498,
        0.85438078453126787, 0.90029480959816899, 0.89611649369699953, 0.86223193940215770,
        0.85215845189584449, 0.89885931832455379, 0.88845925686474325, 0.85689880924916606,
        0.81596592301006066, 0.86223180256362797, 0.85689880912303662, 0.82885106760495830
    };
    const auto [Nx, Ny, Nz] = grid.Dimensions();
    ASSERT_EQ(Nx * Ny * Nz, grid.Values().size());
    const auto discreteMaxIter = std::ranges::max_element(grid.Values());
    const auto discreteMaxGridIndex = std::distance(grid.Values().begin(), discreteMaxIter);
    const unsigned int ix = discreteMaxGridIndex % Nx;
    const unsigned int iy = (discreteMaxGridIndex / Nx) % Ny;
    const unsigned int iz = discreteMaxGridIndex / (Nx * Ny);

    // Act
    const auto localMax = FindLocalMaximumNearScalarGridCell(grid, ix, iy, iz, 1);

    // Assert
    ASSERT_TRUE(localMax.has_value());
}

TEST(ScalarGridLocalMaximumTests3D, DistanceFieldMaximumOnAnother4x4x4Grid)
{
    // Arrange
    ScalarGrid grid(1.0, 
        BoundingBox(Point(-1 + BOX_EPSILON, -1 + BOX_EPSILON, -1 + BOX_EPSILON),
            Point(2 - BOX_EPSILON, 2 - BOX_EPSILON, 2 - BOX_EPSILON)));
    grid.Values() = {
        0.75249668594631991, 0.81023629286903087, 0.85238934794879040, 0.83914933112888346,
        0.78400121110930876, 0.85239291668881567, 0.89974235459188201, 0.88511633497951381,
        0.76232971046297848, 0.83913910333340569, 0.90035563805955754, 0.89886805137362868,
        0.74224524396086133, 0.80352713017607513, 0.85438082927200598, 0.85216718931308788,
        0.79286617148826111, 0.85241283968808357, 0.89975094099893516, 0.90036883229000197,
        0.82550094061990265, 0.89979855218883731, 0.95752867446787882  /* this one's a discrete max */, 0.95256121730760279,
        0.80412204249388619, 0.88509994314652107, 0.95253730050092644, 0.95567193715342036,
        0.78551562098155625, 0.84862831557833784, 0.90029482979158448, 0.89886538570375674,
        0.78592159359624147, 0.83913910333340569, 0.88509994358805388, 0.89886236419191556,
        0.83602171358270672, 0.90035561966136124, 0.95253728894960732, 0.95566664730415296,
        0.82517466842385556, 0.89886214466203296, 0.95566642738693652, 0.94619437949884622,
        0.79122078044506561, 0.84998556810887471, 0.89611624156026271, 0.88845911192517624,
        0.74718598834219774, 0.80352713017607513, 0.84862831731183797, 0.84998591379318611,
        0.79396975654794855, 0.85438078453126787, 0.90029480959816899, 0.89611649369699953,
        0.78280678616200861, 0.85215845189584449, 0.89885931832455379, 0.88845925686474325,
        0.75873701617371925, 0.81596592301006066, 0.86223180256362797, 0.85689880912303662
    };
    const auto [Nx, Ny, Nz] = grid.Dimensions();
    ASSERT_EQ(Nx * Ny * Nz, grid.Values().size());
    const auto discreteMaxIter = std::ranges::max_element(grid.Values());
    const auto discreteMaxGridIndex = std::distance(grid.Values().begin(), discreteMaxIter);
    const unsigned int ix = discreteMaxGridIndex % Nx;
    const unsigned int iy = (discreteMaxGridIndex / Nx) % Ny;
    const unsigned int iz = discreteMaxGridIndex / (Nx * Ny);

    // Act
    const auto localMax = FindLocalMaximumNearScalarGridCell(grid, ix, iy, iz, 1);

    // Assert
    ASSERT_TRUE(localMax.has_value());
    //EXPECT_NEAR((*localMax)[0], /* expected x-coordinate */, 1e-5);
    //EXPECT_NEAR((*localMax)[1], /* expected y-coordinate */, 1e-5);
    //EXPECT_NEAR((*localMax)[2], /* expected z-coordinate */, 1e-5);
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

//#include <fstream>
//
//namespace
//{
//    void ExportToVTI_TODO_REMOVE(const std::string& filename, const Geometry::ScalarGrid& scalarGrid)
//    {
//        if (!scalarGrid.IsValid())
//            throw std::invalid_argument("ExportToVTI: scalarGrid to be exported is invalid!\n");
//        if (filename.empty())
//            throw std::invalid_argument("ExportToVTI: filename cannot be empty!\n");
//
//        std::fstream vti(filename + ".vti", std::fstream::out);
//
//        const auto& dims = scalarGrid.Dimensions();
//        const auto nx = static_cast<unsigned int>(dims.Nx);
//        const auto ny = static_cast<unsigned int>(dims.Ny);
//        const auto nz = static_cast<unsigned int>(dims.Nz);
//        const pmp::Scalar dx = scalarGrid.CellSize();
//
//        const pmp::vec3 min = scalarGrid.Box().min();
//        //const pmp::vec3 max = min + pmp::vec3(static_cast<float>(nx), static_cast<float>(ny),static_cast<float>(nz)) * dx;
//        const pmp::vec3 max = scalarGrid.Box().max();
//
//        vti << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
//        vti << "	<ImageData WholeExtent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1 << "\" Origin=\"" << min[0] + 0.5 * dx << " " << min[1] + 0.5 * dx << " " << min[2] + 0.5 * dx << "\" Spacing=\"" << dx << " " << dx << " " << dx << "\">\n";
//        vti << "		<Piece Extent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1 << "\">\n";
//        vti << "			<PointData Scalars=\"Scalars_\">\n";
//        vti << "				<DataArray type=\"Float32\" Name=\"Scalars_\" format=\"ascii\" RangeMin=\"" << min << "\" RangeMax=\"" << max << "\">\n";
//
//        for (const auto& val : scalarGrid.Values()) {
//            vti << static_cast<float>(val) << "\n";
//        }
//
//        vti << "				</DataArray>\n";
//        vti << "			</PointData>\n";
//        vti << "		<CellData>\n";
//        vti << "		</CellData>\n";
//        vti << "	</Piece>\n";
//        vti << "	</ImageData>\n";
//        vti << "</VTKFile>\n";
//
//        vti.close();
//    }
//}
//
//TEST(ScalarGridApplyPrimitivesTests, ApplyHyperboloidDistanceField)
//{
//	// Arrange
//    constexpr pmp::Scalar roiHalfDim = 5.0;
//    constexpr pmp::Scalar roiDim = 2.0 * roiHalfDim;
//    const auto box = BoundingBox{ vec3{0.0, 0.0, -roiHalfDim}, vec3{roiDim, roiDim, roiHalfDim} };
//    constexpr pmp::Scalar cellSize = 0.1;
//    constexpr double initVal = DEFAULT_SCALAR_GRID_INIT_VAL;
//    ScalarGrid result(cellSize, box, initVal);
//
//    const HyperboloidParams params{
//        box.center(),
//        1.0 / roiHalfDim, 1.5 / roiHalfDim, 2.0 / roiHalfDim,
//        vec3{1, 1, 1} * (roiHalfDim / 4.0),
//        DistanceUnion
//    };
//
//    // Act
//    ApplyHyperboloidDistanceFieldToGrid(result, params);
//
//    ExportToVTI_TODO_REMOVE("C:\\Users\\Martin\\source\\repos\\MCInversionPhDCollection\\output\\HyperboloidDF", result);
//
//    // Assert
//
//}