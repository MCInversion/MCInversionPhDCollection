
#include "gtest/gtest.h"

#include "pmp/Types.h"
#include "pmp/algorithms/DifferentialGeometry.h"
#include "pmp/algorithms/Normals.h"

#include "geometry/GeometryConversionUtils.h"
#include "geometry/GridUtil.h"

#include "core/ConversionUtils.h"

#include <filesystem>

using namespace Geometry;
using namespace pmp;

// set up root directory
const std::filesystem::path fsRootPath = DROOT_DIR;
const auto fsDataDirPath = fsRootPath / "data\\";
const auto fsDataOutPath = fsRootPath / "output\\";
const std::string dataDirPath = fsDataDirPath.string();
const std::string dataOutPath = fsDataOutPath.string();

TEST(ApproximateMeshPatchWithQuadric, ApplyHyperboloidDistanceFieldFromMesh1Ring)
{
    // Arrange
    constexpr float roiHalfDim = 5.0f;
    constexpr float roiDim = 2.0f * roiHalfDim;
    const auto box = BoundingBox{ vec3{0.0f, 0.0f, -roiHalfDim}, vec3{roiDim, roiDim, roiHalfDim} };
    const auto fieldCenter = box.center();
    constexpr float cellSize = 0.1f;
    constexpr double initVal = DEFAULT_SCALAR_GRID_INIT_VAL;
    ScalarGrid grid(cellSize, box, initVal);

    BaseMeshGeometryData saddleBaseMesh;
    saddleBaseMesh.Vertices = std::vector{
        Point{0.035850f, 0.421721f, 0.525021f},
        Point{-0.065786f, 0.509711f, -0.487154f},
        Point{0.549778f, 0.034613f, 0.180541f},
        Point{-0.413034f, 0.047552f, 0.450408f},
        Point{-0.011141f, 0.240242f, 0.022254f}, // central vertex
        Point{0.383242f, 0.051138f, -0.414398f},
        Point{-0.579564f, -0.043537f, -0.144535f},
    };
    saddleBaseMesh.PolyIndices = std::vector{
        std::vector<unsigned int>{0, 2, 4},
        std::vector<unsigned int>{6, 4, 1},
        std::vector<unsigned int>{6, 3, 4},
        std::vector<unsigned int>{3, 0, 4},
        std::vector<unsigned int>{4, 5, 1},
        std::vector<unsigned int>{4, 2, 5}
    };
    const auto saddleSurfaceMesh = ConvertBufferGeomToPMPSurfaceMesh(saddleBaseMesh);
    //saddleSurfaceMesh.write(dataOutPath + "\\core_tests\\saddleSurfaceMesh.vtk");

    const auto curvature = vertex_curvature(saddleSurfaceMesh, Vertex(4));
    const auto normal = Normals::compute_vertex_normal(saddleSurfaceMesh, Vertex(4));
    const auto standardNormal = vec3{ 0, 0, 1 };  // for Z-axis alignment
    const auto reorientToStandardNormal = rotation_matrix(standardNormal, normal);

    const auto standardX = vec3{ 1, 0, 0 }; // for X-axis alignment

    const vec3 arbitraryVec = (std::fabs(normal[0]) > 0.9) ? vec3{ 0.0, 1.0, 0.0 } : vec3{ 1.0, 0.0, 0.0 };
    const auto tangent1 = normalize(cross(normal, arbitraryVec));
    const auto reorientToStandardX = rotation_matrix(standardX, tangent1);

    const auto tangentHyperboloidCenter = saddleBaseMesh.Vertices[4] - (1.0f / curvature.max) * normal;

    const auto translateToHyperboloidCenter = translation_matrix(-tangentHyperboloidCenter);
    const auto translateToFieldCenter = translation_matrix(fieldCenter);
    const auto affineTransformToStandardTangentPoint = translateToFieldCenter * (reorientToStandardNormal * reorientToStandardX) * translateToHyperboloidCenter;

    constexpr Scalar CURVATURE_EPSILON{ 1e-6f };
    QuadricParams params{
        fieldCenter,
        0.5f / std::sqrt(std::max(std::fabs(curvature.max), CURVATURE_EPSILON)),
        0.5f / std::sqrt(std::max(std::fabs(curvature.max), CURVATURE_EPSILON)),
        0.5f / std::sqrt(std::max(std::fabs(curvature.min), CURVATURE_EPSILON)),
        vec3{1, 1, 1} * (roiHalfDim / 2.5f),
        DistanceUnion
    };

    ApplyHyperboloidDistanceFieldToGrid(grid, params);

    SurfaceMesh transformedSaddleMesh(saddleSurfaceMesh);
    transformedSaddleMesh *= affineTransformToStandardTangentPoint;
    transformedSaddleMesh.write(dataOutPath + "\\core_tests\\transformedSaddleMesh.vtk");
    ExportToVTI(dataOutPath + "\\core_tests\\saddleMeshHyperboloid", grid);
}

TEST(ApproximateMeshPatchWithQuadric, ApplyEllipsoidDistanceFieldFromMesh1RingConvex)
{
    // Arrange
    constexpr float roiHalfDim = 5.0f;
    constexpr float roiDim = 2.0f * roiHalfDim;
    const auto box = BoundingBox{ vec3{0.0f, 0.0f, -roiHalfDim}, vec3{roiDim, roiDim, roiHalfDim} };
    const auto fieldCenter = box.center();
    constexpr float cellSize = 0.1f;
    constexpr double initVal = DEFAULT_SCALAR_GRID_INIT_VAL;
    ScalarGrid grid(cellSize, box, initVal);
    
    BaseMeshGeometryData convexBaseMesh;
    convexBaseMesh.Vertices = std::vector{
        Point{0.035850f, 0.248725f, 0.525021f},
        Point{0.549778f, 0.208385f, 0.180541f},
        Point{-0.011141f, 0.539888f, 0.022254f},
        Point{-0.441162f, 0.128599f, -0.144535f},
        Point{-0.065786f, 0.314575f, -0.487154f},
        Point{-0.413034f, 0.047552f, 0.450408f},
        Point{0.383242f, 0.184812f, -0.414398f},
    };
    convexBaseMesh.PolyIndices = std::vector{
        std::vector<unsigned int>{0, 1, 2},
        std::vector<unsigned int>{3, 2, 4},
        std::vector<unsigned int>{3, 5, 2},
        std::vector<unsigned int>{5, 0, 2},
        std::vector<unsigned int>{2, 6, 4},
        std::vector<unsigned int>{2, 1, 6}
    };

    const auto convexSurfaceMesh = ConvertBufferGeomToPMPSurfaceMesh(convexBaseMesh);
    convexSurfaceMesh.write(dataOutPath + "\\core_tests\\convexSurfaceMesh.vtk");

    const auto curvature = vertex_curvature(convexSurfaceMesh, Vertex(4));
    const auto normal = Normals::compute_vertex_normal(convexSurfaceMesh, Vertex(4));
    const auto standardNormal = vec3{ 0, 0, 1 }; // For Z-axis alignment
    const auto reorientToStandardNormal = rotation_matrix(standardNormal, normal);
    
    const vec3 arbitraryVec = (std::fabs(normal[0]) > 0.9) ? vec3{ 0.0, 1.0, 0.0 } : vec3{ 1.0, 0.0, 0.0 };
    const auto tangent1 = normalize(cross(normal, arbitraryVec));
    const auto reorientToStandardX = rotation_matrix(vec3{ 1, 0, 0 }, tangent1); // Reorient to the X-axis

    const auto ellipsoidCenter = convexBaseMesh.Vertices[4] - (1.0f / curvature.max) * normal;
    const auto translateToEllipsoidCenter = translation_matrix(-ellipsoidCenter);
    const auto translateToFieldCenter = translation_matrix(fieldCenter);
    const auto affineTransform = /*translateToFieldCenter * (reorientToStandardNormal * reorientToStandardX) **/ translateToEllipsoidCenter;

    constexpr Scalar CURVATURE_EPSILON{ 1e-6f };
    QuadricParams params{
        fieldCenter,
        1.0f / std::sqrt(std::max(std::fabs(curvature.min), CURVATURE_EPSILON)),
        1.0f / std::sqrt(std::max(std::fabs(curvature.mean), CURVATURE_EPSILON)),
        1.0f / std::sqrt(std::max(std::fabs(curvature.max), CURVATURE_EPSILON)),
        vec3{1, 1, 1} *(roiHalfDim / 2.5f),
        DistanceUnion
    };

    // Act
    ApplyEllipsoidDistanceFieldToGrid(grid, params);

    SurfaceMesh transformedConvexMesh(convexSurfaceMesh);
    transformedConvexMesh *= affineTransform;
    transformedConvexMesh.write(dataOutPath + "\\core_tests\\transformedConvexMesh.vtk");
    ExportToVTI(dataOutPath + "\\core_tests\\convexMeshEllipsoid", grid);
}

TEST(ApproximateMeshPatchWithQuadric, ApplyEllipsoidDistanceFieldFromMesh1RingConcave)
{
    // Arrange
    constexpr float roiHalfDim = 5.0f;
    constexpr float roiDim = 2.0f * roiHalfDim;
    const auto box = BoundingBox{ vec3{0.0f, 0.0f, -roiHalfDim}, vec3{roiDim, roiDim, roiHalfDim} };
    const auto fieldCenter = box.center();
    constexpr float cellSize = 0.1f;
    constexpr double initVal = DEFAULT_SCALAR_GRID_INIT_VAL;
    ScalarGrid grid(cellSize, box, initVal);

    BaseMeshGeometryData concaveBaseMesh;
    concaveBaseMesh.Vertices = std::vector{
        Point{0.035850f, 0.248725f, 0.525021f},
        Point{0.549778f, 0.208385f, 0.180541f},
        Point{-0.011141f, -0.149904f, 0.022254f},
        Point{-0.441162f, 0.128599f, -0.144535f},
        Point{-0.065786f, 0.314575f, -0.487154f},
        Point{-0.413034f, 0.047552f, 0.450408f},
        Point{0.383242f, 0.184812f, -0.414398f},
    };
    concaveBaseMesh.PolyIndices = std::vector{
        std::vector<unsigned int>{0, 1, 2},
        std::vector<unsigned int>{3, 2, 4},
        std::vector<unsigned int>{3, 5, 2},
        std::vector<unsigned int>{5, 0, 2},
        std::vector<unsigned int>{2, 6, 4},
        std::vector<unsigned int>{2, 1, 6}
    };

    const auto concaveSurfaceMesh = ConvertBufferGeomToPMPSurfaceMesh(concaveBaseMesh);
    concaveSurfaceMesh.write(dataOutPath + "\\core_tests\\concaveSurfaceMesh.vtk");

    const auto curvature = vertex_curvature(concaveSurfaceMesh, Vertex(4));
    const auto normal = Normals::compute_vertex_normal(concaveSurfaceMesh, Vertex(4));
    const auto standardNormal = vec3{ 0, 0, 1 };  // For Z-axis alignment
    const auto reorientToStandardNormal = rotation_matrix(standardNormal, normal);
    
    const vec3 arbitraryVec = (std::fabs(normal[0]) > 0.9) ? vec3{ 0.0, 1.0, 0.0 } : vec3{ 1.0, 0.0, 0.0 };
    const auto tangent1 = normalize(cross(normal, arbitraryVec));
    const auto reorientToStandardX = rotation_matrix(vec3{ 1, 0, 0 }, tangent1); // Reorient to the X-axis

    const auto ellipsoidCenter = concaveBaseMesh.Vertices[4] - (1.0f / curvature.max) * normal;
    const auto translateToEllipsoidCenter = translation_matrix(-ellipsoidCenter);
    const auto translateToFieldCenter = translation_matrix(fieldCenter);
    const auto affineTransform = translateToFieldCenter * (reorientToStandardNormal * reorientToStandardX) * translateToEllipsoidCenter;

    constexpr Scalar CURVATURE_EPSILON{ 1e-6f };
    QuadricParams params{
        fieldCenter,
        1.0f / std::sqrt(std::max(std::fabs(curvature.min), CURVATURE_EPSILON)),
        1.0f / std::sqrt(std::max(std::fabs(curvature.mean), CURVATURE_EPSILON)),
        1.0f / std::sqrt(std::max(std::fabs(curvature.max), CURVATURE_EPSILON)),
        vec3{1, 1, 1} *(roiHalfDim / 2.5f),
        DistanceUnion
    };

    // Act
    ApplyEllipsoidDistanceFieldToGrid(grid, params);

    SurfaceMesh transformedConcaveMesh(concaveSurfaceMesh);
    transformedConcaveMesh *= affineTransform;
    transformedConcaveMesh.write(dataOutPath + "\\core_tests\\transformedConcaveMesh.vtk");
    ExportToVTI(dataOutPath + "\\core_tests\\concaveMeshEllipsoid", grid);
}