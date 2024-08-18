#include "gtest/gtest.h"

#include "geometry/MeshAnalysis.h"

#include <pmp/Types.h>
#include <optional>

#include "geometry/CollisionKdTree.h"
#include "geometry/GeometryAdapters.h"
#include "geometry/IcoSphereBuilder.h"
#include "pmp/algorithms/CurveFactory.h"

using namespace Geometry;
using namespace pmp;

// Suite: PMPManifoldCurve2DHasSelfIntersectionsTests

TEST(PMPManifoldCurve2DHasSelfIntersectionsTests, EmptyCurveThrowsException)
{
    // Arrange
    const ManifoldCurve2D emptyCurve;

    // Act & Assert
    EXPECT_THROW(auto ret = PMPManifoldCurve2DHasSelfIntersections(emptyCurve), std::invalid_argument);
}

TEST(PMPManifoldCurve2DHasSelfIntersectionsTests, SimpleClosedCurveHasNoIntersections)
{
    // Arrange
    const ManifoldCurve2D simpleCurve = CurveFactory::circle(Point2(0.0f, 0.0f), 1.0f, 32);

    // Act
    const bool hasIntersections = PMPManifoldCurve2DHasSelfIntersections(simpleCurve);

    // Assert
    EXPECT_FALSE(hasIntersections);
}

TEST(PMPManifoldCurve2DHasSelfIntersectionsTests, SelfIntersectingCurveReturnsTrue)
{
    // Arrange
    ManifoldCurve2D selfIntersectingCurve;
    selfIntersectingCurve.add_vertex(Point2(0.0f, 0.0f));
    selfIntersectingCurve.add_vertex(Point2(1.0f, 0.0f));
    selfIntersectingCurve.add_vertex(Point2(0.5f, 1.0f));
    selfIntersectingCurve.add_vertex(Point2(1.5f, 1.0f));
    selfIntersectingCurve.add_edge(Vertex(0), Vertex(1));
    selfIntersectingCurve.add_edge(Vertex(1), Vertex(2));
    selfIntersectingCurve.add_edge(Vertex(2), Vertex(3));
    selfIntersectingCurve.add_edge(Vertex(3), Vertex(0)); // Self-intersection occurs here

    // Act
    const bool hasIntersections = PMPManifoldCurve2DHasSelfIntersections(selfIntersectingCurve);

    // Assert
    EXPECT_TRUE(hasIntersections);
}

// Suite: IsPointInsidePMPManifoldCurveTests

TEST(IsPointInsidePMPManifoldCurveTests, PointInsideCircleReturnsTrue)
{
    // Arrange
    const ManifoldCurve2D circle = CurveFactory::circle(Point2(0.0f, 0.0f), 1.0f, 32);
    const Point2 pointInside(0.5f, 0.5f);

    // Act
    const bool isInside = IsPointInsidePMPManifoldCurve(pointInside, circle);

    // Assert
    EXPECT_TRUE(isInside);
}

TEST(IsPointInsidePMPManifoldCurveTests, PointOutsideCircleReturnsFalse)
{
    // Arrange
    const ManifoldCurve2D circle = CurveFactory::circle(Point2(0.0f, 0.0f), 1.0f, 32);
    const Point2 pointOutside(2.0f, 2.0f);

    // Act
    const bool isInside = IsPointInsidePMPManifoldCurve(pointOutside, circle);

    // Assert
    EXPECT_FALSE(isInside);
}

TEST(IsPointInsidePMPManifoldCurveTests, PointOnCircleBoundaryReturnsFalse)
{
    // Arrange
    const ManifoldCurve2D circle = CurveFactory::circle(Point2(0.0f, 0.0f), 1.0f, 32);
    const Point2 pointOnBoundary(1.0f, 0.0f);

    // Act
    const bool isInside = IsPointInsidePMPManifoldCurve(pointOnBoundary, circle);

    // Assert
    EXPECT_FALSE(isInside);
}

TEST(IsPointInsidePMPManifoldCurveTests, PointInsideConcaveCurveReturnsTrue)
{
    // Arrange
    ManifoldCurve2D concaveCurve;
    concaveCurve.add_vertex(Point2(0.0f, 0.0f));
    concaveCurve.add_vertex(Point2(2.0f, 0.0f));
    concaveCurve.add_vertex(Point2(1.0f, 1.0f));
    concaveCurve.add_vertex(Point2(1.0f, 2.0f));
    concaveCurve.add_edge(Vertex(0), Vertex(1));
    concaveCurve.add_edge(Vertex(1), Vertex(2));
    concaveCurve.add_edge(Vertex(2), Vertex(3));
    concaveCurve.add_edge(Vertex(3), Vertex(0));
    const Point2 pointInside(1.0f, 0.5f);

    // Act
    const bool isInside = IsPointInsidePMPManifoldCurve(pointInside, concaveCurve);

    // Assert
    EXPECT_TRUE(isInside);
}

TEST(IsPointInsidePMPManifoldCurveTests, PointOutsideConcaveCurveReturnsFalse)
{
    // Arrange
    ManifoldCurve2D concaveCurve;
    concaveCurve.add_vertex(Point2(0.0f, 0.0f));
    concaveCurve.add_vertex(Point2(2.0f, 0.0f));
    concaveCurve.add_vertex(Point2(1.0f, 1.0f));
    concaveCurve.add_vertex(Point2(1.0f, 2.0f));
    concaveCurve.add_edge(Vertex(0), Vertex(1));
    concaveCurve.add_edge(Vertex(1), Vertex(2));
    concaveCurve.add_edge(Vertex(2), Vertex(3));
    concaveCurve.add_edge(Vertex(3), Vertex(0));
    const Point2 pointOutside(1.5f, 1.5f);

    // Act
    const bool isInside = IsPointInsidePMPManifoldCurve(pointOutside, concaveCurve);

    // Assert
    EXPECT_FALSE(isInside);
}


// Suite: IsPointInsidePMPSurfaceMeshTests

TEST(IsPointInsidePMPSurfaceMeshTests, IcoSphere_SamplePtRaysPassingThroughFace)
{
    // Arrange
    constexpr unsigned int icoSphereSubdiv = 2;
    constexpr float icoSphereRadius = 1.0f;
    IcoSphereBuilder icoBuilder({ icoSphereSubdiv, icoSphereRadius });
    icoBuilder.BuildBaseData();
    icoBuilder.BuildPMPSurfaceMesh();
    SurfaceMesh icoSphereMesh = icoBuilder.GetPMPSurfaceMeshResult();
    const auto rotMatrix = rotation_matrix(vec3{ 0.0f, 0.0f, 1.0f }, 0.1f * static_cast<float>(M_PI));
    icoSphereMesh *= rotMatrix;

    const PMPSurfaceMeshAdapter icoSphereAdapter(std::make_shared<SurfaceMesh>(icoSphereMesh));
    const auto surfaceKdTree = std::make_shared<CollisionKdTree>(icoSphereAdapter, CenterSplitFunction);

    const Point insidePoint(0.5f, 0.0f, 0.0f); // Clearly inside the sphere
    const Point onSurfacePoint(1.0f, 0.0f, 0.0f); // On the surface of the sphere
    const Point outsidePoint(1.5f, 0.0f, 0.0f); // Clearly outside the sphere
    const Point outsideNegXPoint(-1.5f, 0.0f, 0.0f); // Clearly outside, but the +X ray intersects the surface twice

    // Act
    const bool isInsidePointInside = IsPointInsidePMPSurfaceMesh(insidePoint, surfaceKdTree);
    const bool isOnSurfacePointInside = IsPointInsidePMPSurfaceMesh(onSurfacePoint, surfaceKdTree);
    const bool isOutsidePointInside = IsPointInsidePMPSurfaceMesh(outsidePoint, surfaceKdTree);
    const bool isOutsideNegXPointInside = IsPointInsidePMPSurfaceMesh(outsideNegXPoint, surfaceKdTree);

    // Assert
    EXPECT_TRUE(isInsidePointInside) << "The inside point should be inside the icosphere.";
    EXPECT_FALSE(isOnSurfacePointInside) << "The point on the surface should not be considered inside.";
    EXPECT_FALSE(isOutsidePointInside) << "The outside point should not be inside the icosphere.";
    EXPECT_FALSE(isOutsideNegXPointInside) << "The outside X<0 point should not be inside the icosphere.";
}

TEST(IsPointInsidePMPSurfaceMeshTests, IcoSphere_SamplePtRaysPassingThroughSharedVertex)
{
    // Arrange
    constexpr unsigned int icoSphereSubdiv = 2;
    constexpr float icoSphereRadius = 1.0f;
    IcoSphereBuilder icoBuilder({ icoSphereSubdiv, icoSphereRadius });
    icoBuilder.BuildBaseData();
    icoBuilder.BuildPMPSurfaceMesh();
    const SurfaceMesh& icoSphereMesh = icoBuilder.GetPMPSurfaceMeshResult();

    const PMPSurfaceMeshAdapter icoSphereAdapter(std::make_shared<SurfaceMesh>(icoSphereMesh));
    const auto surfaceKdTree = std::make_shared<CollisionKdTree>(icoSphereAdapter, CenterSplitFunction);

    const Point insidePoint(0.5f, 0.0f, 0.0f); // Clearly inside the sphere
    const Point onSurfacePoint(1.0f, 0.0f, 0.0f); // On the surface of the sphere
    const Point outsidePoint(1.5f, 0.0f, 0.0f); // Clearly outside the sphere
    const Point outsideNegXPoint(-1.5f, 0.0f, 0.0f); // Clearly outside, but the +X ray intersects the surface twice

    // Act
    const bool isInsidePointInside = IsPointInsidePMPSurfaceMesh(insidePoint, surfaceKdTree);
    const bool isOnSurfacePointInside = IsPointInsidePMPSurfaceMesh(onSurfacePoint, surfaceKdTree);
    const bool isOutsidePointInside = IsPointInsidePMPSurfaceMesh(outsidePoint, surfaceKdTree);
    const bool isOutsideNegXPointInside = IsPointInsidePMPSurfaceMesh(outsideNegXPoint, surfaceKdTree);

    // Assert
    EXPECT_TRUE(isInsidePointInside) << "The inside point should be inside the icosphere.";
    EXPECT_FALSE(isOnSurfacePointInside) << "The point on the surface should not be considered inside.";
    EXPECT_FALSE(isOutsidePointInside) << "The outside point should not be inside the icosphere.";
    EXPECT_FALSE(isOutsideNegXPointInside) << "The outside X<0 point should not be inside the icosphere.";
}

TEST(IsPointInsidePMPSurfaceMeshTests, IcoSphere_SamplePtRaysPassingThroughSharedEdge)
{
	// Arrange
    constexpr unsigned int icoSphereSubdiv = 2;
    constexpr float icoSphereRadius = 1.0f;
    IcoSphereBuilder icoBuilder({ icoSphereSubdiv, icoSphereRadius });
    icoBuilder.BuildBaseData();
    icoBuilder.BuildPMPSurfaceMesh();
    SurfaceMesh icoSphereMesh = icoBuilder.GetPMPSurfaceMeshResult();
    const auto rotMatrix = rotation_matrix(vec3{ 0.0f, 1.0f, 0.0f }, 0.1f * static_cast<float>(M_PI));
    icoSphereMesh *= rotMatrix;

    const PMPSurfaceMeshAdapter icoSphereAdapter(std::make_shared<SurfaceMesh>(icoSphereMesh));
    const auto surfaceKdTree = std::make_shared<CollisionKdTree>(icoSphereAdapter, CenterSplitFunction);

    const Point insidePoint(0.5f, 0.0f, 0.0f); // Clearly inside the sphere
    const Point onSurfacePoint(1.0f, 0.0f, 0.0f); // On the surface of the sphere
    const Point outsidePoint(1.5f, 0.0f, 0.0f); // Clearly outside the sphere
    const Point outsideNegXPoint(-1.5f, 0.0f, 0.0f); // Clearly outside, but the +X ray intersects the surface twice

    // Act
    const bool isInsidePointInside = IsPointInsidePMPSurfaceMesh(insidePoint, surfaceKdTree);
    const bool isOnSurfacePointInside = IsPointInsidePMPSurfaceMesh(onSurfacePoint, surfaceKdTree);
    const bool isOutsidePointInside = IsPointInsidePMPSurfaceMesh(outsidePoint, surfaceKdTree);
    const bool isOutsideNegXPointInside = IsPointInsidePMPSurfaceMesh(outsideNegXPoint, surfaceKdTree);

    // Assert
    EXPECT_TRUE(isInsidePointInside) << "The inside point should be inside the icosphere.";
    EXPECT_FALSE(isOnSurfacePointInside) << "The point on the surface should not be considered inside.";
    EXPECT_FALSE(isOutsidePointInside) << "The outside point should not be inside the icosphere.";
    EXPECT_FALSE(isOutsideNegXPointInside) << "The outside X<0 point should not be inside the icosphere.";
}

TEST(IsPointInsidePMPSurfaceMeshTests, ConcaveDumbbellShapeSamplePtRays)
{
    // Arrange
    BaseMeshGeometryData meshData;
    meshData.Vertices = {
        Point{ 0.955864f,  0.797707f, -0.797707f }, Point{ 0.894073f, -0.838166f, -0.838166f },
        Point{ 0.894073f,  0.838166f,  0.838166f }, Point{ 0.955864f, -0.797707f,  0.797707f },
        Point{-0.894073f,  0.838166f, -0.838166f }, Point{-0.955864f, -0.797707f, -0.797707f },
        Point{-0.955864f,  0.797707f,  0.797707f }, Point{-0.894073f, -0.838166f,  0.838166f },
        Point{-0.058848f, -0.401921f, -0.401921f }, Point{-0.058848f,  0.401921f,  0.401921f },
        Point{ 0.058848f, -0.401921f,  0.401921f }, Point{ 0.058848f,  0.401921f, -0.401921f }
    };
    meshData.PolyIndices = {
        { 4, 9, 11 }, { 9, 7, 10 }, { 6, 5, 7 }, { 1, 10, 8 },
        { 0, 3, 1 }, {11, 1, 8 }, { 4, 8, 5 }, { 8, 7, 5 },
        { 2, 10, 3 }, {11, 2, 0 }, { 4, 6, 9 }, { 9, 6, 7 },
        { 6, 4, 5 }, { 1, 3, 10 }, { 0, 2, 3 }, {11, 0, 1 },
        { 4, 11, 8 }, { 8, 10, 7 }, { 2, 9, 10 }, {11, 9, 2 }
    };
    const BaseMeshAdapter baseAdapter(std::make_shared<BaseMeshGeometryData>(meshData));
    const auto surfaceKdTree = std::make_shared<CollisionKdTree>(baseAdapter, CenterSplitFunction);

    const Point outsidePointTwoTransitions(0.02f, 0.11925f, 0.61f);
    const Point insidePointOneTransition(0.02f, 0.11925f, -0.014083f);
    const Point outsidePointNoTransitions(0.02f, 0.11925f, 0.90422f);
    const Point outsidePointFourTransitions(-1.6135f, -0.029103f, 0.70568f);
    const Point insidePointThreeTransitions(-0.77095f, -0.2763f, 0.53413f);

    // Act
    const bool resultOutsideTwoTransitions = IsPointInsidePMPSurfaceMesh(outsidePointTwoTransitions, surfaceKdTree);
    const bool resultInsideOneTransition = IsPointInsidePMPSurfaceMesh(insidePointOneTransition, surfaceKdTree);
    const bool resultOutsideNoTransitions = IsPointInsidePMPSurfaceMesh(outsidePointNoTransitions, surfaceKdTree);
    const bool resultOutsideFourTransitions = IsPointInsidePMPSurfaceMesh(outsidePointFourTransitions, surfaceKdTree);
    const bool resultInsideThreeTransitions = IsPointInsidePMPSurfaceMesh(insidePointThreeTransitions, surfaceKdTree);

    // Assert
    ASSERT_FALSE(resultOutsideTwoTransitions);
    ASSERT_TRUE(resultInsideOneTransition);
    ASSERT_FALSE(resultOutsideNoTransitions);
    ASSERT_FALSE(resultOutsideFourTransitions);
    ASSERT_TRUE(resultInsideThreeTransitions);
}
