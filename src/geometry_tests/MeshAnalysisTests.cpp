#include "gtest/gtest.h"

#include "geometry/MeshAnalysis.h"

#include <pmp/Types.h>
#include <optional>

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