#include "gtest/gtest.h"
#include "pmp/ManifoldCurve2D.h"
#include "pmp/algorithms/CurveFactory.h"

using namespace pmp;

constexpr float COORD_TOLERANCE = 1e-6f;

// Iterator Tests
TEST(EmptyManifoldCurve2D_Tests, VertexIterator_Empty)
{
    ManifoldCurve2D curve;
    auto it = curve.vertices_begin();
    EXPECT_EQ(it, curve.vertices_end());
}

TEST(EmptyManifoldCurve2D_Tests, VertexIterator_NonEmpty)
{
    ManifoldCurve2D curve;
    curve.add_vertex(Point2(0, 0));
    curve.add_vertex(Point2(1, 0));

    size_t count = 0;
    for (auto v : curve.vertices())
        ++count;

    EXPECT_EQ(count, 2);
}

TEST(EmptyManifoldCurve2D_Tests, EdgeIterator_Empty)
{
    ManifoldCurve2D curve;
    auto it = curve.edges_begin();
    EXPECT_EQ(it, curve.edges_end());
}

TEST(EmptyManifoldCurve2D_Tests, EdgeIterator_NonEmpty)
{
    ManifoldCurve2D curve;
    auto v0 = curve.add_vertex(Point2(0, 0));
    auto v1 = curve.add_vertex(Point2(1, 0));
    curve.add_edge(v0, v1);

    size_t count = 0;
    for (auto e : curve.edges())
        ++count;

    EXPECT_EQ(count, 1);
}

// Property Management Tests
TEST(EmptyManifoldCurve2D_Tests, AddGetVertexProperty)
{
    ManifoldCurve2D curve;
    auto prop = curve.add_vertex_property<int>("v:prop", 0);
    const auto v = curve.add_vertex(Point2(0, 0));
    prop[v] = 42;

    auto retrieved_prop = curve.get_vertex_property<int>("v:prop");
    EXPECT_EQ(retrieved_prop[v], 42);
}
TEST(EmptyManifoldCurve2D_Tests, AddGetEdgeProperty)
{
    ManifoldCurve2D curve;
    auto prop = curve.add_edge_property<int>("e:prop", 0);
    const auto v0 = curve.add_vertex(Point2(0, 0));
    const auto v1 = curve.add_vertex(Point2(1, 0));
    const auto e = curve.add_edge(v0, v1);
    prop[e] = 42;

    auto retrieved_prop = curve.get_edge_property<int>("e:prop");
    EXPECT_EQ(retrieved_prop[e], 42);
}

// Element Allocation Tests
TEST(EmptyManifoldCurve2D_Tests, AllocateVertex)
{
    ManifoldCurve2D curve;
    const auto v = curve.new_vertex();
    EXPECT_TRUE(curve.is_valid(v));
}

TEST(EmptyManifoldCurve2D_Tests, AllocateEdge)
{
    ManifoldCurve2D curve;
    const auto e = curve.new_edge();
    EXPECT_TRUE(curve.is_valid(e));
}

TEST(EmptyManifoldCurve2D_Tests, AllocateEdgeBetweenVertices)
{
    ManifoldCurve2D curve;
    const auto v0 = curve.add_vertex(Point2(0, 0));
    const auto v1 = curve.add_vertex(Point2(1, 0));
    const auto e = curve.new_edge(v0, v1);
    EXPECT_TRUE(curve.is_valid(e));
}

// Deletion and Validity Checks
TEST(EmptyManifoldCurve2D_Tests, DeleteVertex)
{
    ManifoldCurve2D curve;
    const auto v = curve.add_vertex(Point2(0, 0));
    curve.delete_vertex(v);
    EXPECT_TRUE(curve.is_deleted(v));
}

TEST(EmptyManifoldCurve2D_Tests, DeleteEdge)
{
    ManifoldCurve2D curve;
    const auto v0 = curve.add_vertex(Point2(0, 0));
    const auto v1 = curve.add_vertex(Point2(1, 0));
    const auto e = curve.add_edge(v0, v1);
    curve.delete_edge(e);
    EXPECT_TRUE(curve.is_deleted(e));
}

// Connectivity Management Tests
TEST(EmptyManifoldCurve2D_Tests, ConnectivityManagement)
{
    ManifoldCurve2D curve;
    const auto v0 = curve.add_vertex(Point2(0, 0));
    const auto v1 = curve.add_vertex(Point2(1, 0));
    const auto e = curve.add_edge(v0, v1);
    curve.set_edge_from(v0, e);
    curve.set_edge_to(v1, e);

    EXPECT_EQ(curve.edge_from(v0), e);
    EXPECT_EQ(curve.edge_to(v1), e);
    EXPECT_TRUE(curve.is_boundary(v0));
    EXPECT_FALSE(curve.is_isolated(v0));
}

class ManifoldCurve2DTest_ClosedArc : public ::testing::Test
{
protected:
    ManifoldCurve2D curve;

    void SetUp() override
    {
        curve = CurveFactory::circle(Point2(0, 0), 1.0, 32);
    }
};

TEST_F(ManifoldCurve2DTest_ClosedArc, VertexIterator_NonEmpty)
{
    size_t count = 0;
    for (auto v : curve.vertices())
        ++count;

    EXPECT_EQ(count, 32); // 32 segments
}

TEST_F(ManifoldCurve2DTest_ClosedArc, EdgeIterator_NonEmpty)
{
    size_t count = 0;
    for (auto e : curve.edges())
        ++count;

    EXPECT_EQ(count, 32); // 32 segments, no duplicate edge
}

static [[nodiscard]] bool IsContinuousBetweenVertices(const ManifoldCurve2D& curve, Vertex v0, Vertex v1)
{
    if (!curve.is_valid(v0) || !curve.is_valid(v1))
        return false;

    const Edge eFrom = curve.edge_from(v0);
    if (!curve.is_valid(eFrom))
        return false;

    // Iterate through the circle and verify connectivity
    Vertex currentVertex = v0;
    Edge currentEdge = eFrom;
    size_t count = 0;

    do
    {
        // Verify the start and end vertices of the current edge
        if (curve.from_vertex(currentEdge) != currentVertex)
            return false;

        Vertex nextVertex = curve.to_vertex(currentEdge);
        if (currentVertex == nextVertex)
            return false;

        // Move to the next vertex and edge
        currentVertex = nextVertex;

        // Check if the current vertex is a boundary vertex
        if (curve.is_boundary(currentVertex) && currentVertex != v1)
            return false;

        currentEdge = curve.edge_from(currentVertex);
        count++;

        // Check if we have looped back to the starting vertex
        if (currentVertex == v0 && v0 != v1)
            return false;

    } while (count < curve.n_vertices() && currentVertex != v1);

    return currentVertex == v1;
}

static [[nodiscard]] bool IsBackwardsContinuousBetweenVertices(const ManifoldCurve2D& curve, Vertex v0, Vertex v1)
{
    const Edge eTo = curve.edge_to(v0);

    // Iterate through the circle in reverse and verify connectivity
    Vertex currentVertex = v0;
    Edge currentEdge = eTo;
    size_t count = 0;

    do
    {
        // Verify the start and end vertices of the current edge
        if (curve.to_vertex(currentEdge) != currentVertex) 
            return false;
        Vertex prevVertex = curve.from_vertex(currentEdge);
        if (currentVertex == prevVertex)
            return false;

        // Move to the previous vertex and edge
        currentVertex = prevVertex;
        currentEdge = curve.edge_to(currentVertex);
        count++;
    } while (count < curve.n_vertices() && currentVertex != v1);

    if (v0 == v1)
        return currentVertex == v0;

    return currentVertex == v1;
}



TEST_F(ManifoldCurve2DTest_ClosedArc, VerifyConnectivity)
{
    // Get the original vertex and its adjacent edges
    const Vertex originalVertex = *curve.vertices_begin();

    // Verify connectivity
    EXPECT_TRUE(IsContinuousBetweenVertices(curve, originalVertex, originalVertex));
    EXPECT_TRUE(IsBackwardsContinuousBetweenVertices(curve, originalVertex, originalVertex));
}

TEST_F(ManifoldCurve2DTest_ClosedArc, RemoveEdgeAndVerifyBoundary)
{
    // Get an edge to delete
    const Edge edgeToDelete = *curve.edges_begin();

    // Get the start and end vertices of the edge
    const Vertex startVertex = curve.from_vertex(edgeToDelete);
    const Vertex endVertex = curve.to_vertex(edgeToDelete);

    // Delete the edge
    curve.delete_edge(edgeToDelete);
    curve.garbage_collection();

    // Verify that both vertices are now boundary vertices
    EXPECT_TRUE(curve.is_boundary(startVertex));
    EXPECT_TRUE(curve.is_boundary(endVertex));
}

TEST_F(ManifoldCurve2DTest_ClosedArc, RemoveConsecutiveEdgesAndVerifyIsolation)
{
    // Get the first edge to delete
    const Edge firstEdgeToDelete = *curve.edges_begin();
    // Get the second edge to delete by traversing the connectivity
    const Edge secondEdgeToDelete = curve.edge_from(curve.to_vertex(firstEdgeToDelete));

    // Get the vertices of the edges
    const Vertex firstStartVertex = curve.from_vertex(firstEdgeToDelete);
    const Vertex middleVertex = curve.to_vertex(firstEdgeToDelete);
    const Vertex secondEndVertex = curve.to_vertex(secondEdgeToDelete);

    // Delete the first edge
    curve.delete_edge(firstEdgeToDelete);
    // Delete the second edge
    curve.delete_edge(secondEdgeToDelete);
    curve.garbage_collection();

    // Verify that the middle vertex is now isolated
    EXPECT_TRUE(curve.is_isolated(middleVertex));

    // Verify that the start and end vertices are still boundary vertices
    EXPECT_TRUE(curve.is_boundary(firstStartVertex));
    EXPECT_TRUE(curve.is_boundary(secondEndVertex));
}

// Memory Management & Vertex/Edge Deletion

TEST_F(ManifoldCurve2DTest_ClosedArc, GarbageCollectionDeleteVertex)
{
    // Arrange
    const auto v0 = Vertex{ 8 };
    curve.delete_vertex(v0);

    // Act
    curve.garbage_collection();

    // Assert
    EXPECT_EQ(curve.n_vertices(), 31);
    EXPECT_EQ(curve.n_edges(), 31);
    EXPECT_EQ(curve.edge_from(Vertex{ 6 }), curve.edge_to(Vertex{ 7 }));
    EXPECT_EQ(curve.edge_to(Vertex{ 7 }), Edge{ 6 });
    EXPECT_EQ(curve.edge_from(Vertex{ 0 }), Edge{ 0 });
    EXPECT_TRUE(IsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
    EXPECT_TRUE(IsBackwardsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
}

TEST_F(ManifoldCurve2DTest_ClosedArc, GarbageCollectionDeleteEdge)
{
    // Arrange
    const auto v1 = Vertex{ 8 };
    const auto v2 = Vertex{ 9 };
    curve.delete_edge(curve.edge_from(v1));

    // Act
    curve.garbage_collection();

    // Assert
    EXPECT_EQ(curve.n_vertices(), 32);
    EXPECT_EQ(curve.n_edges(), 31);
    EXPECT_EQ(curve.edge_from(Vertex{ 31 }), curve.edge_to(Vertex{ 0 }));
    EXPECT_EQ(curve.edge_to(Vertex{ 31 }), Edge{ 30 });
    EXPECT_EQ(curve.edge_from(Vertex{ 0 }), Edge{ 0 });
    EXPECT_TRUE(curve.is_boundary(v1));
    EXPECT_TRUE(curve.is_boundary(v2));
    EXPECT_TRUE(IsContinuousBetweenVertices(curve, Vertex{ 9 }, Vertex{ 8 }));
    EXPECT_TRUE(IsBackwardsContinuousBetweenVertices(curve, Vertex{ 8 }, Vertex{ 9 }));
    EXPECT_FALSE(IsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
    EXPECT_FALSE(IsBackwardsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
}

TEST_F(ManifoldCurve2DTest_ClosedArc, GarbageCollectionDeleteVertexAndEdge)
{
    // Arrange
    const auto v0 = Vertex{ 0 };
    const auto v1 = Vertex{ 8 };
    const auto v2 = Vertex{ 9 };
    curve.delete_edge(curve.edge_from(v1));
    curve.delete_vertex(v0);

    // Act
    curve.garbage_collection();

    // Assert
    EXPECT_EQ(curve.n_vertices(), 31);
    EXPECT_EQ(curve.n_edges(), 30);
    EXPECT_EQ(curve.edge_from(Vertex{ 30 }), curve.edge_to(Vertex{ 0 }));
    EXPECT_EQ(curve.edge_to(Vertex{ 30 }), Edge{ 29 });
    EXPECT_EQ(curve.edge_from(Vertex{ 0 }), Edge{ 0 });
    EXPECT_TRUE(curve.is_boundary(v1));
    EXPECT_TRUE(curve.is_boundary(v2));
    EXPECT_TRUE(IsContinuousBetweenVertices(curve, Vertex{ 9 }, Vertex{ 8 }));
    EXPECT_TRUE(IsBackwardsContinuousBetweenVertices(curve, Vertex{ 8 }, Vertex{ 9 }));
}

TEST_F(ManifoldCurve2DTest_ClosedArc, GarbageCollectionDeleteConsecutiveVertices)
{
    // Arrange
    const auto v0 = Vertex{ 8 };
    const auto v1 = Vertex{ 9 };
    curve.delete_vertex(v0);
    curve.delete_vertex(v1);

    // Act
    curve.garbage_collection();

    // Assert
    EXPECT_EQ(curve.n_vertices(), 30);
    EXPECT_EQ(curve.n_edges(), 30);
    EXPECT_EQ(curve.edge_from(Vertex{ 7 }), curve.edge_to(Vertex{ 10 }));
    EXPECT_EQ(curve.edge_to(Vertex{ 10 }), Edge{ 7 });
    EXPECT_EQ(curve.edge_from(Vertex{ 0 }), Edge{ 0 });
    EXPECT_TRUE(IsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
    EXPECT_TRUE(IsBackwardsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
}

TEST_F(ManifoldCurve2DTest_ClosedArc, DeleteEdgeAndCheckIsolatedVertex)
{
    // Arrange
    const auto e1 = Edge{ 8 };
    const auto e2 = curve.edge_from(curve.to_vertex(e1));
    const auto vMiddle = curve.to_vertex(e1);
    curve.delete_edge(e1);
    curve.delete_edge(e2);

    // Act
    curve.garbage_collection();

    // Assert
    EXPECT_EQ(curve.n_vertices(), 32);
    EXPECT_EQ(curve.n_edges(), 30);
    EXPECT_TRUE(curve.is_isolated(vMiddle));
    EXPECT_FALSE(IsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
    EXPECT_FALSE(IsBackwardsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
}

// Memory Management & Additional Tessellation-Changing Operations

TEST_F(ManifoldCurve2DTest_ClosedArc, SplitEdge)
{
    // Arrange
    const auto e = Edge{ 8 };
    const auto param = 0.5f;

    // Act
    const auto newVertex = curve.split_edge(e, param);
    // no garbage_collection needed

    // Assert
    EXPECT_EQ(curve.n_vertices(), 33);
    EXPECT_EQ(curve.n_edges(), 33);
    EXPECT_TRUE(curve.is_valid(newVertex));
    EXPECT_EQ(curve.to_vertex(e), newVertex);
    EXPECT_EQ(curve.from_vertex(curve.edge_from(newVertex)), newVertex);
    EXPECT_TRUE(IsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
    EXPECT_TRUE(IsBackwardsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
}

TEST_F(ManifoldCurve2DTest_ClosedArc, CollapseEdgeKeepStart)
{
    // Arrange
    const auto e = Edge{ 8 };
    const auto v0 = curve.from_vertex(e);
    const auto v1 = curve.to_vertex(e);

    // Act
    curve.collapse_edge(e, true);
    curve.garbage_collection();

    // Assert
    EXPECT_EQ(curve.n_vertices(), 31);
    EXPECT_EQ(curve.n_edges(), 31);
    EXPECT_TRUE(curve.is_deleted(v1));
    EXPECT_EQ(curve.edge_from(v0), curve.edge_to(Vertex{ 9 }));
    EXPECT_EQ(curve.edge_to(Vertex{ 9 }), Edge{ 7 });
    EXPECT_EQ(curve.edge_from(Vertex{ 0 }), Edge{ 0 });
    EXPECT_TRUE(IsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
    EXPECT_TRUE(IsBackwardsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
}

TEST_F(ManifoldCurve2DTest_ClosedArc, CollapseEdgeKeepEnd)
{
    // Arrange
    const auto e = Edge{ 8 };
    const auto v0 = curve.from_vertex(e);
    const auto v1 = curve.to_vertex(e);
    const auto ePrev = curve.edge_to(v0);
    const auto eNext = curve.edge_from(v1);

    // Act
    curve.collapse_edge(e, false);
    curve.garbage_collection();

    // Assert
    EXPECT_EQ(curve.n_vertices(), 31);
    EXPECT_EQ(curve.n_edges(), 31);
    EXPECT_TRUE(curve.is_deleted(v0));
    EXPECT_FALSE(curve.is_deleted(v1));
    EXPECT_EQ(curve.from_vertex(ePrev), curve.to_vertex(eNext));
    EXPECT_EQ(curve.to_vertex(eNext), v1);
    EXPECT_TRUE(IsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
    EXPECT_TRUE(IsBackwardsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
}

TEST_F(ManifoldCurve2DTest_ClosedArc, SplitEdgeWithVertex)
{
    // Arrange
    const auto e = Edge{ 8 };
    const Point2 newPos(0.5f, 0.5f);

    // Act
    const auto newVertex = curve.split_edge_with_vertex(e, newPos);
    // No garbage_collection needed

    // Assert
    EXPECT_EQ(curve.n_vertices(), 33);
    EXPECT_EQ(curve.n_edges(), 33);
    EXPECT_TRUE(curve.is_valid(newVertex));
    EXPECT_EQ(curve.position(newVertex), newPos);
    EXPECT_EQ(curve.to_vertex(e), newVertex);
    EXPECT_EQ(curve.from_vertex(curve.edge_from(newVertex)), newVertex);
    EXPECT_TRUE(IsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
    EXPECT_TRUE(IsBackwardsContinuousBetweenVertices(curve, Vertex{ 0 }, Vertex{ 0 }));
}

class ManifoldCurve2DTest_OpenArc : public ::testing::Test
{
protected:
    ManifoldCurve2D curve;

    void SetUp() override
    {
        curve = CurveFactory::circle(Point2(0, 0), 1.0, 16, 0, M_PI); // Half circle
    }
};

TEST_F(ManifoldCurve2DTest_OpenArc, VertexIterator_NonEmpty)
{
    size_t count = 0;
    for (auto v : curve.vertices())
        ++count;

    EXPECT_EQ(count, 16); // 16 segments
}

TEST_F(ManifoldCurve2DTest_OpenArc, EdgeIterator_NonEmpty)
{
    size_t count = 0;
    for (auto e : curve.edges())
        ++count;

    EXPECT_EQ(count, 16); // 16 segments, no duplicate edge
}

// Transformation Tests
TEST_F(ManifoldCurve2DTest_OpenArc, Transformation)
{
    const Vertex v1 = *curve.vertices_begin();
    const Vertex v2 = *(curve.vertices_begin() + 8);
    const mat3 transform = scaling_matrix_2d(2.0f);
    curve *= transform;

    auto pos1 = curve.position(v1);
    EXPECT_NEAR(pos1[0], 2.0f, COORD_TOLERANCE);
    EXPECT_NEAR(pos1[1], 0.0f, COORD_TOLERANCE);
    auto pos2 = curve.position(v2);
    EXPECT_NEAR(pos2[0], 0.0f, COORD_TOLERANCE);
    EXPECT_NEAR(pos2[1], 2.0f, COORD_TOLERANCE);
}

// Tessellation Operations
TEST_F(ManifoldCurve2DTest_OpenArc, CollapseEdge)
{
    auto v0 = curve.add_vertex(Point2(0, 0));
    auto v1 = curve.add_vertex(Point2(1, 0));
    auto e = curve.add_edge(v0, v1);

    curve.collapse_edge(e, true);
    EXPECT_TRUE(curve.is_deleted(v1));
    EXPECT_TRUE(curve.is_deleted(e));
}

TEST_F(ManifoldCurve2DTest_OpenArc, SplitEdge)
{
    auto v0 = curve.add_vertex(Point2(0, 0));
    auto v1 = curve.add_vertex(Point2(1, 0));
    auto e = curve.add_edge(v0, v1);

    auto v = curve.split_edge(e);
    EXPECT_TRUE(curve.is_valid(v));
}

TEST_F(ManifoldCurve2DTest_OpenArc, SplitEdgeWithVertex)
{
    auto v0 = curve.add_vertex(Point2(0, 0));
    auto v1 = curve.add_vertex(Point2(1, 0));
    auto e = curve.add_edge(v0, v1);
    auto newPos = Point2(0.5, 0);

    auto v = curve.split_edge_with_vertex(e, newPos);
    EXPECT_TRUE(curve.is_valid(v));
    EXPECT_EQ(curve.position(v), newPos);
}