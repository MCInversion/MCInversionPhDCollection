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

TEST_F(ManifoldCurve2DTest_ClosedArc, VerifyConnectivity)
{
    // Get the original vertex and its adjacent edges
    const Vertex originalVertex = *curve.vertices_begin();
    const Edge originalEdgeFrom = curve.edge_from(originalVertex);
    const Edge originalEdgeTo = curve.edge_to(originalVertex);

    // Iterate through the circle and verify connectivity
    Vertex currentVertex = originalVertex;
    Edge currentEdge = originalEdgeFrom;
    size_t count = 0;

    do
    {
        // Verify the start and end vertices of the current edge
        EXPECT_EQ(curve.from_vertex(currentEdge), currentVertex);
        Vertex nextVertex = curve.to_vertex(currentEdge);
        EXPECT_NE(currentVertex, nextVertex);

        // Move to the next vertex and edge
        currentVertex = nextVertex;
        currentEdge = curve.edge_from(currentVertex);
        count++;
    } while (count < curve.n_vertices());

    // Ensure we completed the full circle
    EXPECT_EQ(currentVertex, originalVertex);
    EXPECT_EQ(curve.edge_from(currentVertex), originalEdgeFrom);
    EXPECT_EQ(curve.edge_to(currentVertex), originalEdgeTo);
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

// Memory Management
TEST_F(ManifoldCurve2DTest_ClosedArc, GarbageCollection)
{
    const auto v0 = *curve.vertices_begin();
    const auto v0Prev = Vertex(curve.from_vertex(curve.edge_to(v0)).idx() - 1); // there will be one less vertex
    const auto v0Next = curve.to_vertex(curve.edge_from(v0));
    const auto v1 = *(curve.vertices_begin() + 8);
    const auto v2 = *(curve.vertices_begin() + 9);

    curve.delete_vertex(v0);
    //curve.delete_edge(curve.edge_from(v1));

    curve.garbage_collection();
    EXPECT_EQ(curve.n_vertices(), 31);
    //EXPECT_EQ(curve.n_edges(), 30);
    EXPECT_EQ(curve.n_edges(), 31);
    EXPECT_EQ(curve.edge_from(v0Prev), curve.edge_to(v0Next));
    //EXPECT_TRUE(curve.is_boundary(v1));
    //EXPECT_TRUE(curve.is_boundary(v2));
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