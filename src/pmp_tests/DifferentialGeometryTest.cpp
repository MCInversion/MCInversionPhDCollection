// Copyright 2017-2019 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "gtest/gtest.h"

#include <pmp/SurfaceMesh.h>
#include <pmp/algorithms/DifferentialGeometry.h>
#include <pmp/algorithms/SurfaceFactory.h>

#include "Helpers.h"

#include <vector>

using namespace pmp;

class DifferentialGeometryTest : public ::testing::Test
{
public:
    SurfaceMesh mesh;
    static SurfaceMesh sphere;

    Vertex v0, v1, v2, v3;
    Face f0;

    Vertex central_vertex;

    void add_triangle()
    {
        v0 = mesh.add_vertex(Point(0, 0, 0));
        v1 = mesh.add_vertex(Point(1, 0, 0));
        v2 = mesh.add_vertex(Point(0, 1, 0));
        f0 = mesh.add_triangle(v0, v1, v2);
    }

    void one_ring()
    {
        mesh = vertex_onering();
        central_vertex = Vertex(3);             // the central vertex
        mesh.position(central_vertex)[2] = 0.1; // lift central vertex
    }
};

SurfaceMesh DifferentialGeometryTest::sphere = SurfaceFactory::icosphere(5);

TEST_F(DifferentialGeometryTest, triangle_area_points)
{
    add_triangle();
    Scalar area =
        triangle_area(mesh.position(v0), mesh.position(v1), mesh.position(v2));
    EXPECT_EQ(area, 0.5);
}

TEST_F(DifferentialGeometryTest, triangle_area_face)
{
    add_triangle();
    Scalar area = triangle_area(mesh, f0);
    EXPECT_EQ(area, 0.5);
}

TEST_F(DifferentialGeometryTest, barycentric_area)
{
    one_ring();
    Scalar area = barycentric_area(mesh, central_vertex);
    EXPECT_FLOAT_EQ(area, 0.024590395);
}

TEST_F(DifferentialGeometryTest, laplace)
{
    one_ring();
    auto lv = laplace_voronoi(mesh, central_vertex);
    EXPECT_GT(norm(lv), 0);
}

TEST_F(DifferentialGeometryTest, vertex_curvature)
{
    one_ring();
    auto vcurv = vertex_curvature(mesh, central_vertex);
    EXPECT_FLOAT_EQ(vcurv.mean, 6.1538467);
    EXPECT_FLOAT_EQ(vcurv.gauss, 50.860939);
    EXPECT_FLOAT_EQ(vcurv.max, 6.1538467);
    EXPECT_FLOAT_EQ(vcurv.min, 6.1538467);
}

TEST_F(DifferentialGeometryTest, surface_area)
{
    auto area = surface_area(sphere);
    EXPECT_NEAR(area, 12.57, 1.0e-2);
}

TEST_F(DifferentialGeometryTest, volume)
{
    auto v = volume(sphere);
    EXPECT_NEAR(v, 4.18, 1.0e-2);
}

TEST_F(DifferentialGeometryTest, centroid)
{
    auto center = centroid(sphere);
    EXPECT_LT(norm(center), 1e-5);
}

class DifferentialGeometryCurveTest : public ::testing::Test
{
public:
    ManifoldCurve2D curve;
    Vertex v0, v1, v;
    Edge e0, e1;

    void SetUp() override
    {
        v0 = curve.add_vertex(Point2(-1.1f, 0.0f));
        v1 = curve.add_vertex(Point2(1.15f, 0.0f));
        v = curve.add_vertex(Point2(-0.1f, 0.5f));
        e0 = curve.add_edge(v1, v);
        e1 = curve.add_edge(v, v0);
    }
};

TEST_F(DifferentialGeometryCurveTest, laplace_1D_TestCase)
{
    const auto laplaceCentral = laplace_1D(curve, v);
    const auto laplaceBoundary0 = laplace_1D(curve, v0);
    const auto laplaceBoundary1 = laplace_1D(curve, v1);
    const auto vIso = curve.add_vertex(Point2(1.0, 1.0));
    const auto laplaceIsolated = laplace_1D(curve, vIso);

    EXPECT_NEAR(laplaceCentral[0], 0.0276f, 1e-4f);
    EXPECT_NEAR(laplaceCentral[1], -0.6644f, 1e-4f);
    EXPECT_NEAR(laplaceBoundary0[0], 0.0f, 1e-4f);
    EXPECT_NEAR(laplaceBoundary0[1], 0.0f, 1e-4f);
    EXPECT_NEAR(laplaceBoundary1[0], 0.0f, 1e-4f);
    EXPECT_NEAR(laplaceBoundary1[1], 0.0f, 1e-4f);
    EXPECT_NEAR(laplaceIsolated[0], 0.0f, 1e-4f);
    EXPECT_NEAR(laplaceIsolated[1], 0.0f, 1e-4f);
}

TEST_F(DifferentialGeometryCurveTest, laplace_implicit_1D_TestCase)
{
    const auto laplaceWeightsCentral = laplace_implicit_1D(curve, v);

    EXPECT_NEAR(laplaceWeightsCentral.vertexWeights.at(v0), 0.7259f, 1e-4f);
    EXPECT_NEAR(laplaceWeightsCentral.vertexWeights.at(v1), 0.6028f, 1e-4f);
    EXPECT_NEAR(laplaceWeightsCentral.weightSum, 1.3287f, 1e-4f);
}