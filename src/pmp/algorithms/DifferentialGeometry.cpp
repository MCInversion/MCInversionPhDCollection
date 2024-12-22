// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "pmp/algorithms/DifferentialGeometry.h"

#include <cmath>
#include <limits>

namespace pmp {

Scalar triangle_area(const Point& p0, const Point& p1, const Point& p2)
{
    return Scalar(0.5) * norm(cross(p1 - p0, p2 - p0));
}

Scalar triangle_area(const SurfaceMesh& mesh, Face f)
{
    assert(mesh.valence(f) == 3);

    auto fv = mesh.vertices(f);
    const auto& p0 = mesh.position(*fv);
    const auto& p1 = mesh.position(*(++fv));
    const auto& p2 = mesh.position(*(++fv));

    return triangle_area(p0, p1, p2);
}

Scalar surface_area(const SurfaceMesh& mesh)
{
    if (!mesh.is_triangle_mesh())
    {
        throw InvalidInputException("Input is not a pure triangle mesh!");
    }
    Scalar area(0);
    for (auto f : mesh.faces())
    {
        area += triangle_area(mesh, f);
    }
    return area;
}

Scalar volume(const SurfaceMesh& mesh)
{
    if (!mesh.is_triangle_mesh())
    {
        throw InvalidInputException("Input is not a pure triangle mesh!");
    }

    Scalar volume(0);
    for (const auto f : mesh.faces())
    {
        auto fv = mesh.vertices(f);
        const auto& p0 = mesh.position(*fv);
        const auto& p1 = mesh.position(*(++fv));
        const auto& p2 = mesh.position(*(++fv));

        volume += Scalar(1.0) / Scalar(6.0) * dot(cross(p0, p1), p2);
    }

    return std::abs(volume);
}

Point centroid(const SurfaceMesh& mesh, Face f)
{
    Point c(0, 0, 0);
    Scalar n(0);
    for (auto v : mesh.vertices(f))
    {
        c += mesh.position(v);
        ++n;
    }
    c /= n;
    return c;
}

Point2 centroid(const ManifoldCurve2D& curve, Edge e)
{
    const auto [v0, v1] = curve.vertices(e);
    return 0.5 * (curve.position(v0) + curve.position(v1));
}

Point centroid(const SurfaceMesh& mesh)
{
    Point center(0, 0, 0), c;
    Scalar area(0), a;
    for (auto f : mesh.faces())
    {
        a = triangle_area(mesh, f);
        c = centroid(mesh, f);
        area += a;
        center += a * c;
    }
    center /= area;
    return center;
}

void dual(SurfaceMesh& mesh)
{
    // the new dualized mesh
    SurfaceMesh tmp;

    // remember new vertices per face
    auto fvertex = mesh.add_face_property<Vertex>("f:vertex");

    // add centroid for each face
    for (auto f : mesh.faces())
        fvertex[f] = tmp.add_vertex(centroid(mesh, f));

    // add new face for each vertex
    for (auto v : mesh.vertices())
    {
        std::vector<Vertex> vertices;
        for (auto f : mesh.faces(v))
            vertices.push_back(fvertex[f]);

        tmp.add_face(vertices);
    }

    // swap old and new meshes, don't copy properties
    mesh.assign(tmp);
}

double cotan_weight(const SurfaceMesh& mesh, Edge e)
{
    double weight = 0.0;

    const Halfedge h0 = mesh.halfedge(e, 0);
    const Halfedge h1 = mesh.halfedge(e, 1);

    const dvec3 p0 = (dvec3)mesh.position(mesh.to_vertex(h0));
    const dvec3 p1 = (dvec3)mesh.position(mesh.to_vertex(h1));

    if (!mesh.is_boundary(h0))
    {
        const dvec3 p2 =
            (dvec3)mesh.position(mesh.to_vertex(mesh.next_halfedge(h0)));
        const dvec3 d0 = p0 - p2;
        const dvec3 d1 = p1 - p2;
        const double area = norm(cross(d0, d1));
        if (area > std::numeric_limits<double>::min())
        {
            const double cot = dot(d0, d1) / area;
            weight += clamp_cot(cot);
        }
    }

    if (!mesh.is_boundary(h1))
    {
        const dvec3 p2 =
            (dvec3)mesh.position(mesh.to_vertex(mesh.next_halfedge(h1)));
        const dvec3 d0 = p0 - p2;
        const dvec3 d1 = p1 - p2;
        const double area = norm(cross(d0, d1));
        if (area > std::numeric_limits<double>::min())
        {
            const double cot = dot(d0, d1) / area;
            weight += clamp_cot(cot);
        }
    }

    assert(!std::isnan(weight));
    assert(!std::isinf(weight));

    return weight;
}

double voronoi_area(const SurfaceMesh& mesh, Vertex v)
{
    double area(0.0);

    if (!mesh.is_isolated(v))
    {
        Halfedge h0, h1, h2;
        dvec3 p, q, r, pq, qr, pr;
        double dotp, dotq, dotr;
        double cotq, cotr;

        for (auto h : mesh.halfedges(v))
        {
            h0 = h;
            h1 = mesh.next_halfedge(h0);
            h2 = mesh.next_halfedge(h1);

            if (mesh.is_boundary(h0))
                continue;

            // three vertex positions
            p = (dvec3)mesh.position(mesh.to_vertex(h2));
            q = (dvec3)mesh.position(mesh.to_vertex(h0));
            r = (dvec3)mesh.position(mesh.to_vertex(h1));

            // edge vectors
            (pq = q) -= p;
            (qr = r) -= q;
            (pr = r) -= p;

            // compute and check triangle area
            const auto triangle_area = norm(cross(pq, pr));
            if (triangle_area <= std::numeric_limits<double>::min())
                continue;

            // dot products for each corner (of its two emanating edge vectors)
            dotp = dot(pq, pr);
            dotq = -dot(qr, pq);
            dotr = dot(qr, pr);

            // angle at p is obtuse
            if (dotp < 0.0)
            {
                area += 0.25 * triangle_area;
            }

            // angle at q or r obtuse
            else if (dotq < 0.0 || dotr < 0.0)
            {
                area += 0.125 * triangle_area;
            }

            // no obtuse angles
            else
            {
                // cot(angle) = cos(angle)/sin(angle) = dot(A,B)/norm(cross(A,B))
                cotq = dotq / triangle_area;
                cotr = dotr / triangle_area;

                // clamp cot(angle) by clamping angle to [3, 177]
                area += 0.125 * (sqrnorm(pr) * clamp_cot(cotq) +
                                 sqrnorm(pq) * clamp_cot(cotr));
            }
        }
    }

    assert(!std::isnan(area));
    assert(!std::isinf(area));

    return area;
}

double barycentric_area(const SurfaceMesh& mesh, Vertex v)
{
    double area(0.0);

    if (!mesh.is_isolated(v))
    {
        const auto p = mesh.position(v);

        for (auto h : mesh.halfedges(v))
        {
            if (mesh.is_boundary(h))
                continue;

            auto h0 = h;
            auto h1 = mesh.next_halfedge(h0);

            auto pq = mesh.position(mesh.to_vertex(h0));
            pq -= p;

            auto pr = mesh.position(mesh.to_vertex(h1));
            pr -= p;

            area += norm(cross(pq, pr)) / 6.0;
        }
    }

    assert(!std::isnan(area));
    assert(!std::isinf(area));

    return area;
}

double midpoint_covolume_length(const ManifoldCurve2D& curve, Vertex v)
{
    if (curve.is_isolated(v) || curve.is_boundary(v))
        return 0.0;

    const auto [eTo, eFrom] = curve.edges(v);
    const double l0 = curve.edge_length(eTo);    // Length from v_prev to v
    const double l1 = curve.edge_length(eFrom);  // Length from v to v_next
    return 0.5 * (l0 + l1);
}

Point laplace_voronoi(const SurfaceMesh& mesh, Vertex v)
{
    Point laplace(0.0, 0.0, 0.0);

    if (!mesh.is_isolated(v))
    {
        Scalar sum_weights(0.0);

        for (auto h : mesh.halfedges(v))
        {
            const auto weight = cotan_weight(mesh, mesh.edge(h));
            sum_weights += weight;
            laplace += weight * mesh.position(mesh.to_vertex(h));
        }

        laplace -= sum_weights * mesh.position(v);
        laplace /= Scalar(2.0) * voronoi_area(mesh, v);
    }

    return laplace;
}

Point laplace_barycentric(const SurfaceMesh& mesh, Vertex v)
{
    Point laplace(0.0, 0.0, 0.0);

    if (!mesh.is_isolated(v))
    {
        Scalar sum_weights(0.0);

        for (auto h : mesh.halfedges(v))
        {
            const auto weight = cotan_weight(mesh, mesh.edge(h));
            sum_weights += weight;
            laplace += weight * mesh.position(mesh.to_vertex(h));
        }

        laplace -= sum_weights * mesh.position(v);
        laplace /= Scalar(2.0) * barycentric_area(mesh, v);
    }

    return laplace;
}

ImplicitLaplaceInfo laplace_implicit_voronoi(const SurfaceMesh& mesh, Vertex v)
{
    ImplicitLaplaceInfo result{};
    if (mesh.is_isolated(v))
        return result;

    Scalar sum_weights(0.0);

    for (const auto h : mesh.halfedges(v))
    {
        const auto weight = static_cast<Scalar>(cotan_weight(mesh, mesh.edge(h)));
        sum_weights += weight;
        result.vertexWeights[mesh.to_vertex(h)] = weight;
    }

    const auto area = static_cast<Scalar>(voronoi_area(mesh, v));
    const Scalar areaNorm = (2.0 * area);
    for (auto& [v, weight] : result.vertexWeights)
        weight /= areaNorm;

    result.weightSum = sum_weights / areaNorm;

    return result;
}

ImplicitLaplaceInfo laplace_implicit_barycentric(const SurfaceMesh& mesh, Vertex v)
{
    ImplicitLaplaceInfo result{};
    if (mesh.is_isolated(v))
        return result;

    Scalar sum_weights(0.0);

    for (const auto h : mesh.halfedges(v))
    {
        const auto weight = static_cast<Scalar>(cotan_weight(mesh, mesh.edge(h)));
        sum_weights += weight;
        result.vertexWeights[mesh.to_vertex(h)] = weight;
    }

    const auto area = static_cast<Scalar>(barycentric_area(mesh, v));
    const Scalar areaNorm = (2.0 * area);
    for (auto& [v, weight] : result.vertexWeights)
        weight /= areaNorm;

    result.weightSum = sum_weights / areaNorm;

    return result;
}

ImplicitLaplaceInfo laplace_implicit_1D(const ManifoldCurve2D& curve, Vertex v)
{
    ImplicitLaplaceInfo result{};
    if (curve.is_isolated(v) || curve.is_boundary(v))
        return result;

    const auto [eTo, eFrom] = curve.edges(v);
    const Scalar l0 = curve.edge_length(eTo);    // Length from v_prev to v
    const Scalar l1 = curve.edge_length(eFrom);  // Length from v to v_next

    const auto vNext = curve.to_vertex(eFrom);     // Next vertex (v+1)
    const auto vPrev = curve.from_vertex(eTo);     // Previous vertex (v-1)

    const Scalar h = l0 + l1;  // Effective double length of the control volume

    result.vertexWeights[vPrev] = 2.0 / (h * l0);
    result.vertexWeights[vNext] = 2.0 / (h * l1);
    result.weightSum = 2.0 / h * (1.0 / l0 + 1.0 / l1);

    return result;
}

Point2 laplace_1D(const ManifoldCurve2D& curve, Vertex v)
{
    Point2 laplace(0.0, 0.0);
    if (curve.is_isolated(v) || curve.is_boundary(v))
        return laplace;

    const auto [eTo, eFrom] = curve.edges(v);
    const Scalar l0 = curve.edge_length(eTo);    // Length from v_prev to v
    const Scalar l1 = curve.edge_length(eFrom);  // Length from v to v_next

    const auto vNext = curve.to_vertex(eFrom);     // Next vertex (v+1)
    const auto vPrev = curve.from_vertex(eTo);     // Previous vertex (v-1)

    const Scalar h = l0 + l1;  // Effective double length of the control volume

    laplace += (2.0 / h) * ((curve.position(vNext) - curve.position(v)) / l1 - (curve.position(v) - curve.position(vPrev)) / l0);

    return laplace;
}

Scalar angle_sum(const SurfaceMesh& mesh, Vertex v)
{
    Scalar angles(0.0);

    if (!mesh.is_boundary(v))
    {
        const Point& p0 = mesh.position(v);

        for (auto h : mesh.halfedges(v))
        {
            const Point& p1 = mesh.position(mesh.to_vertex(h));
            const Point& p2 =
                mesh.position(mesh.to_vertex(mesh.ccw_rotated_halfedge(h)));

            const Point p01 = normalize(p1 - p0);
            const Point p02 = normalize(p2 - p0);

            Scalar cos_angle = clamp_cos(dot(p01, p02));

            angles += acos(cos_angle);
        }
    }

    return angles;
}

VertexCurvature vertex_curvature(const SurfaceMesh& mesh, Vertex v)
{
    VertexCurvature c;

    const Scalar area = voronoi_area(mesh, v);
    if (area > std::numeric_limits<Scalar>::min())
    {
        c.mean = Scalar(0.5) * norm(laplace_voronoi(mesh, v));
        c.gauss = (2.0 * M_PI - angle_sum(mesh, v)) / area;

        const Scalar s = sqrt(std::max(Scalar(0.0), c.mean * c.mean - c.gauss));
        c.min = c.mean - s;
        c.max = c.mean + s;

        assert(!std::isnan(c.mean));
        assert(!std::isnan(c.gauss));
        assert(!std::isinf(c.mean));
        assert(!std::isinf(c.gauss));

        assert(c.min <= c.mean);
        assert(c.mean <= c.max);
    }

    return c;
}

Scalar vertex_curvature(const ManifoldCurve2D& curve, Vertex v)
{
    if (curve.is_isolated(v) || curve.is_boundary(v))
        return 0.0;

    const auto [eTo, eFrom] = curve.edges(v);
    const auto l0 = curve.edge_length(eTo);
    const auto l1 = curve.edge_length(eFrom);

    const auto [vPrev, vNext] = curve.vertices(v);
    const auto p = curve.position(v);
    const auto pPrev = curve.position(vPrev);
    const auto pNext = curve.position(vNext);

    // Calculate the area of triangle (pPrev, p, pNext) using the determinant method
    const Scalar area = (pPrev[0] * (p[1] - pNext[1]) +
        p[0] * (pNext[1] - pPrev[1]) +
        pNext[0] * (pPrev[1] - p[1])) * 0.5;

    const Scalar l2 = norm(pNext - pPrev); // triangle hypotenuse

    Scalar curvature = 0.0;
    if (l0 > 0 && l1 > 0 && l2 > 0 && std::fabs(area) > 0.0)
    {
        // Calculate circumradius
        const Scalar R = (l0 * l1 * l2) / (4.0 * std::fabs(area));
        curvature = 1.0 / R;

        // Determine the sign of the curvature based on the area sign
        if (area < 0.0)
        {
            curvature = -curvature; // Concave vertex
        }
        else
        {
            curvature = curvature; // Convex vertex (positive curvature)
        }
    }

    return curvature;
}
} // namespace pmp
