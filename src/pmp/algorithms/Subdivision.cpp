// Copyright 2011-2022 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "pmp/algorithms/Subdivision.h"
#include "pmp/algorithms/DifferentialGeometry.h"

namespace
{
    [[nodiscard]] size_t CountBoundaryEdges(const pmp::SurfaceMesh& mesh)
    {
        size_t result = 0;
        for (const auto e : mesh.edges())
        {
            if (!mesh.is_boundary(e))
                continue;

            result++;
        }
        return result;
    }

    [[nodiscard]] std::tuple<size_t, size_t, size_t> EstimateMeshCountsForLoopSubdivision(const pmp::SurfaceMesh& mesh, const size_t& subdivLevel)
    {
        const auto& s = subdivLevel;

        const auto nBdEdges0 = CountBoundaryEdges(mesh);

        if (nBdEdges0 == 0)
        {
            // mesh is watertight
            const auto nEdges0 = mesh.n_edges();
            const auto nVerts0 = mesh.n_vertices();
            const auto nFaces0 = mesh.n_faces();

            const size_t nEdges = (static_cast<size_t>(pow(4, s)) * nEdges0);
            const size_t nVerts = ((nEdges0 * static_cast<size_t>(pow(4, s) - 1) + 3 * nVerts0) / 3);
            const size_t nFaces = (static_cast<size_t>(pow(4, s)) * nFaces0);

            return { nEdges, nVerts, nFaces };
        }

        const auto nEdges0 = mesh.n_edges();
        const auto nVerts0 = mesh.n_vertices();
        const auto nFaces0 = mesh.n_faces();

        const auto nIntEdges0 = nEdges0 - nBdEdges0;

        const size_t nIntEdges = (pow(2, s - 1)) * ((pow(2, s) - 1) * nBdEdges0 + nIntEdges0 * pow(2, s + 1));
        const size_t nBdEdges = (pow(2, s)) * nBdEdges0;
        const size_t nVerts = nBdEdges0 * (pow(4, s) - 4 + 3 * pow(2, s)) / 6 + nIntEdges0 * (pow(4, s) - 1) / 3 + nVerts0;
        const size_t nFaces = (static_cast<size_t>(pow(4, s)) * nFaces0);

        return { nIntEdges + nBdEdges, nVerts, nFaces };
    }
	
} // anonymous namespace

namespace pmp {

Subdivision::Subdivision(SurfaceMesh& mesh) : mesh_(mesh)
{
    points_ = mesh_.vertex_property<Point>("v:point");
    vfeature_ = mesh_.get_vertex_property<bool>("v:feature");
    efeature_ = mesh_.get_edge_property<bool>("e:feature");
}

void Subdivision::catmull_clark()
{
    // reserve memory
    size_t nv = mesh_.n_vertices();
    size_t ne = mesh_.n_edges();
    size_t nf = mesh_.n_faces();
    mesh_.reserve(nv + ne + nf, 2 * ne + 4 * nf, 4 * nf);

    // get properties
    auto vpoint = mesh_.add_vertex_property<Point>("catmull:vpoint");
    auto epoint = mesh_.add_edge_property<Point>("catmull:epoint");
    auto fpoint = mesh_.add_face_property<Point>("catmull:fpoint");

    // compute face vertices
    for (auto f : mesh_.faces())
    {
        fpoint[f] = centroid(mesh_, f);
    }

    // compute edge vertices
    for (auto e : mesh_.edges())
    {
        // boundary or feature edge?
        if (mesh_.is_boundary(e) || (efeature_ && efeature_[e]))
        {
            epoint[e] = 0.5 * (points_[mesh_.vertex(e, 0)] +
                                points_[mesh_.vertex(e, 1)]);
        }

        // interior edge
        else
        {
            Point p(0, 0, 0);
            p += points_[mesh_.vertex(e, 0)];
            p += points_[mesh_.vertex(e, 1)];
            p += fpoint[mesh_.face(e, 0)];
            p += fpoint[mesh_.face(e, 1)];
            p *= 0.25;
            epoint[e] = p;
        }
    }

    // compute new positions for old vertices
    for (auto v : mesh_.vertices())
    {
        // isolated vertex?
        if (mesh_.is_isolated(v))
        {
            vpoint[v] = points_[v];
        }

        // boundary vertex?
        else if (mesh_.is_boundary(v))
        {
            auto h1 = mesh_.halfedge(v);
            auto h0 = mesh_.prev_halfedge(h1);

            Point p = points_[v];
            p *= 6.0;
            p += points_[mesh_.to_vertex(h1)];
            p += points_[mesh_.from_vertex(h0)];
            p *= 0.125;

            vpoint[v] = p;
        }

        // interior feature vertex?
        else if (vfeature_ && vfeature_[v])
        {
            Point p = points_[v];
            p *= 6.0;
            int count(0);

            for (auto h : mesh_.halfedges(v))
            {
                if (efeature_[mesh_.edge(h)])
                {
                    p += points_[mesh_.to_vertex(h)];
                    ++count;
                }
            }

            if (count == 2) // vertex is on feature edge
            {
                p *= 0.125;
                vpoint[v] = p;
            }
            else // keep fixed
            {
                vpoint[v] = points_[v];
            }
        }

        // interior vertex
        else
        {
            // weights from SIGGRAPH paper "Subdivision Surfaces in Character Animation"

            const Scalar k = mesh_.valence(v);
            Point p(0, 0, 0);

            for (auto vv : mesh_.vertices(v))
                p += points_[vv];

            for (auto f : mesh_.faces(v))
                p += fpoint[f];

            p /= (k * k);

            p += ((k - 2.0) / k) * points_[v];

            vpoint[v] = p;
        }
    }

    // assign new positions to old vertices
    for (auto v : mesh_.vertices())
    {
        points_[v] = vpoint[v];
    }

    // split edges
    for (auto e : mesh_.edges())
    {
        // feature edge?
        if (efeature_ && efeature_[e])
        {
            auto h = mesh_.insert_vertex(e, epoint[e]);
            auto v = mesh_.to_vertex(h);
            auto e0 = mesh_.edge(h);
            auto e1 = mesh_.edge(mesh_.next_halfedge(h));

            vfeature_[v] = true;
            efeature_[e0] = true;
            efeature_[e1] = true;
        }

        // normal edge
        else
        {
            mesh_.insert_vertex(e, epoint[e]);
        }
    }

    // split faces
    for (auto f : mesh_.faces())
    {
        auto h0 = mesh_.halfedge(f);
        mesh_.insert_edge(h0, mesh_.next_halfedge(mesh_.next_halfedge(h0)));

        auto h1 = mesh_.next_halfedge(h0);
        mesh_.insert_vertex(mesh_.edge(h1), fpoint[f]);

        auto h =
            mesh_.next_halfedge(mesh_.next_halfedge(mesh_.next_halfedge(h1)));
        while (h != h0)
        {
            mesh_.insert_edge(h1, h);
            h = mesh_.next_halfedge(
                mesh_.next_halfedge(mesh_.next_halfedge(h1)));
        }
    }

    // clean-up properties
    mesh_.remove_vertex_property(vpoint);
    mesh_.remove_edge_property(epoint);
    mesh_.remove_face_property(fpoint);
}

void Subdivision::loop()
{
    if (!mesh_.is_triangle_mesh())
    {
        auto what = "Subdivision: Not a triangle mesh.";
        throw InvalidInputException(what);
    }

    // reserve memory
    size_t nv = mesh_.n_vertices();
    size_t ne = mesh_.n_edges();
    size_t nf = mesh_.n_faces();
    mesh_.reserve(nv + ne, 2 * ne + 3 * nf, 4 * nf);

    // add properties
    auto vpoint = mesh_.add_vertex_property<Point>("loop:vpoint");
    auto epoint = mesh_.add_edge_property<Point>("loop:epoint");

    // compute vertex positions
    for (auto v : mesh_.vertices())
    {
        // isolated vertex?
        if (mesh_.is_isolated(v))
        {
            vpoint[v] = points_[v];
        }

        // boundary vertex?
        else if (mesh_.is_boundary(v))
        {
            auto h1 = mesh_.halfedge(v);
            auto h0 = mesh_.prev_halfedge(h1);

            Point p = points_[v];
            p *= 6.0;
            p += points_[mesh_.to_vertex(h1)];
            p += points_[mesh_.from_vertex(h0)];
            p *= 0.125;
            vpoint[v] = p;
        }

        // interior feature vertex?
        else if (vfeature_ && vfeature_[v])
        {
            Point p = points_[v];
            p *= 6.0;
            int count(0);

            for (auto h : mesh_.halfedges(v))
            {
                if (efeature_[mesh_.edge(h)])
                {
                    p += points_[mesh_.to_vertex(h)];
                    ++count;
                }
            }

            if (count == 2) // vertex is on feature edge
            {
                p *= 0.125;
                vpoint[v] = p;
            }
            else // keep fixed
            {
                vpoint[v] = points_[v];
            }
        }

        // interior vertex
        else
        {
            Point p(0, 0, 0);
            Scalar k(0);

            for (auto vv : mesh_.vertices(v))
            {
                p += points_[vv];
                ++k;
            }
            p /= k;

            Scalar beta =
                (0.625 - pow(0.375 + 0.25 * std::cos(2.0 * M_PI / k), 2.0));

            vpoint[v] = points_[v] * (Scalar)(1.0 - beta) + beta * p;
        }
    }

    // compute edge positions
    for (auto e : mesh_.edges())
    {
        // boundary or feature edge?
        if (mesh_.is_boundary(e) || (efeature_ && efeature_[e]))
        {
            epoint[e] =
                (points_[mesh_.vertex(e, 0)] + points_[mesh_.vertex(e, 1)]) *
                Scalar(0.5);
        }

        // interior edge
        else
        {
            auto h0 = mesh_.halfedge(e, 0);
            auto h1 = mesh_.halfedge(e, 1);
            Point p = points_[mesh_.to_vertex(h0)];
            p += points_[mesh_.to_vertex(h1)];
            p *= 3.0;
            p += points_[mesh_.to_vertex(mesh_.next_halfedge(h0))];
            p += points_[mesh_.to_vertex(mesh_.next_halfedge(h1))];
            p *= 0.125;
            epoint[e] = p;
        }
    }

    // set new vertex positions
    for (auto v : mesh_.vertices())
    {
        points_[v] = vpoint[v];
    }

    // insert new vertices on edges
    for (auto e : mesh_.edges())
    {
        // feature edge?
        if (efeature_ && efeature_[e])
        {
            auto h = mesh_.insert_vertex(e, epoint[e]);
            auto v = mesh_.to_vertex(h);
            auto e0 = mesh_.edge(h);
            auto e1 = mesh_.edge(mesh_.next_halfedge(h));

            vfeature_[v] = true;
            efeature_[e0] = true;
            efeature_[e1] = true;
        }

        // normal edge
        else
        {
            mesh_.insert_vertex(e, epoint[e]);
        }
    }

    // split faces
    Halfedge h;
    for (auto f : mesh_.faces())
    {
        h = mesh_.halfedge(f);
        mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
        h = mesh_.next_halfedge(h);
        mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
        h = mesh_.next_halfedge(h);
        mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
    }

    // clean-up properties
    mesh_.remove_vertex_property(vpoint);
    mesh_.remove_edge_property(epoint);
}

void Subdivision::loop_prealloc(const size_t& steps)
{
    if (!mesh_.is_triangle_mesh())
    {
        auto what = "Subdivision: Not a triangle mesh.";
        throw InvalidInputException(what);
    }

    // return { nEdges, nVerts, nFaces };
    const auto [ne_final, nv_final, nf_final] = EstimateMeshCountsForLoopSubdivision(mesh_, steps);
    mesh_.reserve(nv_final, ne_final, nf_final);

    for (size_t s = 0; s < steps; s++)
    {
        // add properties
        auto vpoint = mesh_.add_vertex_property<Point>("loop:vpoint");
        auto epoint = mesh_.add_edge_property<Point>("loop:epoint");

        // compute vertex positions
        for (auto v : mesh_.vertices())
        {
            // isolated vertex?
            if (mesh_.is_isolated(v))
            {
                vpoint[v] = points_[v];
            }

            // boundary vertex?
            else if (mesh_.is_boundary(v))
            {
                auto h1 = mesh_.halfedge(v);
                auto h0 = mesh_.prev_halfedge(h1);

                Point p = points_[v];
                p *= 6.0;
                p += points_[mesh_.to_vertex(h1)];
                p += points_[mesh_.from_vertex(h0)];
                p *= 0.125;
                vpoint[v] = p;
            }

            // interior feature vertex?
            else if (vfeature_ && vfeature_[v])
            {
                Point p = points_[v];
                p *= 6.0;
                int count(0);

                for (auto h : mesh_.halfedges(v))
                {
                    if (efeature_[mesh_.edge(h)])
                    {
                        p += points_[mesh_.to_vertex(h)];
                        ++count;
                    }
                }

                if (count == 2) // vertex is on feature edge
                {
                    p *= 0.125;
                    vpoint[v] = p;
                }
                else // keep fixed
                {
                    vpoint[v] = points_[v];
                }
            }

            // interior vertex
            else
            {
                Point p(0, 0, 0);
                Scalar k(0);

                for (auto vv : mesh_.vertices(v))
                {
                    p += points_[vv];
                    ++k;
                }
                p /= k;

                Scalar beta =
                    (0.625 - pow(0.375 + 0.25 * std::cos(2.0 * M_PI / k), 2.0));

                vpoint[v] = points_[v] * (Scalar)(1.0 - beta) + beta * p;
            }
        }

        // compute edge positions
        for (auto e : mesh_.edges())
        {
            // boundary or feature edge?
            if (mesh_.is_boundary(e) || (efeature_ && efeature_[e]))
            {
                epoint[e] =
                    (points_[mesh_.vertex(e, 0)] + points_[mesh_.vertex(e, 1)]) *
                    Scalar(0.5);
            }

            // interior edge
            else
            {
                auto h0 = mesh_.halfedge(e, 0);
                auto h1 = mesh_.halfedge(e, 1);
                Point p = points_[mesh_.to_vertex(h0)];
                p += points_[mesh_.to_vertex(h1)];
                p *= 3.0;
                p += points_[mesh_.to_vertex(mesh_.next_halfedge(h0))];
                p += points_[mesh_.to_vertex(mesh_.next_halfedge(h1))];
                p *= 0.125;
                epoint[e] = p;
            }
        }

        // set new vertex positions
        for (auto v : mesh_.vertices())
        {
            points_[v] = vpoint[v];
        }

        // insert new vertices on edges
        for (auto e : mesh_.edges())
        {
            // feature edge?
            if (efeature_ && efeature_[e])
            {
                auto h = mesh_.insert_vertex(e, epoint[e]);
                auto v = mesh_.to_vertex(h);
                auto e0 = mesh_.edge(h);
                auto e1 = mesh_.edge(mesh_.next_halfedge(h));

                vfeature_[v] = true;
                efeature_[e0] = true;
                efeature_[e1] = true;
            }

            // normal edge
            else
            {
                mesh_.insert_vertex(e, epoint[e]);
            }
        }

        // split faces
        Halfedge h;
        for (auto f : mesh_.faces())
        {
            h = mesh_.halfedge(f);
            mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
            h = mesh_.next_halfedge(h);
            mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
            h = mesh_.next_halfedge(h);
            mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
        }

        // clean-up properties
        mesh_.remove_vertex_property(vpoint);
        mesh_.remove_edge_property(epoint);
    }
}

void Subdivision::quad_tri()
{
    // split each edge evenly into two parts
    for (auto e : mesh_.edges())
    {
        mesh_.insert_vertex(e, 0.5 * (points_[mesh_.vertex(e, 0)] +
                                       points_[mesh_.vertex(e, 1)]));
    }

    // subdivide faces without repositioning
    for (auto f : mesh_.faces())
    {
        size_t f_val = mesh_.valence(f) / 2;
        if (f_val == 3)
        {
            // face was a triangle
            Halfedge h0 = mesh_.halfedge(f);
            Halfedge h1 = mesh_.next_halfedge(mesh_.next_halfedge(h0));
            mesh_.insert_edge(h0, h1);

            h0 = mesh_.next_halfedge(h0);
            h1 = mesh_.next_halfedge(mesh_.next_halfedge(h0));
            mesh_.insert_edge(h0, h1);

            h0 = mesh_.next_halfedge(h0);
            h1 = mesh_.next_halfedge(mesh_.next_halfedge(h0));
            mesh_.insert_edge(h0, h1);
        }
        else
        {
            // quadrangulate the rest
            Halfedge h0 = mesh_.halfedge(f);
            Halfedge h1 = mesh_.next_halfedge(mesh_.next_halfedge(h0));
            h1 = mesh_.insert_edge(h0, h1);
            mesh_.insert_vertex(mesh_.edge(h1), centroid(mesh_, f));

            auto h = mesh_.next_halfedge(
                mesh_.next_halfedge(mesh_.next_halfedge(h1)));
            while (h != h0)
            {
                mesh_.insert_edge(h1, h);
                h = mesh_.next_halfedge(
                    mesh_.next_halfedge(mesh_.next_halfedge(h1)));
            }
        }
    }

    auto new_pos =
        mesh_.add_vertex_property<Point>("quad_tri:new_position", Point(0));

    for (auto v : mesh_.vertices())
    {
        if (mesh_.is_boundary(v))
        {
            new_pos[v] = 0.5 * points_[v];

            // add neighbouring vertices on boundary
            for (auto vv : mesh_.vertices(v))
            {
                if (mesh_.is_boundary(vv))
                {
                    new_pos[v] += 0.25 * points_[vv];
                }
            }
        }
        else
        {
            // count the number of faces and quads surrounding the vertex
            int n_faces = 0;
            int n_quads = 0;
            for (auto f : mesh_.faces(v))
            {
                n_faces++;
                if (mesh_.valence(f) == 4)
                    n_quads++;
            }

            if (n_quads == 0)
            {
                // vertex is surrounded only by triangles
                double a =
                    2.0 * pow(3.0 / 8.0 +
                                  (std::cos(2.0 * M_PI / n_faces) - 1.0) / 4.0,
                              2.0);
                double b = (1.0 - a) / n_faces;

                new_pos[v] = a * points_[v];
                for (auto vv : mesh_.vertices(v))
                {
                    new_pos[v] += b * points_[vv];
                }
            }
            else if (n_quads == n_faces)
            {
                // vertex is surrounded only by quads
                double c = (n_faces - 3.0) / n_faces;
                double d = 2.0 / pow(n_faces, 2.0);
                double e = 1.0 / pow(n_faces, 2.0);

                new_pos[v] = c * points_[v];
                for (auto h : mesh_.halfedges(v))
                {
                    new_pos[v] += d * points_[mesh_.to_vertex(h)];
                    new_pos[v] +=
                        e * points_[mesh_.to_vertex(mesh_.next_halfedge(h))];
                }
            }
            else
            {
                // vertex is surrounded by triangles and quads
                double alpha = 1.0 / (1.0 + 0.5 * n_faces + 0.25 * n_quads);
                double beta = 0.5 * alpha;
                double gamma = 0.25 * alpha;

                new_pos[v] = alpha * points_[v];
                for (auto h : mesh_.halfedges(v))
                {
                    new_pos[v] += beta * points_[mesh_.to_vertex(h)];
                    if (mesh_.valence(mesh_.face(h)) == 4)
                    {
                        new_pos[v] +=
                            gamma *
                            points_[mesh_.to_vertex(mesh_.next_halfedge(h))];
                    }
                }
            }
        }
    }

    // apply new positions to the mesh
    for (auto v : mesh_.vertices())
    {
        points_[v] = new_pos[v];
    }

    mesh_.remove_vertex_property(new_pos);
}

} // namespace pmp
