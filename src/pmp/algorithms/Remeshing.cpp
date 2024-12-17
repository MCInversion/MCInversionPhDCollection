// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "pmp/algorithms/Remeshing.h"

#include <cmath>

#include <algorithm>
#include <stdexcept>

#include "pmp/algorithms/TriangleKdTree.h"
#include "pmp/algorithms/Curvature.h"
#include "pmp/algorithms/Normals.h"
#include "pmp/algorithms/BarycentricCoordinates.h"
#include "pmp/algorithms/DifferentialGeometry.h"
#include "pmp/Timer.h"

namespace pmp {

Remeshing::Remeshing(SurfaceMesh& mesh)
    : mesh_(mesh), refmesh_(nullptr), kd_tree_(nullptr)
{
    if (!mesh_.is_triangle_mesh())
        throw InvalidInputException("Input is not a pure triangle mesh!");

    points_ = mesh_.vertex_property<Point>("v:point");

    Normals::compute_vertex_normals(mesh_);
    vnormal_ = mesh_.vertex_property<Point>("v:normal");

    has_feature_vertices_ = mesh_.has_vertex_property("v:feature");
    has_feature_edges_ = mesh_.has_edge_property("e:feature");
}

void Remeshing::uniform_remeshing(Scalar edge_length, unsigned int iterations,
                                  bool use_projection)
{
    uniform_ = true;
    use_projection_ = use_projection;
    target_edge_length_ = edge_length;

    preprocessing();

    for (unsigned int i = 0; i < iterations; ++i)
    {
        split_long_edges();

        Normals::compute_vertex_normals(mesh_);

        collapse_short_edges();

        flip_edges();

        tangential_smoothing(5);
    }

    remove_caps();

    postprocessing();
}

void Remeshing::adaptive_remeshing(const AdaptiveRemeshingSettings& settings, const bool& is_convex_hull)
{
    uniform_ = false;
    min_edge_length_ = settings.MinEdgeLength;
    max_edge_length_ = settings.MaxEdgeLength;
    approx_error_ = settings.ApproxError;
    use_projection_ = settings.UseProjection;

    preprocessing();

    for (unsigned int i = 0; i < settings.NRemeshingIterations; ++i)
    {
        split_long_edges();

        Normals::compute_vertex_normals(mesh_);

        collapse_short_edges();

        flip_edges();

        tangential_smoothing(settings.NTangentialSmoothingIters);
    }

    remove_caps();

    if (!is_convex_hull)
		postprocessing();
}

void Remeshing::convex_hull_adaptive_remeshing(const AdaptiveRemeshingSettings& settings)
{
    max_edge_length_ = settings.MaxEdgeLength; // needed for the initial split_long_edges
    convex_hull_preprocessing();
    split_long_edges(3);

    // Now proceed with general adaptive remeshing
    adaptive_remeshing(settings, true);

    // marks all locked vertices as feature
    convex_hull_postprocessing();
}


void Remeshing::preprocessing()
{
    // properties
    if (!vfeature_)
		vfeature_ = mesh_.vertex_property<bool>("v:feature", false);
    if (!efeature_)
		efeature_ = mesh_.edge_property<bool>("e:feature", false);
    if (!vlocked_)
		vlocked_ = mesh_.add_vertex_property<bool>("v:locked", false);
    if (!elocked_)
		elocked_ = mesh_.add_edge_property<bool>("e:locked", false);
    if (!vsizing_)
		vsizing_ = mesh_.add_vertex_property<Scalar>("v:sizing");

    // lock unselected vertices if some vertices are selected
    auto vselected = mesh_.get_vertex_property<bool>("v:selected");
    if (vselected)
    {
        bool has_selection = false;
        for (auto v : mesh_.vertices())
        {
            if (vselected[v])
            {
                has_selection = true;
                break;
            }
        }

        if (has_selection)
        {
            for (auto v : mesh_.vertices())
            {
                vlocked_[v] = !vselected[v];
            }

            // lock an edge if one of its vertices is locked
            for (auto e : mesh_.edges())
            {
                elocked_[e] = (vlocked_[mesh_.vertex(e, 0)] ||
                               vlocked_[mesh_.vertex(e, 1)]);
            }
        }
    }

    // lock feature corners
    for (auto v : mesh_.vertices())
    {
        if (vfeature_[v])
        {
            int c = 0;
            for (auto h : mesh_.halfedges(v))
                if (efeature_[mesh_.edge(h)])
                    ++c;

            if (c != 2)
                vlocked_[v] = true;
        }
    }

    // compute sizing field
    if (uniform_)
    {
        for (auto v : mesh_.vertices())
        {
            vsizing_[v] = target_edge_length_;
        }
    }
    else
    {
        // compute curvature for all mesh vertices, using cotan or Cohen-Steiner
        // don't use two-ring neighborhood, since we otherwise compute
        // curvature over sharp features edges, leading to high curvatures.
        // prefer tensor analysis over cotan-Laplace, since the former is more
        // robust and gives better results on the boundary.
        Curvature curv(mesh_);
        curv.analyze_tensor(1);

        // use vsizing_ to store/smooth curvatures to avoid another vertex property

        // curvature values for feature vertices and boundary vertices
        // are not meaningful. mark them as negative values.
        for (auto v : mesh_.vertices())
        {
            if (mesh_.is_boundary(v) || (vfeature_ && vfeature_[v]))
                vsizing_[v] = -1.0;
            else
                vsizing_[v] = curv.max_abs_curvature(v);
        }

        // curvature values might be noisy. smooth them.
        // don't consider feature vertices' curvatures.
        // don't consider boundary vertices' curvatures.
        // do this for two iterations, to propagate curvatures
        // from non-feature regions to feature vertices.
        for (int iters = 0; iters < 2; ++iters)
        {
            for (auto v : mesh_.vertices())
            {
                Scalar w, ww = 0.0;
                Scalar c, cc = 0.0;

                for (auto h : mesh_.halfedges(v))
                {
                    c = vsizing_[mesh_.to_vertex(h)];
                    if (c > 0.0)
                    {
                        w = std::max(0.0, cotan_weight(mesh_, mesh_.edge(h)));
                        ww += w;
                        cc += w * c;
                    }
                }

                if (ww)
                    cc /= ww;
                vsizing_[v] = cc;
            }
        }

        // now convert per-vertex curvature into target edge length
        for (auto v : mesh_.vertices())
        {
            Scalar c = vsizing_[v];

            // get edge length from curvature
            const Scalar r = 1.0 / c;
            const Scalar e = approx_error_;
            Scalar h;
            if (e < r)
            {
                // see mathworld: "circle segment" and "equilateral triangle"
                //h = sqrt(2.0*r*e-e*e) * 3.0 / sqrt(3.0);
                h = sqrt(6.0 * e * r - 3.0 * e * e); // simplified...
            }
            else
            {
                // this does not really make sense
                h = e * 3.0 / sqrt(3.0);
            }

            // clamp to min. and max. edge length
            if (h < min_edge_length_)
                h = min_edge_length_;
            else if (h > max_edge_length_)
                h = max_edge_length_;

            // store target edge length
            vsizing_[v] = h;
        }
    }

    if (use_projection_)
    {
        // build reference mesh
        refmesh_ = std::make_shared<SurfaceMesh>();
        refmesh_->assign(mesh_);
        Normals::compute_vertex_normals(*refmesh_);
        refpoints_ = refmesh_->vertex_property<Point>("v:point");
        refnormals_ = refmesh_->vertex_property<Point>("v:normal");

        // copy sizing field from mesh_
        refsizing_ = refmesh_->add_vertex_property<Scalar>("v:sizing");
        for (auto v : refmesh_->vertices())
        {
            refsizing_[v] = vsizing_[v];
        }

        // build kd-tree
        kd_tree_ = std::make_unique<TriangleKdTree>(refmesh_, 0);
    }
}

void Remeshing::postprocessing()
{
    // remove properties
    mesh_.remove_vertex_property(vlocked_);
    mesh_.remove_edge_property(elocked_);
    mesh_.remove_vertex_property(vsizing_);

    if (!has_feature_vertices_)
    {
        mesh_.remove_vertex_property(vfeature_);
    }
    if (!has_feature_edges_)
    {
        mesh_.remove_edge_property(efeature_);
    }
}

void Remeshing::convex_hull_preprocessing()
{
    // Initialize and lock vertex features
    vfeature_ = mesh_.vertex_property<bool>("v:feature", false);
    efeature_ = mesh_.edge_property<bool>("e:feature", false);
    vlocked_ = mesh_.vertex_property<bool>("v:locked", false);
    elocked_ = mesh_.edge_property<bool>("e:locked", false);
    vsizing_ = mesh_.vertex_property<Scalar>("v:sizing", max_edge_length_);

    // Compute vertex normals for the mesh
    Normals::compute_vertex_normals(mesh_);

    // Manage edge features based on length
    for (const auto e : mesh_.edges())
    {
        if (mesh_.edge_length(e) > max_edge_length_) 
            continue;
        // DISCLAIMER: splitting a feature edge will propagate the feature property to its split vertex which should be subject to further remeshing.
        elocked_[e] = true;
        const auto v0 = mesh_.vertex(e, 0);
        const auto v1 = mesh_.vertex(e, 1);
        vlocked_[v0] = true;
        vlocked_[v1] = true;
    }

    if (use_projection_)
    {
        // build reference mesh
        refmesh_ = std::make_shared<SurfaceMesh>();
        refmesh_->assign(mesh_);
        Normals::compute_vertex_normals(*refmesh_);
        refpoints_ = refmesh_->vertex_property<Point>("v:point");
        refnormals_ = refmesh_->vertex_property<Point>("v:normal");

        // copy sizing field from mesh_
        refsizing_ = refmesh_->add_vertex_property<Scalar>("v:sizing");
        for (auto v : refmesh_->vertices())
        {
            refsizing_[v] = vsizing_[v];
        }

        // build kd-tree
        kd_tree_ = std::make_unique<TriangleKdTree>(refmesh_, 0);
    }
}

void Remeshing::convex_hull_postprocessing()
{
    for (const auto v : mesh_.vertices())
    {
        if (!vlocked_[v])
            continue;

        vfeature_[v] = true;
    }

    has_feature_vertices_ = true;
    postprocessing();
}


void Remeshing::project_to_reference(Vertex v)
{
    if (!use_projection_)
    {
        return;
    }

    // find closest triangle of reference mesh
    auto nn = kd_tree_->nearest(points_[v]);
    const Point p = nn.nearest;
    const Face f = nn.face;

    // get face data
    auto fvIt = refmesh_->vertices(f);
    const Point p0 = refpoints_[*fvIt];
    const Point n0 = refnormals_[*fvIt];
    const Scalar s0 = refsizing_[*fvIt];
    ++fvIt;
    const Point p1 = refpoints_[*fvIt];
    const Point n1 = refnormals_[*fvIt];
    const Scalar s1 = refsizing_[*fvIt];
    ++fvIt;
    const Point p2 = refpoints_[*fvIt];
    const Point n2 = refnormals_[*fvIt];
    const Scalar s2 = refsizing_[*fvIt];

    // get barycentric coordinates
    Point b = barycentric_coordinates(p, p0, p1, p2);

    // interpolate normal
    Point n;
    n = (n0 * b[0]);
    n += (n1 * b[1]);
    n += (n2 * b[2]);
    n.normalize();
    assert(!std::isnan(n[0]));

    // interpolate sizing field
    Scalar s;
    s = (s0 * b[0]);
    s += (s1 * b[1]);
    s += (s2 * b[2]);

    // set result
    points_[v] = p;
    vnormal_[v] = n;
    vsizing_[v] = s;
}

void Remeshing::split_long_edges(unsigned int nIterations)
{
    Vertex vnew, v0, v1;
    Edge enew;
    bool ok, is_feature, is_boundary;
    unsigned int i;

    for (ok = false, i = 0; !ok && i < nIterations; ++i)
    {
        ok = true;

        for (auto e : mesh_.edges())
        {
            v0 = mesh_.vertex(e, 0);
            v1 = mesh_.vertex(e, 1);

            if (!elocked_[e] && is_too_long(v0, v1))
            {
                const Point& p0 = points_[v0];
                const Point& p1 = points_[v1];

                is_feature = efeature_[e];
                is_boundary = mesh_.is_boundary(e);

                vnew = mesh_.add_vertex((p0 + p1) * 0.5f);
                mesh_.split(e, vnew);

                // need normal or sizing for adaptive refinement
                vnormal_[vnew] = Normals::compute_vertex_normal(mesh_, vnew);
                vsizing_[vnew] = 0.5f * (vsizing_[v0] + vsizing_[v1]);

                if (is_feature)
                {
                    enew = is_boundary ? Edge(mesh_.n_edges() - 2)
                                       : Edge(mesh_.n_edges() - 3);
                    efeature_[enew] = true;
                    vfeature_[vnew] = true;
                }
                else
                {
                    project_to_reference(vnew);
                }

                ok = false;
            }
        }
    }
}

void Remeshing::collapse_short_edges()
{
    Vertex v0, v1;
    Halfedge h0, h1, h01, h10;
    bool ok, b0, b1, l0, l1, f0, f1;
    int i;
    bool hcol01, hcol10;

    for (ok = false, i = 0; !ok && i < 10; ++i)
    {
        ok = true;

        for (auto e : mesh_.edges())
        {
            if (!mesh_.is_deleted(e) && !elocked_[e])
            {
                h10 = mesh_.halfedge(e, 0);
                h01 = mesh_.halfedge(e, 1);
                v0 = mesh_.to_vertex(h10);
                v1 = mesh_.to_vertex(h01);

                if (is_too_short(v0, v1))
                {
                    // get status
                    b0 = mesh_.is_boundary(v0);
                    b1 = mesh_.is_boundary(v1);
                    l0 = vlocked_[v0];
                    l1 = vlocked_[v1];
                    f0 = vfeature_[v0];
                    f1 = vfeature_[v1];
                    hcol01 = hcol10 = true;

                    // boundary rules
                    if (b0 && b1)
                    {
                        if (!mesh_.is_boundary(e))
                            continue;
                    }
                    else if (b0)
                        hcol01 = false;
                    else if (b1)
                        hcol10 = false;

                    // locked rules
                    if (l0 && l1)
                        continue;
                    else if (l0)
                        hcol01 = false;
                    else if (l1)
                        hcol10 = false;

                    // feature rules
                    if (f0 && f1)
                    {
                        // edge must be feature
                        if (!efeature_[e])
                            continue;

                        // the other two edges removed by collapse must not be features
                        h0 = mesh_.prev_halfedge(h01);
                        h1 = mesh_.next_halfedge(h10);
                        if (efeature_[mesh_.edge(h0)] ||
                            efeature_[mesh_.edge(h1)])
                            hcol01 = false;
                        // the other two edges removed by collapse must not be features
                        h0 = mesh_.prev_halfedge(h10);
                        h1 = mesh_.next_halfedge(h01);
                        if (efeature_[mesh_.edge(h0)] ||
                            efeature_[mesh_.edge(h1)])
                            hcol10 = false;
                    }
                    else if (f0)
                        hcol01 = false;
                    else if (f1)
                        hcol10 = false;

                    // topological rules
                    bool collapse_ok = mesh_.is_collapse_ok(h01);

                    if (hcol01)
                        hcol01 = collapse_ok;
                    if (hcol10)
                        hcol10 = collapse_ok;

                    // both collapses possible: collapse into vertex w/ higher valence
                    if (hcol01 && hcol10)
                    {
                        if (mesh_.valence(v0) < mesh_.valence(v1))
                            hcol10 = false;
                        else
                            hcol01 = false;
                    }

                    // try v1 -> v0
                    if (hcol10)
                    {
                        // don't create too long edges
                        for (auto vv : mesh_.vertices(v1))
                        {
                            if (is_too_long(v0, vv))
                            {
                                hcol10 = false;
                                break;
                            }
                        }

                        if (hcol10)
                        {
                            mesh_.collapse(h10);
                            ok = false;
                        }
                    }

                    // try v0 -> v1
                    else if (hcol01)
                    {
                        // don't create too long edges
                        for (auto vv : mesh_.vertices(v0))
                        {
                            if (is_too_long(v1, vv))
                            {
                                hcol01 = false;
                                break;
                            }
                        }

                        if (hcol01)
                        {
                            mesh_.collapse(h01);
                            ok = false;
                        }
                    }
                }
            }
        }
    }

    mesh_.garbage_collection();
}

void Remeshing::flip_edges()
{
    Vertex v0, v1, v2, v3;
    Halfedge h;
    int val0, val1, val2, val3;
    int val_opt0, val_opt1, val_opt2, val_opt3;
    int ve0, ve1, ve2, ve3, ve_before, ve_after;
    bool ok;
    int i;

    // precompute valences
    auto valence = mesh_.add_vertex_property<int>("valence");
    for (auto v : mesh_.vertices())
    {
        valence[v] = mesh_.valence(v);
    }

    for (ok = false, i = 0; !ok && i < 10; ++i)
    {
        ok = true;

        for (auto e : mesh_.edges())
        {
            if (!elocked_[e] && !efeature_[e])
            {
                h = mesh_.halfedge(e, 0);
                v0 = mesh_.to_vertex(h);
                v2 = mesh_.to_vertex(mesh_.next_halfedge(h));
                h = mesh_.halfedge(e, 1);
                v1 = mesh_.to_vertex(h);
                v3 = mesh_.to_vertex(mesh_.next_halfedge(h));

                if (!vlocked_[v0] && !vlocked_[v1] && !vlocked_[v2] &&
                    !vlocked_[v3])
                {
                    val0 = valence[v0];
                    val1 = valence[v1];
                    val2 = valence[v2];
                    val3 = valence[v3];

                    val_opt0 = (mesh_.is_boundary(v0) ? 4 : 6);
                    val_opt1 = (mesh_.is_boundary(v1) ? 4 : 6);
                    val_opt2 = (mesh_.is_boundary(v2) ? 4 : 6);
                    val_opt3 = (mesh_.is_boundary(v3) ? 4 : 6);

                    ve0 = (val0 - val_opt0);
                    ve1 = (val1 - val_opt1);
                    ve2 = (val2 - val_opt2);
                    ve3 = (val3 - val_opt3);

                    ve0 *= ve0;
                    ve1 *= ve1;
                    ve2 *= ve2;
                    ve3 *= ve3;

                    ve_before = ve0 + ve1 + ve2 + ve3;

                    --val0;
                    --val1;
                    ++val2;
                    ++val3;

                    ve0 = (val0 - val_opt0);
                    ve1 = (val1 - val_opt1);
                    ve2 = (val2 - val_opt2);
                    ve3 = (val3 - val_opt3);

                    ve0 *= ve0;
                    ve1 *= ve1;
                    ve2 *= ve2;
                    ve3 *= ve3;

                    ve_after = ve0 + ve1 + ve2 + ve3;

                    if (ve_before > ve_after && mesh_.is_flip_ok(e))
                    {
                        mesh_.flip(e);
                        --valence[v0];
                        --valence[v1];
                        ++valence[v2];
                        ++valence[v3];
                        ok = false;
                    }
                }
            }
        }
    }

    mesh_.remove_vertex_property(valence);
}

void Remeshing::tangential_smoothing(unsigned int iterations)
{
    Vertex v1, v2, v3, vv;
    Edge e;
    Scalar w, ww;
    Point u, n, t, b;

    // add property
    auto update = mesh_.add_vertex_property<Point>("v:update");

    // project at the beginning to get valid sizing values and normal vectors
    // for vertices introduced by splitting
    if (use_projection_)
    {
        for (auto v : mesh_.vertices())
        {
            if (!mesh_.is_boundary(v) && !vlocked_[v])
            {
                project_to_reference(v);
            }
        }
    }

    for (unsigned int iters = 0; iters < iterations; ++iters)
    {
        for (auto v : mesh_.vertices())
        {
            if (!mesh_.is_boundary(v) && !vlocked_[v])
            {
                if (vfeature_[v])
                {
                    u = Point(0.0);
                    t = Point(0.0);
                    ww = 0;
                    int c = 0;

                    for (auto h : mesh_.halfedges(v))
                    {
                        if (efeature_[mesh_.edge(h)])
                        {
                            vv = mesh_.to_vertex(h);

                            b = points_[v];
                            b += points_[vv];
                            b *= 0.5;

                            w = distance(points_[v], points_[vv]) /
                                (0.5 * (vsizing_[v] + vsizing_[vv]));
                            ww += w;
                            u += w * b;

                            if (c == 0)
                            {
                                t += normalize(points_[vv] - points_[v]);
                                ++c;
                            }
                            else
                            {
                                ++c;
                                t -= normalize(points_[vv] - points_[v]);
                            }
                        }
                    }

                    assert(c == 2);

                    u *= (1.0 / ww);
                    u -= points_[v];
                    t = normalize(t);
                    u = t * dot(u, t);

                    update[v] = u;
                }
                else
                {
                    Point p(0);
                    try
                    {
                        p = minimize_squared_areas(v);
                    }
                    catch (SolverException&)
                    {
                        p = weighted_centroid(v);
                    }
                    u = p - mesh_.position(v);

                    n = vnormal_[v];
                    u -= n * dot(u, n);

                    update[v] = u;
                }
            }
        }

        // update vertex positions
        for (auto v : mesh_.vertices())
        {
            if (!mesh_.is_boundary(v) && !vlocked_[v])
            {
                points_[v] += update[v];
            }
        }

        // update normal vectors (if not done so through projection)
        Normals::compute_vertex_normals(mesh_);
    }

    // project at the end
    if (use_projection_)
    {
        for (auto v : mesh_.vertices())
        {
            if (!mesh_.is_boundary(v) && !vlocked_[v])
            {
                project_to_reference(v);
            }
        }
    }

    // remove property
    mesh_.remove_vertex_property(update);
}

void Remeshing::remove_caps()
{
    Halfedge h;
    Vertex v, vb, vd;
    Face fb, fd;
    Scalar a0, a1, amin, aa(::cos(170.0 * M_PI / 180.0));
    Point a, b, c, d;

    for (auto e : mesh_.edges())
    {
        if (!elocked_[e] && mesh_.is_flip_ok(e))
        {
            h = mesh_.halfedge(e, 0);
            a = points_[mesh_.to_vertex(h)];

            h = mesh_.next_halfedge(h);
            b = points_[vb = mesh_.to_vertex(h)];

            h = mesh_.halfedge(e, 1);
            c = points_[mesh_.to_vertex(h)];

            h = mesh_.next_halfedge(h);
            d = points_[vd = mesh_.to_vertex(h)];

            a0 = dot(normalize(a - b), normalize(c - b));
            a1 = dot(normalize(a - d), normalize(c - d));

            if (a0 < a1)
            {
                amin = a0;
                v = vb;
            }
            else
            {
                amin = a1;
                v = vd;
            }

            // is it a cap?
            if (amin < aa)
            {
                // feature edge and feature vertex -> seems to be intended
                if (efeature_[e] && vfeature_[v])
                    continue;

                // project v onto feature edge
                if (efeature_[e])
                    points_[v] = (a + c) * 0.5f;

                // flip
                mesh_.flip(e);
            }
        }
    }
}

Point Remeshing::minimize_squared_areas(Vertex v)
{
    dmat3 A(0);
    dvec3 b(0), x;

    for (auto h : mesh_.halfedges(v))
    {
        assert(!mesh_.is_boundary(h));

        // get edge opposite to vertex v
        auto v0 = mesh_.to_vertex(h);
        auto v1 = mesh_.to_vertex(mesh_.next_halfedge(h));
        auto p = (dvec3)points_[v0];
        auto q = (dvec3)points_[v1];
        auto d = q - p;
        auto w = 1.0 / norm(d);

        // build squared cross-product-with-d matrix
        dmat3 D;
        D(0, 0) = d[1] * d[1] + d[2] * d[2];
        D(1, 1) = d[0] * d[0] + d[2] * d[2];
        D(2, 2) = d[0] * d[0] + d[1] * d[1];
        D(1, 0) = D(0, 1) = -d[0] * d[1];
        D(2, 0) = D(0, 2) = -d[0] * d[2];
        D(1, 2) = D(2, 1) = -d[1] * d[2];
        A += w * D;

        // build right-hand side
        b += w * D * p;
    }

    // compute minimizer
    try
    {
        x = inverse(A) * b;
    }
    catch (...)
    {
        auto what = "Remeshing: Matrix not invertible.";
        throw SolverException(what);
    }

    return Point(x);
}

Point Remeshing::weighted_centroid(Vertex v)
{
    auto p = Point(0);
    double ww = 0;

    for (auto h : mesh_.halfedges(v))
    {
        auto v1 = v;
        auto v2 = mesh_.to_vertex(h);
        auto v3 = mesh_.to_vertex(mesh_.next_halfedge(h));

        auto b = points_[v1];
        b += points_[v2];
        b += points_[v3];
        b *= (1.0 / 3.0);

        double area =
            norm(cross(points_[v2] - points_[v1], points_[v3] - points_[v1]));

        // take care of degenerate faces to avoid all zero weights and division
        // by zero later on
        if (area == 0)
            area = 1.0;

        double w =
            area / pow((vsizing_[v1] + vsizing_[v2] + vsizing_[v3]) / 3.0, 2.0);

        p += w * b;
        ww += w;
    }

    p /= ww;

    return p;
}

CurveRemeshing::CurveRemeshing(ManifoldCurve2D& curve, const std::shared_ptr<EvolvingArcLengthCalculator>& arcLengthCalc)
    : curve_(curve), refcurve_(nullptr), kd_tree_(nullptr), m_LengthCalculator(arcLengthCalc)
{
    if (curve_.is_empty())
        throw InvalidInputException("CurveRemeshing::CurveRemeshing: Input has to be a non-empty manifold curve!");

    points_ = curve_.vertex_property<Point2>("v:point");

    Normals2::compute_vertex_normals(curve_);
    vnormal_ = curve_.vertex_property<Point2>("v:normal");

    has_feature_vertices_ = curve_.has_vertex_property("v:feature");
}

void CurveRemeshing::uniform_remeshing(Scalar edge_length, unsigned int iterations, bool use_projection)
{
    uniform_ = true;
    use_projection_ = use_projection;
    target_edge_length_ = edge_length;

    preprocessing();

    for (unsigned int i = 0; i < iterations; ++i)
    {
        split_long_edges();

        Normals2::compute_vertex_normals(curve_);

        collapse_short_edges();

        tangential_smoothing(5);
    }

    postprocessing();
}

void CurveRemeshing::adaptive_remeshing(const AdaptiveRemeshingSettings& settings)
{
    uniform_ = false;
    min_edge_length_ = settings.MinEdgeLength;
    max_edge_length_ = settings.MaxEdgeLength;
    approx_error_ = settings.ApproxError;
    use_projection_ = settings.UseProjection;

    preprocessing();

    for (unsigned int i = 0; i < settings.NRemeshingIterations; ++i)
    {
        split_long_edges();

        Normals2::compute_vertex_normals(curve_);

        collapse_short_edges();

        tangential_smoothing(settings.NTangentialSmoothingIters);
    }

    postprocessing();
}

void CurveRemeshing::preprocessing()
{
    // properties
    if (!vfeature_)
        vfeature_ = curve_.vertex_property<bool>("v:feature", false);
    if (!vlocked_)
        vlocked_ = curve_.add_vertex_property<bool>("v:locked", false);
    if (!vsizing_)
        vsizing_ = curve_.add_vertex_property<Scalar>("v:sizing");

    // lock unselected vertices if some vertices are selected
    auto vselected = curve_.get_vertex_property<bool>("v:selected");
    if (vselected)
    {
        bool has_selection = false;
        for (auto v : curve_.vertices())
        {
            if (vselected[v])
            {
                has_selection = true;
                break;
            }
        }

        if (has_selection)
        {
            for (auto v : curve_.vertices())
            {
                vlocked_[v] = !vselected[v];
            }

            // lock an edge if one of its vertices is locked
            for (auto e : curve_.edges())
            {
                vlocked_[curve_.from_vertex(e)] = vlocked_[curve_.from_vertex(e)] || vlocked_[curve_.to_vertex(e)];
                vlocked_[curve_.to_vertex(e)] = vlocked_[curve_.from_vertex(e)] || vlocked_[curve_.to_vertex(e)];
            }
        }
    }

    // lock feature corners
    for (auto v : curve_.vertices())
    {
        if (vfeature_[v])
        {
            int c = 0;
            auto [e1, e2] = curve_.edges(v);
            if (curve_.is_valid(e1) && curve_.is_valid(e2))
            {
                if (vfeature_[curve_.from_vertex(e1)] && vfeature_[curve_.to_vertex(e1)])
                    ++c;
                if (vfeature_[curve_.from_vertex(e2)] && vfeature_[curve_.to_vertex(e2)])
                    ++c;
            }

            if (c != 2)
                vlocked_[v] = true;
        }
    }

    // compute sizing field
    if (uniform_)
    {
        for (auto v : curve_.vertices())
        {
            vsizing_[v] = target_edge_length_;
        }
    }
    else
    {
        // compute curvature for all curve vertices
        Curvature1D curv(curve_);
        curv.analyze();

        // use vsizing_ to store/smooth curvatures to avoid another vertex property

        // curvature values for feature vertices are not meaningful. mark them as negative values.
        for (auto v : curve_.vertices())
        {
            if (vfeature_ && vfeature_[v])
                vsizing_[v] = -1.0;
            else
                vsizing_[v] = curv.curvature(v);
        }

        // curvature values might be noisy. smooth them.
        // don't consider feature vertices' curvatures.
        for (int iters = 0; iters < 2; ++iters)
        {
            for (auto v : curve_.vertices())
            {
                Scalar w, ww = 0.0;
                Scalar c, cc = 0.0;

                auto [e1, e2] = curve_.edges(v);
                if (curve_.is_valid(e1))
                {
                    auto tv = curve_.to_vertex(e1);
                    c = vsizing_[tv];
                    if (c > 0.0)
                    {
                        w = 1.0; // simple uniform weight for 2D curves
                        ww += w;
                        cc += w * c;
                    }
                }
                if (curve_.is_valid(e2))
                {
                    auto tv = curve_.from_vertex(e2);
                    c = vsizing_[tv];
                    if (c > 0.0)
                    {
                        w = 1.0; // simple uniform weight for 2D curves
                        ww += w;
                        cc += w * c;
                    }
                }

                if (ww)
                    cc /= ww;
                vsizing_[v] = cc;
            }
        }

        // now convert per-vertex curvature into target edge length
        for (auto v : curve_.vertices())
        {
            Scalar c = vsizing_[v];

            // get edge length from curvature
            const Scalar r = 1.0 / c;
            const Scalar e = approx_error_;
            Scalar h;
            if (e < r)
            {
                // see mathworld: "circle segment" and "equilateral triangle"
                // h = sqrt(2.0*r*e-e*e) * 3.0 / sqrt(3.0);
                h = sqrt(6.0 * e * r - 3.0 * e * e); // simplified...
            }
            else
            {
                // this does not really make sense
                h = e * 3.0 / sqrt(3.0);
            }

            // clamp to min. and max. edge length
            if (h < min_edge_length_)
                h = min_edge_length_;
            else if (h > max_edge_length_)
                h = max_edge_length_;

            // store target edge length
            vsizing_[v] = h;
        }
    }

    if (use_projection_)
    {
        // build reference curve
        refcurve_ = std::make_shared<ManifoldCurve2D>();
        refcurve_->assign(curve_);
        Normals2::compute_vertex_normals(*refcurve_);
        refpoints_ = refcurve_->vertex_property<Point2>("v:point");
        refnormals_ = refcurve_->vertex_property<Normal2>("v:normal");

        // copy sizing field from curve_
        refsizing_ = refcurve_->add_vertex_property<Scalar>("v:sizing");
        for (auto v : refcurve_->vertices())
        {
            refsizing_[v] = vsizing_[v];
        }

        // build kd-tree
        kd_tree_ = std::make_unique<EdgeKdTree>(refcurve_, 0);
    }

    if (m_LengthCalculator)
    {
        m_LengthCalculator->RecordPrevRefEdgePositions();
    }
}


void CurveRemeshing::postprocessing()
{
    // remove properties
    curve_.remove_vertex_property(vlocked_);
    curve_.remove_vertex_property(vsizing_);

    if (!has_feature_vertices_)
    {
        curve_.remove_vertex_property(vfeature_);
    }

    if (m_LengthCalculator)
    {
        m_LengthCalculator->UpdateRefEdge();
    }
}

void CurveRemeshing::split_long_edges(unsigned int nIterations)
{
    Vertex vnew, v0, v1;
    bool ok;
    unsigned int i;

    for (ok = false, i = 0; !ok && i < nIterations; ++i)
    {
        ok = true;

        for (auto e : curve_.edges())
        {
            std::tie(v0, v1) = curve_.vertices(e);

            if (!(vlocked_[v0] || vlocked_[v1]) && is_too_long(v0, v1))
            {
                const Point2& p0 = curve_.position(v0);
                const Point2& p1 = curve_.position(v1);

                vnew = curve_.split_edge_with_vertex(e, (p0 + p1) * 0.5f);

                // need normal or sizing for adaptive refinement
                vnormal_[vnew] = Normals2::compute_vertex_normal(curve_, vnew);
                vsizing_[vnew] = 0.5f * (vsizing_[v0] + vsizing_[v1]);

                // Check if either vertex is a feature vertex
                if (vfeature_[v0] || vfeature_[v1])
                {
                    vfeature_[vnew] = true;
                }
                else
                {
                    project_to_reference(vnew);
                }

                ok = false;
            }
        }
    }
}

void CurveRemeshing::collapse_short_edges()
{
    Vertex v0, v1;
    bool ok, b0, b1, l0, l1, f0, f1;
    int i;

    for (ok = false, i = 0; !ok && i < 10; ++i)
    {
        ok = true;

        for (auto e : curve_.edges())
        {
            if (!curve_.is_deleted(e))
            {
                v0 = curve_.from_vertex(e);
                v1 = curve_.to_vertex(e);

                if (is_too_short(v0, v1))
                {
                    // get status
                    b0 = curve_.is_boundary(v0);
                    b1 = curve_.is_boundary(v1);
                    l0 = vlocked_[v0];
                    l1 = vlocked_[v1];
                    f0 = vfeature_[v0];
                    f1 = vfeature_[v1];

                    // Locked rules
                    if (l0 && l1)
                        continue;

                    // Feature rules
                    if (f0 && f1)
                        continue;

                    // Collapse v1 to v0 if possible
                    if (!l0 && (!f0 || f1) && !(b0 && !b1))
                    {
                        curve_.collapse_edge(v0, v1, true);
                        ok = false;
                    }
                    // Collapse v0 to v1 if possible
                    else if (!l1 && (!f1 || f0) && !(b1 && !b0))
                    {
                        curve_.collapse_edge(v0, v1, false);
                        ok = false;
                    }
                }
            }
        }
    }

    curve_.garbage_collection();
}


void CurveRemeshing::tangential_smoothing(unsigned int iterations)
{
    Vertex v_prev, v_next;
    Point2 u, t, b;

    // add property
    auto update = curve_.add_vertex_property<Point2>("v:update");

    // project at the beginning to get valid sizing values and normal vectors
    // for vertices introduced by splitting
    if (use_projection_)
    {
        for (const auto v : curve_.vertices())
        {
            if (!curve_.is_boundary(v) && !vlocked_[v])
            {
                project_to_reference(v);
            }
        }
    }

    for (unsigned int iters = 0; iters < iterations; ++iters)
    {
        for (const auto v : curve_.vertices())
        {
            std::tie(v_prev, v_next) = curve_.vertices(v);

            if (!curve_.is_boundary(v) && !vlocked_[v])
            {
                if (vfeature_[v])
                {
                    u = Point2(0.0);
                    t = Point2(0.0);
                    int c = 0;

                    if (curve_.is_valid(v_prev) && vfeature_[v_prev])
                    {
                        b = (curve_.position(v) + curve_.position(v_prev)) * 0.5;

                        u += b;

                        if (c == 0)
                        {
                            t = normalize(curve_.position(v_prev) - curve_.position(v));
                            ++c;
                        }
                        else
                        {
                            t -= normalize(curve_.position(v_prev) - curve_.position(v));
                        }
                    }

                    if (curve_.is_valid(v_next) && vfeature_[v_next])
                    {
                        b = (curve_.position(v) + curve_.position(v_next)) * 0.5;

                        u += b;

                        if (c == 0)
                        {
                            t = normalize(curve_.position(v_next) - curve_.position(v));
                            ++c;
                        }
                        else
                        {
                            t -= normalize(curve_.position(v_next) - curve_.position(v));
                        }
                    }

                    assert(c == 2);

                    u *= 0.5; // average of the two midpoints
                    u -= curve_.position(v);
                    u = t * dot(u, t);

                    update[v] = u;
                }
                else
                {
                    u = (curve_.position(v_prev) + curve_.position(v_next)) * 0.5 - curve_.position(v);
                    t = normalize(curve_.position(v_next) - curve_.position(v_prev));
                    u = t * dot(u, t);

                    update[v] = u;
                }
            }
        }

        // update vertex positions
        for (auto v : curve_.vertices())
        {
            if (!curve_.is_boundary(v) && !vlocked_[v])
            {
                points_[v] += update[v];
            }
        }

        // update normal vectors (if not done so through projection)
        Normals2::compute_vertex_normals(curve_);
    }

    // project at the end
    if (use_projection_)
    {
        for (auto v : curve_.vertices())
        {
            if (!curve_.is_boundary(v) && !vlocked_[v])
            {
                project_to_reference(v);
            }
        }
    }

    // remove property
    curve_.remove_vertex_property(update);
}


void CurveRemeshing::project_to_reference(Vertex v)
{
    if (!use_projection_)
    {
        return;
    }

    // find closest edge of reference curve
    const auto nn = kd_tree_->nearest(points_[v]);
    const Point2 p = nn.nearest;
    const Edge e = nn.edge;

    // Get edge data
    const Vertex v0 = refcurve_->from_vertex(e);
    const Vertex v1 = refcurve_->to_vertex(e);

    const Point2 p0 = refpoints_[v0];
    const Point2 n0 = refnormals_[v0];
    const Scalar s0 = refsizing_[v0];

    const Point2 p1 = refpoints_[v1];
    const Point2 n1 = refnormals_[v1];
    const Scalar s1 = refsizing_[v1];

    // Compute the barycentric coordinate t along the edge
    const Scalar edge_length = distance(p0, p1);
    const Scalar t = distance(p0, p) / edge_length;

    // Interpolate normal
    Point2 n = (1 - t) * n0 + t * n1;
    n.normalize();
    assert(!std::isnan(n[0]));

    // Interpolate sizing field
    const Scalar s = (1 - t) * s0 + t * s1;

    // Set result
    curve_.position(v) = p;
    vnormal_[v] = n;
    vsizing_[v] = s;
}

} // namespace pmp
