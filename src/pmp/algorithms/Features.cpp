// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "pmp/algorithms/Features.h"

#include "Curvature.h"
#include "pmp/algorithms/Normals.h"

namespace pmp {

Features::Features(SurfaceMesh& mesh) : mesh_(mesh)
{
    vfeature_ = mesh_.vertex_property("v:feature", false);
    efeature_ = mesh_.edge_property("e:feature", false);
}

void Features::clear()
{
    for (auto v : mesh_.vertices())
        vfeature_[v] = false;

    for (auto e : mesh_.edges())
        efeature_[e] = false;
}

size_t Features::detect_boundary()
{
    for (auto v : mesh_.vertices())
        if (mesh_.is_boundary(v))
            vfeature_[v] = true;

    size_t n_edges = 0;
    for (auto e : mesh_.edges())
        if (mesh_.is_boundary(e))
        {
            efeature_[e] = true;
            n_edges++;
        }
    return n_edges;
}

size_t Features::detect_angle(Scalar angle)
{
    const Scalar feature_cosine = cos(angle / 180.0 * M_PI);
    size_t n_edges = 0;
    for (auto e : mesh_.edges())
    {
        if (!mesh_.is_boundary(e))
        {
            const auto f0 = mesh_.face(mesh_.halfedge(e, 0));
            const auto f1 = mesh_.face(mesh_.halfedge(e, 1));

            const Normal n0 = Normals::compute_face_normal(mesh_, f0);
            const Normal n1 = Normals::compute_face_normal(mesh_, f1);

            if (dot(n0, n1) < feature_cosine)
            {
                efeature_[e] = true;
                vfeature_[mesh_.vertex(e, 0)] = true;
                vfeature_[mesh_.vertex(e, 1)] = true;
                n_edges++;
            }
        }
    }
    return n_edges;
}

size_t Features::detect_angle_within_bounds(Scalar minAngle, Scalar maxAngle)
{
    const Scalar feature_cosine_min = cos(minAngle / 180.0 * M_PI);
    const Scalar feature_cosine_max = cos(maxAngle / 180.0 * M_PI);
    size_t n_edges = 0;
    for (auto e : mesh_.edges())
    {
        if (!mesh_.is_boundary(e))
        {
            const auto f0 = mesh_.face(mesh_.halfedge(e, 0));
            const auto f1 = mesh_.face(mesh_.halfedge(e, 1));

            const Normal n0 = Normals::compute_face_normal(mesh_, f0);
            const Normal n1 = Normals::compute_face_normal(mesh_, f1);

            const auto dotNormals = dot(n0, n1);
            if (dotNormals < feature_cosine_min && dotNormals > feature_cosine_max)
            {
                efeature_[e] = true;
                vfeature_[mesh_.vertex(e, 0)] = true;
                vfeature_[mesh_.vertex(e, 1)] = true;
                n_edges++;
            }
        }
    }
    return n_edges;
}

size_t Features::detect_vertices_with_curvatures_imbalance(const Scalar& principalCurvatureFactor, const bool& excludeEdgesWithoutTwoFeatureVerts)
{
    assert(principalCurvatureFactor >= 1.0f);
    Curvature curvAlg{ mesh_ };
    curvAlg.analyze_tensor(1);

    size_t n_edges = 0;
    for (const auto e : mesh_.edges())
    {
        if (mesh_.is_boundary(e))
            continue;

        const auto v0 = mesh_.vertex(e, 0);
        auto vMinCurvature = curvAlg.min_curvature(v0);
        auto vMaxCurvature = curvAlg.max_curvature(v0);
        //if (vMinCurvature > 0.0 && vMaxCurvature < 0.0)
        //   continue;

	    auto vAbsMinCurvature = std::fabs(vMinCurvature);
        auto vAbsMaxCurvature = std::fabs(vMaxCurvature);
        if (vAbsMaxCurvature > principalCurvatureFactor * vAbsMinCurvature)
            vfeature_[v0] = true;

        const auto v1 = mesh_.vertex(e, 1);
        vMinCurvature = curvAlg.min_curvature(v1);
        vMaxCurvature = curvAlg.max_curvature(v1);
        //if (vMinCurvature > 0.0 && vMaxCurvature < 0.0)
        //    continue;

        vAbsMinCurvature = std::fabs(vMinCurvature);
        vAbsMaxCurvature = std::fabs(vMaxCurvature);

        if (vAbsMaxCurvature > principalCurvatureFactor * vAbsMinCurvature)
            vfeature_[v1] = true;

        if (excludeEdgesWithoutTwoFeatureVerts && (vfeature_[v0] && vfeature_[v1]))
        {
            efeature_[e] = true;
            n_edges++;
            continue;
        }

        if (vfeature_[v0] || vfeature_[v1])
        {
            efeature_[e] = true;
            n_edges++;
        }
    }

    return n_edges;
}

} // namespace pmp
