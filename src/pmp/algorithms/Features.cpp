// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "pmp/algorithms/Features.h"
#include "pmp/algorithms/Normals.h"
#include "DifferentialGeometry.h"

#include "Curvature.h"

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
    //const Scalar feature_cosine_min = cos(minAngle / 180.0 * M_PI);
    //const Scalar feature_cosine_max = cos(maxAngle / 180.0 * M_PI);
    size_t n_edges = 0;
    for (auto e : mesh_.edges())
    {
        if (!mesh_.is_boundary(e))
        {
            const auto f0 = mesh_.face(mesh_.halfedge(e, 0));
            const auto f1 = mesh_.face(mesh_.halfedge(e, 1));

            const Normal n0 = Normals::compute_face_normal(mesh_, f0);
            const Normal n1 = Normals::compute_face_normal(mesh_, f1);

            const auto angleBetweenNormals = angle(n0, n1);
            if (angleBetweenNormals < minAngle && angleBetweenNormals > maxAngle)
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
    assert(principalCurvatureFactor >= 1.0);
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
    	const bool isV0Saddle = (vMinCurvature < 0.0 && vMaxCurvature > 0.0) || (vMinCurvature > 0.0 && vMaxCurvature < 0.0);
        const bool isV0Concave = vMinCurvature < 0.0 && vMaxCurvature < 0.0;
    	auto vAbsMinCurvature = std::fabs(vMinCurvature);
        auto vAbsMaxCurvature = std::fabs(vMaxCurvature);
        if (!isV0Saddle && vAbsMaxCurvature > principalCurvatureFactor * vAbsMinCurvature)
            vfeature_[v0] = true;
        //if (!isV0Saddle && !isV0Concave && vAbsMaxCurvature > principalCurvatureFactor * vAbsMinCurvature)
        //   vfeature_[v0] = true;

        const auto v1 = mesh_.vertex(e, 1);
        vMinCurvature = curvAlg.min_curvature(v1);
        vMaxCurvature = curvAlg.max_curvature(v1);
        const bool isV1Saddle = (vMinCurvature < 0.0 && vMaxCurvature > 0.0) || (vMinCurvature > 0.0 && vMaxCurvature < 0.0);
        const bool isV1Concave = vMinCurvature < 0.0 && vMaxCurvature < 0.0;
        vAbsMinCurvature = std::fabs(vMinCurvature);
        vAbsMaxCurvature = std::fabs(vMaxCurvature);
        if (!isV1Saddle && vAbsMaxCurvature > principalCurvatureFactor * vAbsMinCurvature)
            vfeature_[v1] = true;
        //if (!isV1Saddle && !isV1Concave && vAbsMaxCurvature > principalCurvatureFactor * vAbsMinCurvature)
        //    vfeature_[v1] = true;

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

/// \brief computes mean length of an outgoing edge from a given vertex.
[[nodiscard]] Scalar ComputeMeanArcLengthAtVertex(const SurfaceMesh& mesh, const Vertex& v)
{
    Scalar edgeLength = 0.0;
    size_t valence = 0;
	for (const auto w : mesh.vertices(v))
	{
        edgeLength += norm(mesh.position(v) - mesh.position(w));
        valence++;
	}
    return 2.0 * edgeLength / static_cast<Scalar>(valence);
}

bool IsConvexDominantSaddle(const Scalar& vMinCurvature, const Scalar& vMaxCurvature, const Scalar& curvatureFactor)
{
    const auto vAbsMinCurvature = std::fabs(vMinCurvature);
    const auto vAbsMaxCurvature = std::fabs(vMaxCurvature);

    const bool isSaddle = (vMinCurvature < 0.0 && vMaxCurvature > 0.0) || (vMinCurvature > 0.0 && vMaxCurvature < 0.0);

    return (vAbsMaxCurvature < curvatureFactor * vAbsMinCurvature) && isSaddle;
}

size_t Features::detect_vertices_with_high_curvature(const Scalar& curvatureAngle, const Scalar& principalCurvatureFactor, const bool& excludeEdgesWithoutTwoFeatureVerts)
{
    assert(curvatureAngle >= 0.0);
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
        auto valence = mesh_.valence(v0);
        Scalar vMeanArcLength, vMeanCurvature, vMeanCurvatureAngle;
        if (!IsConvexDominantSaddle(vMinCurvature, vMaxCurvature, principalCurvatureFactor) && valence < 7)
        {
            // only vertices with low principal curvature imbalance are tested for mean curvature angle
			vMeanCurvature = vMinCurvature + vMaxCurvature;
	        vMeanArcLength = ComputeMeanArcLengthAtVertex(mesh_, v0);
	        vMeanCurvatureAngle = vMeanCurvature * vMeanArcLength;
	        if (vMeanCurvatureAngle < curvatureAngle)
	            vfeature_[v0] = true;
        }

        const auto v1 = mesh_.vertex(e, 1);
        vMinCurvature = curvAlg.min_curvature(v1);
        vMaxCurvature = curvAlg.max_curvature(v1);
        valence = mesh_.valence(v1);
        if (!IsConvexDominantSaddle(vMinCurvature, vMaxCurvature, principalCurvatureFactor) && valence < 7)
        {
            // only vertices with low principal curvature imbalance are tested for mean curvature angle
	        vMeanCurvature = vMinCurvature + vMaxCurvature;
	        vMeanArcLength = ComputeMeanArcLengthAtVertex(mesh_, v1);
	        vMeanCurvatureAngle = vMeanCurvature * vMeanArcLength;
	        if (vMeanCurvatureAngle < curvatureAngle)
	            vfeature_[v1] = true;	        
        }

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
