// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#pragma once

#include "pmp/SurfaceMesh.h"
#include "pmp/ManifoldCurve2D.h"

namespace pmp {

//! \brief Compute per-vertex curvature (min,max,mean,Gaussian).
//! \details Curvature values for boundary vertices are interpolated from their
//! interior neighbors. Curvature values can be smoothed. See
//! \cite meyer_2003_discrete and \cite cohen-steiner_2003_restricted for
//! details.
//! \ingroup algorithms
class Curvature
{
public:
    //! construct with mesh to be analyzed
    Curvature(SurfaceMesh& mesh);

    //! destructor
    ~Curvature();

    //! compute curvature information for each vertex, optionally followed
    //! by some smoothing iterations of the curvature values
    void analyze(unsigned int post_smoothing_steps = 0);

    //! compute curvature information for each vertex, optionally followed
    //! by some smoothing iterations of the curvature values
    void analyze_tensor(unsigned int post_smoothing_steps = 0,
                        bool two_ring_neighborhood = false);

    //! return mean curvature
    Scalar mean_curvature(Vertex v) const
    {
        return Scalar(0.5) * (min_curvature_[v] + max_curvature_[v]);
    }

    //! return Gaussian curvature
    Scalar gauss_curvature(Vertex v) const
    {
        return min_curvature_[v] * max_curvature_[v];
    }

    //! return minimum (signed) curvature
    Scalar min_curvature(Vertex v) const { return min_curvature_[v]; }

    //! return maximum (signed) curvature
    Scalar max_curvature(Vertex v) const { return max_curvature_[v]; }

    //! return maximum absolute curvature
    Scalar max_abs_curvature(Vertex v) const
    {
        return std::max(fabs(min_curvature_[v]), fabs(max_curvature_[v]));
    }

    //! convert (precomputed) mean curvature to 1D texture coordinates
    void mean_curvature_to_texture_coordinates() const;

    //! convert (precomputed) Gauss curvature to 1D texture coordinates
    void gauss_curvature_to_texture_coordinates() const;

    //! convert (precomputed) max. abs. curvature to 1D texture coordinates
    void max_curvature_to_texture_coordinates() const;

private:
    // smooth curvature values
    void smooth_curvatures(unsigned int iterations);

    // convert curvature values ("v:curv") to 1D texture coordinates
    void curvature_to_texture_coordinates() const;

    SurfaceMesh& mesh_;
    VertexProperty<Scalar> min_curvature_;
    VertexProperty<Scalar> max_curvature_;
};

//! \brief Compute per-vertex curvature for a 1D curve.
//! \details Curvature values are computed for each vertex.
//! \ingroup algorithms
class Curvature1D
{
public:
    //! construct with curve to be analyzed
    Curvature1D(ManifoldCurve2D& curve);

    //! destructor
    ~Curvature1D();

    //! compute curvature information for each vertex
    void analyze();

    //! return curvature
    Scalar curvature(Vertex v) const { return curvature_[v]; }

private:
    ManifoldCurve2D& curve_;
    VertexProperty<Scalar> curvature_;
};

} // namespace pmp
