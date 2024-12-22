// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#pragma once

#include "pmp/SurfaceMesh.h"
#include "pmp/ManifoldCurve2D.h"
#include "pmp/algorithms/ArcLengthCalculator.h"
#include "pmp/algorithms/TriangleKdTree.h"
#include "pmp/algorithms/EdgeKdTree.h"

#include <memory>

namespace pmp {
/**
 * \brief A wrapper for adaptive remeshing settings.
 * \struct AdaptiveRemeshingSettings
 */
struct AdaptiveRemeshingSettings
{
    Scalar MinEdgeLength{ 1.0 };
    Scalar MaxEdgeLength{ 2.0 };
    Scalar ApproxError{ 1.5 };
    unsigned int NRemeshingIterations{ 10 };
    unsigned int NTangentialSmoothingIters{ 6 };
    bool UseProjection{ true };
};

//! \brief A class for uniform and adaptive surface remeshing.
//! \details The algorithm implemented here performs incremental remeshing based
//! on edge collapse, split, flip, and tangential relaxation.
//! See \cite botsch_2004_remeshing and \cite dunyach_2013_adaptive for a more
//! detailed description.
//! \ingroup algorithms
class Remeshing
{
public:
    //! \brief Construct with mesh to be remeshed.
    //! \pre Input mesh needs to be a pure triangle mesh.
    //! \throw InvalidInputException if the input precondition is violated.
    Remeshing(SurfaceMesh& mesh);

    //! \brief Perform uniform remeshing.
    //! \param edge_length the target edge length.
    //! \param iterations the number of iterations
    //! \param use_projection use back-projection to the input surface
    void uniform_remeshing(Scalar edge_length, unsigned int iterations = 10,
                           bool use_projection = true);

    //! \brief Perform adaptive remeshing.
    //! \param settings         input settings.
    //! \param is_convex_hull   if true post-processing will be postponed for later.
    void adaptive_remeshing(const AdaptiveRemeshingSettings& settings, const bool& is_convex_hull = false);

    //! \brief Perform adaptive remeshing on a convex hull starting surface. Prioritizes splitting of elongated triangles.
    //! \param settings      input settings.
    void convex_hull_adaptive_remeshing(const AdaptiveRemeshingSettings& settings);

private:
    void preprocessing();
    void postprocessing();
    void convex_hull_preprocessing();
    void convex_hull_postprocessing();

    void split_long_edges(unsigned int nIterations = 10);
    void collapse_short_edges();
    void flip_edges();
    void tangential_smoothing(unsigned int iterations);
    void remove_caps();

    Point minimize_squared_areas(Vertex v);
    Point weighted_centroid(Vertex v);

    void project_to_reference(Vertex v);

    bool is_too_long(Vertex v0, Vertex v1) const
    {
        return distance(points_[v0], points_[v1]) >
               4.0 / 3.0 * std::min(vsizing_[v0], vsizing_[v1]);
    }
    bool is_too_short(Vertex v0, Vertex v1) const
    {
        return distance(points_[v0], points_[v1]) <
               4.0 / 5.0 * std::min(vsizing_[v0], vsizing_[v1]);
    }

    SurfaceMesh& mesh_;
    std::shared_ptr<SurfaceMesh> refmesh_;

    bool use_projection_;
    std::unique_ptr<TriangleKdTree> kd_tree_;

    bool uniform_;
    Scalar target_edge_length_;
    Scalar min_edge_length_;
    Scalar max_edge_length_;
    Scalar approx_error_;

    bool has_feature_vertices_{false};
    bool has_feature_edges_{false};
    VertexProperty<Point> points_;
    VertexProperty<Point> vnormal_;
    VertexProperty<bool> vfeature_;
    EdgeProperty<bool> efeature_;
    VertexProperty<bool> vlocked_;
    EdgeProperty<bool> elocked_;
    VertexProperty<Scalar> vsizing_;

    VertexProperty<Point> refpoints_;
    VertexProperty<Point> refnormals_;
    VertexProperty<Scalar> refsizing_;
};

//! \brief A class for uniform and adaptive 2D curve remeshing.
//! \details The algorithm implemented here performs incremental remeshing based
//! on edge collapse, split, and tangential relaxation, with optional reprojection.
//! \ingroup algorithms
class CurveRemeshing
{
public:
    //! \brief Construct with curve to be remeshed.
    //! \pre Input curve needs to be a pure manifold curve.
    //! \param arcLengthCalc       an optional arc length calculator which needs to be updated before garbage collection.
    //! \throw InvalidInputException if the input precondition is violated.
    CurveRemeshing(ManifoldCurve2D& curve, const std::shared_ptr<EvolvingArcLengthCalculator>& arcLengthCalc = nullptr);

    //! \brief Perform uniform remeshing.
    //! \param edge_length the target edge length.
    //! \param iterations the number of iterations.
    //! \param use_projection use back-projection to the input curve.
    void uniform_remeshing(Scalar edge_length, unsigned int iterations = 10,
        bool use_projection = true);

    //! \brief Perform adaptive remeshing.
    //! \param settings input settings.
    void adaptive_remeshing(const AdaptiveRemeshingSettings& settings);

private:
    void preprocessing();
    void postprocessing();

    void split_long_edges(unsigned int nIterations = 10);
    void collapse_short_edges();
    void tangential_smoothing(unsigned int iterations);
    void project_to_reference(Vertex v);

    bool is_too_long(Vertex v0, Vertex v1) const
    {
        return distance(points_[v0], points_[v1]) >
            4.0 / 3.0 * std::min(vsizing_[v0], vsizing_[v1]);
    }
    bool is_too_short(Vertex v0, Vertex v1) const
    {
        return distance(points_[v0], points_[v1]) <
            4.0 / 5.0 * std::min(vsizing_[v0], vsizing_[v1]);
    }

    ManifoldCurve2D& curve_;
    std::shared_ptr<ManifoldCurve2D> refcurve_;

    std::unique_ptr<EdgeKdTree> kd_tree_;

    bool uniform_;
    Scalar target_edge_length_;
    Scalar min_edge_length_;
    Scalar max_edge_length_;
    Scalar approx_error_;

    bool use_projection_;
    bool has_feature_vertices_{ false };
    VertexProperty<Point2> points_;
    VertexProperty<Point2> vnormal_;
    VertexProperty<bool> vfeature_;
    VertexProperty<bool> vlocked_;
    VertexProperty<Scalar> vsizing_;

    VertexProperty<Point2> refpoints_;
    VertexProperty<Point2> refnormals_;
    VertexProperty<Scalar> refsizing_;

    std::shared_ptr<EvolvingArcLengthCalculator> m_LengthCalculator{ nullptr };
};

} // namespace pmp
