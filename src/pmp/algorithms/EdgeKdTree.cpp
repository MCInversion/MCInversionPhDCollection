#include "pmp/algorithms/EdgeKdTree.h"

#include "pmp/algorithms/DistanceUtils.h"

namespace pmp {

    EdgeKdTree::EdgeKdTree(const ManifoldCurve2D& curve, unsigned int max_edges, unsigned int max_depth)
    {
        // init
        root_ = new Node();
        root_->edges = new Edges();

        // collect edges and points
        root_->edges->reserve(curve.n_edges());
        edge_points_.reserve(curve.n_edges());
        auto points = curve.get_vertex_property<Point2>("v:point");

        for (const auto& e : curve.edges())
        {
            root_->edges->push_back(e);

            auto [v0, v1] = curve.vertices(e);
            const auto& p0 = points[v0];
            const auto& p1 = points[v1];
            edge_points_.push_back({ p0, p1 });
        }

        // call recursive helper
        build_recurse(root_, max_edges, max_depth);
    }

    void EdgeKdTree::build_recurse(Node* node, unsigned int max_edges, unsigned int depth)
    {
        // should we stop at this level?
        if ((depth == 0) || (node->edges->size() <= max_edges))
            return;

        // compute bounding box
        BoundingBox2 bbox;
        for (const auto& e : *node->edges)
        {
            const auto idx = e.idx();
            bbox += edge_points_[idx][0];
            bbox += edge_points_[idx][1];
        }

        // split longest side of bounding box
        Point2 bb = bbox.max() - bbox.min();
        Scalar length = bb[0];
        int axis = 0;
        if (bb[1] > length)
            length = bb[(axis = 1)];

        // split in the middle
        Scalar split = bbox.center()[axis];

        // create children
        auto* left = new Node();
        left->edges = new Edges();
        left->edges->reserve(node->edges->size() / 2);
        auto* right = new Node();
        right->edges = new Edges;
        right->edges->reserve(node->edges->size() / 2);

        // partition for left and right child
        for (const auto& e : *node->edges)
        {
            bool l = false, r = false;

            const auto& pos = edge_points_[e.idx()];

            if (pos[0][axis] <= split)
                l = true;
            else
                r = true;
            if (pos[1][axis] <= split)
                l = true;
            else
                r = true;

            if (l)
            {
                left->edges->push_back(e);
            }

            if (r)
            {
                right->edges->push_back(e);
            }
        }

        // stop here?
        if (left->edges->size() == node->edges->size() ||
            right->edges->size() == node->edges->size())
        {
            // compact my memory
            node->edges->shrink_to_fit();

            // delete new nodes
            delete left;
            delete right;

            // stop recursion
            return;
        }

        // or recurse further?
        else
        {
            // free my memory
            delete node->edges;
            node->edges = nullptr;

            // store internal data
            node->axis = axis;
            node->split = split;
            node->left_child = left;
            node->right_child = right;

            // recurse to children
            build_recurse(node->left_child, max_edges, depth - 1);
            build_recurse(node->right_child, max_edges, depth - 1);
        }
    }

    EdgeKdTree::NearestNeighbor EdgeKdTree::nearest(const Point2& p) const
    {
        NearestNeighbor data;
        data.dist = std::numeric_limits<Scalar>::max();
        nearest_recurse(root_, p, data);
        return data;
    }

    void EdgeKdTree::nearest_recurse(Node* node, const Point2& point, NearestNeighbor& data) const
    {
        // terminal node?
        if (!node->left_child)
        {
            for (const auto& e : *node->edges)
            {
                Point2 n;
                const auto& pos = edge_points_[e.idx()];
                auto d = dist_point_line_segment(point, pos[0], pos[1], n);
                if (d < data.dist)
                {
                    data.dist = d;
                    data.edge = e;
                    data.nearest = n;
                }
            }
        }

        // non-terminal node
        else
        {
            Scalar dist = point[node->axis] - node->split;

            if (dist <= 0.0)
            {
                nearest_recurse(node->left_child, point, data);
                if (fabs(dist) < data.dist)
                    nearest_recurse(node->right_child, point, data);
            }
            else
            {
                nearest_recurse(node->right_child, point, data);
                if (fabs(dist) < data.dist)
                    nearest_recurse(node->left_child, point, data);
            }
        }
    }

} // namespace pmp
