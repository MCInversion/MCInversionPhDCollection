#pragma once

#include <pmp/ManifoldCurve2D.h>
#include <pmp/BoundingBox.h>
#include <vector>
#include <array>
#include <limits>

namespace pmp {

    class EdgeKdTree
    {
    public:
        //! Construct with curve.
        EdgeKdTree(const ManifoldCurve2D& curve, unsigned int max_edges = 10, unsigned int max_depth = 30);

        //! Construct with curve.
        EdgeKdTree(std::shared_ptr<const ManifoldCurve2D> curve, unsigned int max_edges = 10, unsigned int max_depth = 30)
            : EdgeKdTree(*curve, max_edges, max_depth)
        {
        }

        //! Destructor
        ~EdgeKdTree() { delete root_; }

        //! Nearest neighbor information
        struct NearestNeighbor
        {
            Scalar dist;
            Edge edge;
            Point2 nearest;
        };

        //! Return handle of the nearest neighbor
        NearestNeighbor nearest(const Point2& p) const;

    private:
        // Vector of Edges
        using Edges = std::vector<Edge>;

        // Node of the tree: contains parent, children and splitting plane
        struct Node
        {
            Node() = default;

            ~Node()
            {
                delete edges;
                delete left_child;
                delete right_child;
            }

            unsigned char axis;
            Scalar split;
            Edges* edges{ nullptr };
            Node* left_child{ nullptr };
            Node* right_child{ nullptr };
        };

        // Recursive part of build()
        void build_recurse(Node* node, unsigned int max_edges, unsigned int depth);

        // Recursive part of nearest()
        void nearest_recurse(Node* node, const Point2& point, NearestNeighbor& data) const;

        Node* root_;

        std::vector<std::array<Point2, 2>> edge_points_;
    };

} // namespace pmp

