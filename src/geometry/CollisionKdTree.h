#pragma once

#include "GeometryAdapters.h"

#include <vector>

namespace Geometry
{
	// forward decl
	struct Ray;
	class CollisionKdTree;
	class Collision2DTree;

	//! \brief the amount by which boxes of kd-tree nodes are inflated to account for round-off errors.
	constexpr pmp::Scalar BOX_INFLATION = 1e-6;

	/**
	 * \brief helper data item for splitting box.
	 */
	struct BoxSplitData
	{
		CollisionKdTree* kdTree{ nullptr };
		const pmp::BoundingBox* box{ nullptr };
		unsigned int axis{0};
	};

	/**
	 * \brief primitive data item for mesh triangle containing indices
	 */
	struct Triangle
	{
		unsigned int v0Id;
		unsigned int v1Id;
		unsigned int v2Id;

		unsigned int t0Id{ UINT_MAX }; // a neighbor sharing edge between v0 and v1
		unsigned int t1Id{ UINT_MAX }; // a neighbor sharing edge between v1 and v2
		unsigned int t2Id{ UINT_MAX }; // a neighbor sharing edge between v2 and v0
	};

	using Triangles = std::vector<Triangle>;

	// function used to find split position and fill child primitive buffers.
	using SplitFunction = std::function<pmp::Scalar(
		const BoxSplitData&,               // data of box to be split
		const std::vector<unsigned int>&,  // input face ids 
		std::vector<unsigned int>&,        // output left face ids
		std::vector<unsigned int>&)>;      // output right face ids

	// ===================== Split functions ================================

	[[nodiscard]] pmp::Scalar CenterSplitFunction(const BoxSplitData& splitData,
		const std::vector<unsigned int>& facesIn, std::vector<unsigned int>& leftFacesOut, std::vector<unsigned int>& rightFacesOut);

	[[nodiscard]] pmp::Scalar AdaptiveSplitFunction(const BoxSplitData& splitData,
		const std::vector<unsigned int>& facesIn, std::vector<unsigned int>& leftFacesOut, std::vector<unsigned int>& rightFacesOut);

	// ======================================================================

	//! \brief A k-d tree for collision detection with triangles
	class CollisionKdTree
	{
	public:
        //! Construct with mesh.
        CollisionKdTree(const MeshAdapter& meshAdapter, const SplitFunction& spltFunc);

		// getters
		[[nodiscard]] const std::vector<pmp::vec3>& VertexPositions() const
		{
			return m_VertexPositions;
		}

		[[nodiscard]] const Triangles& TriVertexIds() const
		{
			return m_Triangles;
		}

		/**
		 * \brief Fills a buffer of indices of triangles intersecting a given box.
		 * \param box                  a box whose contents are to be queried.
		 * \param foundTriangleIds     buffer to be filled.
		 */
		void GetTrianglesInABox(const pmp::BoundingBox& box, std::vector<unsigned int>& foundTriangleIds) const;

		/**
		 * \brief A stackless approach to fill a buffer of indices of triangles intersecting a given box.
		 * \param box                  a box whose contents are to be queried.
		 * \param foundTriangleIds     buffer to be filled.
		 */
		void GetTrianglesInABox_Stackless(const pmp::BoundingBox& box, std::vector<unsigned int>& foundTriangleIds) const;

		/**
		 * \brief Performs an intersection test between a given box and any triangle.
		 * \param box    a box to be queried.
		 * \return true if an intersection was detected.
		 */
		[[nodiscard]] bool BoxIntersectsATriangle(const pmp::BoundingBox& box) const;

		/**
		 * \brief A stackless approach for performing an intersection test between a given box and any triangle.
		 * \param box    a box to be queried.
		 * \return true if an intersection was detected.
		 */
		[[nodiscard]] bool BoxIntersectsATriangle_Stackless(const pmp::BoundingBox& box) const;

        /**
         * \brief Performs an intersection test between a given ray and any triangle.
         * \param ray    ray to intersect a triangle.
         * \return true if an intersection was detected.
         */
        [[nodiscard]] bool RayIntersectsATriangle(Geometry::Ray& ray) const;

        /**
         * \brief Counts the number of inward/outward transitions between a given ray and triangles in this kd-tree.
         * \param ray    intersecting ray.
         * \return true if the number of transitions is negative.
         */
        [[nodiscard]] bool IsRayStartPointInsideTriangleMesh(Geometry::Ray& ray) const;

	private:

		/// \brief a node object of this tree.
		struct Node
		{
			Node() = default;

            ~Node()
            {
                delete left_child;
                delete right_child;
            }

			[[nodiscard]] bool IsALeaf() const
            {
				return (left_child == nullptr && right_child == nullptr);
            }

			unsigned int axis{0};
			pmp::Scalar splitPosition{ 0.0 };
			pmp::BoundingBox box{};
			std::vector<unsigned int> triangleIds{};
            Node* left_child{ nullptr };
            Node* right_child{ nullptr };
		};

        /**
		 * \brief The recursive part of building this KD tree.
		 * \param node             node to be initialized.
		 * \param box              bounding box of the node to be initialized.
		 * \param triangleIds      indices of triangles to construct node from.
		 * \param remainingDepth   depth remaining for node construction.
		 */
        void BuildRecurse(Node* node, const pmp::BoundingBox& box, const std::vector<unsigned int>& triangleIds, unsigned int remainingDepth);

        Node* m_Root{ nullptr };

        size_t m_NodeCount{0};

		SplitFunction m_FindSplitAndClassify{};
		std::vector<pmp::vec3> m_VertexPositions{};
		std::vector<Triangle> m_Triangles{};
	};

	// ===========================================================================================

	/**
	 * \brief helper data item for splitting box.
	 */
	struct BoxSplitData2D
	{
		Collision2DTree* kdTree{ nullptr };
		const pmp::BoundingBox2* box{ nullptr };
		unsigned int axis{ 0 };
	};

	/**
	 * \brief primitive data item for an edge of a curve containing vertex indices
	 */
	struct Edge
	{
		unsigned int v0Id;
		unsigned int v1Id;
	};

	using Edges = std::vector<Edge>;

	// function used to find split position and fill child primitive buffers.
	using SplitFunction2D = std::function<pmp::Scalar(
		const BoxSplitData2D&,             // data of box to be split
		const std::vector<unsigned int>&,  // input edge ids 
		std::vector<unsigned int>&,        // output left edge ids
		std::vector<unsigned int>&)>;      // output right edge ids

	// ===================== Split functions ================================

	[[nodiscard]] pmp::Scalar CenterSplitFunction2D(const BoxSplitData2D& splitData,
		const std::vector<unsigned int>& edgesIn, std::vector<unsigned int>& leftEdgesOut, std::vector<unsigned int>& rightEdgesOut);

	//[[nodiscard]] pmp::Scalar AdaptiveSplitFunction2D(const BoxSplitData2D& splitData,
	//	const std::vector<unsigned int>& edgesIn, std::vector<unsigned int>& leftEdgesOut, std::vector<unsigned int>& rightEdgesOut);

	// ======================================================================

	//! \brief A 2D binary tree for collision detection with curve edges (line segments)
	class Collision2DTree
	{
	public:
		//! Construct with mesh.
		Collision2DTree(const CurveAdapter& curveAdapter, const SplitFunction2D& spltFunc);

		// getters
		[[nodiscard]] const std::vector<pmp::Point2>& VertexPositions() const
		{
			return m_VertexPositions;
		}

		[[nodiscard]] const Edges& EdgeVertexIds() const
		{
			return m_Edges;
		}

		/**
		 * \brief Fills a buffer of indices of edge 2D line segments intersecting a given 2D box.
		 * \param box               a box whose contents are to be queried.
		 * \param foundEdgeIds      buffer to be filled.
		 */
		void GetEdgesInABox(const pmp::BoundingBox2& box, std::vector<unsigned int>& foundEdgeIds) const;

		/**
		 * \brief Performs an intersection test between a given box and any 2D line segment logged as a primitive.
		 * \param box    a box to be queried.
		 * \return true if an intersection was detected.
		 */
		[[nodiscard]] bool BoxIntersectsAnEdge(const pmp::BoundingBox2& box) const;


	private:

		/// \brief a node object of this tree.
		struct Node
		{
			Node() = default;

			~Node()
			{
				delete left_child;
				delete right_child;
			}

			[[nodiscard]] bool IsALeaf() const
			{
				return (left_child == nullptr && right_child == nullptr);
			}

			unsigned int axis{ 0 };
			pmp::Scalar splitPosition{ 0.0 };
			pmp::BoundingBox2 box{};
			std::vector<unsigned int> edgeIds{};
			Node* left_child{ nullptr };
			Node* right_child{ nullptr };
		};

		/**
		 * \brief The recursive part of building this 2D tree.
		 * \param node             node to be initialized.
		 * \param box              bounding box of the node to be initialized.
		 * \param edgeIds          indices of edges to construct node from.
		 * \param remainingDepth   depth remaining for node construction.
		 */
		void BuildRecurse(Node* node, const pmp::BoundingBox2& box, const std::vector<unsigned int>& edgeIds, unsigned int remainingDepth);

		Node* m_Root{ nullptr };

		size_t m_NodeCount{ 0 };

		SplitFunction2D m_FindSplitAndClassify{};
		std::vector<pmp::Point2> m_VertexPositions{};
		Edges m_Edges{};

	};
}
