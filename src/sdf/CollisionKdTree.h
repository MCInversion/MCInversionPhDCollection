#pragma once

#include "pmp/SurfaceMesh.h"
#include <vector>

// forward decl
namespace Geometry
{
	struct Ray;
}

namespace SDF
{
	// forward decl
	class CollisionKdTree;

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
	};

	using Triangles = std::vector<Triangle>;

	// function used to find split position and fill child primitive buffers.
	using SplitFunction = std::function<float(
		const BoxSplitData&,               // data of box to be split
		const std::vector<unsigned int>&,  // input face ids 
		std::vector<unsigned int>&,        // output left face ids
		std::vector<unsigned int>&)>;      // output right face ids

	// ===================== Split functions ================================

	[[nodiscard]] float CenterSplitFunction(const BoxSplitData& splitData,
		const std::vector<unsigned int>& facesIn, std::vector<unsigned int>& leftFacesOut, std::vector<unsigned int>& rightFacesOut);

	[[nodiscard]] float AdaptiveSplitFunction(const BoxSplitData& splitData,
		const std::vector<unsigned int>& facesIn, std::vector<unsigned int>& leftFacesOut, std::vector<unsigned int>& rightFacesOut);

	// ======================================================================

	//! \brief A k-d tree for collision detection with pmp::SurfaceMesh triangles
	class CollisionKdTree
	{
	public:
        //! Construct with mesh.
        CollisionKdTree(const pmp::SurfaceMesh& mesh, const SplitFunction& spltFunc);

		// getters
		[[nodiscard]] const std::vector<pmp::vec3>& VertexPositions() const
		{
			return m_VertexPositions;
		}

		[[nodiscard]] const std::vector<Triangle>& TriVertexIds() const
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
         * \brief Counts the number of intersections between a given ray and triangles in this kd-tree.
         * \param ray    intersecting ray.
         * \return the number of intersections.
         */
        [[nodiscard]] unsigned int GetRayTriangleIntersectionCount(Geometry::Ray& ray) const;

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
			float splitPosition{ 0.0f };
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
}
