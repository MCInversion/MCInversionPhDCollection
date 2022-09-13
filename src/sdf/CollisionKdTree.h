#pragma once

#include "pmp/SurfaceMesh.h"
#include <vector>

namespace SDF
{
	/**
	 * \brief Preference for split axis during KD tree construction.
	 */
	enum class [[nodiscard]] SplitAxisPreference
	{
		XAxis = 0, //>! the node's box is longest in the x direction, it should be split along the x-axis.
		YAxis = 1, //>! the node's box is longest in the y direction, it should be split along the y-axis.
		ZAxis = 2  //>! the node's box is longest in the z direction, it should be split along the z-axis.
	};

	// forward decl
	class CollisionKdTree;

	/**
	 * \brief helper data item for splitting box.
	 */
	struct BoxSplitData
	{
		CollisionKdTree* kdTree;
		const pmp::BoundingBox* box;
		SplitAxisPreference axis;
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

			SplitAxisPreference axis{0};

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
