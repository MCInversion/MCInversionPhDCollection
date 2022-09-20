#pragma once

#include "CollisionKdTree.h"

#include "pmp/BoundingBox.h"

namespace SDF
{
	//! \brief maximum number of octree children.
	constexpr size_t MAX_OCTREE_CHILDREN = 8;

	//! \brief maximum octree depth.
	constexpr size_t MAX_OCTREE_DEPTH = 50;

	//! \brief An Octree for generating a "voxel outline" of a mesh using its KD-tree.
	class OctreeVoxelizer
	{
	public:
		/**
		 * \brief Constructor. Initialize from a kd-tree, starting box, and a targetLeafSize.
		 * \param kdTree            a non-null collision kd-tree to be used for intersection queries.
		 * \param startBox          starting box, will be converted to a cube centered at the original box center.
		 * \param targetLeafSize    the leaf size preference for this octree. This value will become the voxel size.
		 */
		OctreeVoxelizer(const CollisionKdTree& kdTree, const pmp::BoundingBox& startBox, const float& targetLeafSize);

		/**
		 * \brief Collects leaf nodes' (cube) boxes and the distance values stored in the leaf nodes themselves.
		 * \param boxBuffer         buffer for cube boxes.
		 * \param valueBuffer       buffer for distance values.
		 */
		void GetLeafBoxesAndValues(std::vector<pmp::BoundingBox*>& boxBuffer, std::vector<double>& valueBuffer) const;

		/// \brief root box getter.
		[[nodiscard]] const pmp::BoundingBox& GetRootBox() const
		{
			assert(m_Root != nullptr);
			return m_Root->cubeBox;
		}

		~OctreeVoxelizer()
		{
			delete m_Root;
		}

	private:
		/// \brief a node object of this tree.
		struct Node
		{
			explicit Node(OctreeVoxelizer* tree)
				: octreeEnvironment(tree)
			{
			}

			~Node()
			{
				for (const auto& ch : children)
					delete ch;
			}

			[[nodiscard]] bool IsALeaf() const
			{
				const float size = cubeBox.max()[0] - cubeBox.min()[0];
				return (children.empty() && octreeEnvironment->m_LeafSize <= size);
			}

			pmp::BoundingBox cubeBox{};
			std::vector<Node*> children{};
			double distanceVal{ DBL_MAX };
			OctreeVoxelizer* octreeEnvironment{ nullptr };
		};

		/**
		 * \brief The recursive part of building this Octree.
		 * \param node             node to be initialized.
		 * \param box              bounding box of the node to be initialized.
		 * \param remainingDepth   depth remaining for node construction.
		 */
		void BuildRecurse(Node* node, const pmp::BoundingBox& box, unsigned int remainingDepth);

		// =====================================================================

		const CollisionKdTree& m_KdTree;
		float m_LeafSize{ 1.0 };
		size_t m_NodeCount{ 0 };

		Node* m_Root{ nullptr };

	};

} // namespace SDF