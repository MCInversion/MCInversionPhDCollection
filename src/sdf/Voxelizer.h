#pragma once

#include "pmp/Types.h"
#include "geometry/CollisionKdTree.h"

#include "pmp/BoundingBox.h"

namespace SDF
{
	//! \brief maximum number of octree children.
	constexpr size_t MAX_OCTREE_CHILDREN = 8;

	//! \brief maximum octree/quadtree depth.
	constexpr size_t MAX_VOXELIZER_DEPTH = 50;

	//! \brief tolerance for leaf size verification.
	constexpr pmp::Scalar LEAF_SIZE_EPSILON = 1e-5;

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
		OctreeVoxelizer(const Geometry::CollisionKdTree& kdTree, const pmp::BoundingBox& startBox, const pmp::Scalar& targetLeafSize);

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
				const pmp::Scalar size = cubeBox.max()[0] - cubeBox.min()[0];
				return (children.empty() && octreeEnvironment->m_LeafSize < size + LEAF_SIZE_EPSILON);
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

		const Geometry::CollisionKdTree& m_KdTree;
		pmp::Scalar m_LeafSize{ 1.0 };
		size_t m_NodeCount{ 0 };

		Node* m_Root{ nullptr };

	};

	//! \brief The maximum number of children in a quadtree.
	constexpr size_t MAX_QUADTREE_CHILDREN = 4;

	//! \brief A Quadtree for generating a "pixel outline" of a mesh using its KD-tree.
	class QuadtreeVoxelizer
	{
	public:
		/**
		 * \brief Constructor. Initialize from a 2D-tree, starting box, and a targetLeafSize.
		 * \param kdTree            a non-null collision 2D-tree to be used for intersection queries.
		 * \param startBox          starting box, will be converted to a square centered at the original box center.
		 * \param targetLeafSize    the leaf size preference for this quadtree. This value will become the pixel size.
		 */
		QuadtreeVoxelizer(const Geometry::Collision2DTree& kdTree, const pmp::BoundingBox2& startBox, const pmp::Scalar& targetLeafSize);

		/**
		 * \brief Collects leaf nodes' (square) boxes and the distance values stored in the leaf nodes themselves.
		 * \param boxBuffer         buffer for square boxes.
		 * \param valueBuffer       buffer for distance values.
		 */
		void GetLeafBoxesAndValues(std::vector<pmp::BoundingBox2*>& boxBuffer, std::vector<double>& valueBuffer) const;

		/// \brief root box getter.
		[[nodiscard]] const pmp::BoundingBox2& GetRootBox() const
		{
			assert(m_Root != nullptr);
			return m_Root->squareBox;
		}

		~QuadtreeVoxelizer()
		{
			delete m_Root;
		}

	private:
		/// \brief a node object of this tree.
		struct Node
		{
			explicit Node(QuadtreeVoxelizer* tree)
				: quadtreeEnvironment(tree)
			{
			}

			~Node()
			{
				for (const auto& ch : children)
					delete ch;
			}

			[[nodiscard]] bool IsALeaf() const
			{
				const pmp::Scalar size = squareBox.max()[0] - squareBox.min()[0];
				return (children.empty() && quadtreeEnvironment->m_LeafSize < size + LEAF_SIZE_EPSILON);
			}

			pmp::BoundingBox2  squareBox{};
			std::vector<Node*> children{};
			double distanceVal{ DBL_MAX };
			QuadtreeVoxelizer* quadtreeEnvironment{ nullptr };
		};

		/**
		 * \brief The recursive part of building this Quadtree.
		 * \param node             node to be initialized.
		 * \param box              bounding box of the node to be initialized.
		 * \param remainingDepth   depth remaining for node construction.
		 */
		void BuildRecurse(Node* node, const pmp::BoundingBox2& box, unsigned int remainingDepth);

		// =====================================================================

		const Geometry::Collision2DTree& m_KdTree;
		pmp::Scalar m_LeafSize{ 1.0 };
		size_t m_NodeCount{ 0 };

		Node* m_Root{ nullptr };

	};

} // namespace SDF