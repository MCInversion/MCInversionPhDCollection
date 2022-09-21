#include "OctreeVoxelizer.h"

#include <stack>

#include "geometry/GeometryUtil.h"


namespace SDF
{
	/**
	 * \brief computes the "cube-box" from a given starting box, centered at the midpoint of the starting box center voxel.
	 * \param startBox         bounding box to start from.
	 * \param targetLeafSize   preferred leaf size for the octree voxelizer.
	 * \return cube bounding box for octree voxelizer.
	 * \throw std::logic error if targetLeafSize <= 0.0f.
	 */
	[[nodiscard]] pmp::BoundingBox ComputeCubeBoxFromTargetLeafSize(const pmp::BoundingBox& startBox, const float& targetLeafSize)
	{
		if (targetLeafSize <= 0.0f)
		{
			throw std::logic_error("ComputeCubeBoxFromTargetLeafSize: targetLeafSize <= 0.0f!\n");
		}

		const auto startBoxSize = startBox.max() - startBox.min();
		const float maxDim = std::max({ startBoxSize[0], startBoxSize[1], startBoxSize[2] });
		const unsigned int expectedDepth = std::floor(log2(maxDim / targetLeafSize));
		const float cubeBoxHalfDim = pow(2.0f, expectedDepth) * targetLeafSize;

		// TODO: this math does not make sense!

		const auto startBoxCenter = startBox.center();
		const pmp::vec3 cubeBoxCenter{
			std::round(startBoxCenter[0] / targetLeafSize) * targetLeafSize,
			std::round(startBoxCenter[1] / targetLeafSize) * targetLeafSize,
			std::round(startBoxCenter[2] / targetLeafSize) * targetLeafSize,
		};

		const pmp::vec3 cubeBoxMin{
			cubeBoxCenter[0] - cubeBoxHalfDim,
			cubeBoxCenter[1] - cubeBoxHalfDim,
			cubeBoxCenter[2] - cubeBoxHalfDim
		};
		const pmp::vec3 cubeBoxMax{
			cubeBoxCenter[0] + cubeBoxHalfDim,
			cubeBoxCenter[1] + cubeBoxHalfDim,
			cubeBoxCenter[2] + cubeBoxHalfDim
		};

		// >>>>>>>>> Debugging In progress >>>>>>>>>>>>>>>>>
		std::cout << "OctreeVoxelizer:\n";
		std::cout << "       expectedDepth = " << expectedDepth << "\n";
		std::cout << "       cubeBoxCenter = " << cubeBoxCenter << "\n";
		std::cout << "       cubeBoxHalfDim = " << cubeBoxHalfDim << "\n";
		// <<<<<<<<<< Remove Afterwards <<<<<<<<<<<<<<<<<<<<

		return { cubeBoxMin, cubeBoxMax };
	}

	OctreeVoxelizer::OctreeVoxelizer(const CollisionKdTree& kdTree, const pmp::BoundingBox& startBox, const float& targetLeafSize)
		: m_KdTree(kdTree), m_LeafSize(targetLeafSize)
	{
		m_Root = new Node(this);
		m_NodeCount++;

		const auto cubeBox = ComputeCubeBoxFromTargetLeafSize(startBox, targetLeafSize);
		BuildRecurse(m_Root, cubeBox, MAX_OCTREE_DEPTH);
	}

	// ======================== Helper macros for Octree subdivision ======================

	#define SET_BOX_MIN_COORD(b, B, i, j, k, size)							        \
			b.min()[0] = B.min()[0] + (static_cast<float>(i) / 2.0f) * (size);		\
			b.min()[1] = B.min()[1] + (static_cast<float>(j) / 2.0f) * (size);		\
			b.min()[2] = B.min()[2] + (static_cast<float>(k) / 2.0f) * (size)

	#define SET_BOX_MAX_COORD(b, B, i, j, k, size)							        \
			b.max()[0] = B.min()[0] + (static_cast<float>(i + 1) / 2.0f) * (size);	\
			b.max()[1] = B.min()[1] + (static_cast<float>(j + 1) / 2.0f) * (size);	\
			b.max()[2] = B.min()[2] + (static_cast<float>(k + 1) / 2.0f) * (size)

	void OctreeVoxelizer::BuildRecurse(Node* node, const pmp::BoundingBox& box, unsigned int remainingDepth)
	{
		assert(!box.is_empty());
		node->cubeBox = box;

		const float boxSize = box.max()[0] - box.min()[0];

		if (remainingDepth == 0 || boxSize < m_LeafSize + LEAF_SIZE_EPSILON)
		{
			// ============  process a leaf node ==================

			std::vector<unsigned int> voxelTriangleIds{};
			m_KdTree.GetTrianglesInABox(box, voxelTriangleIds);
			const auto center = box.center();

			assert(!voxelTriangleIds.empty());

			const auto& vertexPositions = m_KdTree.VertexPositions();
			const auto& triangles = m_KdTree.TriVertexIds();

			double distToTriSq = DBL_MAX;
			std::vector triangle{ pmp::vec3(), pmp::vec3(), pmp::vec3() };
			for (const auto& triId : voxelTriangleIds)
			{
				triangle[0] = vertexPositions[triangles[triId].v0Id];
				triangle[1] = vertexPositions[triangles[triId].v1Id];
				triangle[2] = vertexPositions[triangles[triId].v2Id];

				const double currentTriDistSq = Geometry::GetDistanceToTriangleSq(triangle, center);

				if (currentTriDistSq < distToTriSq)
					distToTriSq = currentTriDistSq;
			}
			assert(distToTriSq < DBL_MAX);
			node->distanceVal = sqrt(distToTriSq);
			return;
		}

		// ============ subdivide a non-leaf node ===============

		node->children.reserve(MAX_OCTREE_CHILDREN);
		pmp::BoundingBox childBox{};
		for (unsigned int i = 0; i < 2; i++) 
		{
			for (unsigned int j = 0; j < 2; j++)
			{
				for (unsigned int k = 0; k < 2; k++)
				{
					SET_BOX_MIN_COORD(childBox, box, i, j, k, boxSize);
					SET_BOX_MAX_COORD(childBox, box, i, j, k, boxSize);

					if (!m_KdTree.BoxIntersectsATriangle(childBox))
						continue;

					node->children.emplace_back(new Node(this));
					m_NodeCount++;
					BuildRecurse(node->children[node->children.size() - 1], childBox, remainingDepth - 1);
				}
			}
		}
		node->children.shrink_to_fit();
	}

	void OctreeVoxelizer::GetLeafBoxesAndValues(std::vector<pmp::BoundingBox*>& boxBuffer, std::vector<double>& valueBuffer) const
	{
		assert(boxBuffer.empty() && valueBuffer.empty());

		boxBuffer.reserve(m_NodeCount);
		valueBuffer.reserve(m_NodeCount);

		std::stack<Node*> nodeStack{};
		nodeStack.push(m_Root);

		while (!nodeStack.empty())
		{
			Node* currentNode = nodeStack.top();
			nodeStack.pop();

			if (currentNode->IsALeaf())
			{
				boxBuffer.push_back(&currentNode->cubeBox);
				valueBuffer.push_back(currentNode->distanceVal);
				continue;
			}

			for (auto& child : currentNode->children)
				nodeStack.push(child);
		}

		boxBuffer.shrink_to_fit();
		valueBuffer.shrink_to_fit();
	}

} // namespace SDF