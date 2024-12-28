#include "Voxelizer.h"

#include <stack>

#include "geometry/GeometryUtil.h"


namespace SDF
{
	/**
	 * \brief computes the "cube-box" from a given starting box, centered at the midpoint of the starting box center voxel.
	 * \param startBox         bounding box to start from.
	 * \param targetLeafSize   preferred leaf size for the octree voxelizer.
	 * \return cube bounding box for octree voxelizer.
	 * \throw std::logic error if targetLeafSize <= 0.0.
	 */
	[[nodiscard]] pmp::BoundingBox ComputeCubeBoxFromTargetLeafSize(const pmp::BoundingBox& startBox, const pmp::Scalar& targetLeafSize)
	{
		if (targetLeafSize <= 0.0)
		{
			throw std::logic_error("ComputeCubeBoxFromTargetLeafSize: targetLeafSize <= 0.0!\n");
		}

		const auto startBoxSize = startBox.max() - startBox.min();
		const pmp::Scalar maxDim = std::max({ startBoxSize[0], startBoxSize[1], startBoxSize[2] });
		const unsigned int expectedDepth = std::floor(log2(maxDim / targetLeafSize));
		const float cubeBoxHalfDim = pow(2.0, expectedDepth) * targetLeafSize;

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

		return { cubeBoxMin, cubeBoxMax };
	}

	OctreeVoxelizer::OctreeVoxelizer(const Geometry::CollisionKdTree& kdTree, const pmp::BoundingBox& startBox, const pmp::Scalar& targetLeafSize)
		: m_KdTree(kdTree), m_LeafSize(targetLeafSize)
	{
		m_Root = new Node(this);
		m_NodeCount++;

		const auto cubeBox = ComputeCubeBoxFromTargetLeafSize(startBox, targetLeafSize);
		BuildRecurse(m_Root, cubeBox, MAX_VOXELIZER_DEPTH);
	}

	// ======================== Helper macros for Octree subdivision ======================

	#define SET_BOX_MIN_COORD(b, B, i, j, k, size)							        \
			b.min()[0] = B.min()[0] + (static_cast<pmp::Scalar>(i) / 2.0) * (size);		\
			b.min()[1] = B.min()[1] + (static_cast<pmp::Scalar>(j) / 2.0) * (size);		\
			b.min()[2] = B.min()[2] + (static_cast<pmp::Scalar>(k) / 2.0) * (size); \
			 do {} while (0)

	#define SET_BOX_MAX_COORD(b, B, i, j, k, size)							        \
			b.max()[0] = B.min()[0] + (static_cast<pmp::Scalar>(i + 1) / 2.0) * (size);	\
			b.max()[1] = B.min()[1] + (static_cast<pmp::Scalar>(j + 1) / 2.0) * (size);	\
			b.max()[2] = B.min()[2] + (static_cast<pmp::Scalar>(k + 1) / 2.0) * (size); \
			 do {} while (0)

	void OctreeVoxelizer::BuildRecurse(Node* node, const pmp::BoundingBox& box, unsigned int remainingDepth)
	{
		assert(!box.is_empty());
		node->cubeBox = box;

		const pmp::Scalar boxSize = box.max()[0] - box.min()[0];

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

	//
	// ================================ 2D Quadtree Implementations ========================
	//

	/**
	 * \brief computes the "square-box" from a given starting box, centered at the midpoint of the starting box center pixel.
	 * \param startBox         bounding box to start from.
	 * \param targetLeafSize   preferred leaf size for the quadtree voxelizer.
	 * \return cube bounding box for quadtree voxelizer.
	 * \throw std::logic error if targetLeafSize <= 0.0.
	 */
	[[nodiscard]] pmp::BoundingBox2 ComputeSquareBoxFromTargetLeafSize(const pmp::BoundingBox2& startBox, const pmp::Scalar& targetLeafSize)
	{
		if (targetLeafSize <= 0.0)
		{
			throw std::logic_error("ComputeSquareBoxFromTargetLeafSize: targetLeafSize <= 0.0!\n");
		}

		const auto startBoxSize = startBox.max() - startBox.min();
		const pmp::Scalar maxDim = std::max(startBoxSize[0], startBoxSize[1] );
		const unsigned int expectedDepth = std::floor(log2(maxDim / targetLeafSize));
		const pmp::Scalar squareBoxHalfDim = pow(2.0, expectedDepth) * targetLeafSize;

		const auto startBoxCenter = startBox.center();
		const pmp::vec2 squareBoxCenter{
			std::round(startBoxCenter[0] / targetLeafSize) * targetLeafSize,
			std::round(startBoxCenter[1] / targetLeafSize) * targetLeafSize
		};

		const pmp::vec2 squareBoxMin{
			squareBoxCenter[0] - squareBoxHalfDim,
			squareBoxCenter[1] - squareBoxHalfDim
		};
		const pmp::vec2 squareBoxMax{
			squareBoxCenter[0] + squareBoxHalfDim,
			squareBoxCenter[1] + squareBoxHalfDim
		};

		return { squareBoxMin, squareBoxMax };
	}

	QuadtreeVoxelizer::QuadtreeVoxelizer(const Geometry::Collision2DTree& kdTree, const pmp::BoundingBox2& startBox, const pmp::Scalar& targetLeafSize)
		: m_KdTree(kdTree), m_LeafSize(targetLeafSize)
	{
		m_Root = new Node(this);
		m_NodeCount++;

		const auto squareBox = ComputeSquareBoxFromTargetLeafSize(startBox, targetLeafSize);
		BuildRecurse(m_Root, squareBox, MAX_VOXELIZER_DEPTH);
	}

	// ======================== Helper macros for Quadtree subdivision ======================

#define SET_BOX_MIN_COORD_2D(b, B, i, j, size)							            \
			b.min()[0] = B.min()[0] + (static_cast<pmp::Scalar>(i) / 2.0) * (size);		\
			b.min()[1] = B.min()[1] + (static_cast<pmp::Scalar>(j) / 2.0) * (size); \
			do {} while (0)

#define SET_BOX_MAX_COORD_2D(b, B, i, j, size)							            \
			b.max()[0] = B.min()[0] + (static_cast<pmp::Scalar>(i + 1) / 2.0) * (size);	\
			b.max()[1] = B.min()[1] + (static_cast<pmp::Scalar>(j + 1) / 2.0) * (size); \
			do {} while (0)

	void QuadtreeVoxelizer::BuildRecurse(Node* node, const pmp::BoundingBox2& box, unsigned int remainingDepth)
	{
		assert(!box.is_empty());
		node->squareBox = box;

		const pmp::Scalar boxSize = box.max()[0] - box.min()[0];

		if (remainingDepth == 0 || boxSize < m_LeafSize + LEAF_SIZE_EPSILON)
		{
			// ============  process a leaf node ==================

			std::vector<unsigned int> pixelEdgeIds{};
			m_KdTree.GetEdgesInABox(box, pixelEdgeIds);
			const auto center = box.center();

			assert(!pixelEdgeIds.empty());

			const auto& vertexPositions = m_KdTree.VertexPositions();
			const auto& edges = m_KdTree.EdgeVertexIds();

			double distToEdgeSq = DBL_MAX;
			std::vector line{ pmp::vec2(), pmp::vec2() };
			for (const auto& edgeId : pixelEdgeIds)
			{
				line[0] = vertexPositions[edges[edgeId].v0Id];
				line[1] = vertexPositions[edges[edgeId].v1Id];

				const double currentEdgeDistSq = Geometry::GetDistanceToLine2DSq(line, center);

				if (currentEdgeDistSq < distToEdgeSq)
					distToEdgeSq = currentEdgeDistSq;
			}
			assert(distToEdgeSq < DBL_MAX);
			node->distanceVal = sqrt(distToEdgeSq);
			return;
		}

		// ============ subdivide a non-leaf node ===============

		node->children.reserve(MAX_QUADTREE_CHILDREN);
		pmp::BoundingBox2 childBox{};
		for (unsigned int i = 0; i < 2; i++)
		{
			for (unsigned int j = 0; j < 2; j++)
			{
				SET_BOX_MIN_COORD_2D(childBox, box, i, j, boxSize);
				SET_BOX_MAX_COORD_2D(childBox, box, i, j, boxSize);

				if (!m_KdTree.BoxIntersectsAnEdge(childBox))
					continue;

				node->children.emplace_back(new Node(this));
				m_NodeCount++;
				BuildRecurse(node->children[node->children.size() - 1], childBox, remainingDepth - 1);
			}
		}
		node->children.shrink_to_fit();
	}

	void QuadtreeVoxelizer::GetLeafBoxesAndValues(std::vector<pmp::BoundingBox2*>& boxBuffer, std::vector<double>& valueBuffer) const
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
				boxBuffer.push_back(&currentNode->squareBox);
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