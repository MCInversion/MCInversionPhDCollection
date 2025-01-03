#include "CollisionKdTree.h"

#include "geometry/GeometryUtil.h"

#include "pmp/algorithms/Triangulation.h"

#include <numeric>
#include <stack>

#include <utility>
#include <functional>

// Custom hash function for std::pair
struct PairHash
{
	template <typename T1, typename T2>
	std::size_t operator()(const std::pair<T1, T2>& pair) const
	{
		// Hash the first and second element of the pair and combine them
		auto hash1 = std::hash<T1>{}(pair.first);
		auto hash2 = std::hash<T2>{}(pair.second);
		return hash1 ^ (hash2 << 1); // Combine the two hashes
	}
};

namespace Geometry
{
	//! \brief maximum allowed depth of the CollisionKdTree.
	constexpr unsigned int MAX_DEPTH = 20;

	/**
	 * \brief Converts polygons from a mesh given by its adapter to triangles.
	 * \param meshAdapter     Input mesh adapter.
	 */
	void TriangulateMesh(MeshAdapter& meshAdapter)
	{
		if (!meshAdapter.IsTriangle())
		{
			if (const auto pmpAdapter = dynamic_cast<PMPSurfaceMeshAdapter*>(&meshAdapter)) {
				pmp::SurfaceMesh& mesh = pmpAdapter->GetMesh();
				/*pmp::Triangulation tri(result);
				tri.triangulate();*/
				// TODO: Use Poly2Tri
				for (const auto f : mesh.faces()) {
					if (mesh.valence(f) == 3) continue;
					const auto vBegin = *mesh.vertices(f).begin();
					mesh.split(f, vBegin);
				}
			}
			// TODO: If BaseMeshAdapter needs triangulation, implement that logic here
			// Example:
			// else if (auto baseAdapter = dynamic_cast<BaseMeshAdapter*>(&adapter)) {
			//     TriangulateBaseMesh(baseAdapter->GetBaseMesh());
			// }
			throw std::runtime_error("Geometry::TriangulateMesh: meshAdapter not supported!\n");
		}
	}

	CollisionKdTree::CollisionKdTree(const MeshAdapter& meshAdapter, const SplitFunction& spltFunc)
	{
		auto bbox = meshAdapter.GetBounds();
		const auto bboxSize = bbox.max() - bbox.min();
		bbox.expand(BOX_INFLATION * bboxSize[0], BOX_INFLATION * bboxSize[1], BOX_INFLATION * bboxSize[1]);
		const auto triMesh = meshAdapter.Clone();
		TriangulateMesh(*triMesh);

		// extract vertex positions
		m_VertexPositions = triMesh->GetVertices();

		// extract triangle ids
		const auto triIds = triMesh->GetPolyIndices();
		const size_t nTriangles = triIds.size();
		m_Triangles.reserve(nTriangles);
		for (const auto& tri : triIds)
		{
			m_Triangles.emplace_back(Triangle{ tri[0], tri[1], tri[2] });
		}

		// Find and set neighboring triangles
		std::unordered_map<std::pair<unsigned int, unsigned int>, unsigned int, PairHash> edgeToTriangleMap;
		for (unsigned int i = 0; i < nTriangles; ++i)
		{
			auto& triangle = m_Triangles[i];

			// Define the three edges of the triangle
			std::pair<unsigned int, unsigned int> edges[3] = {
				{triangle.v0Id, triangle.v1Id},
				{triangle.v1Id, triangle.v2Id},
				{triangle.v2Id, triangle.v0Id}
			};

			// Ensure the smaller index comes first in each edge
			for (auto& edge : edges)
			{
				if (edge.first > edge.second)
				{
					std::swap(edge.first, edge.second);
				}
			}

			// Check for neighbors
			for (int j = 0; j < 3; ++j)
			{
				auto& edge = edges[j];
				auto it = edgeToTriangleMap.find(edge);
				if (it != edgeToTriangleMap.end())
				{
					unsigned int neighborIndex = it->second;
					m_Triangles[neighborIndex].t0Id = i; // Link this triangle to its neighbor
					triangle.t0Id = neighborIndex;
					edgeToTriangleMap.erase(it); // Remove edge since it's now matched
				}
				else
				{
					edgeToTriangleMap[edge] = i; // Add this edge to the map
				}
			}
		}

		m_FindSplitAndClassify = spltFunc;

		// generate index array to primitives
		std::vector<unsigned int> triangleIds(m_Triangles.size());
		std::iota(triangleIds.begin(), triangleIds.end(), 0);

		m_Root = new Node();
		m_NodeCount++;
		BuildRecurse(m_Root, bbox, triangleIds, MAX_DEPTH);
	}

	/**
	 * \brief Computes the preference for split axis during KD tree construction.
	 * \param box    evaluated bounding box.
	 * \return index of the preferred axis.
	 */
	[[nodiscard]] unsigned int GetSplitAxisPreference(const pmp::BoundingBox& box)
	{
		if ((box.max()[0] - box.min()[0]) > (box.max()[1] - box.min()[1]) &&
			(box.max()[0] - box.min()[0]) > (box.max()[2] - box.min()[2]))
		{
			return 0; // X-axis;
		}

		if ((box.max()[1] - box.min()[1]) > (box.max()[2] - box.min()[2]))
		{
			return 1; //  Y-axis;
		}

		return 2; // Z-axis;
	}

	//
	// =================== Split functions ==================================================================
	//

	[[nodiscard]] pmp::Scalar TriangleMin(const Triangle& tri, const std::vector<pmp::vec3>& vertices, const unsigned int& axisId)
	{
		pmp::Scalar min = FLT_MAX;

		if (vertices[tri.v0Id][axisId] < min)
			min = vertices[tri.v0Id][axisId];

		if (vertices[tri.v1Id][axisId] < min)
			min = vertices[tri.v1Id][axisId];

		if (vertices[tri.v2Id][axisId] < min)
			min = vertices[tri.v2Id][axisId];

		return min;
	}

	[[nodiscard]] pmp::Scalar TriangleMax(const Triangle& tri, const std::vector<pmp::vec3>& vertices, const unsigned int& axisId)
	{
		pmp::Scalar max = -FLT_MAX;

		if (vertices[tri.v0Id][axisId] > max)
			max = vertices[tri.v0Id][axisId];

		if (vertices[tri.v1Id][axisId] > max)
			max = vertices[tri.v1Id][axisId];

		if (vertices[tri.v2Id][axisId] > max)
			max = vertices[tri.v2Id][axisId];

		return max;
	}

	// simple center split function
	pmp::Scalar CenterSplitFunction(const BoxSplitData& splitData,
		const std::vector<unsigned int>& facesIn, std::vector<unsigned int>& leftFacesOut, std::vector<unsigned int>& rightFacesOut)
	{
		const size_t nFaces = facesIn.size();

		const auto center = splitData.box->center();
		const unsigned int axisId = splitData.axis;

		const auto& vertices = splitData.kdTree->VertexPositions();
		const auto& triVertexIds = splitData.kdTree->TriVertexIds();

		const pmp::Scalar splitPos = center[axisId];
		leftFacesOut.reserve(nFaces);
		rightFacesOut.reserve(nFaces);
		for (const auto& fId : facesIn)
		{
			const pmp::Scalar min = TriangleMin(triVertexIds[fId], vertices, axisId);
			const pmp::Scalar max = TriangleMax(triVertexIds[fId], vertices, axisId);

			if (min <= splitPos)
				leftFacesOut.emplace_back(fId);

			if (max >= splitPos)
				rightFacesOut.emplace_back(fId);
		}
		leftFacesOut.shrink_to_fit();
		rightFacesOut.shrink_to_fit();

		return splitPos;
	}

	/// \brief default number of box sampling points
	constexpr unsigned int BOX_CUTS = 4;

	// TODO: fix this & add a version with intrinsics
	// Fast kd-tree Construction with an Adaptive Error-Bounded Heuristic (Hunt, Mark, Stoll)
	//
	pmp::Scalar AdaptiveSplitFunction(const BoxSplitData& splitData, const std::vector<unsigned int>& facesIn, std::vector<unsigned int>& leftFacesOut, std::vector<unsigned int>& rightFacesOut)
	{
		const auto nFaces = static_cast<unsigned int>(facesIn.size());
		const auto& box = *splitData.box;
		const unsigned int axisId = splitData.axis;
		const auto& vertices = splitData.kdTree->VertexPositions();
		const auto& triVertexIds = splitData.kdTree->TriVertexIds();

		const pmp::Scalar tot_BoxArea =
			2.0 * (box.max()[0] - box.min()[0]) +
			2.0 * (box.max()[1] - box.min()[1]) +
			2.0 * (box.max()[2] - box.min()[2]);

		unsigned int i, j;

		// === Stage 1: Initial sampling of C_L(x) and C_R(x) ========================================================

		// split interval limits:
		const pmp::Scalar a = box.min()[axisId];
		const pmp::Scalar b = box.max()[axisId];

		// splits:
		std::vector<pmp::Scalar> bCutPos = std::vector<pmp::Scalar>(BOX_CUTS + 2);
		for (i = 0; i <= BOX_CUTS + 1; i++) 
			bCutPos[i] = a * (1.0 - (static_cast<pmp::Scalar>(i) / static_cast<pmp::Scalar>(BOX_CUTS + 1))) + b * (static_cast<pmp::Scalar>(i) / static_cast<pmp::Scalar>(BOX_CUTS + 1));

		// set C_L(x) = 0, C_R(x) = 0
		auto C_L = std::vector<unsigned int>(BOX_CUTS + 2);
		auto C_R = std::vector<unsigned int>(BOX_CUTS + 2);

		auto faceMins = std::vector<pmp::Scalar>(nFaces);
		auto faceMaxes = std::vector<pmp::Scalar>(nFaces);

		for (i = 0; i < nFaces; i++)
		{
			faceMins[i] = TriangleMin(triVertexIds[i], vertices, axisId);
			faceMaxes[i] = TriangleMax(triVertexIds[i], vertices, axisId);

			for (j = 1; j <= BOX_CUTS; j++) 
			{
				C_L[j] += (faceMins[i] < bCutPos[j] ? 1 : 0);
				C_R[j] += (faceMaxes[i] > bCutPos[j] ? 1 : 0);
			}
		}
		C_L[5] = C_R[0] = nFaces;
		std::ranges::reverse(C_R); // store in reverse since C_R(x) is non-increasing

		// ===== Stage 2: Sample range [0, nFaces] uniformly & count the number of samples within each segment ======

		auto S_L = std::vector<pmp::Scalar>(BOX_CUTS + 1);
		auto S_R = std::vector<pmp::Scalar>(BOX_CUTS + 1);

		pmp::Scalar ran_s;

		for (i = 0; i < BOX_CUTS; i++) 
		{
			ran_s = static_cast<pmp::Scalar>(i + 1) / static_cast<pmp::Scalar>(BOX_CUTS + 1) * static_cast<pmp::Scalar>(nFaces);

			for (j = 0; j <= BOX_CUTS; j++) 
			{
				S_L[j] += (ran_s > static_cast<pmp::Scalar>(C_L[j]) && ran_s < static_cast<pmp::Scalar>(C_L[j + 1]) ? 1 : 0);
				S_R[j] += (ran_s > static_cast<pmp::Scalar>(C_R[j]) && ran_s < static_cast<pmp::Scalar>(C_R[j + 1]) ? 1 : 0);
			}
		}
		std::ranges::reverse(S_R);

		// ==== Stage 3: add more sampling positions to subdivided segments ===========================================

		auto all_splt_L = std::vector<pmp::Scalar>(2 * BOX_CUTS);
		auto all_splt_R = std::vector<pmp::Scalar>(2 * BOX_CUTS);
		pmp::Scalar segLen = (b - a) / static_cast<pmp::Scalar>(BOX_CUTS + 1);
		unsigned int nSeg_L = 0, nSeg_R = 0;

		for (i = 0; i <= BOX_CUTS; i++)
		{
			if (i > 0) 
			{
				all_splt_L[nSeg_L++] = bCutPos[i];
				all_splt_R[nSeg_R++] = bCutPos[i];
			}

			for (j = 0; j < S_L[i]; j++) 
			{
				all_splt_L[nSeg_L++] = bCutPos[i] + static_cast<pmp::Scalar>(j + 1) / (S_L[i] + 1) * segLen;
			}
			for (j = 0; j < S_R[i]; j++) 
			{
				all_splt_R[nSeg_R++] = bCutPos[i] + static_cast<pmp::Scalar>(j + 1) / (S_R[i] + 1) * segLen;
			}
		}

		// Compute surface area heuristic SAH:
		// remaining two dimensions of the child box candidates
		const pmp::Scalar boxDim0 = box.max()[(axisId + 1) % 3] - box.min()[(axisId + 1) % 3];
		const pmp::Scalar boxDim1 = box.max()[(axisId + 2) % 3] - box.min()[(axisId + 2) % 3];

		auto SA_L = std::vector<pmp::Scalar>(2 * BOX_CUTS);
		auto SA_R = std::vector<pmp::Scalar>(2 * BOX_CUTS);

		// SA_L(x) = (boxDim_L(x) + boxDim0 + boxDim1) * 2.0 / tot_BoxArea
		// SA_R(x) = (boxDim_R(x) + boxDim0 + boxDim1) * 2.0 / tot_BoxArea

		for (i = 0; i < 2 * BOX_CUTS; i++) 
		{
			SA_L[i] = ((all_splt_L[i] - a) + boxDim0 + boxDim1) * 2.0 / tot_BoxArea;
			SA_R[i] = ((b - all_splt_R[i]) + boxDim0 + boxDim1) * 2.0 / tot_BoxArea;
		}

		// ==== Stage 4: RESAMPLE C_L and C_R on all sample points & construct an approximation of cost(x) to minimize
		pmp::Scalar min, max;
		auto cost = std::vector<pmp::Scalar>(2 * BOX_CUTS);

		// cost(x) = C_L(x) * SA_L(x) + C_R(x) * SA_R(x):
		for (i = 0; i < nFaces; i++)
		{
			min = faceMins[i];
			max = faceMaxes[i];

			for (j = 0; j < 2 * BOX_CUTS; j++) 
			{
				cost[j] += static_cast<pmp::Scalar>(min < all_splt_L[j] ? 1 : 0) * SA_L[j] + static_cast<pmp::Scalar>(max > all_splt_R[j] ? 1 : 0) * SA_R[j];
			}
		}

		// ==== Stage 5: Minimize cost(x) & classify primitives  =====================================================

		pmp::Scalar bestSplit = 0.5 * (a + b); // if this loop fails to initialize bestSplit, set it to middle
		pmp::Scalar minCost = FLT_MAX;

		for (i = 1; i < 2 * BOX_CUTS; i++) 
		{
			if (cost[i] >= minCost)
				continue;

			minCost = cost[i];
			bestSplit = all_splt_L[i];
		}

		// fill left and right arrays now that best split position is known:
		leftFacesOut.reserve(nFaces);
		rightFacesOut.reserve(nFaces);
		for (i = 0; i < nFaces; i++)
		{
			min = faceMins[i];
			max = faceMaxes[i];

			if (min <= bestSplit) 
			{
				leftFacesOut.push_back(facesIn[i]);
			}
			if (max >= bestSplit) 
			{
				rightFacesOut.push_back(facesIn[i]);
			}
		}
		leftFacesOut.shrink_to_fit();
		rightFacesOut.shrink_to_fit();

		return bestSplit;
	}

	//
	// ====================================================================================================
	//

	/**
	 * \brief primitive count verification after split.
	 * \param nLeft    number of primitives in the left child node candidate.
	 * \param nRight   number of primitives in the right child node candidate.
	 * \param n        number of primitives in the current node.
	 * \return true if the imbalance is too high.
	 */
	[[nodiscard]] bool ShouldStopBranching(size_t nLeft, size_t nRight, size_t n)
	{
		return static_cast<double>(nLeft + nRight) >= 1.5 * static_cast<double>(n);
	}

	/**
	 * \brief Filters the triangles intersecting a given box.
	 * \param box                    box in question.
	 * \param triangleIds            indices of relevant triangles.
	 * \param triangles              index triples for triangle vertices.
	 * \param vertexPositions        actual vertex positions.
	 * \return true if the imbalance is too high.
	 */
	[[nodiscard]] std::vector<unsigned int> FilterTriangles(
		const pmp::BoundingBox& box, const std::vector<unsigned int>& triangleIds, const Triangles& triangles, const std::vector<pmp::vec3>& vertexPositions)
	{
		const pmp::vec3 boxCenter = box.center();
		const pmp::vec3 boxHalfSize = (box.max() - box.min()) * 0.5;

		std::vector<unsigned int> result{};
		result.reserve(triangleIds.size());
		for (const auto& triId : triangleIds)
		{
			const std::vector triVertices{
				vertexPositions[triangles[triId].v0Id],
				vertexPositions[triangles[triId].v1Id],
				vertexPositions[triangles[triId].v2Id]
			};

			if (!Geometry::TriangleIntersectsBox(triVertices, boxCenter, boxHalfSize))
				continue;

			result.emplace_back(triId);
		}

		result.shrink_to_fit();
		return result;
	}

	[[nodiscard]] pmp::BoundingBox GetChildBox(
		const pmp::BoundingBox& parentBox, const pmp::Scalar& splitPos, const unsigned int& axisId, const bool& isLeft)
	{
		pmp::BoundingBox childBox(parentBox);
		if (isLeft)
			childBox.max()[axisId] = splitPos;
		else
			childBox.min()[axisId] = splitPos;

		const auto chBoxSize = childBox.max() - childBox.min();
		childBox.expand(BOX_INFLATION * chBoxSize[0], BOX_INFLATION * chBoxSize[1], BOX_INFLATION * chBoxSize[2]);

		return childBox;
	}

	//! minimum number of primitives for node
	constexpr size_t MIN_NODE_PRIMITIVE_COUNT = 2;

	void CollisionKdTree::BuildRecurse(Node* node, const pmp::BoundingBox& box, const std::vector<unsigned int>& triangleIds, unsigned int remainingDepth)
	{
		assert(node != nullptr);
		assert(!box.is_empty());
		/*if (box.is_empty())
		{
		    // TODO: fix adaptive resampling!
			std::cout << "CollisionKdTree::BuildRecurse: empty box!\n";
		}*/

		node->box = box;

		if (remainingDepth == 0 || triangleIds.size() <= MIN_NODE_PRIMITIVE_COUNT)
		{
			node->triangleIds = FilterTriangles(box, triangleIds, m_Triangles, m_VertexPositions);
			return;
		}

		const auto axisPreference = GetSplitAxisPreference(box);

		std::vector<unsigned int> leftTriangleIds{};
		std::vector<unsigned int> rightTriangleIds{};
		node->splitPosition = m_FindSplitAndClassify({ this, &box, axisPreference }, triangleIds, leftTriangleIds, rightTriangleIds);

		if (ShouldStopBranching(leftTriangleIds.size(), rightTriangleIds.size(), triangleIds.size()))
		{
			node->triangleIds = FilterTriangles(box, triangleIds, m_Triangles, m_VertexPositions);
			return;
		}

		if (!leftTriangleIds.empty())
		{
			const auto leftBox = GetChildBox(box, node->splitPosition, axisPreference, true);
			node->left_child = new Node();
			m_NodeCount++;
			BuildRecurse(node->left_child, leftBox, leftTriangleIds, remainingDepth - 1);

			if (node->left_child->triangleIds.empty() &&
				node->left_child->left_child == nullptr && node->left_child->right_child == nullptr)
			{
				node->left_child = nullptr;
				m_NodeCount--;
			}
		}

		if (!rightTriangleIds.empty())
		{
			const auto rightBox = GetChildBox(box, node->splitPosition, axisPreference, false);
			node->right_child = new Node();
			m_NodeCount++;
			BuildRecurse(node->right_child, rightBox, rightTriangleIds, remainingDepth - 1);

			if (node->right_child->triangleIds.empty() &&
				node->right_child->left_child == nullptr && node->right_child->right_child == nullptr)
			{
				node->right_child = nullptr;
				m_NodeCount--;
			}
		}
	}

	/**
	 * \brief Computes the average stack height for a binary search tree from node count.
	 * \param nNodes     number of nodes in a binary tree.
	 * \return average stack height.
	 *
	 * Flajolet, Odlyzko - The Average Height of Binary Trees and Other Simple Trees, Journal of Computer and System Sciences, 25, 171-213 (1982).
	 */
	[[nodiscard]] size_t GetAverageStackHeight(const size_t& nNodes)
	{
		return std::round(2.0 * sqrt(M_PI * nNodes));
	}

	/**/
	void CollisionKdTree::GetTrianglesInABox_Stackless(const pmp::BoundingBox& box, std::vector<unsigned int>& foundTriangleIds) const
	{
		assert(foundTriangleIds.empty());

		foundTriangleIds.reserve(m_Triangles.size());
		const size_t expectedStackHeight = GetAverageStackHeight(m_NodeCount);
		std::vector<Node*> nodeStack{};
		nodeStack.reserve(expectedStackHeight);
		nodeStack.emplace_back(m_Root);

		while(!nodeStack.empty())
		{
			Node* currentNode = nodeStack[nodeStack.size() - 1];
			nodeStack.erase(std::prev(nodeStack.end()));

			if (currentNode->IsALeaf())
			{
				for (const auto& triId : currentNode->triangleIds)
					foundTriangleIds.emplace_back(triId);

				continue;
			}

			if (currentNode->left_child)
			{
				if (box.Intersects(currentNode->left_child->box))
					nodeStack.emplace_back(currentNode->left_child);
			}

			if (currentNode->right_child)
			{
				if (box.Intersects(currentNode->right_child->box))
					nodeStack.emplace_back(currentNode->right_child);
			}
		}

		foundTriangleIds.shrink_to_fit();
	}

	void CollisionKdTree::GetTrianglesInABox(const pmp::BoundingBox& box, std::vector<unsigned int>& foundTriangleIds) const
	{
		assert(foundTriangleIds.empty());

		foundTriangleIds.reserve(m_Triangles.size());
		std::stack<Node*> nodeStack{};
		nodeStack.push(m_Root);

		while (!nodeStack.empty())
		{
			const Node* currentNode = nodeStack.top();
			nodeStack.pop();

			if (currentNode->IsALeaf())
			{
				for (const auto& triId : currentNode->triangleIds)
					foundTriangleIds.emplace_back(triId);

				continue;
			}

			if (currentNode->left_child)
			{
				assert(!currentNode->left_child->box.is_empty());
				if (box.Intersects(currentNode->left_child->box))
				{
					nodeStack.push(currentNode->left_child);
				}
			}

			if (currentNode->right_child)
			{
				assert(!currentNode->right_child->box.is_empty());
				if (box.Intersects(currentNode->right_child->box))
				{
					nodeStack.push(currentNode->right_child);	
				}
			}
		}

		foundTriangleIds.shrink_to_fit();
	}

	bool CollisionKdTree::BoxIntersectsATriangle(const pmp::BoundingBox& box) const
	{
		const auto center = box.center();
		const pmp::vec3 halfSize{
			(pmp::Scalar)0.5 * (box.max()[0] - box.min()[0]),
			(pmp::Scalar)0.5 * (box.max()[1] - box.min()[1]),
			(pmp::Scalar)0.5 * (box.max()[2] - box.min()[2])
		};
		std::vector triVerts{ pmp::vec3(), pmp::vec3(), pmp::vec3() };
		std::stack<Node*> nodeStack{};
		nodeStack.push(m_Root);

		while (!nodeStack.empty())
		{
			const Node* currentNode = nodeStack.top();
			nodeStack.pop();

			if (currentNode->IsALeaf())
			{
				for (const auto& triId : currentNode->triangleIds)
				{
					triVerts[0] = m_VertexPositions[m_Triangles[triId].v0Id];
					triVerts[1] = m_VertexPositions[m_Triangles[triId].v1Id];
					triVerts[2] = m_VertexPositions[m_Triangles[triId].v2Id];
					if (Geometry::TriangleIntersectsBox(triVerts, center, halfSize))
						return true;
				}

				continue;
			}

			if (currentNode->left_child)
			{
				assert(!currentNode->left_child->box.is_empty());
				if (box.Intersects(currentNode->left_child->box))
				{
					nodeStack.push(currentNode->left_child);
				}
			}

			if (currentNode->right_child)
			{
				assert(!currentNode->right_child->box.is_empty());
				if (box.Intersects(currentNode->right_child->box))
				{
					nodeStack.push(currentNode->right_child);
				}
			}
		}

		return false;
	}

	bool CollisionKdTree::BoxIntersectsATriangle_Stackless(const pmp::BoundingBox& box) const
	{
		const auto center = box.center();
		const pmp::vec3 halfSize{
			(pmp::Scalar)0.5 * (box.max()[0] - box.min()[0]),
			(pmp::Scalar)0.5 * (box.max()[1] - box.min()[1]),
			(pmp::Scalar)0.5 * (box.max()[2] - box.min()[2])
		};
		std::vector triVerts{ pmp::vec3(), pmp::vec3(), pmp::vec3() };
		const size_t expectedStackHeight = GetAverageStackHeight(m_NodeCount);
		std::vector<Node*> nodeStack{};
		nodeStack.reserve(expectedStackHeight);
		nodeStack.emplace_back(m_Root);

		while (!nodeStack.empty())
		{
			Node* currentNode = nodeStack[nodeStack.size() - 1];
			nodeStack.erase(std::prev(nodeStack.end()));

			if (currentNode->IsALeaf())
			{
				for (const auto& triId : currentNode->triangleIds)
				{
					triVerts[0] = m_VertexPositions[m_Triangles[triId].v0Id];
					triVerts[1] = m_VertexPositions[m_Triangles[triId].v1Id];
					triVerts[2] = m_VertexPositions[m_Triangles[triId].v2Id];
					if (Geometry::TriangleIntersectsBox(triVerts, center, halfSize))
						return true;
				}

				continue;
			}

			if (currentNode->left_child)
			{
				assert(!currentNode->left_child->box.is_empty());
				if (box.Intersects(currentNode->left_child->box))
				{
					nodeStack.emplace_back(currentNode->left_child);
				}
			}

			if (currentNode->right_child)
			{
				assert(!currentNode->right_child->box.is_empty());
				if (box.Intersects(currentNode->right_child->box))
				{
					nodeStack.emplace_back(currentNode->right_child);
				}
			}
		}

		return false;
	}

	bool CollisionKdTree::RayIntersectsATriangle(Geometry::Ray& ray) const
	{
		Node* nearNode = nullptr;
		Node* farNode = nullptr;
		float t_split;
		bool leftIsNear, force_near, force_far, t_splitAtInfinity;
		std::vector triVerts{ pmp::vec3(), pmp::vec3(), pmp::vec3() };

		std::stack<Node*> stack = {};
		stack.push(m_Root);

		while (!stack.empty()) 
		{
			Node* currentNode = stack.top();
			stack.pop();

			if (currentNode->IsALeaf())
			{
				for (const auto& triId : currentNode->triangleIds)
				{
					triVerts[0] = m_VertexPositions[m_Triangles[triId].v0Id];
					triVerts[1] = m_VertexPositions[m_Triangles[triId].v1Id];
					triVerts[2] = m_VertexPositions[m_Triangles[triId].v2Id];
					if (Geometry::RayIntersectsTriangle(ray, triVerts))
						return true;
				}

				continue;
			}

			leftIsNear = ray.StartPt[currentNode->axis] < currentNode->splitPosition;
			nearNode = currentNode->left_child;
			farNode = currentNode->right_child;
			if (!leftIsNear)
			{
				nearNode = currentNode->right_child;
				farNode = currentNode->left_child;
			}

			t_split = (currentNode->splitPosition - ray.StartPt[currentNode->axis]) * ray.InvDirection[currentNode->axis];
			force_near = false;
			force_far = false;
			t_splitAtInfinity = t_split > FLT_MAX || t_split < -FLT_MAX;
			if (t_splitAtInfinity) 
			{
				if (ray.StartPt[currentNode->axis] <= currentNode->splitPosition &&
					ray.StartPt[currentNode->axis] >= currentNode->box.min()[currentNode->axis])
				{
					if (leftIsNear) force_near = true;
					else force_far = true;
				}
				if (ray.StartPt[currentNode->axis] >= currentNode->splitPosition &&
					ray.StartPt[currentNode->axis] <= currentNode->box.max()[currentNode->axis])
				{
					if (leftIsNear) force_near = true;
					else force_far = true;
				}
			}

			if (farNode && ((!t_splitAtInfinity && ray.ParamMax >= t_split) || force_far))
			{
				stack.push(farNode);
			}
			if (nearNode && ((!t_splitAtInfinity && ray.ParamMin <= t_split) || force_near))
			{
				stack.push(nearNode);
			}
		}
		return false;
	}

	static pmp::vec3 ComputeTriangleNormal(const std::vector<pmp::vec3>& triVerts)
	{
		// Calculate the vectors representing two edges of the triangle
		const pmp::vec3 edge1 = triVerts[1] - triVerts[0];
		const pmp::vec3 edge2 = triVerts[2] - triVerts[0];

		// Compute the cross product of the two edges to get the normal vector
		const pmp::vec3 normal = cross(edge1, edge2);

		// Normalize the normal vector to ensure it's a unit vector
		return normalize(normal);
	}

	bool CollisionKdTree::IsRayStartPointInsideTriangleMesh(Ray& ray) const
	{
		int transitionCount = 0;
		std::unordered_set<unsigned int> visitedTriangles;

		if (!RayIntersectsABox(ray, m_Root->box))
			return false;

		std::vector triVerts{ pmp::vec3(), pmp::vec3(), pmp::vec3() };
		std::stack<Node*> stack = {};
		stack.push(m_Root);

		while (!stack.empty())
		{
			const Node* currentNode = stack.top();
			stack.pop();

			if (currentNode->IsALeaf())
			{
				for (const auto& triId : currentNode->triangleIds)
				{
					if (visitedTriangles.contains(triId))
						continue; // Skip if already processed

					triVerts[0] = m_VertexPositions[m_Triangles[triId].v0Id];
					triVerts[1] = m_VertexPositions[m_Triangles[triId].v1Id];
					triVerts[2] = m_VertexPositions[m_Triangles[triId].v2Id];

					if (RayIntersectsTriangle(ray, triVerts))
					{
						// Determine if the ray is entering or exiting the mesh
						pmp::vec3 normal = ComputeTriangleNormal(triVerts);

						if (dot(ray.Direction, normal) < 0)
						{
							// entering
							transitionCount++;
						}
						else
						{
							// exiting
							transitionCount--;
						}

						// Mark the triangle and its neighbors as processed
						visitedTriangles.insert(triId);
						const auto& triangle = m_Triangles[triId];
						if (triangle.t0Id != UINT_MAX) visitedTriangles.insert(triangle.t0Id);
						if (triangle.t1Id != UINT_MAX) visitedTriangles.insert(triangle.t1Id);
						if (triangle.t2Id != UINT_MAX) visitedTriangles.insert(triangle.t2Id);
					}
				}
				continue;
			}

			// currentNode is not a leaf
			if (currentNode->left_child && RayIntersectsABox(ray, currentNode->left_child->box))
			{
				stack.push(currentNode->left_child);
			}

			if (currentNode->right_child && RayIntersectsABox(ray, currentNode->right_child->box))
			{
				stack.push(currentNode->right_child);
			}
		}

		// Final determination: if transitionCount > 0, the point is inside
		return transitionCount < 0;
	}


	//
	// ========================================================================================================================
	//

	/**
	 * \brief Filters the edges intersecting a given box.
	 * \param box                    box in question.
	 * \param edgeIds                indices of relevant edges.
	 * \param edges                  index pairs for edge vertices.
	 * \param vertexPositions        actual vertex positions.
	 * \return true if the imbalance is too high.
	 */
	[[nodiscard]] std::vector<unsigned int> FilterEdges(
		const pmp::BoundingBox2& box, const std::vector<unsigned int>& edgeIds, const Edges& edges, const std::vector<pmp::Point2>& vertexPositions)
	{
		const pmp::Point2 boxCenter = box.center();
		const pmp::vec2 boxHalfSize = (box.max() - box.min()) * 0.5;

		std::vector<unsigned int> result{};
		result.reserve(edgeIds.size());
		for (const auto& edgeId : edgeIds)
		{
			const std::vector edgeVertices{
				vertexPositions[edges[edgeId].v0Id],
				vertexPositions[edges[edgeId].v1Id]
			};

			if (!Geometry::Line2DIntersectsBox(edgeVertices, boxCenter, boxHalfSize))
				continue;

			result.emplace_back(edgeId);
		}

		result.shrink_to_fit();
		return result;
	}

	[[nodiscard]] pmp::BoundingBox2 GetChildBox(
		const pmp::BoundingBox2& parentBox, const pmp::Scalar& splitPos, const unsigned int& axisId, const bool& isLeft)
	{
		pmp::BoundingBox2 childBox(parentBox);
		if (isLeft)
			childBox.max()[axisId] = splitPos;
		else
			childBox.min()[axisId] = splitPos;

		const auto chBoxSize = childBox.max() - childBox.min();
		childBox.expand(BOX_INFLATION * chBoxSize[0], BOX_INFLATION * chBoxSize[1]);

		return childBox;
	}

	Collision2DTree::Collision2DTree(const CurveAdapter& curveAdapter, const SplitFunction2D& spltFunc)
	{
		auto bbox = curveAdapter.GetBounds();
		const auto bboxSize = bbox.max() - bbox.min();
		bbox.expand(BOX_INFLATION * bboxSize[0], BOX_INFLATION * bboxSize[1]);
		const auto curve = curveAdapter.Clone();

		// extract vertex positions
		m_VertexPositions = curve->GetVertices();

		// extract edge ids
		const auto eIds = curve->GetEdgeIndices();
		const size_t nEdges = eIds.size();
		m_Edges.reserve(nEdges);
		for (const auto& e : eIds)
		{
			m_Edges.emplace_back(Edge{ e.first, e.second});
		}

		m_FindSplitAndClassify = spltFunc;

		// generate index array to primitives
		std::vector<unsigned int> edgeIds(m_Edges.size());
		std::iota(edgeIds.begin(), edgeIds.end(), 0);

		m_Root = new Node();
		m_NodeCount++;
		BuildRecurse(m_Root, bbox, edgeIds, MAX_DEPTH);
	}

	void Collision2DTree::GetEdgesInABox(const pmp::BoundingBox2& box, std::vector<unsigned int>& foundEdgeIds) const
	{
		assert(foundEdgeIds.empty());

		foundEdgeIds.reserve(m_Edges.size());
		std::stack<Node*> nodeStack{};
		nodeStack.push(m_Root);

		while (!nodeStack.empty())
		{
			const Node* currentNode = nodeStack.top();
			nodeStack.pop();

			if (currentNode->IsALeaf())
			{
				for (const auto& eId : currentNode->edgeIds)
					foundEdgeIds.emplace_back(eId);

				continue;
			}

			if (currentNode->left_child)
			{
				assert(!currentNode->left_child->box.is_empty());
				if (box.Intersects(currentNode->left_child->box))
				{
					nodeStack.push(currentNode->left_child);
				}
			}

			if (currentNode->right_child)
			{
				assert(!currentNode->right_child->box.is_empty());
				if (box.Intersects(currentNode->right_child->box))
				{
					nodeStack.push(currentNode->right_child);
				}
			}
		}

		foundEdgeIds.shrink_to_fit();
	}

	bool Collision2DTree::BoxIntersectsAnEdge(const pmp::BoundingBox2& box) const
	{
		const auto center = box.center();
		const pmp::vec2 halfSize{
			(pmp::Scalar)0.5 * (box.max()[0] - box.min()[0]),
			(pmp::Scalar)0.5 * (box.max()[1] - box.min()[1])
		};
		std::vector eVerts{ pmp::vec2(), pmp::vec2() };
		std::stack<Node*> nodeStack{};
		nodeStack.push(m_Root);

		while (!nodeStack.empty())
		{
			const Node* currentNode = nodeStack.top();
			nodeStack.pop();

			if (currentNode->IsALeaf())
			{
				for (const auto& eId : currentNode->edgeIds)
				{
					eVerts[0] = m_VertexPositions[m_Edges[eId].v0Id];
					eVerts[1] = m_VertexPositions[m_Edges[eId].v1Id];
					if (Geometry::Line2DIntersectsBox(eVerts, center, halfSize))
						return true;
				}

				continue;
			}

			if (currentNode->left_child)
			{
				assert(!currentNode->left_child->box.is_empty());
				if (box.Intersects(currentNode->left_child->box))
				{
					nodeStack.push(currentNode->left_child);
				}
			}

			if (currentNode->right_child)
			{
				assert(!currentNode->right_child->box.is_empty());
				if (box.Intersects(currentNode->right_child->box))
				{
					nodeStack.push(currentNode->right_child);
				}
			}
		}

		return false;
	}

	// ====================================================

	/**
	 * \brief Computes the preference for split axis during 2D tree construction.
	 * \param box    evaluated bounding box.
	 * \return index of the preferred axis.
	 */
	[[nodiscard]] unsigned int GetSplitAxisPreference(const pmp::BoundingBox2& box)
	{
		if ((box.max()[0] - box.min()[0]) > (box.max()[1] - box.min()[1]))
		{
			return 0; // X-axis;
		}

		return 1; //  Y-axis;
	}

	[[nodiscard]] pmp::Scalar EdgeMin(const Edge& e, const std::vector<pmp::vec2>& vertices, const unsigned int& axisId)
	{
		pmp::Scalar min = FLT_MAX;

		if (vertices[e.v0Id][axisId] < min)
			min = vertices[e.v0Id][axisId];

		if (vertices[e.v1Id][axisId] < min)
			min = vertices[e.v1Id][axisId];

		return min;
	}

	[[nodiscard]] pmp::Scalar EdgeMax(const Edge& e, const std::vector<pmp::vec2>& vertices, const unsigned int& axisId)
	{
		pmp::Scalar max = -FLT_MAX;

		if (vertices[e.v0Id][axisId] > max)
			max = vertices[e.v0Id][axisId];

		if (vertices[e.v1Id][axisId] > max)
			max = vertices[e.v1Id][axisId];

		return max;
	}

	pmp::Scalar CenterSplitFunction2D(const BoxSplitData2D& splitData, const std::vector<unsigned int>& edgesIn, std::vector<unsigned int>& leftEdgesOut, std::vector<unsigned int>& rightEdgesOut)
	{
		const size_t nEdges = edgesIn.size();

		const auto center = splitData.box->center();
		const unsigned int axisId = splitData.axis;

		const auto& vertices = splitData.kdTree->VertexPositions();
		const auto& edgeVertexIds = splitData.kdTree->EdgeVertexIds();

		const pmp::Scalar splitPos = center[axisId];
		leftEdgesOut.reserve(nEdges);
		rightEdgesOut.reserve(nEdges);
		for (const auto& eId : edgesIn)
		{
			const pmp::Scalar min = EdgeMin(edgeVertexIds[eId], vertices, axisId);
			const pmp::Scalar max = EdgeMax(edgeVertexIds[eId], vertices, axisId);

			if (min <= splitPos)
				leftEdgesOut.emplace_back(eId);

			if (max >= splitPos)
				rightEdgesOut.emplace_back(eId);
		}
		leftEdgesOut.shrink_to_fit();
		rightEdgesOut.shrink_to_fit();

		return splitPos;
	}

	//pmp::Scalar AdaptiveSplitFunction2D(const BoxSplitData2D& splitData, const std::vector<unsigned int>& edgesIn, std::vector<unsigned int>& leftEdgesOut, std::vector<unsigned int>& rightEdgesOut)
	//{
	//	return 0.0;
	//}

	void Collision2DTree::BuildRecurse(Node* node, const pmp::BoundingBox2& box, const std::vector<unsigned int>& edgeIds, unsigned int remainingDepth)
	{
		assert(node != nullptr);
		assert(!box.is_empty());

		node->box = box;

		if (remainingDepth == 0 || edgeIds.size() <= MIN_NODE_PRIMITIVE_COUNT)
		{
			node->edgeIds = FilterEdges(box, edgeIds, m_Edges, m_VertexPositions);
			return;
		}

		const auto axisPreference = GetSplitAxisPreference(box);

		std::vector<unsigned int> leftEdgeIds{};
		std::vector<unsigned int> rightEdgeIds{};
		node->splitPosition = m_FindSplitAndClassify({ this, &box, axisPreference }, edgeIds, leftEdgeIds, rightEdgeIds);

		if (ShouldStopBranching(leftEdgeIds.size(), rightEdgeIds.size(), edgeIds.size()))
		{
			node->edgeIds = FilterEdges(box, edgeIds, m_Edges, m_VertexPositions);
			return;
		}

		if (!leftEdgeIds.empty())
		{
			const auto leftBox = GetChildBox(box, node->splitPosition, axisPreference, true);
			node->left_child = new Node();
			m_NodeCount++;
			BuildRecurse(node->left_child, leftBox, leftEdgeIds, remainingDepth - 1);

			if (node->left_child->edgeIds.empty() &&
				node->left_child->left_child == nullptr && node->left_child->right_child == nullptr)
			{
				node->left_child = nullptr;
				m_NodeCount--;
			}
		}

		if (!rightEdgeIds.empty())
		{
			const auto rightBox = GetChildBox(box, node->splitPosition, axisPreference, false);
			node->right_child = new Node();
			m_NodeCount++;
			BuildRecurse(node->right_child, rightBox, rightEdgeIds, remainingDepth - 1);

			if (node->right_child->edgeIds.empty() &&
				node->right_child->left_child == nullptr && node->right_child->right_child == nullptr)
			{
				node->right_child = nullptr;
				m_NodeCount--;
			}
		}
	}

} // namespace SDF