#include "CollisionKdTree.h"

#include "geometry/GeometryUtil.h"

#include "pmp/algorithms/Triangulation.h"

#include <numeric>
#include <stack>


namespace Geometry
{
	//! \brief the amount by which boxes of kd-tree nodes are inflated to account for round-off errors.
	constexpr float BOX_INFLATION = 1e-6f;

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
		bbox.expand(BOX_INFLATION, BOX_INFLATION, BOX_INFLATION);
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

	[[nodiscard]] float TriangleMin(const Triangle& tri, const std::vector<pmp::vec3>& vertices, const unsigned int& axisId)
	{
		float min = FLT_MAX;

		if (vertices[tri.v0Id][axisId] < min)
			min = vertices[tri.v0Id][axisId];

		if (vertices[tri.v1Id][axisId] < min)
			min = vertices[tri.v1Id][axisId];

		if (vertices[tri.v2Id][axisId] < min)
			min = vertices[tri.v2Id][axisId];

		return min;
	}

	[[nodiscard]] float TriangleMax(const Triangle& tri, const std::vector<pmp::vec3>& vertices, const unsigned int& axisId)
	{
		float max = -FLT_MAX;

		if (vertices[tri.v0Id][axisId] > max)
			max = vertices[tri.v0Id][axisId];

		if (vertices[tri.v1Id][axisId] > max)
			max = vertices[tri.v1Id][axisId];

		if (vertices[tri.v2Id][axisId] > max)
			max = vertices[tri.v2Id][axisId];

		return max;
	}

	// simple center split function
	float CenterSplitFunction(const BoxSplitData& splitData, 
		const std::vector<unsigned int>& facesIn, std::vector<unsigned int>& leftFacesOut, std::vector<unsigned int>& rightFacesOut)
	{
		const size_t nFaces = facesIn.size();

		const auto center = splitData.box->center();
		const unsigned int axisId = splitData.axis;

		const auto& vertices = splitData.kdTree->VertexPositions();
		const auto& triVertexIds = splitData.kdTree->TriVertexIds();

		const float splitPos = center[axisId];
		leftFacesOut.reserve(nFaces);
		rightFacesOut.reserve(nFaces);
		for (const auto& fId : facesIn)
		{
			const float min = TriangleMin(triVertexIds[fId], vertices, axisId);
			const float max = TriangleMax(triVertexIds[fId], vertices, axisId);

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
	float AdaptiveSplitFunction(const BoxSplitData& splitData, const std::vector<unsigned int>& facesIn, std::vector<unsigned int>& leftFacesOut, std::vector<unsigned int>& rightFacesOut)
	{
		const auto nFaces = static_cast<unsigned int>(facesIn.size());
		const auto& box = *splitData.box;
		const unsigned int axisId = splitData.axis;
		const auto& vertices = splitData.kdTree->VertexPositions();
		const auto& triVertexIds = splitData.kdTree->TriVertexIds();

		const float tot_BoxArea = 
			2.0f * (box.max()[0] - box.min()[0]) +
			2.0f * (box.max()[1] - box.min()[1]) +
			2.0f * (box.max()[2] - box.min()[2]);

		unsigned int i, j;

		// === Stage 1: Initial sampling of C_L(x) and C_R(x) ========================================================

		// split interval limits:
		const float a = box.min()[axisId];
		const float b = box.max()[axisId];

		// splits:
		std::vector<float> bCutPos = std::vector<float>(BOX_CUTS + 2);
		for (i = 0; i <= BOX_CUTS + 1; i++) 
			bCutPos[i] = a * (1.0f - (static_cast<float>(i) / static_cast<float>(BOX_CUTS + 1))) + b * (static_cast<float>(i) / static_cast<float>(BOX_CUTS + 1));

		// set C_L(x) = 0, C_R(x) = 0
		auto C_L = std::vector<unsigned int>(BOX_CUTS + 2);
		auto C_R = std::vector<unsigned int>(BOX_CUTS + 2);

		auto faceMins = std::vector<float>(nFaces);
		auto faceMaxes = std::vector<float>(nFaces);

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

		auto S_L = std::vector<float>(BOX_CUTS + 1);
		auto S_R = std::vector<float>(BOX_CUTS + 1);

		float ran_s;

		for (i = 0; i < BOX_CUTS; i++) 
		{
			ran_s = static_cast<float>(i + 1) / static_cast<float>(BOX_CUTS + 1) * static_cast<float>(nFaces);

			for (j = 0; j <= BOX_CUTS; j++) 
			{
				S_L[j] += (ran_s > static_cast<float>(C_L[j]) && ran_s < static_cast<float>(C_L[j + 1]) ? 1 : 0);
				S_R[j] += (ran_s > static_cast<float>(C_R[j]) && ran_s < static_cast<float>(C_R[j + 1]) ? 1 : 0);
			}
		}
		std::ranges::reverse(S_R);

		// ==== Stage 3: add more sampling positions to subdivided segments ===========================================

		auto all_splt_L = std::vector<float>(2 * BOX_CUTS);
		auto all_splt_R = std::vector<float>(2 * BOX_CUTS);
		float segLen = (b - a) / static_cast<float>(BOX_CUTS + 1);
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
				all_splt_L[nSeg_L++] = bCutPos[i] + static_cast<float>(j + 1) / (S_L[i] + 1) * segLen;
			}
			for (j = 0; j < S_R[i]; j++) 
			{
				all_splt_R[nSeg_R++] = bCutPos[i] + static_cast<float>(j + 1) / (S_R[i] + 1) * segLen;
			}
		}

		// Compute surface area heuristic SAH:
		// remaining two dimensions of the child box candidates
		const float boxDim0 = box.max()[(axisId + 1) % 3] - box.min()[(axisId + 1) % 3];
		const float boxDim1 = box.max()[(axisId + 2) % 3] - box.min()[(axisId + 2) % 3];

		auto SA_L = std::vector<float>(2 * BOX_CUTS);
		auto SA_R = std::vector<float>(2 * BOX_CUTS);

		// SA_L(x) = (boxDim_L(x) + boxDim0 + boxDim1) * 2.0 / tot_BoxArea
		// SA_R(x) = (boxDim_R(x) + boxDim0 + boxDim1) * 2.0 / tot_BoxArea

		for (i = 0; i < 2 * BOX_CUTS; i++) 
		{
			SA_L[i] = ((all_splt_L[i] - a) + boxDim0 + boxDim1) * 2.0f / tot_BoxArea;
			SA_R[i] = ((b - all_splt_R[i]) + boxDim0 + boxDim1) * 2.0f / tot_BoxArea;
		}

		// ==== Stage 4: RESAMPLE C_L and C_R on all sample points & construct an approximation of cost(x) to minimize
		float min, max;
		auto cost = std::vector<float>(2 * BOX_CUTS);

		// cost(x) = C_L(x) * SA_L(x) + C_R(x) * SA_R(x):
		for (i = 0; i < nFaces; i++)
		{
			min = faceMins[i];
			max = faceMaxes[i];

			for (j = 0; j < 2 * BOX_CUTS; j++) 
			{
				cost[j] += static_cast<float>(min < all_splt_L[j] ? 1 : 0) * SA_L[j] + static_cast<float>(max > all_splt_R[j] ? 1 : 0) * SA_R[j];
			}
		}

		// ==== Stage 5: Minimize cost(x) & classify primitives  =====================================================

		float bestSplit = 0.5f * (a + b); // if this loop fails to initialize bestSplit, set it to middle
		float minCost = FLT_MAX;

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
		const pmp::BoundingBox& parentBox, const float& splitPos, const unsigned int& axisId, const bool& isLeft)
	{
		pmp::BoundingBox childBox(parentBox);
		if (isLeft)
			childBox.max()[axisId] = splitPos;
		else
			childBox.min()[axisId] = splitPos;

		childBox.expand(BOX_INFLATION, BOX_INFLATION, BOX_INFLATION);

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
			0.5f * (box.max()[0] - box.min()[0]),
			0.5f * (box.max()[1] - box.min()[1]),
			0.5f * (box.max()[2] - box.min()[2])
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
			0.5f * (box.max()[0] - box.min()[0]),
			0.5f * (box.max()[1] - box.min()[1]),
			0.5f * (box.max()[2] - box.min()[2])
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

	unsigned int CollisionKdTree::GetRayTriangleIntersectionCount(Geometry::Ray& ray) const
	{
		unsigned int hitCount = 0;
		if (!Geometry::RayIntersectsABox(ray, m_Root->box))
			return hitCount;

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
					triVerts[0] = m_VertexPositions[m_Triangles[triId].v0Id];
					triVerts[1] = m_VertexPositions[m_Triangles[triId].v1Id];
					triVerts[2] = m_VertexPositions[m_Triangles[triId].v2Id];
					if (Geometry::RayIntersectsTriangle(ray, triVerts))
					{
						hitCount++;
					}
				}
				continue;
			}

			// currentNode is not a leaf

			if (currentNode->left_child)
			{
				if (Geometry::RayIntersectsABox(ray, currentNode->left_child->box))
				{
					stack.push(currentNode->left_child);
				}				
			}

			if (currentNode->right_child)
			{
				if (Geometry::RayIntersectsABox(ray, currentNode->right_child->box))
				{
					stack.push(currentNode->right_child);
				}				
			}
		}

		return hitCount;
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
		const pmp::BoundingBox2& box, const std::vector<unsigned int>& edgeIds, const Triangles& edges, const std::vector<pmp::Point2>& vertexPositions)
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

			if (!Geometry::LineIntersectsBox(edgeVertices, boxCenter, boxHalfSize))
				continue;

			result.emplace_back(edgeId);
		}

		result.shrink_to_fit();
		return result;
	}

	[[nodiscard]] pmp::BoundingBox2 GetChildBox(
		const pmp::BoundingBox2& parentBox, const float& splitPos, const unsigned int& axisId, const bool& isLeft)
	{
		pmp::BoundingBox2 childBox(parentBox);
		if (isLeft)
			childBox.max()[axisId] = splitPos;
		else
			childBox.min()[axisId] = splitPos;

		childBox.expand(BOX_INFLATION, BOX_INFLATION);

		return childBox;
	}

	float CenterSplitFunction2D(const BoxSplitData2D& splitData, const std::vector<unsigned int>& edgesIn, std::vector<unsigned int>& leftEdgesOut, std::vector<unsigned int>& rightEdgesOut)
	{
		return 0.0f;
	}

	float AdaptiveSplitFunction2D(const BoxSplitData2D& splitData, const std::vector<unsigned int>& edgesIn, std::vector<unsigned int>& leftEdgesOut, std::vector<unsigned int>& rightEdgesOut)
	{
		return 0.0f;
	}

	void Collision2DTree::BuildRecurse(Node* node, const pmp::BoundingBox2& box, const std::vector<unsigned int>& edgeIds, unsigned int remainingDepth)
	{
	}

} // namespace SDF