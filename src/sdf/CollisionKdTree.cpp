#include "CollisionKdTree.h"

#include "geometry/GeometryUtil.h"

#include <numeric>
#include <stack>

namespace SDF
{
	constexpr float BOX_INFLATION = 1e-6f;
	constexpr unsigned int MAX_DEPTH = 20;

	/**
	 * \brief Converts polygons from pmp::SurfaceMesh to triangles.
	 * \param mesh     Input surface mesh.
	 * \return a triangulated pmp::SurfaceMesh.
	 */
	pmp::SurfaceMesh GetTriangulatedSurfaceMesh(const pmp::SurfaceMesh& mesh)
	{
		pmp::SurfaceMesh result(mesh);

		// TODO: Use Poly2Tri

		for (const auto f : result.faces())
		{
			if (result.valence(f) == 3)
				continue;

			const auto vBegin = *result.vertices(f).begin();
			result.split(f, vBegin);
		}

		return result;
	}

	CollisionKdTree::CollisionKdTree(const pmp::SurfaceMesh& mesh, const SplitFunction& spltFunc)
	{
		auto bbox = mesh.bounds();
		bbox.expand(BOX_INFLATION, BOX_INFLATION, BOX_INFLATION);
		const auto triMesh = mesh.is_triangle_mesh() ? mesh : GetTriangulatedSurfaceMesh(mesh);

		// extract vertex positions
		const auto vpointPropData = mesh.get_vertex_property<pmp::Point>("v:point").data();
		const size_t nVerts = mesh.n_vertices();
		m_VertexPositions.reserve(nVerts);
		for (size_t i = 0; i < nVerts; i++)
			m_VertexPositions.emplace_back(vpointPropData[i]);

		// extract triangle ids
		const size_t nTriangles = mesh.n_faces();
		m_Triangles.reserve(nTriangles);
		for (size_t i = 0; i < nTriangles; i++)
		{
			const pmp::Face f(i);
			std::vector<unsigned int> vIds;
			vIds.reserve(3);
			for (const auto v : mesh.vertices(f))
				vIds.emplace_back(v.idx());
			m_Triangles.emplace_back(Triangle{ vIds[0], vIds[1], vIds[2] });
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
	 */
	[[nodiscard]] SplitAxisPreference GetSplitAxisPreference(const pmp::BoundingBox& box)
	{
		if ((box.max()[0] - box.min()[0]) > (box.max()[1] - box.min()[1]) &&
			(box.max()[0] - box.min()[0]) > (box.max()[2] - box.min()[2]))
		{
			return SplitAxisPreference::XAxis;
		}

		if ((box.max()[1] - box.min()[1]) > (box.max()[2] - box.min()[2]))
		{
			return SplitAxisPreference::YAxis;
		}

		return SplitAxisPreference::ZAxis;
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
		const unsigned int axisId = (splitData.axis == SplitAxisPreference::XAxis ? 0 : (splitData.axis == SplitAxisPreference::YAxis ? 1 : 2));

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

	// Fast kd-tree Construction with an Adaptive Error-Bounded Heuristic (Hunt, Mark, Stoll)
	//
	float AdaptiveSplitFunction(const BoxSplitData& splitData, const std::vector<unsigned int>& facesIn, std::vector<unsigned int>& leftFacesOut, std::vector<unsigned int>& rightFacesOut)
	{
		const unsigned int nFaces = facesIn.size();
		const auto& box = *splitData.box;
		const unsigned int axisId = (splitData.axis == SplitAxisPreference::XAxis ? 0 : (splitData.axis == SplitAxisPreference::YAxis ? 1 : 2));
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
		for (i = 0; i <= BOX_CUTS + 1; i++) bCutPos[i] = a * (1.0 - ((float)i / (float)(BOX_CUTS + 1))) + b * ((float)i / (float)(BOX_CUTS + 1));

		// set C_L(x) = 0, C_R(x) = 0
		auto C_L = std::vector<unsigned int>(BOX_CUTS + 2);
		auto C_R = std::vector<unsigned int>(BOX_CUTS + 2);

		auto faceMins = std::vector<float>(nFaces);
		auto faceMaxes = std::vector<float>(nFaces);

		for (const auto& fId : facesIn)
		{
			faceMins[i] = TriangleMin(triVertexIds[fId], vertices, axisId);
			faceMaxes[i] = TriangleMax(triVertexIds[fId], vertices, axisId);

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
			ran_s = (float)(i + 1) / (float)(BOX_CUTS + 1) * nFaces;

			for (j = 0; j <= BOX_CUTS; j++) 
			{
				S_L[j] += (ran_s > C_L[j] && ran_s < C_L[j + 1] ? 1 : 0);
				S_R[j] += (ran_s > C_R[j] && ran_s < C_R[j + 1] ? 1 : 0);
			}
		}
		std::reverse(S_R.begin(), S_R.end());

		// ==== Stage 3: add more sampling positions to subdivided segments ===========================================

		auto all_splt_L = std::vector<float>(2 * BOX_CUTS);
		auto all_splt_R = std::vector<float>(2 * BOX_CUTS);
		float segLen = (float)(b - a) / (float)(BOX_CUTS + 1);
		unsigned int nSeg_L = 0, nSeg_R = 0;

		for (i = 0; i <= BOX_CUTS; i++) {
			if (i > 0) {
				all_splt_L[nSeg_L++] = bCutPos[i];
				all_splt_R[nSeg_R++] = bCutPos[i];
			}
			for (j = 0; j < S_L[i]; j++) {
				all_splt_L[nSeg_L++] = bCutPos[i] + (j + 1) / (S_L[i] + 1) * segLen;
			}
			for (j = 0; j < S_R[i]; j++) {
				all_splt_R[nSeg_R++] = bCutPos[i] + (j + 1) / (S_R[i] + 1) * segLen;
			}
		}

		// Compute surface area heuristic SAH:
		// remaining two dimensions of the child box candidates
		// TODO: verify with SurfaceEvolver
		const float boxDim0 = box.max()[(axisId + 1) % 3] - box.min()[(axisId + 1) % 3];
		const float boxDim1 = box.max()[(axisId + 2) % 3] - box.min()[(axisId + 2) % 3];

		auto SA_L = std::vector<float>(2 * BOX_CUTS);
		auto SA_R = std::vector<float>(2 * BOX_CUTS);

		// SA_L(x) = (boxDim_L(x) + boxDim0 + boxDim1) * 2.0 / tot_BoxArea
		// SA_R(x) = (boxDim_R(x) + boxDim0 + boxDim1) * 2.0 / tot_BoxArea

		for (i = 0; i < 2 * BOX_CUTS; i++) 
		{
			SA_L[i] = ((all_splt_L[i] - a) + boxDim0 + boxDim1) * 2.0 / tot_BoxArea;
			SA_R[i] = ((b - all_splt_R[i]) + boxDim0 + boxDim1) * 2.0 / tot_BoxArea;
		}

		// ==== Stage 4: RESAMPLE C_L and C_R on all sample points & construct an approximation of cost(x) to minimize
		float min, max;
		auto cost = std::vector<float>(2 * BOX_CUTS);

		// cost(x) = C_L(x) * SA_L(x) + C_R(x) * SA_R(x):
		for (const auto& fId : facesIn)
		{
			min = faceMins[fId];
			max = faceMaxes[fId];

			for (j = 0; j < 2 * BOX_CUTS; j++) {
				cost[j] += (min < all_splt_L[j] ? 1 : 0) * SA_L[j] + (max > all_splt_R[j] ? 1 : 0) * SA_R[j];
			}
		}

		// ==== Stage 5: Minimize cost(x) & classify primitives  =====================================================

		float bestSplit = 0.5f * (a + b); // if this loop fails to initialize bestSplit, set it to middle
		float minCost = FLT_MAX;

		for (i = 1; i < 2 * BOX_CUTS; i++) {
			if (cost[i] < minCost) {
				minCost = cost[i];
				bestSplit = all_splt_L[i];
			}
		}

		// fill left and right arrays now that best split position is known:
		leftFacesOut.reserve(nFaces);
		rightFacesOut.reserve(nFaces);
		for (const auto& fId : facesIn)
		{
			min = faceMins[fId];
			max = faceMaxes[fId];

			if (min <= bestSplit) 
			{
				leftFacesOut.push_back(fId);
			}
			if (max >= bestSplit) 
			{
				rightFacesOut.push_back(fId);
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
	 * \brief Triangle count verification after split.
	 * \param nLeftTris    number of triangles in the left child node candidate.
	 * \param nRightTris   number of triangles in the right child node candidate.
	 * \param nTris        number of triangles in the current node.
	 * \return true if the imbalance is too high.
	 */
	[[nodiscard]] bool ShouldStopBranching(size_t nLeftTris, size_t nRightTris, size_t nTris)
	{
		return static_cast<double>(nLeftTris + nRightTris) >= 1.5 * static_cast<double>(nTris);
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
		result.reserve(triangles.size());
		for (const auto& tri : triangleIds)
		{
			const std::vector triVertices{
				vertexPositions[triangles[tri].v0Id],
				vertexPositions[triangles[tri].v1Id],
				vertexPositions[triangles[tri].v2Id]
			};

			if (!Geometry::TriangleIntersectsBox(triVertices, boxCenter, boxHalfSize))
				continue;

			result.emplace_back(tri);
		}

		result.shrink_to_fit();
		return result;
	}

	[[nodiscard]] pmp::BoundingBox GetChildBox(
		const pmp::BoundingBox& parentBox, const float& splitPos, const SplitAxisPreference& axisPreference, const bool& isLeft)
	{
		const unsigned int axisId = (axisPreference == SplitAxisPreference::XAxis ? 0 : (axisPreference == SplitAxisPreference::YAxis ? 1 : 2));
		pmp::BoundingBox childBox(parentBox);
		if (isLeft)
			childBox.max()[axisId] = splitPos;
		else
			childBox.min()[axisId] = splitPos;

		childBox.expand(BOX_INFLATION, BOX_INFLATION, BOX_INFLATION);

		return childBox;
	}

	//! minimum number of triangles for node
	constexpr size_t MIN_NODE_TRIANGLE_COUNT = 2;

	void CollisionKdTree::BuildRecurse(Node* node, const pmp::BoundingBox& box, const std::vector<unsigned int>& triangleIds, unsigned int remainingDepth)
	{
		assert(node != nullptr);

		if (remainingDepth == 0 || triangleIds.size() <= MIN_NODE_TRIANGLE_COUNT)
		{
			node->triangleIds = FilterTriangles(box, triangleIds, m_Triangles, m_VertexPositions);
			return;
		}

		const auto axisPreference = GetSplitAxisPreference(box);

		std::vector<unsigned int> leftTriangleIds{};
		std::vector<unsigned int> rightTriangleIds{};
		const float splitPosition = m_FindSplitAndClassify({ this, &box, axisPreference }, triangleIds, leftTriangleIds, rightTriangleIds);

		if (ShouldStopBranching(leftTriangleIds.size(), rightTriangleIds.size(), triangleIds.size()))
		{
			node->triangleIds = FilterTriangles(box, triangleIds, m_Triangles, m_VertexPositions);
			return;
		}

		if (!leftTriangleIds.empty())
		{
			const auto leftBox = GetChildBox(box, splitPosition, axisPreference, true);
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
			const auto rightBox = GetChildBox(box, splitPosition, axisPreference, false);
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

	/*void CollisionKdTree::GetTrianglesInABox(const pmp::BoundingBox& box, std::vector<unsigned int>& foundTriangleIds) const
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
	}*/

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
				if (box.Intersects(currentNode->left_child->box))
					nodeStack.push(currentNode->left_child);
			}

			if (currentNode->right_child)
			{
				if (box.Intersects(currentNode->right_child->box))
					nodeStack.push(currentNode->right_child);
			}
		}

		foundTriangleIds.shrink_to_fit();
	}

} // namespace SDF