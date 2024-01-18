#include "IcoSphereBuilder.h"

#include <map>
#include <unordered_map>

namespace IcoSphere
{
	using VertexList = std::vector<pmp::vec3>;
	using Triangle = std::vector<unsigned int>;
	using TriangleList = std::vector<Triangle>;

	/// \brief golden ratio.
	const float phi = (1.0f + sqrt(5.0f)) / 2.0f;

	/// \brief golden ratio ico-sphere radius.
	const float norm = sqrt(1.0f + phi * phi);

	/// \brief vertices of an icosahedron.
	inline const VertexList ICOSAHEDRON_BASE_VERTICES{
		/* v0 */ pmp::vec3{-1.0f / norm, phi / norm, 0.0f},  /* v1 */ pmp::vec3{1.0f / norm, phi / norm, 0.0f},
		/* v2 */ pmp::vec3{-1.0f / norm, -phi / norm, 0.0f}, /* v3 */ pmp::vec3{1.0f / norm, -phi / norm, 0.0f},
		/* v4 */ pmp::vec3{0.0f, -1.0f / norm, phi / norm},  /* v5 */ pmp::vec3{0.0f, 1.0f / norm, phi / norm},
		/* v6 */ pmp::vec3{0.0f, -1.0f / norm, -phi / norm}, /* v7 */ pmp::vec3{0.0f, 1.0f / norm, -phi / norm},
		/* v8 */ pmp::vec3{phi / norm, 0.0f, -1.0f / norm},  /* v9 */ pmp::vec3{phi / norm, 0.0f, 1.0f / norm},
		/* v10 */pmp::vec3{-phi / norm, 0.0f, -1.0f / norm},/* v11 */ pmp::vec3{-phi / norm, 0.0f, 1.0f / norm}
	};

	/// \brief vertex index triples for icosahedron triangle faces.
	inline const TriangleList ICOSAHEDRON_BASE_VERTEX_INDICES{
		{0, 11, 5},    {0, 5, 1},    {0, 1, 7},    {0, 7, 10},    {0, 10, 11},
		{1, 5, 9},     {5, 11, 4},   {11, 10, 2},  {10, 7, 6},    {7, 1, 8},
		{3, 9, 4},     {3, 4, 2},    {3, 2, 6},    {3, 6, 8},     {3, 8, 9},
		{4, 9, 5},     {2, 4, 11},   {6, 2, 10},   {8, 6, 7},     {9, 8, 1}
	};
	
} // namespace IcoSphere

using Lookup = std::map<std::pair<unsigned int, unsigned int>, unsigned int>;
using LookupMulti = std::unordered_map<size_t, std::vector<unsigned int>>;

namespace 
{

	//
	// ===============================================================================
	//

	size_t createKey(unsigned int a, unsigned int b)
	{
		static_assert(sizeof(size_t) >= 2 * sizeof(unsigned int),
			"size_t must be at least twice the size of unsigned int to safely store two unsigned int values");
		size_t key = static_cast<size_t>(a);
		key <<= 32;
		key |= b;
		return key;
	}

	//std::pair<unsigned int, unsigned int> extractKey(size_t key)
	//{
	//	return { static_cast<unsigned int>(key >> 32), static_cast<unsigned int>(key) };
	//}

	void CalculateEdgePoints(
		LookupMulti& lookup,
		IcoSphere::VertexList& vertices,
		const unsigned int& startPtId,
		const unsigned int& endPtId,
		const size_t& nInteriorEdgePts,
		std::vector<unsigned int>& ids)
	{
		unsigned int keyVal0 = startPtId;
		unsigned int keyVal1 = endPtId;
		const LookupMulti::key_type key = createKey(keyVal0, keyVal1);

		if (keyVal0 > keyVal1)
		{
			std::swap(keyVal0, keyVal1);
		}

		const auto it = lookup.find(key);
		if (it != lookup.end())
		{
			// Already calculated for this edge, just return the precomputed IDs
			// but reverse the order, because of opposite half-edge orientation
			ids = it->second;
			std::reverse(ids.begin(), ids.end());
			return;
		}

		ids.reserve(nInteriorEdgePts + 2);
		ids.push_back(startPtId);
		const pmp::vec3 start = vertices[startPtId];
		const pmp::vec3 end = vertices[endPtId];

		const auto nTotalEdgeSegments = nInteriorEdgePts + 1;
		for (size_t i = 1; i <= nInteriorEdgePts; ++i)
		{
			const auto t = static_cast<float>(i) / static_cast<float>(nTotalEdgeSegments);
			pmp::vec3 newPoint = normalize(start + t * (end - start));

			ids.push_back(static_cast<unsigned int>(vertices.size()));
			vertices.push_back(newPoint);
		}
		ids.push_back(endPtId);

		// Save the calculated ids in the lookup map
		lookup[key] = ids;
	}

	void CalculateInteriorPoints(
		IcoSphere::VertexList& vertices,
		const std::vector<std::vector<unsigned int>>& edgePoints,
		const size_t& nInteriorPts,
		std::vector<std::vector<unsigned int>>& interiorPoints)
	{
		if (nInteriorPts == 0)
		{
			return;
		}

		if (edgePoints.size() < 3)
		{
			throw std::logic_error("CalculateInteriorPoints: edgePoints.size() < 3!\n");
		}

		const auto nMaxPtsInRow = static_cast<size_t>(sqrt(1 + 8 * nInteriorPts) - 1) / 2;
		interiorPoints.resize(nMaxPtsInRow); // Initialize each row

		// interpolate between the triangle pts with parameter values starting at p3
		const pmp::vec3& p1 = vertices[edgePoints[0][0]];
		const pmp::vec3& p2 = vertices[edgePoints[1][0]];
		const pmp::vec3& p3 = vertices[edgePoints[2][0]];

		const auto nEdgePts = edgePoints[0].size();

		for (size_t i = 1; i < nEdgePts - 1; i++)
		{
			const float iParam = static_cast<float>(i) / static_cast<float>(nEdgePts - 1);
			interiorPoints.reserve(nMaxPtsInRow - i + 1);
			for (size_t j = 1; j < nEdgePts - i - 1; j++)
			{
				const float jParam = static_cast<float>(j) / static_cast<float>(nEdgePts - 1);
				const pmp::vec3 newPoint = normalize(p1 * (1.0f - iParam - jParam) + p2 * jParam + p3 * iParam);

				interiorPoints[i - 1].push_back(static_cast<unsigned int>(vertices.size()));
				vertices.push_back(newPoint);
			}
		}
	}

	void GenerateTriangles(
		const std::vector<std::vector<unsigned int>>& edgePoints,
		const std::vector<std::vector<unsigned int>>& interiorPoints,
		const size_t& subdivLevel,
		IcoSphere::TriangleList& newTriangles)
	{
		if (edgePoints.size() < 3)
		{
			throw std::logic_error("GenerateTriangles: No points to generate triangles from!\n");
		}

		if (edgePoints[0].size() < 3 || edgePoints[1].size() < 3 || edgePoints[2].size() < 3)
		{
			throw std::logic_error("GenerateTriangles: Incorrect number of points to generate triangles from!\n");
		}

		//newTriangles.reserve(static_cast<size_t>(pow(4, subdivLevel))); // already resesrved

		// Special case: if there are no interior points (s = 1), generate four triangles
		if (interiorPoints.empty())
		{
			newTriangles.push_back({ edgePoints[0][0], edgePoints[0][1], edgePoints[2][1] });
			newTriangles.push_back({ edgePoints[0][1], edgePoints[0][2], edgePoints[1][1] });
			newTriangles.push_back({ edgePoints[2][1], edgePoints[1][1], edgePoints[1][2] });
			newTriangles.push_back({ edgePoints[0][1], edgePoints[1][1], edgePoints[2][1] });

			return;
		}

		const size_t nIntPtsRowSize = interiorPoints.size(); // should be the same as interiorPoints[0].size()
		// because the edge blocks of triangles need to complement each other to avoid overlap,
		// we stop 2 triangles earlier than we normally would to make space for the first two triangles of the
		// consecutive edge
		const size_t nTrisPerEdge = static_cast<size_t>(pow(2, subdivLevel));
		const size_t nCutOffTrisPerEdge = nTrisPerEdge - 1;

		// triangulate along edge 0
		for (size_t i = 0; i < nCutOffTrisPerEdge; i++)
		{
			// triangle pointing from edge 0
			{
				const auto v0Id = edgePoints[0][i];
				const auto v1Id = edgePoints[0][i + 1];
				const auto v2Id = (i > 0 ? interiorPoints[0][i - 1] : edgePoints[2][nTrisPerEdge - 1]);

				newTriangles.push_back({ v0Id, v1Id, v2Id });
			}

			// triangle pointing towards edge 0
			if (i < nCutOffTrisPerEdge - 1)
			{
				const auto v0Id = edgePoints[0][i + 1];
				const auto v1Id = interiorPoints[0][i];
				const auto v2Id = (i > 0 ? interiorPoints[0][i - 1] : edgePoints[2][nTrisPerEdge - 1]);

				newTriangles.push_back({ v0Id, v1Id, v2Id });
			}
		}

		// triangulate along edge 1
		for (size_t i = 0; i < nCutOffTrisPerEdge; i++)
		{
			// iterating diagonally across interiorPoints

			// triangle pointing from edge 1
			{
				const auto v0Id = edgePoints[1][i];
				const auto v1Id = edgePoints[1][i + 1];
				const auto v2Id = (i > 0 ? interiorPoints[i - 1][nIntPtsRowSize - i] : edgePoints[0][nTrisPerEdge - 1]);

				newTriangles.push_back({ v0Id, v1Id, v2Id });
			}

			// triangle pointing towards edge 1
			if (i < nCutOffTrisPerEdge - 1)
			{
				const auto v0Id = edgePoints[1][i + 1];
				const auto v1Id = interiorPoints[i][nIntPtsRowSize - i - 1];
				const auto v2Id = (i > 0 ? interiorPoints[i - 1][nIntPtsRowSize - i] : edgePoints[0][nTrisPerEdge - 1]);

				newTriangles.push_back({ v0Id, v1Id, v2Id });
			}
		}

		// triangulate along edge 2
		for (size_t i = 0; i < nCutOffTrisPerEdge; i++)
		{
			// triangle pointing from edge 2
			{
				const auto v0Id = edgePoints[2][i];
				const auto v1Id = edgePoints[2][i + 1];
				const auto v2Id = (i > 0 ? interiorPoints[nIntPtsRowSize - i][0] : edgePoints[1][nTrisPerEdge - 1]);

				newTriangles.push_back({ v0Id, v1Id, v2Id });
			}

			// triangle pointing towards edge 2
			if (i < nCutOffTrisPerEdge - 1)
			{
				const auto v0Id = edgePoints[2][i + 1];
				const auto v1Id = interiorPoints[nIntPtsRowSize - i - 1][0];
				const auto v2Id = (i > 0 ? interiorPoints[nIntPtsRowSize - i][0] : edgePoints[1][nTrisPerEdge - 1]);

				newTriangles.push_back({ v0Id, v1Id, v2Id });
			}
		}

		// triangulate interior points
		for (size_t i = 0; i < nIntPtsRowSize - 1; i++)
		{
			for (size_t j = 0; j < nIntPtsRowSize - i - 1; j++)
			{
				// upward pointing triangles
				{
					const auto v0Id = interiorPoints[i][j];
					const auto v1Id = interiorPoints[i][j + 1];
					const auto v2Id = interiorPoints[i + 1][j];

					newTriangles.push_back({ v0Id, v1Id, v2Id });
				}

				// downward pointing triangles
				if (j < nIntPtsRowSize - i - 2)
				{
					const auto v0Id = interiorPoints[i][j + 1];
					const auto v1Id = interiorPoints[i + 1][j + 1];
					const auto v2Id = interiorPoints[i + 1][j];

					newTriangles.push_back({ v0Id, v1Id, v2Id });
				}
			}
		}
	}

	void SubdivideSingleTriangle(LookupMulti& lookup, IcoSphere::VertexList& vertices, const IcoSphere::Triangle& triangle, const unsigned int& subdivLevel, IcoSphere::TriangleList& newTriangles)
	{
		const auto nInteriorEdgePts = (static_cast<size_t>(pow(2, subdivLevel)) - 1);
		std::vector<std::vector<unsigned int>> edgePoints(3);
		for (size_t e = 0; e < 3; ++e)
		{
			CalculateEdgePoints(lookup, vertices, triangle[e], triangle[(e + 1) % 3], nInteriorEdgePts, edgePoints[e]);
		}

		const auto nInteriorTrianglePts = static_cast<size_t>(pow(4, subdivLevel) - 3 * pow(2, subdivLevel) + 2) / 2;
		std::vector<std::vector<unsigned int>> interiorPoints;
		CalculateInteriorPoints(vertices, edgePoints, nInteriorTrianglePts, interiorPoints);

		GenerateTriangles(edgePoints, interiorPoints, subdivLevel, newTriangles);
	}

	[[nodiscard]] std::pair<size_t, size_t> IcoSpherePreallocationCapacities(const size_t& subdivLvl, const size_t& nFaces0, const size_t& nVerts0)
	{
		const size_t nEdges0 = nVerts0 + nFaces0 - 2; // from Euler characteristic
		const size_t nVerts = ((nEdges0 * static_cast<size_t>(pow(4, subdivLvl) - 1) + 3 * nVerts0) / 3);
		const size_t nFaces = (static_cast<size_t>(pow(4, subdivLvl)) * nFaces0);

		return { nVerts, nFaces };
	}
	
} // anonymous namespace

namespace Geometry
{
	/**
	 * \brief Inserts a midpoint back-projected onto unit sphere and logs it into the lookup table.
	 * \param lookup        Lookup table for logging edge -> inserted vertex mapping.
	 * \param vertices      mesh vertices buffer.
	 * \param startPtId     index of edge start point.
	 * \param endPtId       index of edge end point.
	 * \return index of newly inserted back-projected midpoint.
	 */
	[[nodiscard]] unsigned int InsertMidpoint(Lookup& lookup, IcoSphere::VertexList& vertices,
		const unsigned int& startPtId, const unsigned int& endPtId)
	{
		Lookup::key_type key(startPtId, endPtId);

		if (key.first > key.second) 
		{
			std::swap(key.first, key.second);
		}

		const auto& [insertedIter, insertionSuccessful] = lookup.insert({ key, vertices.size() });
		if (insertionSuccessful)
		{
			const auto& edge0 = vertices[startPtId];
			const auto& edge1 = vertices[endPtId];
			const auto midPt = normalize(edge0 + edge1);
			vertices.push_back(midPt);
		}

		return insertedIter->second;
	};

	/**
	 * \brief Subdivides a triangle mesh representation using 4:1 subdivision.
	 * \param vertices      modifiable list of mesh vertices.
	 * \param triangles     triangle vertex indices before subdivision.
	 * \return triangle vertex indices after subdivision.
	 */
	[[nodiscard]] IcoSphere::TriangleList Subdivide(IcoSphere::VertexList& vertices, const IcoSphere::TriangleList& triangles)
	{
		Lookup lookup;
		IcoSphere::TriangleList result;
		result.reserve(4 * triangles.size());

		for (const auto& triVertId : triangles)
		{
			std::array midPtIds{UINT_MAX, UINT_MAX, UINT_MAX};
			for (unsigned int edgeId = 0; edgeId < 3; edgeId++)
			{
				midPtIds[edgeId] = InsertMidpoint(lookup, vertices, triVertId[edgeId], triVertId[(edgeId + 1) % 3]);
			}

			result.push_back({ triVertId[0], midPtIds[0], midPtIds[2] });
			result.push_back({ triVertId[1], midPtIds[1], midPtIds[0] });
			result.push_back({ triVertId[2], midPtIds[2], midPtIds[1] });
			result.push_back({ midPtIds[0], midPtIds[1], midPtIds[2] });
		}

		return result;
	}

	//
	// ===========================================================
	//

	void IcoSphereBuilder::BuildBaseData()
	{
		IcoSphere::VertexList vertices = IcoSphere::ICOSAHEDRON_BASE_VERTICES;
		IcoSphere::TriangleList triangles = IcoSphere::ICOSAHEDRON_BASE_VERTEX_INDICES;

		if (m_SubdivisionLevel > 0)
		{
			// base icosahedron triangles need to be subdivided.
			if (m_UseRecursiveStrategy)
			{
				// old "recursive" strategy
				for (unsigned int i = 0; i < m_SubdivisionLevel; i++) 
				{
					triangles = Subdivide(vertices, triangles);
				}			
			}
			else
			{
				std::cerr << "!!!--------------------------------------------------------------------------------------\n";
				std::cerr << "IcoSphereBuilder::BuildBaseData [WARNING]: Using an experimental predictive construction with std::unordered_map lookup. Throws an exception for subdiv = 7 and higher.\n";
				std::cerr << "!!!--------------------------------------------------------------------------------------\n";
				// new "predictive" strategy:
				// generates new points uniformly across each triangle's spherical projection
				IcoSphere::TriangleList newTriangles{};
				IcoSphere::VertexList newVertices{};
				// reserve memory
				const auto [vertexCapacity, faceCapacity] = IcoSpherePreallocationCapacities(m_SubdivisionLevel, triangles.size(), vertices.size());
				newVertices.reserve(vertexCapacity);
				newTriangles.reserve(faceCapacity);

				LookupMulti lookup;
				lookup.reserve(vertexCapacity + faceCapacity - 2);
				newVertices = vertices;

				for (const auto& triangle : triangles)
				{
					SubdivideSingleTriangle(lookup, newVertices, triangle, m_SubdivisionLevel, newTriangles);
				}

				vertices = newVertices;
				triangles = newTriangles;
			}			
		}

		m_BaseResult = std::make_unique<BaseMeshGeometryData>();
		m_BaseResult->PolyIndices = triangles;

		auto& resultVertices = m_BaseResult->Vertices;

		if (m_ComputeNormals)
		{
			auto& resultVertexNormals = m_BaseResult->VertexNormals;
			resultVertexNormals.reserve(vertices.size());
			for (const auto& vert : vertices)
			{
				resultVertexNormals.push_back(vert);
			}
		}

		resultVertices.reserve(vertices.size());
		for (const auto& vert : vertices)
		{
			resultVertices.push_back(vert * m_Radius);
		}
	}

} // namespace Geometry