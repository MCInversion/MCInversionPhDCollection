#include "IcoSphereBuilder.h"

#include <map>

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

namespace Geometry
{
	using Lookup = std::map<std::pair<unsigned int, unsigned int>, unsigned int>;
	using LookupMulti = std::map<std::pair<unsigned int, unsigned int>, std::vector<unsigned int>>;

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
	// ===============================================================================
	//

	std::vector<unsigned int> CalculateEdgePoints(
		LookupMulti& lookup,
		IcoSphere::VertexList& vertices,
		unsigned int startPtId,
		unsigned int endPtId,
		const size_t& nInteriorEdgePts)
	{
		std::vector<unsigned int> ids;
		LookupMulti::key_type key(startPtId, endPtId);

		if (key.first > key.second)
		{
			std::swap(key.first, key.second);
		}

		const auto it = lookup.find(key);
		if (it != lookup.end())
		{
			// Already calculated for this edge, just return the precomputed IDs
			// but reverse the order, because of opposite half-edge orientation
			ids = it->second;
			std::reverse(ids.begin(), ids.end());
			return ids;
		}

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

		return ids;
	}

	std::vector<std::vector<unsigned int>> CalculateInteriorPoints(
		IcoSphere::VertexList& vertices,
		const std::vector<std::vector<unsigned int>>& edgePoints,
		const size_t& nInteriorPts)
	{
		if (nInteriorPts == 0)
		{
			return {};
		}
		const auto nMaxPtsInRow = static_cast<size_t>(sqrt(1 + 8 * nInteriorPts) - 1) / 2;
		std::vector<std::vector<unsigned int>> interiorPoints;
		interiorPoints.resize(nMaxPtsInRow); // Initialize each row

		size_t pointCount = 0;

		for (size_t i = 0; i < nMaxPtsInRow; ++i)
		{
			interiorPoints[i].resize(nMaxPtsInRow - i); // Initialize each column for each row
			for (size_t j = 0; j < nMaxPtsInRow - i; ++j)
			{
				const pmp::vec3 p1 = vertices[edgePoints[0][i]];
				const pmp::vec3 p2 = vertices[edgePoints[1][j]];
				const pmp::vec3 p3 = vertices[edgePoints[2][nMaxPtsInRow - i - j]];

				const pmp::vec3 newPoint = normalize((p1 + p2 + p3) / 3.0f);

				interiorPoints[i][j] = static_cast<unsigned int>(vertices.size());
				vertices.push_back(newPoint);

				++pointCount;

				if (pointCount >= nInteriorPts)
					return interiorPoints;
			}
		}

		return interiorPoints;
	}

	IcoSphere::TriangleList GenerateTriangles(
		const std::vector<std::vector<unsigned int>>& edgePoints,
		const std::vector<std::vector<unsigned int>>& interiorPoints,
		const size_t& nTrisPerEdge)
	{
		if (edgePoints.size() < 3)
		{
			throw std::logic_error("GenerateTriangles: No points to generate triangles from!\n");
		}

		if (edgePoints[0].size() < 3 || edgePoints[1].size() < 3 || edgePoints[2].size() < 3)
		{
			throw std::logic_error("GenerateTriangles: Incorrect number of points to generate triangles from!\n");
		}

		IcoSphere::TriangleList newTriangles;

		// Special case: if there are no interior points (s = 1), generate four triangles
		if (interiorPoints.empty())
		{
			const IcoSphere::Triangle tri1 = { edgePoints[0][0], edgePoints[0][1], edgePoints[2][1] };
			const IcoSphere::Triangle tri2 = { edgePoints[0][1], edgePoints[0][2], edgePoints[1][1] };
			const IcoSphere::Triangle tri3 = { edgePoints[2][1], edgePoints[1][1], edgePoints[1][2] };
			const IcoSphere::Triangle tri4 = { edgePoints[0][1], edgePoints[1][1], edgePoints[2][1] };

			newTriangles.push_back(tri1);
			newTriangles.push_back(tri2);
			newTriangles.push_back(tri3);
			newTriangles.push_back(tri4);

			return newTriangles;
		}

		// First row of triangles attached to the first edge
		for (size_t i = 0; i < nTrisPerEdge; ++i)
		{
			IcoSphere::Triangle tri = { edgePoints[0][i], edgePoints[0][i + 1], edgePoints[1][nTrisPerEdge - 1 - i] };
			newTriangles.push_back(tri);
		}

		// Generate triangles for the interior rows
		size_t rowStartIdx = 0; // Keeps track of where each row starts in the edgePoints
		for (size_t i = 1; i <= nTrisPerEdge; ++i)
		{
			for (size_t j = 0; j < nTrisPerEdge - i; ++j)
			{
				// Determine the vertices for this set of triangles
				const unsigned int v1 = edgePoints[0][rowStartIdx + j];
				const unsigned int v2 = edgePoints[0][rowStartIdx + j + 1];
				const unsigned int v3 = edgePoints[1][rowStartIdx + j + 1];
				const unsigned int v4 = edgePoints[1][rowStartIdx + j];
				const unsigned int v5 = interiorPoints[i - 1][j];

				// Generate the four triangles using these vertices
				IcoSphere::Triangle tri1 = { v1, v2, v5 };
				IcoSphere::Triangle tri2 = { v2, v3, v5 };
				IcoSphere::Triangle tri3 = { v3, v4, v5 };
				IcoSphere::Triangle tri4 = { v4, v1, v5 };

				// Add the new triangles to the list
				newTriangles.push_back(tri1);
				newTriangles.push_back(tri2);
				newTriangles.push_back(tri3);
				newTriangles.push_back(tri4);
			}
			// Update the row start index for edgePoints
			rowStartIdx += nTrisPerEdge - i + 1;
		}

		return newTriangles;
	}

	IcoSphere::TriangleList SubdivideSingleTriangle(LookupMulti& lookup, IcoSphere::VertexList& vertices, const IcoSphere::Triangle& triangle, const unsigned int& subdivLevel)
	{
		const auto nInteriorEdgePts = (static_cast<size_t>(pow(2, subdivLevel)) - 1);
		std::vector<std::vector<unsigned int>> edgePoints(3);
		for (size_t e = 0; e < 3; ++e)
		{
			edgePoints[e] = CalculateEdgePoints(lookup, vertices, triangle[e], triangle[(e + 1) % 3], nInteriorEdgePts);
		}

		const auto nInteriorTrianglePts = static_cast<size_t>(pow(4, subdivLevel) - 3 * pow(2, subdivLevel) + 2) / 2;
		const auto interiorPoints = CalculateInteriorPoints(vertices, edgePoints, nInteriorTrianglePts);

		IcoSphere::TriangleList newTriangles = GenerateTriangles(edgePoints, interiorPoints, nInteriorEdgePts + 1);

		return newTriangles;
	}

	[[nodiscard]] std::pair<size_t, size_t> IcoSpherePreallocationCapacities(const size_t& subdivLvl, const size_t& nFaces0, const size_t& nVerts0)
	{
		const size_t nEdges0 = nVerts0 + nFaces0 - 2; // from Euler characteristic
		const size_t nVerts = ((nEdges0 * static_cast<size_t>(pow(4, subdivLvl) - 1) + 3 * nVerts0) / 3);
		const size_t nFaces = (static_cast<size_t>(pow(4, subdivLvl)) * nFaces0);

		return { nVerts, nFaces };
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
				// new "predictive" strategy:
				// generates new points uniformly across each triangle's spherical projection
				IcoSphere::TriangleList newTriangles{};
				IcoSphere::VertexList newVertices{};
				// reserve memory
				const auto [vertexCapacity, faceCapacity] = IcoSpherePreallocationCapacities(m_SubdivisionLevel, triangles.size(), vertices.size());
				newVertices.reserve(vertexCapacity);
				newTriangles.reserve(faceCapacity);

				LookupMulti lookup;
				newVertices = vertices;

				for (const auto& triangle : triangles)
				{
					// Subdivide this triangle
					IcoSphere::TriangleList subdividedTriangles = SubdivideSingleTriangle(lookup, newVertices, triangle, m_SubdivisionLevel);

					// Append these to the list of newTriangles
					newTriangles.insert(newTriangles.end(), subdividedTriangles.begin(), subdividedTriangles.end());
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