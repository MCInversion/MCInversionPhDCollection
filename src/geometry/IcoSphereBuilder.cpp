#include "IcoSphereBuilder.h"

#include <map>

namespace
{
	using VertexList = std::vector<pmp::vec3>;
	using TriangleList = std::vector<std::vector<unsigned int>>;

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
	
} // namespace

namespace Geometry
{
	using Lookup = std::map<std::pair<unsigned int, unsigned int>, unsigned int>;

	/**
	 * \brief Inserts a midpoint back-projected onto unit sphere and logs it into the lookup table.
	 * \param lookup        Lookup table for logging edge -> inserted vertex mapping.
	 * \param vertices      mesh vertices buffer.
	 * \param startPtId     index of edge start point.
	 * \param endPtId       index of edge end point.
	 * \return index of newly inserted back-projected midpoint.
	 */
	[[nodiscard]] unsigned int InsertMidpoint(Lookup& lookup, VertexList& vertices, 
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
	[[nodiscard]] TriangleList Subdivide(VertexList& vertices, const TriangleList& triangles)
	{
		Lookup lookup;
		TriangleList result;
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

	void IcoSphereBuilder::BuildBaseData()
	{
		VertexList vertices = ICOSAHEDRON_BASE_VERTICES;
		TriangleList triangles = ICOSAHEDRON_BASE_VERTEX_INDICES;

		for (unsigned int i = 0; i < m_SubdivisionLevel; i++) 
		{
			triangles = Subdivide(vertices, triangles);
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