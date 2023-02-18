#include "GeometryConversionUtils.h"

#include <set>

namespace Geometry
{
	pmp::SurfaceMesh ConvertBufferGeomToPMPSurfaceMesh(const BaseMeshGeometryData& geomData)
	{
		pmp::SurfaceMesh result;

		// count edges
		std::set<std::pair<unsigned int, unsigned int>> edgeIdsSet;
		for (const auto& indexTuple : geomData.PolyIndices)
		{
			for (unsigned int i = 0; i < indexTuple.size(); i++)
			{
				unsigned int vertId0 = indexTuple[i];
				unsigned int vertId1 = indexTuple[(static_cast<size_t>(i) + 1) % indexTuple.size()];

				if (vertId0 > vertId1) std::swap(vertId0, vertId1);

				edgeIdsSet.insert({ vertId0, vertId1 });
			}
		}

		result.reserve(geomData.Vertices.size(), edgeIdsSet.size(), geomData.PolyIndices.size());
		for (const auto& v : geomData.Vertices)
			result.add_vertex(pmp::Point(v[0], v[1], v[2]));

		if (!geomData.VertexNormals.empty())
		{
			auto vNormal = result.vertex_property<pmp::Normal>("v:normal");
			for (auto v : result.vertices())
				vNormal[v] = geomData.VertexNormals[v.idx()];
		}

        for (const auto& indexTuple : geomData.PolyIndices)
        {
			std::vector<pmp::Vertex> vertices;
			vertices.reserve(indexTuple.size());
			for (const auto& vId : indexTuple)
				vertices.emplace_back(pmp::Vertex(vId));

			result.add_face(vertices);	        
        }

		return result;
	}

	pmp::SurfaceMesh ConvertMCMeshToPMPSurfaceMesh(const MC_Mesh& mcMesh)
	{
		pmp::SurfaceMesh result;

		// count edges
		std::set<std::pair<unsigned int, unsigned int>> edgeIdsSet;
		for (unsigned int i = 0; i < mcMesh.faceCount * 3; i += 3)
		{
			for (unsigned int j = 0; j < 3; j++)
			{
				unsigned int vertId0 = mcMesh.faces[i + j];
				unsigned int vertId1 = mcMesh.faces[i + (j + 1) % 3];

				if (vertId0 > vertId1) std::swap(vertId0, vertId1);

				edgeIdsSet.insert({ vertId0, vertId1 });
			}
		}

		result.reserve(mcMesh.vertexCount, edgeIdsSet.size(), mcMesh.faceCount);
		// MC produces normals by default
		auto vNormal = result.vertex_property<pmp::Normal>("v:normal");

		for (unsigned int i = 0; i < mcMesh.vertexCount; i++)
		{
			result.add_vertex(
			   pmp::Point(
				   mcMesh.vertices[i][0],
				   mcMesh.vertices[i][1],
				   mcMesh.vertices[i][2]
			));
			vNormal[pmp::Vertex(i)] = pmp::Normal{
				mcMesh.normals[i][0],
				mcMesh.normals[i][1],
				mcMesh.normals[i][2]
			};
		}

		for (unsigned int i = 0; i < mcMesh.faceCount * 3; i += 3)
		{
			std::vector<pmp::Vertex> vertices;
			vertices.reserve(3);
			for (unsigned int j = 0; j < 3; j++)
			{
				vertices.emplace_back(pmp::Vertex(mcMesh.faces[i + j]));
			}

			result.add_face(vertices);
		}

		return result;
	}

} // namespace Geometry