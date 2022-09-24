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

		// TODO: normals

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
} // namespace Geometry