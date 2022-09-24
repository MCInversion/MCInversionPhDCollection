#pragma once

#include "pmp/SurfaceMesh.h"

namespace Geometry
{
	/**
	 * \brief A simple data structure for mesh geometry containing only vertices, index triples, triangulation indices and, optionally, vertex normal coordinate triples.
	 * \struct BaseMeshGeometryData
	*/
	struct BaseMeshGeometryData
	{
		std::vector<pmp::vec3> Vertices{};
		std::vector<std::vector<unsigned int>> PolyIndices{};
		std::vector<pmp::vec3> VertexNormals{};
	};

	/**
	 * \brief Converts given BaseMeshGeometryData to pmp::SurfaceMesh.
	 * \param geomData     input base mesh geometry data.
	 * \return pmp::SurfaceMesh result.
	 */
	[[nodiscard]] pmp::SurfaceMesh ConvertBufferGeomToPMPSurfaceMesh(const BaseMeshGeometryData& geomData);
	
} // namespace Geometry