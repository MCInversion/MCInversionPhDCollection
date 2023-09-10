#pragma once

#include "pmp/SurfaceMesh.h"
#include "MarchingCubes.h"

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

	/**
	 * \brief Converts given MC_Mesh to pmp::SurfaceMesh.
	 * \param mcMesh       Marching cubes mesh to be converted.
	 * \return pmp::SurfaceMesh result.
	 */
	[[nodiscard]] pmp::SurfaceMesh ConvertMCMeshToPMPSurfaceMesh(const MC_Mesh& mcMesh);

	/**
	 * \brief For testing out the BaseMeshGeometryData by exporting it to a Wavefront OBJ file.
	 * \param geomData       input geom data.
	 * \param absFileName    absolute file path for the created file.
	 * \return if true, the export was successful.
	 */
	[[nodiscard]] bool ExportBaseMeshGeometryDataToOBJ(const BaseMeshGeometryData& geomData, const std::string& absFileName);
	
} // namespace Geometry