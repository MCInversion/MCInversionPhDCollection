#pragma once

#include <optional>

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

	/**
	 * \brief For importing very large OBJ mesh files with option for parallel.
	 * \param absFileName                absolute file path for the opened file.
	 * \param importInParallel           if true, a parallel version of the importer will be used.
	 * \param chunkIdsVertexPropPtrOpt   an optional ptr to a vector of "chunk" ids (1 chunk = 1 thread).
	 * \return optional BaseMeshGeometryData.
	 */
	[[nodiscard]] std::optional<BaseMeshGeometryData> ImportOBJMeshGeometryData(const std::string& absFileName, const bool& importInParallel = false, std::optional<std::vector<float>*> chunkIdsVertexPropPtrOpt = std::nullopt);

	/**
	 * \brief For importing PLY point cloud files with option for parallel.
	 * \param absFileName        absolute file path for the opened file.
	 * \param importInParallel   if true, a parallel version of the importer will be used.
	 * \return optional vector of points (pmp::vec3).
	 */
	[[nodiscard]] std::optional<std::vector<pmp::vec3>> ImportPLYPointCloudData(const std::string& absFileName, const bool& importInParallel = false);

	/**
	 * \brief For importing PLY point cloud files.
	 * \param absFileName        absolute file path for the opened file.
	 * \return optional vector of points (pmp::vec3).
	 */
	[[nodiscard]] std::optional<std::vector<pmp::vec3>> ImportPLYPointCloudDataMainThread(const std::string& absFileName);

	/**
	 * \brief At this point this is a demo function for MeshLab's ball-pivoting.
	 * \param meshData       input geom data.
	 * \param nVerts         the number of vertices to be sampled from meshData.
	 * \param absFileName    absolute file path for the opened file.
	 * \return if true, the export was successful.
	 */
	[[nodiscard]] bool ExportSampledVerticesToPLY(const BaseMeshGeometryData& meshData, size_t nVerts, const std::string& absFileName);

	/**
	 * \brief A utility for exporting polylines as Wavefront OBJ file.
	 * \param polylines      vector of polylines to be exported.
	 * \param absFileName    absolute file path for the opened file.
	 */
	[[nodiscard]] bool ExportPolylinesToOBJ(const std::vector<std::vector<pmp::vec3>>& polylines, const std::string& absFileName);

	/**
	 * \brief Computes the convex hull of an input point cloud.
	 * \param points           input point cloud.
	 * \return optional resulting pmp::SurfaceMesh if the computation is successful.
	 */
	[[nodiscard]] std::optional<pmp::SurfaceMesh> ComputePMPConvexHullFromPoints(const std::vector<pmp::Point>& points);

	/**
	 * \brief Computes the convex hull of an input point cloud.
	 * \param points           input point cloud.
	 * \return optional resulting BaseMeshGeometryData if the computation is successful.
	 */
	[[nodiscard]] std::optional<BaseMeshGeometryData> ComputeConvexHullFromPoints(const std::vector<pmp::Point>& points);

} // namespace Geometry