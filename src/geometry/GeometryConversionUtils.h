#pragma once

#include "pmp/SurfaceMesh.h"
#include "MarchingCubes.h"
#include <optional>

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
	 * \brief Converts given pmp::SurfaceMesh to BaseMeshGeometryData.
	 * \param pmpMesh     input pmp::SurfaceMesh
	 * \return BaseMeshGeometryData result.
	 */
	[[nodiscard]] BaseMeshGeometryData ConvertPMPSurfaceMeshToBaseMeshGeometryData(const pmp::SurfaceMesh& pmpMesh);

	/**
	 * \brief Converts given MC_Mesh to pmp::SurfaceMesh.
	 * \param mcMesh       Marching cubes mesh to be converted.
	 * \return pmp::SurfaceMesh result.
	 */
	[[nodiscard]] pmp::SurfaceMesh ConvertMCMeshToPMPSurfaceMesh(const MarchingCubes::MC_Mesh& mcMesh);

	/**
	 * \brief For testing out the BaseMeshGeometryData by exporting it to a Wavefront OBJ file.
	 * \param geomData       input geom data.
	 * \param absFileName    absolute file path for the created file.
	 * \return if true, the export was successful.
	 */
	[[nodiscard]] bool ExportBaseMeshGeometryDataToOBJ(const BaseMeshGeometryData& geomData, const std::string& absFileName);

	/**
	 * \brief For testing out the BaseMeshGeometryData by exporting it to VTK polydata file.
	 * \param geomData       input geom data.
	 * \param absFileName    absolute file path for the created file.
	 * \return if true, the export was successful.
	 */
	[[nodiscard]] bool ExportBaseMeshGeometryDataToVTK(const BaseMeshGeometryData& geomData, const std::string& absFileName);

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
	 * \brief Randomly sample vertices from given mesh data.
	 * \param meshData       input geom data.
	 * \param nVerts         the number of vertices to be sampled from meshData.
	 * \param seed           seed value for random index generation.
	 * \return The resulting point cloud
	 */
	[[nodiscard]] std::vector<pmp::Point> SamplePointsFromMeshData(const BaseMeshGeometryData& meshData, size_t nVerts, const std::optional<unsigned int>& seed = std::nullopt);

	/**
	 * \brief Randomly sample vertices and export them as *.ply point cloud.
	 * \param meshData       input geom data.
	 * \param nVerts         the number of vertices to be sampled from meshData.
	 * \param absFileName    absolute file path for the opened file.
	 * \param seed           seed value for random index generation.
	 * \return if true, the export was successful.
	 */
	[[nodiscard]] bool ExportSampledVerticesToPLY(const BaseMeshGeometryData& meshData, size_t nVerts, const std::string& absFileName, const std::optional<unsigned int>& seed = std::nullopt);

	/**
	 * \brief Randomly sample vertices with normals and export them as *.vtk vector data.
	 * \param meshData       input geom data.
	 * \param nVerts         the number of vertex normals to be sampled from meshData.
	 * \param absFileName    absolute file path for the opened file.
	 * \param seed           seed value for random index generation.
	 * \return if true, the export was successful.
	 */
	[[nodiscard]] bool ExportSampledVerticesWithNormalsToVTK(const BaseMeshGeometryData& meshData, size_t nVerts, const std::string& absFileName, const std::optional<unsigned int>& seed = std::nullopt);

	/**
	 * \brief Randomly sample vertices with normals and export them as *.ply oriented point cloud data.
	 * \param meshData       input geom data.
	 * \param nVerts         the number of vertex normals to be sampled from meshData.
	 * \param absFileName    absolute file path for the opened file.
	 * \param seed           seed value for random index generation.
	 * \return if true, the export was successful.
	 */
	[[nodiscard]] bool ExportSampledVerticesWithNormalsToPLY(const BaseMeshGeometryData& meshData, size_t nVerts, const std::string& absFileName, const std::optional<unsigned int>& seed = std::nullopt);

	/**
	 * \brief Export mesh data vertices as *.ply point cloud.
	 * \param meshData       input geom data.
	 * \param absFileName    absolute file path for the opened file.
	 * \return if true, the export was successful.
	 */
	[[nodiscard]] bool ExportPointsToPLY(const BaseMeshGeometryData& meshData, const std::string& absFileName);

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

	/// \brief Returns a bounding sphere with a center and a radius combined in a pair.
	///	\throw std::invalid_argument if the mesh is empty.
	[[nodiscard]] std::pair<pmp::Point, pmp::Scalar> ComputeMeshBoundingSphere(const pmp::SurfaceMesh& mesh);

	/// \brief Returns a bounding sphere with a center and a radius combined in a pair.
	///	\throw std::invalid_argument if the point cloud is empty.
	[[nodiscard]] std::pair<pmp::Point, pmp::Scalar> ComputePointCloudBoundingSphere(const std::vector<pmp::Point>& points);

	/**
	 * \brief Computes a mesh from the given point cloud using the Ball-Pivoting algorithm.
	 *        F. Bernardini, J. Mittleman, H. Rushmeier, C. Silva and G. Taubin, "The ball-pivoting algorithm for surface reconstruction," in IEEE Transactions on Visualization and Computer Graphics, vol. 5, no. 4, pp. 349-359, Oct.-Dec. 1999
	 * \param points                                 input point cloud.
	 * \param ballRadius                             ball radius parameter.
	 * \param clusteringPercentageOfBallRadius       this percentage of ballRadius will be used for clustering.
	 * \param angleThreshold                         angle threshold.
	 * \return optional resulting BaseMeshGeometryData if the computation is successful.
	 */
	[[nodiscard]] std::optional<BaseMeshGeometryData> ComputeBallPivotingMeshFromPoints(const std::vector<pmp::Point>& points, 
		const pmp::Scalar& ballRadius, 
		const pmp::Scalar& clusteringPercentageOfBallRadius = 20, 
		const pmp::Scalar& angleThreshold = 90.0);

	/// \brief Computes the minimum distance between points in the input point cloud.
	[[nodiscard]] pmp::Scalar ComputeMinInterVertexDistance(const std::vector<pmp::Point>& points);

	/// \brief Computes the average distance between points in the input point cloud.
	[[nodiscard]] pmp::Scalar ComputeNearestNeighborMeanInterVertexDistance(const std::vector<pmp::Point>& points, const size_t& nNeighbors = 6);

	/// \brief Computes the minimum distance between points in the input point cloud.
	[[nodiscard]] pmp::Scalar ComputeMinInterVertexDistanceBruteForce(const std::vector<pmp::Point>& points);

	/// \brief Computes the maximum distance between points in the input point cloud.
	[[nodiscard]] pmp::Scalar ComputeMaxInterVertexDistanceBruteForce(const std::vector<pmp::Point>& points);

	/// \brief Computes the average distance between points in the input point cloud.
	[[nodiscard]] pmp::Scalar ComputeMeanInterVertexDistanceBruteForce(const std::vector<pmp::Point>& points);

	/// \brief Computes the index of the closest point
	[[nodiscard]] int GetClosestPointIndex(const std::vector<pmp::Point>& points, const pmp::Point& sampledPoint);

	/// \brief Computes the index of the closest point in the x,y plane.
	[[nodiscard]] int GetClosestPointIndex2D(const std::vector<pmp::Point>& points, const pmp::Point& sampledPoint);

	/**
	 * \brief Obtains a 2D slice of a given 3D point cloud from given slicing plane and distance tolerance.
	 * \param points              Processed point cloud.
	 * \param planePt             Reference point of the slicing plane.
	 * \param planeNormal         Unit normal vector of the slicing plane.
	 * \param distTolerance       Distance tolerance for slicing (points which will be closer to the slicing plane than this distance are projected and pushed to the resulting 2D slice).
	 * \return The sliced point cloud.
	 */
	[[nodiscard]] std::vector<pmp::Point2> GetSliceOfThePointCloud(const std::vector<pmp::Point>& points, const pmp::Point& planePt, const pmp::vec3& planeNormal, const float& distTolerance);

} // namespace Geometry