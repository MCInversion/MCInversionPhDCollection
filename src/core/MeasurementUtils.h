#pragma once

#include "pmp/Types.h"

#include "geometry/Grid.h"
#include "geometry/GeometryConversionUtils.h"

namespace Geometry
{
	/// \brief Computes an interpolated point cloud to mesh vertex distance histogram by computing the distance field.
	/// \param[in] mesh      input mesh.
	/// \param[in] ptCloud   input point cloud. 
	/// \param[in] nBins     the number of histogram bins fitting [maxDistVal, minDistVal] interval.
	/// \return optional histogram output as pair { [maxDistVal, minDistVal] range, bin counts }.
	[[nodiscard]] std::optional<std::pair<std::pair<pmp::Scalar, pmp::Scalar>, std::vector<unsigned int>>>
		ComputeMeshDistanceToPointCloudPerVertexHistogram(const pmp::SurfaceMesh& mesh, const std::vector<pmp::vec3>& ptCloud, const unsigned int& nBins);

	/// \brief Computes Hausdorff distance between mesh and a point cloud.
	/// \param[in] mesh                      input mesh.
	/// \param[in] ptCloud                   input point cloud.
	/// \param[in] nVoxelsPerMinDimension    the key parameter to compute distance voxel field resolution: sampling per minimum bbox dimension.
	/// \return optional evaluated Hausdorff distance dH(X, Y) = max(sup d(x, Y), sup d(X, y)).
	[[nodiscard]] std::optional<double> ComputeMeshToPointCloudHausdorffDistance(
		const pmp::SurfaceMesh& mesh,
		const std::vector<pmp::Point>& ptCloud,
		const unsigned int& nVoxelsPerMinDimension);

	/// \brief Computes Hausdorff distance between mesh and a point cloud.
	/// \param[in] mesh                      input mesh.
	/// \param[in] ptCloud                   input point cloud.
	/// \param[in] ptCloudDf                 input point cloud's pre-computed distance field.
	/// \param[in] nVoxelsPerMinDimension    the key parameter to compute distance voxel field resolution: sampling per minimum bbox dimension.
	/// \return optional evaluated Hausdorff distance dH(X, Y) = max(sup d(x, Y), sup d(X, y)).
	[[nodiscard]] std::optional<double> ComputeMeshToPointCloudHausdorffDistance(
		const pmp::SurfaceMesh& mesh,
		const std::vector<pmp::Point>& ptCloud,
		const ScalarGrid& ptCloudDf, // Use the existing ScalarGrid
		const unsigned int& nVoxelsPerMinDimension);

	/// \brief Computes Hausdorff distance between a mesh and a reference mesh.
	/// \param[in] mesh                      input mesh.
	/// \param[in] refMesh                   input reference mesh.
	/// \param[in] refMeshDf                 input reference mesh's pre-computed distance field.
	/// \param[in] nVoxelsPerMinDimension    the key parameter to compute distance voxel field resolution: sampling per minimum bbox dimension.
	/// \return optional evaluated Hausdorff distance dH(X, Y) = max(sup d(x, Y), sup d(X, y)).
	[[nodiscard]] std::optional<double> ComputeMeshToMeshHausdorffDistance(
		const pmp::SurfaceMesh& mesh,
		const pmp::SurfaceMesh& refMesh,
		const ScalarGrid& refMeshDf, // Use the existing ScalarGrid
		const unsigned int& nVoxelsPerMinDimension);

	/// \brief Computes Hausdorff distance between a mesh and a reference mesh.
	/// \param[in] mesh                      input mesh.
	/// \param[in] refMesh                   input reference mesh.
	/// \param[in] refMeshDf                 input reference mesh's pre-computed distance field.
	/// \param[in] nVoxelsPerMinDimension    the key parameter to compute distance voxel field resolution: sampling per minimum bbox dimension.
	/// \return optional evaluated Hausdorff distance dH(X, Y) = max(sup d(x, Y), sup d(X, y)).
	[[nodiscard]] std::optional<double> ComputeMeshToMeshHausdorffDistance(
		const BaseMeshGeometryData& mesh,
		const BaseMeshGeometryData& refMesh,
		const ScalarGrid& refMeshDf, // Use the existing ScalarGrid
		const unsigned int& nVoxelsPerMinDimension);

} // namespace Geometry