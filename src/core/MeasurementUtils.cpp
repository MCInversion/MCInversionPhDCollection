
#include "geometry/GridUtil.h"
#include "geometry/GeometryAdapters.h"

#include "sdf/SDF.h"

#include "core/MeasurementUtils.h"

namespace Geometry
{
	constexpr unsigned int DEFAULT_N_VOXELS_PER_MIN_DIMENSION = 40;

	std::optional<std::pair<std::pair<float, float>, std::vector<unsigned int>>> ComputeMeshDistanceToPointCloudPerVertexHistogram(const pmp::SurfaceMesh& mesh, const std::vector<pmp::vec3>& ptCloud, const unsigned int& nBins)
	{
		// Check for empty mesh or point cloud
		if (mesh.n_vertices() == 0 || ptCloud.empty())
		{
			return {};
		}

		// Compute distance field for the point cloud
		const pmp::BoundingBox ptCloudBBox(ptCloud);
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float cellSize = minSize / DEFAULT_N_VOXELS_PER_MIN_DIMENSION;
		constexpr float volExpansionFactor = 1.0f;
		const SDF::PointCloudDistanceFieldSettings dfSettings{
			cellSize,
			volExpansionFactor,
			Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::BlurPostprocessingType::None
		};

		auto df = SDF::PointCloudDistanceFieldGenerator::Generate(ptCloud, dfSettings);

		// Prepare for histogram computation
		std::vector<unsigned int> histogram(nBins, 0);
		float maxDistVal = std::numeric_limits<float>::lowest();
		float minDistVal = std::numeric_limits<float>::max();

		// Iterate through vertices to compute distances and populate the histogram
		for (const auto& v : mesh.vertices())
		{
			const auto vPos = mesh.position(v);
			const double vDistanceToTarget = TrilinearInterpolateScalarValue(vPos, df);
			maxDistVal = std::max(maxDistVal, static_cast<float>(vDistanceToTarget));
			minDistVal = std::min(minDistVal, static_cast<float>(vDistanceToTarget));
		}

		// Calculate histogram bin size and populate histogram
		const float binSize = (maxDistVal - minDistVal) / nBins;
		for (const auto& v : mesh.vertices())
		{
			const auto vPos = mesh.position(v);
			const double vDistanceToTarget = TrilinearInterpolateScalarValue(vPos, df);
			int binIndex = static_cast<int>((vDistanceToTarget - minDistVal) / binSize);
			binIndex = std::min(binIndex, static_cast<int>(nBins) - 1); // Ensure binIndex is within range
			histogram[binIndex]++;
		}

		return { {{minDistVal, maxDistVal}, histogram} };
	}

	std::optional<double> ComputeMeshToPointCloudHausdorffDistance(const pmp::SurfaceMesh& mesh, const std::vector<pmp::Point>& ptCloud, const unsigned int& nVoxelsPerMinDimension)
	{
		if (mesh.n_vertices() == 0 || ptCloud.empty())
		{
			return {};
		}

		// Compute distance field for the point cloud
		const pmp::BoundingBox ptCloudBBox(ptCloud);
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const float ptCloudMinSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float ptCloudCellSize = ptCloudMinSize / static_cast<float>(nVoxelsPerMinDimension);
		const SDF::PointCloudDistanceFieldSettings ptCloudDfSettings{
			ptCloudCellSize,
			1.0f, // volExpansionFactor
			DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::BlurPostprocessingType::None
		};
		const auto ptCloudDf = SDF::PointCloudDistanceFieldGenerator::Generate(ptCloud, ptCloudDfSettings);

		// Compute distance field for the mesh
		const auto meshBBox = mesh.bounds();
		const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
		const float meshMinSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		const float meshCellSize = meshMinSize / static_cast<float>(nVoxelsPerMinDimension);
		const SDF::DistanceFieldSettings meshDfSettings{
			meshCellSize,
			1.0f, // volExpansionFactor
			DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::KDTreeSplitType::Center,
			SDF::SignComputation::None, // Unsigned distance field
			SDF::BlurPostprocessingType::None,
			SDF::PreprocessingType::Octree
		};
		const PMPSurfaceMeshAdapter meshAdapter(std::make_shared<pmp::SurfaceMesh>(mesh));
		const auto meshDf = SDF::DistanceFieldGenerator::Generate(meshAdapter, meshDfSettings);

		double maxDistMeshToPointCloud = std::numeric_limits<double>::lowest();
		double maxDistPointCloudToMesh = std::numeric_limits<double>::lowest();

		// Mesh to Point Cloud: Compute max distance using the distance field
		for (const auto v : mesh.vertices())
		{
			const auto vPos = mesh.position(v);
			const double vDistanceToPointCloud = TrilinearInterpolateScalarValue(vPos, ptCloudDf);
			maxDistMeshToPointCloud = std::max(maxDistMeshToPointCloud, vDistanceToPointCloud);
		}

		// Point Cloud to Mesh: Compute max distance using the distance field
		for (const auto& p : ptCloud)
		{
			const double pDistanceToMesh = TrilinearInterpolateScalarValue(p, meshDf);
			maxDistPointCloudToMesh = std::max(maxDistPointCloudToMesh, pDistanceToMesh);
		}

		// Compute Hausdorff Distance as the maximum of these two distances
		return std::max(maxDistMeshToPointCloud, maxDistPointCloudToMesh);
	}

	std::optional<double> ComputeMeshToPointCloudHausdorffDistance(const pmp::SurfaceMesh& mesh, const std::vector<pmp::Point>& ptCloud, const ScalarGrid& ptCloudDf, const unsigned int& nVoxelsPerMinDimension)
	{
		if (mesh.n_vertices() == 0)
		{
			return {};
		}

		// Compute distance field for the mesh
		const auto meshBBox = mesh.bounds();
		const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
		const float meshMinSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		const float meshCellSize = meshMinSize / static_cast<float>(nVoxelsPerMinDimension);
		const SDF::DistanceFieldSettings meshDfSettings{
			meshCellSize,
			1.0f, // volExpansionFactor
			DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::KDTreeSplitType::Center,
			SDF::SignComputation::None, // Unsigned distance field
			SDF::BlurPostprocessingType::None,
			SDF::PreprocessingType::Octree
		};
		const PMPSurfaceMeshAdapter meshAdapter(std::make_shared<pmp::SurfaceMesh>(mesh));
		const auto meshDf = SDF::DistanceFieldGenerator::Generate(meshAdapter, meshDfSettings);

		double maxDistMeshToPointCloud = std::numeric_limits<double>::lowest();
		double maxDistPointCloudToMesh = std::numeric_limits<double>::lowest();

		// Mesh to Point Cloud: Compute max distance using the distance field
		for (const auto v : mesh.vertices())
		{
			const auto vPos = mesh.position(v);
			const double vDistanceToPointCloud = TrilinearInterpolateScalarValue(vPos, ptCloudDf);
			maxDistMeshToPointCloud = std::max(maxDistMeshToPointCloud, vDistanceToPointCloud);
		}

		// Point Cloud to Mesh: Compute max distance using the distance field
		for (const auto& p : ptCloud)
		{
			const double pDistanceToMesh = TrilinearInterpolateScalarValue(p, meshDf);
			maxDistPointCloudToMesh = std::max(maxDistPointCloudToMesh, pDistanceToMesh);
		}

		// Compute Hausdorff Distance as the maximum of these two distances
		return std::max(maxDistMeshToPointCloud, maxDistPointCloudToMesh);
	}

	std::optional<double> ComputeMeshToMeshHausdorffDistance(const pmp::SurfaceMesh& mesh, const pmp::SurfaceMesh& refMesh, const ScalarGrid& refMeshDf, const unsigned int& nVoxelsPerMinDimension)
	{
		if (refMesh.is_empty() || mesh.is_empty())
		{
			return {};
		}

		// Compute distance field for the mesh
		const auto meshBBox = mesh.bounds();
		const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
		const float meshMinSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		const float meshCellSize = meshMinSize / static_cast<float>(nVoxelsPerMinDimension);
		const SDF::DistanceFieldSettings meshDfSettings{
			meshCellSize,
			1.0f, // volExpansionFactor
			DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::KDTreeSplitType::Center,
			SDF::SignComputation::None, // Unsigned distance field
			SDF::BlurPostprocessingType::None,
			SDF::PreprocessingType::Octree
		};
		const PMPSurfaceMeshAdapter meshAdapter(std::make_shared<pmp::SurfaceMesh>(mesh));
		const auto meshDf = SDF::DistanceFieldGenerator::Generate(meshAdapter, meshDfSettings);

		double maxDistMeshToRefMesh = std::numeric_limits<double>::lowest();
		double maxDistRefMeshToMesh = std::numeric_limits<double>::lowest();

		// Mesh to Ref Mesh: Compute max distance using the distance field
		for (const auto v : mesh.vertices())
		{
			const auto vPos = mesh.position(v);
			const double vDistanceToRefMesh = TrilinearInterpolateScalarValue(vPos, refMeshDf);
			maxDistMeshToRefMesh = std::max(maxDistMeshToRefMesh, vDistanceToRefMesh);
		}

		// Ref Mesh to Mesh: Compute max distance using the distance field
		for (const auto v : refMesh.vertices())
		{
			const auto vPos = refMesh.position(v);
			const double pDistanceToMesh = TrilinearInterpolateScalarValue(vPos, meshDf);
			maxDistRefMeshToMesh = std::max(maxDistRefMeshToMesh, pDistanceToMesh);
		}

		// Compute Hausdorff Distance as the maximum of these two distances
		return std::max(maxDistMeshToRefMesh, maxDistRefMeshToMesh);
	}

	std::optional<double> ComputeMeshToMeshHausdorffDistance(const BaseMeshGeometryData& mesh, const BaseMeshGeometryData& refMesh, const ScalarGrid& refMeshDf, const unsigned int& nVoxelsPerMinDimension)
	{
		if (refMesh.Vertices.empty() || mesh.Vertices.empty())
		{
			return {};
		}

		// Compute distance field for the mesh
		const BaseMeshAdapter meshAdapter(std::make_shared<BaseMeshGeometryData>(mesh));
		const auto meshBBox = meshAdapter.GetBounds();
		const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
		const float meshMinSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		const float meshCellSize = meshMinSize / static_cast<float>(nVoxelsPerMinDimension);
		const SDF::DistanceFieldSettings meshDfSettings{
			meshCellSize,
			1.0f, // volExpansionFactor
			DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::KDTreeSplitType::Center,
			SDF::SignComputation::None, // Unsigned distance field
			SDF::BlurPostprocessingType::None,
			SDF::PreprocessingType::Octree
		};
		const auto meshDf = SDF::DistanceFieldGenerator::Generate(meshAdapter, meshDfSettings);

		double maxDistMeshToRefMesh = std::numeric_limits<double>::lowest();
		double maxDistRefMeshToMesh = std::numeric_limits<double>::lowest();

		// Mesh to Ref Mesh: Compute max distance using the distance field
		for (const auto& v : mesh.Vertices)
		{
			const double vDistanceToRefMesh = TrilinearInterpolateScalarValue(v, refMeshDf);
			maxDistMeshToRefMesh = std::max(maxDistMeshToRefMesh, vDistanceToRefMesh);
		}

		// Ref Mesh to Mesh: Compute max distance using the distance field
		for (const auto& v : refMesh.Vertices)
		{
			const double pDistanceToMesh = TrilinearInterpolateScalarValue(v, meshDf);
			maxDistRefMeshToMesh = std::max(maxDistRefMeshToMesh, pDistanceToMesh);
		}

		// Compute Hausdorff Distance as the maximum of these two distances
		return std::max(maxDistMeshToRefMesh, maxDistRefMeshToMesh);
	}
	
} // namespace Geometry
