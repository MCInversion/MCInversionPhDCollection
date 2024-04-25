#include "PointCloudMeshingStrategies.h"

#include <iostream>

#include "geometry/GeometryConversionUtils.h"

namespace IMB
{
	void PointCloudMeshingStrategy::Process(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds)
	{
		if (ioPoints.empty())
		{
			std::cerr << "PointCloudMeshingStrategy::Process: No points to triangulate.\n";
			return;
		}
		// TODO: this is the right place for the new dynamic mesh data structure which will be the novel research contribution of this project.
		resultPolyIds.clear(); // there is a pre-existing triangulation from previous steps, so clear it
		ProcessImpl(ioPoints, resultPolyIds);
	}

	void EmptyMeshingStrategy::ProcessImpl(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds)
	{
		std::cerr << "EmptyMeshingStrategy::ProcessImpl: attempting to NOT triangulate a mesh with " << ioPoints.size() << " vertices.\n";
	}

	constexpr float MAGIC_RADIUS_MULTIPLIER = 1.9f;

	void BallPivotingMeshingStrategy::ProcessImpl(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds)
	{
		std::cout << "BallPivotingMeshingStrategy::ProcessImpl: attempting to triangulate a mesh with " << ioPoints.size() << " vertices.\n";

		// Compute an appropriate radius based on point distribution
		const auto meanDistance = Geometry::ComputeNearestNeighborMeanInterVertexDistance(ioPoints, 6);
		const auto ballRadius = meanDistance * MAGIC_RADIUS_MULTIPLIER;
		constexpr auto clusteringPercentage = (1.0f / MAGIC_RADIUS_MULTIPLIER) * 100.0f;

		const auto meshDataOpt = Geometry::ComputeBallPivotingMeshFromPoints(ioPoints, ballRadius, clusteringPercentage);
		if (!meshDataOpt.has_value())
		{
			std::cerr << "BallPivotingMeshingStrategy::ProcessImpl: Internal algorithm error!\n";
			return;
		}
		const auto& meshData = meshDataOpt.value();
		resultPolyIds = meshData.PolyIndices;
	}

	void PoissonMeshingStrategy::ProcessImpl(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds)
	{
		std::cerr << "PoissonMeshingStrategy::ProcessImpl: Not implemented!\n";
		std::cerr << "PoissonMeshingStrategy::ProcessImpl: attempting to triangulate a mesh with " << ioPoints.size() << " vertices.\n";
	}

	void MarchingCubesMeshingStrategy::ProcessImpl(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds)
	{
		std::cerr << "MarchingCubesMeshingStrategy::ProcessImpl: Not implemented!\n";
		std::cerr << "MarchingCubesMeshingStrategy::ProcessImpl: attempting to triangulate a mesh with " << ioPoints.size() << " vertices.\n";
	}

	void LagrangianShrinkWrappingMeshingStrategy::ProcessImpl(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds)
	{
		std::cerr << "LagrangianShrinkWrappingMeshingStrategy::ProcessImpl: Not implemented!\n";
		std::cerr << "LagrangianShrinkWrappingMeshingStrategy::ProcessImpl: attempting to triangulate a mesh with " << ioPoints.size() << " vertices.\n";
	}
} // namespace IMB
