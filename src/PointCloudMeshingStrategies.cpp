#include "pmp/Types.h"

#include "PointCloudMeshingStrategies.h"

#include <iostream>

namespace IMB
{
	void BallPivotingMeshingStrategy::Process(const std::vector<pmp::Point>& points, std::vector<std::vector<unsigned int>>& result)
	{
		std::cerr << "BallPivotingMeshingStrategy: Not implemented!\n";
		std::cerr << "BallPivotingMeshingStrategy: attempting to triangulate a mesh with " << points.size() << " vertices.\n";
	}

	void PoissonMeshingStrategy::Process(const std::vector<pmp::Point>& points, std::vector<std::vector<unsigned int>>& result)
	{
		std::cerr << "PoissonMeshingStrategy: Not implemented!\n";
		std::cerr << "PoissonMeshingStrategy: attempting to triangulate a mesh with " << points.size() << " vertices.\n";
	}

	void MarchingCubesMeshingStrategy::Process(const std::vector<pmp::Point>& points, std::vector<std::vector<unsigned int>>& result)
	{
		std::cerr << "MarchingCubesMeshingStrategy: Not implemented!\n";
		std::cerr << "MarchingCubesMeshingStrategy: attempting to triangulate a mesh with " << points.size() << " vertices.\n";
	}

	void LagrangianShrinkWrappingMeshingStrategy::Process(const std::vector<pmp::Point>& points, std::vector<std::vector<unsigned int>>& result)
	{
		std::cerr << "LagrangianShrinkWrappingMeshingStrategy: Not implemented!\n";
		std::cerr << "LagrangianShrinkWrappingMeshingStrategy: attempting to triangulate a mesh with " << points.size() << " vertices.\n";
	}
} // namespace IMB
