#include "pmp/Types.h"

#include "PointCloudMeshingStrategies.h"

#include <iostream>

namespace IMB
{
	void BallPivotingMeshingStrategy::Process(const std::vector<pmp::Point>& points, std::vector<std::vector<unsigned int>>& result)
	{
		std::cerr << "BallPivotingMeshingStrategy Not implemented\n";
	}

	void PoissonMeshingStrategy::Process(const std::vector<pmp::Point>& points, std::vector<std::vector<unsigned int>>& result)
	{
		std::cerr << "PoissonMeshingStrategy Not implemented\n";
	}

	void MarchingCubesMeshingStrategy::Process(const std::vector<pmp::Point>& points, std::vector<std::vector<unsigned int>>& result)
	{
		std::cerr << "MarchingCubesMeshingStrategy Not implemented\n";
	}

	void LagrangianShrinkWrappingMeshingStrategy::Process(const std::vector<pmp::Point>& points, std::vector<std::vector<unsigned int>>& result)
	{
		std::cerr << "LagrangianShrinkWrappingMeshingStrategy Not implemented\n";
	}
} // namespace IMB
