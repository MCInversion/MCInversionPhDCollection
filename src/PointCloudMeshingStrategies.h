#pragma once

#include <vector>

namespace IMB
{
	class PointCloudMeshingStrategy
	{
	public:
		virtual ~PointCloudMeshingStrategy() = default;
		virtual void Process(const std::vector<pmp::Point>& points, std::vector<std::vector<unsigned int>>& result) = 0;
	};

	class BallPivotingMeshingStrategy : public PointCloudMeshingStrategy
	{
	public:
		void Process(const std::vector<pmp::Point>& points, std::vector<std::vector<unsigned int>>& result) override;
	};

	class PoissonMeshingStrategy : public PointCloudMeshingStrategy
	{
	public:
		void Process(const std::vector<pmp::Point>& points, std::vector<std::vector<unsigned int>>& result) override;
	};

	class MarchingCubesMeshingStrategy : public PointCloudMeshingStrategy
	{
	public:
		void Process(const std::vector<pmp::Point>& points, std::vector<std::vector<unsigned int>>& result) override;
	};

	class LagrangianShrinkWrappingMeshingStrategy : public PointCloudMeshingStrategy
	{
	public:
		void Process(const std::vector<pmp::Point>& points, std::vector<std::vector<unsigned int>>& result) override;
	};
	
} // namespace IMB