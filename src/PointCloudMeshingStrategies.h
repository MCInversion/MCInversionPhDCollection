#pragma once

#include <vector>

namespace IMB
{
	class PointCloudMeshingStrategy
	{
	public:
		virtual ~PointCloudMeshingStrategy() = default;

		/// =====================================================================================================
		/// \brief Process the input points (with checks and general preprocessing) and generate a mesh.
		/// \param[in,out] ioPoints       The input/output points. DISCLAIMER: Some strategies may replace the input pt list with a different point list.
		/// \param[out] resultPolyIds     The output mesh indexing. Each element is a list of point indices that form a polygon.
		/// =====================================================================================================
		virtual void Process(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds);

	private:
		/// =====================================================================================================
		/// \brief Process the input points and generate a mesh.
		/// \param[in,out] ioPoints       The input/output points. DISCLAIMER: Some strategies may replace the input pt list with a different point list.
		/// \param[out] resultPolyIds     The output mesh indexing. Each element is a list of point indices that form a polygon.
		/// =====================================================================================================
		virtual void ProcessImpl(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds) = 0;
	};

	class BallPivotingMeshingStrategy : public PointCloudMeshingStrategy
	{
	private:
		/// =====================================================================================================
		/// \brief Process the input points and generate a mesh using the ball pivoting algorithm.
		/// \param[in,out] ioPoints       The input/output points. DISCLAIMER: This strategy does not modify the input point list.
		/// \param[out] resultPolyIds     The output mesh indexing. Each element is a list of point indices that form a polygon.
		/// =====================================================================================================
		void ProcessImpl(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds) override;
	};

	class PoissonMeshingStrategy : public PointCloudMeshingStrategy
	{
	private:
		/// =====================================================================================================
		/// \brief Process the input points and generate a mesh using the Poisson reconstruction algorithm.
		/// \param[in,out] ioPoints       The input/output points. DISCLAIMER: This strategy does not modify the input point list.
		/// \param[out] resultPolyIds     The output mesh indexing. Each element is a list of point indices that form a polygon.
		/// =====================================================================================================
		void ProcessImpl(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds) override;
	};

	class MarchingCubesMeshingStrategy : public PointCloudMeshingStrategy
	{
	private:
		/// =====================================================================================================
		/// \brief Process the input points and generate a mesh using the marching cubes algorithm of the point cloud distance field.
		/// \param[in,out] ioPoints       The input/output points. DISCLAIMER: This strategy DOES modify the input point list.
		/// \param[out] resultPolyIds     The output mesh indexing. Each element is a list of point indices that form a polygon.
		/// =====================================================================================================
		void ProcessImpl(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds) override;
	};

	class LagrangianShrinkWrappingMeshingStrategy : public PointCloudMeshingStrategy
	{
	private:
		/// =====================================================================================================
		/// \brief Process the input points and generate a mesh using the Lagrangian shrink wrapping algorithm.
		/// \param[in,out] ioPoints       The input/output points. DISCLAIMER: This strategy DOES modify the input point list.
		/// \param[out] resultPolyIds     The output mesh indexing. Each element is a list of point indices that form a polygon.
		/// =====================================================================================================
		void ProcessImpl(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds) override;
	};
	
} // namespace IMB