#include "PointCloudMeshingStrategies.h"

#include "EvolverUtilsCommon.h"
#include "SurfaceEvolver.h"

#include "geometry/GeometryConversionUtils.h"

#include "sdf/SDF.h"

#include <iostream>

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

	constexpr pmp::Scalar MAGIC_RADIUS_MULTIPLIER = 2.9;

	void BallPivotingMeshingStrategy::ProcessImpl(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds)
	{
		std::cout << "BallPivotingMeshingStrategy::ProcessImpl: attempting to triangulate a mesh with " << ioPoints.size() << " vertices.\n";

		// Compute an appropriate radius based on point distribution
		const auto meanDistance = Geometry::ComputeNearestNeighborMeanInterVertexDistance(ioPoints, 6);
		const auto ballRadius = meanDistance * MAGIC_RADIUS_MULTIPLIER;
		constexpr auto clusteringPercentage = 40.0; // (1.0 / MAGIC_RADIUS_MULTIPLIER) * 100.0;

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
		std::cout << "PoissonMeshingStrategy::ProcessImpl: attempting to triangulate a mesh with " << ioPoints.size() << " vertices.\n";


	}

	void MarchingCubesMeshingStrategy::ProcessImpl(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds)
	{
		std::cerr << "MarchingCubesMeshingStrategy::ProcessImpl: Not implemented!\n";
		std::cerr << "MarchingCubesMeshingStrategy::ProcessImpl: attempting to triangulate a mesh with " << ioPoints.size() << " vertices.\n";
	}

	void LagrangianShrinkWrappingMeshingStrategy::ProcessImpl(std::vector<pmp::Point>& ioPoints, std::vector<std::vector<unsigned int>>& resultPolyIds)
	{
		std::cerr << "LagrangianShrinkWrappingMeshingStrategy::ProcessImpl: attempting to triangulate a mesh with " << ioPoints.size() << " vertices.\n";	

		constexpr size_t NTimeSteps = 40;
		constexpr unsigned int nVoxelsPerMinDimension = 40;
		constexpr double tau = 0.1;
		constexpr double defaultOffsetFactor = 1.0;

		const pmp::BoundingBox ptCloudBBox(ioPoints);
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
		constexpr pmp::Scalar volExpansionFactor = 0.5;
		const SDF::PointCloudDistanceFieldSettings dfSettings{
			cellSize,
				volExpansionFactor,
				Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
				SDF::BlurPostprocessingType::None
		};
		auto sdf = SDF::PointCloudDistanceFieldGenerator::Generate(ioPoints, dfSettings);

		MeshTopologySettings topoParams;
		topoParams.FixSelfIntersections = true;
		topoParams.MinEdgeMultiplier = 0.14;
		topoParams.UseBackProjection = false;
		topoParams.PrincipalCurvatureFactor = 3.2;
		topoParams.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);
		topoParams.EdgeLengthDecayFactor = 0.7;
		topoParams.ExcludeEdgesWithoutBothFeaturePts = true;
		topoParams.FeatureType = FeatureDetectionType::MeanCurvature;

		AdvectionDiffusionParameters adParams{
			1.0, 1.0,
			2.0, 1.0
		};
		const double fieldIsoLevel = defaultOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);
		SurfaceEvolutionSettings seSettings{
			"IncrementalLSWMesh",
			NTimeSteps,
			tau,
			fieldIsoLevel,
			2, // IcoSphereSubdivisionLevel
			adParams,
			topoParams,
			minSize, maxSize,
			ptCloudBBox.center(),
			false, false,
			"",
			MeshLaplacian::Voronoi,
			{"equilateralJacobianCondition"},
			0.05,
			true,
			false
		};
		//ReportInput(seSettings, std::cout);
		SurfaceEvolver evolver(sdf, volExpansionFactor, seSettings);

		try
		{
			evolver.Evolve();
		}
		catch (...)
		{
			std::cerr << "> > > > > > > > > > > > > > SurfaceEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
		}

		const auto& resultMesh = evolver.GetResultSurface();
		const auto resultBaseMesh = Geometry::ConvertPMPSurfaceMeshToBaseMeshGeometryData(resultMesh);
		ioPoints = resultBaseMesh.Vertices;
		resultPolyIds = resultBaseMesh.PolyIndices;
	}
} // namespace IMB
