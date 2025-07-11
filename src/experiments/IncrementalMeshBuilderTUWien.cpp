// ------------------------- IncrementalMeshBuilderTUWien ------------------------------
// Late Spring / Summer of 2025. TU Wien June 6th Talk.
// ..................................................................................
/*
 * This talk builds on the previously introduced Incremental Mesh Builder (IMB) framework for progressive 
 * Level-of-Detail (LOD) mesh reconstruction from large 3D datasets. The current focus is on evaluating 
 * the viability of a geometry-aware stochastic vertex sampling process as a foundation for adaptive mesh 
 * reconstruction — balancing simplification speed with feature preservation.
 *
 * The latest developments test this sampling approach on both structured and noisy 3D datasets, using 
 * simplified artificial heightmaps and real-world scan-derived meshes. Performance results include 
 * comparisons between uniform and softmax-based sampling, along with LOD convergence times and visual 
 * fidelity metrics.
 *
 * The talk will also outline the ongoing effort to shape this framework into a lightweight, open-source 
 * application tailored to experts in 3D scanning, offering an efficient and scalable way to browse 
 * massive raw mesh files without full preprocessing.
 *
 * DISCLAIMER: This is work-in-progress!
 */
 // --------------------------------------------------------------------------------

#include "utils/TimingUtils.h"
#include "utils/FileMappingWrapper.h"

#include "geometry/GeometryIOUtils.h"

#include "core/IMB_ShrinkWrapper.h"
#include "core/IncrementalProgressUtils.h"
#include "core/IncrementalMeshFileHandler.h"
#include "core/PointCloudMeshingStrategies.h"
#include "core/VertexSamplingStrategies.h"

#include "IOEnvironment.h"

#include "../Experiments.h"

void PointCloudClustering()
{
	const auto ptCloudName = "spheres";
	const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataDirPath + ptCloudName + ".ply", true);
	if (!ptCloudOpt.has_value())
	{
		std::cerr << "PointCloudClustering: ptCloudOpt == nullopt!\n";
		return;
	}
	const auto& pts = *ptCloudOpt;

	const pmp::Scalar criticalRadius = Geometry::ComputeNearestNeighborMeanInterVertexDistance(pts, 10) * 2.5;
	std::cout << "PointCloudClustering: criticalRadius evaluated to: " << criticalRadius << ".\n";
	const auto ptClusters = Geometry::GetPointClusters(pts, criticalRadius);
	constexpr size_t nExpected = 4;
	std::cout << "PointCloudClustering: Found " << ptClusters.size() << (ptClusters.size() == nExpected ? " == " : " != ") << std::to_string(nExpected) << " clusters in " << pts.size() << " points.\n";

	for (int i = 0; i < ptClusters.size(); ++i)
	{
		std::cout << "PointCloudClustering: cluster " << std::to_string(i) << ": " << ptClusters[i].size() << " pts.\n";
		const std::string outputFileName = dataOutPath + ptCloudName + "_Cluster" + std::to_string(i) + ".ply";
		if (!Geometry::Export3DPointCloudToPLY(ptClusters[i], outputFileName))
		{
			std::cerr << "PointCloudClustering: Failed to export point cloud." << "\n";
			return;
		}
	}
}

void PointCloudClusteringPipeline()
{
	const auto ptCloudName = "spheres";
	const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataDirPath + ptCloudName + ".ply", true);
	if (!ptCloudOpt.has_value())
	{
		std::cerr << "PointCloudClusteringPipeline: ptCloudOpt == nullopt!\n";
		return;
	}
	const auto& pts = *ptCloudOpt;

	const auto pt3DIndex = Geometry::Get3DPointSearchIndex(pts);
	if (!pt3DIndex)
	{
		std::cerr << "PointCloudClusteringPipeline: pt3DIndex == nullptr!\n";
		return;
	}

	const auto& ptCloud = pt3DIndex->cloud;
	auto& kdTree = pt3DIndex->tree;

	const pmp::Scalar criticalRadius = Geometry::ComputeNearestNeighborMeanInterVertexDistance(ptCloud, kdTree, 10) * 2.5;
	std::cout << "PointCloudClusteringPipeline: criticalRadius evaluated to: " << criticalRadius << ".\n";
	const auto ptClusters = Geometry::GetPointClusters(ptCloud, kdTree, criticalRadius);
	constexpr size_t nExpected = 4;
	std::cout << "PointCloudClusteringPipeline: Found " << ptClusters.size() << (ptClusters.size() == nExpected ? " == " : " != ") << std::to_string(nExpected) << " clusters in " << pts.size() << " points.\n";

	for (int i = 0; i < ptClusters.size(); ++i)
	{
		std::cout << "PointCloudClusteringPipeline: cluster " << std::to_string(i) << ": " << ptClusters[i].size() << " pts.\n";
		const std::string outputFileName = dataOutPath + ptCloudName + "_Cluster" + std::to_string(i) + ".ply";
		if (!Geometry::Export3DPointCloudToPLY(ptClusters[i], outputFileName))
		{
			std::cerr << "PointCloudClusteringPipeline: Failed to export point cloud." << "\n";
			return;
		}
	}
}

void PointCloudNormalsVCG()
{
	const auto ptCloudName = "bunnyPts_3";
	const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataDirPath + ptCloudName + ".ply", true);
	if (!ptCloudOpt.has_value())
	{
		std::cerr << "PointCloudNormalsVCG: ptCloudOpt == nullopt!\n";
		return;
	}
	const auto& pts = *ptCloudOpt;

	constexpr size_t nNeighbors = 10;
	constexpr size_t smoothingIters = 1;
	const pmp::Point viewPoint{ 0, 0, 0 };
	constexpr bool useViewPoint = false;
	const auto normals = Geometry::EstimatePointCloudNormalsVCG(pts, nNeighbors, smoothingIters, viewPoint, useViewPoint);

	Geometry::BaseMeshGeometryData orientedPts;
	orientedPts.Vertices = pts;
	orientedPts.VertexNormals = normals;

	const std::string outputFileName = dataOutPath + ptCloudName + "_WithVCGNormals.ply";
	if (!Geometry::ExportPointsToPLY(orientedPts, outputFileName))
	{
		std::cerr << "PointCloudNormalsVCG: Failed to export point cloud." << "\n";
		return;
	}
}

void TestPoissonMeshingStrategy()
{
	const auto ptCloudName = "bunnyPts_3";
	const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataDirPath + ptCloudName + ".ply", true);
	if (!ptCloudOpt.has_value())
	{
		std::cerr << "ptCloudOpt == nullopt!\n";
		return;
	}
	const auto& pts = *ptCloudOpt;

	const auto meshingStrategy = IMB::GetReconstructionStrategy(IMB::ReconstructionFunctionType::Poisson);
	Geometry::BaseMeshGeometryData mesh;
	mesh.Vertices = pts;
	meshingStrategy->Process(mesh.Vertices, mesh.PolyIndices);

	const std::string outputFileName = dataOutPath + ptCloudName + "_PoissonRecon.vtk";
	if (!Geometry::ExportBaseMeshGeometryDataToVTK(mesh, outputFileName))
	{
		std::cerr << "Failed to export mesh data." << "\n";
		return;
	}
}

void TestPoissonMeshingWithClustering()
{
	const auto ptCloudName = "spheres";
	const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataDirPath + ptCloudName + ".ply", true);
	if (!ptCloudOpt.has_value())
	{
		std::cerr << "PointCloudClusteringPipeline: ptCloudOpt == nullopt!\n";
		return;
	}
	const auto& pts = *ptCloudOpt;

	const auto meshingStrategy = IMB::GetReconstructionStrategy(IMB::ReconstructionFunctionType::Poisson);
	Geometry::BaseMeshGeometryData mesh;
	mesh.Vertices = pts;
	meshingStrategy->Process(mesh.Vertices, mesh.PolyIndices);

	const std::string outputFileName = dataOutPath + ptCloudName + "_PoissonRecon.vtk";
	if (!Geometry::ExportBaseMeshGeometryDataToVTK(mesh, outputFileName))
	{
		std::cerr << "PointCloudClusteringPipeline: Failed to export mesh data." << "\n";
		return;
	}
}

void TestIMBShrinkWrapper()
{
	const auto ptCloudName = "bunnyPts_3";
	const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataDirPath + ptCloudName + ".ply", true);
	if (!ptCloudOpt.has_value())
	{
		std::cerr << "TestIMBShrinkWrapper: ptCloudOpt == nullopt!\n";
		return;
	}
	const auto& pts = *ptCloudOpt;

	const auto pt3DIndex = Geometry::Get3DPointSearchIndex(pts);
	if (!pt3DIndex)
	{
		std::cerr << "TestIMBShrinkWrapper: pt3DIndex == nullptr!\n";
		return;
	}

	const auto& ptCloud = pt3DIndex->cloud;
	auto& kdTree = pt3DIndex->tree;

	const pmp::Scalar criticalRadius = Geometry::ComputeNearestNeighborMeanInterVertexDistance(ptCloud, kdTree, 10) * 2.5;

	// ~~~
	
	IMB_ShrinkWrapperSettings swSettings;
	swSettings.ProcedureName = ptCloudName;

	swSettings.LevelOfDetail = 3;

	swSettings.TimeStep = 0.05;

	swSettings.TangentialVelocityWeight = 0.05;
	swSettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
	swSettings.RemeshingSettings.UseBackProjection = true;

	swSettings.FieldSettings.NVoxelsPerMinDimension = 20;
	swSettings.FieldSettings.FieldIsoLevel = 1.3;

	swSettings.PreStep = [&](unsigned int step, const SurfaceInOrigScaleFunction& getSurface) {
		const auto surface = getSurface();
		if (!surface)
			return;
		//std::cout << "Exporting " << swSettings.ProcedureName << ", step " << step << " ... ";
		surface->write(dataOutPath + swSettings.ProcedureName + "_Evol_" + std::to_string(step) + ".vtk");
		//std::cout << "done.\n";
	};

	//const double minDistancePercentageEpsilon = 0.1;
	//const double minDistancePercentageEta = 0.1;
	//swSettings.Epsilon.Bind(minDistancePercentageEpsilon, [](double distance)
	//{
	//	return 0.5 * (1.0 - exp(-distance * distance / 1.0));
	//});
	//swSettings.Eta.Bind(minDistancePercentageEta, [](double distance, double negGradDotNormal)
	//{
	//	return -1.0 * distance * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
	//});

	swSettings.Epsilon = [](double distance)
	{
		return 0.5 * (1.0 - exp(-distance * distance / 1.0));
	};
	swSettings.Eta = [](double distance, double negGradDotNormal)
	{
		return -1.0 * distance * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
	};

	swSettings.PointActivationRadius = criticalRadius * 4.0;
	swSettings.PointActivationAlignmentAngle = M_PI_2 * 0.4;
	swSettings.ActivatedPointPercentageThreshold = 1.0; //0.98;

	IMB_ShrinkWrapper sw{ swSettings, pts };

	auto result = sw.Perform();
	if (!result)
	{
		std::cerr << "TestIMBShrinkWrapper: internal error!\n";
		return;
	}

	// ~~~

	const std::string outputFileName = dataOutPath + ptCloudName + "_IMBSW.vtk";
	result->write(outputFileName);
}

void TestDeadlinedIMBShrinkWrapper()
{
	const auto ptCloudName = "bunnyPts_3";
	const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataDirPath + ptCloudName + ".ply", true);
	if (!ptCloudOpt.has_value())
	{
		std::cerr << "TestDeadlinedIMBShrinkWrapper: ptCloudOpt == nullopt!\n";
		return;
	}
	const auto& pts = *ptCloudOpt;

	const auto pt3DIndex = Geometry::Get3DPointSearchIndex(pts);
	if (!pt3DIndex)
	{
		std::cerr << "TestDeadlinedIMBShrinkWrapper: pt3DIndex == nullptr!\n";
		return;
	}

	const auto& ptCloud = pt3DIndex->cloud;
	auto& kdTree = pt3DIndex->tree;

	const pmp::Scalar criticalRadius = Geometry::ComputeNearestNeighborMeanInterVertexDistance(ptCloud, kdTree, 10) * 2.5;

	std::optional<pmp::SurfaceMesh> result;
	if (!Utils::DeadlinedScope::Run(
		[&]() {
			IMB_ShrinkWrapperSettings swSettings;
			swSettings.ProcedureName = ptCloudName;

			swSettings.LevelOfDetail = 3;

			swSettings.TimeStep = 0.05;

			swSettings.TangentialVelocityWeight = 0.05;
			swSettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
			swSettings.RemeshingSettings.UseBackProjection = true;

			swSettings.FieldSettings.NVoxelsPerMinDimension = 40;
			swSettings.FieldSettings.FieldIsoLevel = 0.5;

			swSettings.Epsilon = [](double distance)
			{
				return 0.5 * (1.0 - exp(-distance * distance / 1.0));
			};
			swSettings.Eta = [](double distance, double negGradDotNormal)
			{
				return -1.0 * distance * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};

			swSettings.PointActivationRadius = criticalRadius * 0.95;
			swSettings.ActivatedPointPercentageThreshold = 0.8;

			IMB_ShrinkWrapper sw{ swSettings, pts };

			auto r = sw.Perform();
			if (!r)
			{
				std::cerr << "TestDeadlinedIMBShrinkWrapper: internal error!\n";
				return;
			}
			result = std::move(*r);
		},
		std::chrono::seconds(2) /* TODO: do a proper analysis of the require time limit */))
	{
		std::cerr << "TestDeadlinedIMBShrinkWrapper: Time for IMB_ShrinkWrapper ran out!\n";
		return;
	}

	const std::string outputFileName = dataOutPath + ptCloudName + "_IMBSW.vtk";
	result->write(outputFileName);
}

void TestIMBShrinkWrapperNormalEstimation()
{
	const auto ptCloudName = "bunnyPts_3";
	const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataDirPath + ptCloudName + ".ply", true);
	if (!ptCloudOpt.has_value())
	{
		std::cerr << "ptCloudOpt == nullopt!\n";
		return;
	}
	const auto& pts = *ptCloudOpt;

	const auto meshingStrategy = IMB::GetReconstructionStrategy(IMB::ReconstructionFunctionType::LagrangianShrinkWrapping);
	Geometry::BaseMeshGeometryData mesh;
	mesh.Vertices = pts;
	meshingStrategy->Process(mesh.Vertices, mesh.PolyIndices);

	const std::string outputFileName = dataOutPath + ptCloudName + "_LSWRecon.vtk";
	if (!Geometry::ExportBaseMeshGeometryDataToVTK(mesh, outputFileName))
	{
		std::cerr << "Failed to export mesh data." << "\n";
		return;
	}
}

void SingleThreadSoftMaxUniformStrategy()
{
	const auto meshName = "sampledIco";
	const auto fileName = dataDirPath + meshName + ".ply";
	const auto fileMapping = std::make_unique<Utils::FileMappingWrapper>(fileName);
	if (!fileMapping->IsValid())
	{
		std::cerr << "SingleThreadSoftMaxUniformStrategy: Error during initialization of FileMappingWrapper!\n";
		return;
	}
	const char* fileStart = fileMapping->GetFileMemory();
	const char* fileEnd = fileStart + fileMapping->GetFileSize();
	const auto fileHandler = IMB::CreateMeshFileHandler(IMB::FileHandlerParams{ fileName, fileStart, fileEnd });
	if (!fileHandler)
	{
		std::cerr << "SingleThreadSoftMaxUniformStrategy: Error during initialization of IncrementalMeshFileHandler!\n";
		return;
	}

	constexpr unsigned int seed = 4999;
	constexpr unsigned int nSamples = 6;
	constexpr unsigned int maxVertexCount = 5000;
	const auto vertexSamplingStrategy = IMB::GetVertexSelectionStrategy(
		IMB::VertexSelectionType::SoftMaxUniform,
		nSamples, maxVertexCount, fileHandler);
	// to compare against
	const auto vertexSamplingStrategy2 = IMB::GetVertexSelectionStrategy(
		IMB::VertexSelectionType::UniformRandom,
		nSamples, maxVertexCount, fileHandler);
	if (!vertexSamplingStrategy || !vertexSamplingStrategy2)
	{
		std::cerr << "SingleThreadSoftMaxUniformStrategy: Error during initialization of VertexSamplingStrategy!\n";
		return;
	}

	IMB::GeometricSamplingParams params;
	params.ErrorMetricGradientMultiplier = 1.0;
	params.DistributionMetricGradientMultiplier = 10'000.0; //0.000'1;
	params.TargetVertexDensity = 5.0;
	params.AvgNeighborhoodDisplacementMultiplier = 1.0;
	params.SelectionSharpness = 0.001;
	params.ConfidenceGM1Imax = 0.05;
	params.ConfidenceNormal = 0.8;

	dynamic_cast<IMB::SoftmaxUniformVertexSamplingStrategy*>(vertexSamplingStrategy.get())->Params = params;

	std::vector<pmp::Point> result;
	std::vector<pmp::Point> result2;
	unsigned int lodIndex = 0;
	unsigned int lodIndex2 = 0;
	const auto lodStepCallback = [&] {
		/*const */std::string outputFileName = dataOutPath + meshName + "IMB_LOD" + std::to_string(lodIndex) + ".ply";
		//++lodIndex;
		if (!Geometry::Export3DPointCloudToPLY(result, outputFileName))
		{
			std::cout << "SingleThreadSoftMaxUniformStrategy: Failed to export sampled point data." << "\n";
			return;
		}

		const auto meshingStrategy = IMB::GetReconstructionStrategy(IMB::ReconstructionFunctionType::BallPivoting);
		Geometry::BaseMeshGeometryData mesh;
		mesh.Vertices = result;
		meshingStrategy->Process(mesh.Vertices, mesh.PolyIndices);

		/*const std::string*/ outputFileName = dataOutPath + meshName + "IMB_LOD_BallPivoting" + std::to_string(lodIndex) + ".vtk";
		++lodIndex;
		if (!Geometry::ExportBaseMeshGeometryDataToVTK(mesh, outputFileName))
		{
			std::cerr << "SingleThreadSoftMaxUniformStrategy: Failed to export mesh data." << "\n";
			return;
		}
	};
	const auto lodStepCallback2 = [&] {
		/*const */std::string outputFileName = dataOutPath + meshName + "IMB2_LOD" + std::to_string(lodIndex2) + ".ply";
		//++lodIndex2;
		if (!Geometry::Export3DPointCloudToPLY(result2, outputFileName))
		{
			std::cout << "SingleThreadSoftMaxUniformStrategy: Failed to export sampled point data." << "\n";
			return;
		}

		const auto meshingStrategy = IMB::GetReconstructionStrategy(IMB::ReconstructionFunctionType::BallPivoting);
		Geometry::BaseMeshGeometryData mesh;
		mesh.Vertices = result2;
		meshingStrategy->Process(mesh.Vertices, mesh.PolyIndices);

		/*const std::string */outputFileName = dataOutPath + meshName + "IMB2_LOD_BallPivoting" + std::to_string(lodIndex2) + ".vtk";
		++lodIndex2;
		if (!Geometry::ExportBaseMeshGeometryDataToVTK(mesh, outputFileName))
		{
			std::cerr << "SingleThreadSoftMaxUniformStrategy: Failed to export mesh data." << "\n";
			return;
		}
	};
	const auto progressTracker = std::make_unique<IMB::IncrementalProgressTracker>(
		vertexSamplingStrategy->GetVertexCountEstimate(), nSamples,
		vertexSamplingStrategy->GetVertexCap(), vertexSamplingStrategy->GetMinVertexCount(), 
		lodStepCallback, [](){});
	const auto progressTracker2 = std::make_unique<IMB::IncrementalProgressTracker>(
		vertexSamplingStrategy2->GetVertexCountEstimate(), nSamples,
		vertexSamplingStrategy2->GetVertexCap(), vertexSamplingStrategy2->GetMinVertexCount(),
		lodStepCallback2, []() {});

	const char* verticesStart = fileHandler->GetMemoryStart();
	const char* verticesEnd = fileHandler->GetMemoryEnd();
	for (unsigned int i = 0; i < nSamples; ++i)
	{
		vertexSamplingStrategy2->Sample(verticesStart, verticesEnd, result2, seed, *progressTracker2);
		vertexSamplingStrategy->Sample(verticesStart, verticesEnd, result, seed, *progressTracker);
	}
}

void SingleThreadPoissonDiscSamplingStrategy()
{
	//const auto meshName = "sampledIco";
	const auto meshName = "Apollon_ArtecEva";
	const auto fileName = dataDirPath + meshName + ".ply";
	const auto fileMapping = std::make_unique<Utils::FileMappingWrapper>(fileName);
	if (!fileMapping->IsValid())
	{
		std::cerr << "SingleThreadPoissonDiscSamplingStrategy: Error during initialization of FileMappingWrapper!\n";
		return;
	}
	const char* fileStart = fileMapping->GetFileMemory();
	const char* fileEnd = fileStart + fileMapping->GetFileSize();
	const auto fileHandler = IMB::CreateMeshFileHandler(IMB::FileHandlerParams{ fileName, fileStart, fileEnd });
	if (!fileHandler)
	{
		std::cerr << "SingleThreadPoissonDiscSamplingStrategy: Error during initialization of IncrementalMeshFileHandler!\n";
		return;
	}

	constexpr unsigned int seed = 4999;
	constexpr unsigned int nSamples = 6;
	constexpr unsigned int maxVertexCount = 5000;
	const auto vertexSamplingStrategy = IMB::GetVertexSelectionStrategy(
		IMB::VertexSelectionType::PoissonDisc,
		nSamples, maxVertexCount, fileHandler);
	// to compare against
	const auto vertexSamplingStrategy2 = IMB::GetVertexSelectionStrategy(
		IMB::VertexSelectionType::UniformRandom,
		nSamples, maxVertexCount, fileHandler);
	if (!vertexSamplingStrategy || !vertexSamplingStrategy2)
	{
		std::cerr << "SingleThreadPoissonDiscSamplingStrategy: Error during initialization of VertexSamplingStrategy!\n";
		return;
	}

	IMB::PoissonDiscSamplingParams params;
	params.MinDiscRadius = 5.0; // 1.25; 2.5;
	params.MaxDiscRadius = 7.0; // 3.5;
	params.ExpectedPtsMultiplier = 50.0;

	dynamic_cast<IMB::PoissonDiscVertexSamplingStrategy*>(vertexSamplingStrategy.get())->Params = params;

	std::vector<pmp::Point> result;
	std::vector<pmp::Point> result2;
	unsigned int lodIndex = 0;
	unsigned int lodIndex2 = 0;
	const auto lodStepCallback = [&] {
		/*const */std::string outputFileName = dataOutPath + meshName + "IMB_LOD" + std::to_string(lodIndex) + ".ply";
		//++lodIndex;
		if (!Geometry::Export3DPointCloudToPLY(result, outputFileName))
		{
			std::cout << "SingleThreadPoissonDiscSamplingStrategy: Failed to export sampled point data." << "\n";
			return;
		}

		const auto meshingStrategy = IMB::GetReconstructionStrategy(IMB::ReconstructionFunctionType::BallPivoting);
		Geometry::BaseMeshGeometryData mesh;
		mesh.Vertices = result;
		meshingStrategy->Process(mesh.Vertices, mesh.PolyIndices);

		/*const std::string*/ outputFileName = dataOutPath + meshName + "IMB_LOD_LSW" + std::to_string(lodIndex) + ".vtk";
		++lodIndex;
		if (!Geometry::ExportBaseMeshGeometryDataToVTK(mesh, outputFileName))
		{
			std::cerr << "SingleThreadPoissonDiscSamplingStrategy: Failed to export mesh data." << "\n";
			return;
		}
	};
	const auto lodStepCallback2 = [&] {
		/*const */std::string outputFileName = dataOutPath + meshName + "IMB2_LOD" + std::to_string(lodIndex2) + ".ply";
		//++lodIndex2;
		if (!Geometry::Export3DPointCloudToPLY(result2, outputFileName))
		{
			std::cout << "SingleThreadPoissonDiscSamplingStrategy: Failed to export sampled point data." << "\n";
			return;
		}

		const auto meshingStrategy = IMB::GetReconstructionStrategy(IMB::ReconstructionFunctionType::BallPivoting);
		Geometry::BaseMeshGeometryData mesh;
		mesh.Vertices = result2;
		meshingStrategy->Process(mesh.Vertices, mesh.PolyIndices);

		/*const std::string */outputFileName = dataOutPath + meshName + "IMB2_LOD_BallPivoting" + std::to_string(lodIndex2) + ".vtk";
		++lodIndex2;
		if (!Geometry::ExportBaseMeshGeometryDataToVTK(mesh, outputFileName))
		{
			std::cerr << "SingleThreadPoissonDiscSamplingStrategy: Failed to export mesh data." << "\n";
			return;
		}
	};
	const auto progressTracker = std::make_unique<IMB::IncrementalProgressTracker>(
		vertexSamplingStrategy->GetVertexCountEstimate(), nSamples,
		vertexSamplingStrategy->GetVertexCap(), vertexSamplingStrategy->GetMinVertexCount(),
		lodStepCallback, []() {});
	const auto progressTracker2 = std::make_unique<IMB::IncrementalProgressTracker>(
		vertexSamplingStrategy2->GetVertexCountEstimate(), nSamples,
		vertexSamplingStrategy2->GetVertexCap(), vertexSamplingStrategy2->GetMinVertexCount(),
		lodStepCallback2, []() {});

	const char* verticesStart = fileHandler->GetMemoryStart();
	const char* verticesEnd = fileHandler->GetMemoryEnd();
	for (unsigned int i = 0; i < nSamples; ++i)
	{
		vertexSamplingStrategy2->Sample(verticesStart, verticesEnd, result2, seed, *progressTracker2);
		vertexSamplingStrategy->Sample(verticesStart, verticesEnd, result, seed, *progressTracker);
	}
}
