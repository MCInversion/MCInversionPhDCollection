// ------------------------- IncrementalMeshBuilderCESCG ----------------------------
// Spring of 2024. CESCG 2024 doctoral colloquium.
// ..................................................................................
/*
 * This was a chance to revisit an idea from the project of dissertation proposal about
 * progressively sampling very large mesh files. This lead to the development of the
 * IMB (Incremental Mesh Builder) framework. The results of these experiments were
 * presented on CESCG 2024 doctoral colloquium.
 */
 // --------------------------------------------------------------------------------

#include "utils/TimingUtils.h"

#include "geometry/GeometryConversionUtils.h"
#include "geometry/MeshAnalysis.h"

#include "sdf/SDF.h"

#include "core/ConversionUtils.h"
#include "core/IncrementalMeshBuilder.h"
#include "core/MeasurementUtils.h"

#include "IOEnvironment.h"
#include "../Experiments.h"

void BPATest()
{
	const std::vector<std::string> meshForPtCloudNames{
		"armadillo",
		"blub",
		"bunny",
		"maxPlanck",
		"nefertiti",
		"ogre",
		"spot"
	};

	constexpr size_t samplingLevel = 3;
	constexpr size_t nSamplings = 10;
	constexpr size_t minVerts = 9; // Minimum number of vertices to sample

	constexpr unsigned int seed = 5000; // seed for the pt cloud sampling RNG
	constexpr pmp::Scalar radiusPercentageOfMinDim = 0.05;

	for (const auto& meshName : meshForPtCloudNames)
	{
		std::cout << "==================================================================\n";
		std::cout << "Mesh To Pt Cloud: " << meshName << ".obj -> " << meshName << "Pts_" << samplingLevel << ".ply\n";
		std::cout << "------------------------------------------------------------------\n";
		const auto baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", false);
		if (!baseDataOpt.has_value())
		{
			std::cerr << "baseDataOpt == nullopt!\n";
			break;
		}
		std::cout << meshName << ".obj" << " imported as BaseMeshGeometryData.\n";
		const auto& baseData = baseDataOpt.value();
		const size_t maxVerts = baseData.Vertices.size(); // Maximum number of vertices available
		size_t nVerts = minVerts + (maxVerts - minVerts) * samplingLevel / (nSamplings - 1);
		nVerts = std::max(minVerts, std::min(nVerts, maxVerts));

		std::cout << "Sampling " << nVerts << "/" << maxVerts << " vertices...\n";

		// Export sampled vertices to PLY
		std::string filename = dataOutPath + meshName + "Pts_" + std::to_string(samplingLevel) + ".ply";
		if (!ExportSampledVerticesToPLY(baseData, nVerts, filename, seed))
		{
			std::cerr << "ExportSampledVerticesToPLY failed!\n";
			break;
		}

		const auto ptCloudName = meshName + "Pts_" + std::to_string(samplingLevel);
		const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataOutPath + ptCloudName + ".ply", true);
		if (!ptCloudOpt.has_value())
		{
			std::cerr << "ptCloudOpt == nullopt!\n";
			break;
		}

		const auto& ptCloud = ptCloudOpt.value();

		const pmp::BoundingBox ptCloudBBox(ptCloud);
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar bpaRadius = radiusPercentageOfMinDim * minSize;
		std::cout << "BPA radius: " << bpaRadius << " units, i.e.: " << radiusPercentageOfMinDim * 100 << " % of min dim.\n";

		const auto bpaMeshOpt = Geometry::ComputeBallPivotingMeshFromPoints(ptCloud, bpaRadius);
		if (!bpaMeshOpt.has_value())
		{
			std::cerr << "bpaMeshOpt == nullopt!\n";
			break;
		}

		const auto& bpaMesh = bpaMeshOpt.value();

		const auto nUnvisitedVerts = CountUnreferencedVertices(bpaMesh);
		std::cout << "We detected " << nUnvisitedVerts << "/" << bpaMesh.Vertices.size() << " unvisited vertices!\n";

		if (!Geometry::ExportBaseMeshGeometryDataToOBJ(bpaMesh, dataOutPath + meshName + "Pts_" + std::to_string(samplingLevel) + "BPAResult.obj"))
		{
			std::cerr << "ExportBaseMeshGeometryDataToOBJ failed!\n";
			break;
		}
	}
}

void IncrementalMeshBuilderTests()
{
	// *.ply format:
	const std::vector<std::string> meshForPtCloudNames{
		//"Apollon_ArtecEva",
		//"armadillo",
		//"bunny",
		//"CaesarBust",
		//"maxPlanck",
		"nefertiti"
	};

	for (const auto& meshName : meshForPtCloudNames)
	{
		constexpr size_t nUpdates = 10;
		unsigned int lodIndex = 0;

		//const IMB::MeshRenderFunction exportToOBJ = [&lodIndex, &meshName](const Geometry::BaseMeshGeometryData& meshData) {
		//	const std::string outputFileName = dataOutPath + "IncrementalMeshBuilder_" + meshName + "/" + meshName + "_IMB_LOD" + std::to_string(lodIndex) + ".obj";
		//	if (!Geometry::ExportBaseMeshGeometryDataToOBJ(meshData, outputFileName))
		//	{
		//		std::cout << "Failed to export mesh data." << "\n";
		//		return;
		//	}
		//	std::cout << "Mesh data exported successfully to " << outputFileName << "\n";
		//	++lodIndex;
		//};
		//const IMB::MeshRenderFunction exportPtsToPLY = [&lodIndex, &meshName](const Geometry::BaseMeshGeometryData& meshData) {
		//	const std::string outputFileName = dataOutPath + "IncrementalMeshBuilder_" + meshName + "/" + meshName + "_IMB_LOD" + std::to_string(lodIndex) + ".ply";
		//	if (!Geometry::ExportPointsToPLY(meshData, outputFileName))
		//	{
		//		std::cout << "Failed to export mesh data." << "\n";
		//		return;
		//	}
		//	std::cout << "Mesh data exported successfully to " << outputFileName << "\n";
		//	++lodIndex;
		//};
		
		std::vector timeTicks = { std::chrono::high_resolution_clock::now() };
		const IMB::MeshRenderFunction exportToVTK = [&lodIndex, &meshName, &timeTicks](const Geometry::BaseMeshGeometryData& meshData) {
			timeTicks.push_back(std::chrono::high_resolution_clock::now());
			const std::string outputFileName = dataOutPath + "IncrementalMeshBuilder_" + meshName + "/" + meshName + "_IMB_LOD" + std::to_string(lodIndex) + ".vtk";
			if (!Geometry::ExportBaseMeshGeometryDataToVTK(meshData, outputFileName))
			{
				std::cout << "Failed to export mesh data." << "\n";
				return;
			}
			std::cout << "Mesh data with " << meshData.Vertices.size() << " vertices exported successfully to " << outputFileName << "\n";

			if (!ExportTimeVectorInSeconds(timeTicks, dataOutPath + "IncrementalMeshBuilder_" + meshName + "/" + meshName + "_Timings.txt"))
				throw std::logic_error("File not exported!");
			++lodIndex;
		};
		auto& meshBuilder = IMB::IncrementalMeshBuilder::GetInstance();
		meshBuilder.Init(
			dataDirPath + meshName + ".ply",
			nUpdates,
			IMB::ReconstructionFunctionType::LagrangianShrinkWrapping,
			//IMB::ReconstructionFunctionType::BallPivoting,
			//IMB::ReconstructionFunctionType::None,
			IMB::VertexSelectionType::UniformRandom,
			//IMB::VertexSelectionType::Sequential,
			//exportPtsToPLY,
			//exportToOBJ,
			exportToVTK,
			40000
		);
		constexpr unsigned int seed = 4999;
		constexpr unsigned int nThreads = 4;
		meshBuilder.DispatchAndSyncWorkers(seed, nThreads);
	}
}

void TwoGBApollonMeshBuilderTest()
{
	// WARNING: This thing is huge
	const std::string inputFileName = "C:/Users/Martin/source/testMeshes/Apollon/Apollon_50MPx_el1-2-3-4-5-6-7_parcial_FS_158201601_111M_scaled.ply";
	const std::string meshName = "Apollon_111M";

	constexpr size_t nUpdates = 10;
	unsigned int lodIndex = 0;
	//const IMB::MeshRenderFunction exportToOBJ = [&lodIndex, &meshName](const Geometry::BaseMeshGeometryData& meshData) {
	//	const std::string outputFileName = dataOutPath + "IncrementalMeshBuilder_" + meshName + "/" + meshName + "_IMB_LOD" + std::to_string(lodIndex) + ".obj";
	//	if (!Geometry::ExportBaseMeshGeometryDataToOBJ(meshData, outputFileName))
	//	{
	//		std::cout << "Failed to export mesh data." << "\n";
	//		return;
	//	}
	//	std::cout << "Mesh data exported successfully to " << outputFileName << "\n";
	//	++lodIndex;
	//};
	//const IMB::MeshRenderFunction exportPtsToPLY = [&lodIndex, &meshName](const Geometry::BaseMeshGeometryData& meshData) {
	//	const std::string outputFileName = dataOutPath + "IncrementalMeshBuilder_" + meshName + "/" + meshName + "_IMB_LOD" + std::to_string(lodIndex) + ".ply";
	//	if (!Geometry::ExportPointsToPLY(meshData, outputFileName))
	//	{
	//		std::cout << "Failed to export mesh data." << "\n";
	//		return;
	//	}
	//	std::cout << "Mesh data exported successfully to " << outputFileName << "\n";
	//	++lodIndex;
	//};
	std::vector<std::pair<std::chrono::high_resolution_clock::time_point, size_t>> apollonTimeTicks = { {std::chrono::high_resolution_clock::now(), 0} };
	const IMB::MeshRenderFunction exportToVTK = [&lodIndex, &meshName, &apollonTimeTicks](const Geometry::BaseMeshGeometryData& meshData) {
		apollonTimeTicks.push_back({ std::chrono::high_resolution_clock::now(), meshData.Vertices.size() });
		const std::string outputFileName = dataOutPath + "IncrementalMeshBuilder_" + meshName + "/" + meshName + "_IMB_LOD" + std::to_string(lodIndex) + ".vtk";
		if (!Geometry::ExportBaseMeshGeometryDataToVTK(meshData, outputFileName))
		{
			std::cout << "Failed to export mesh data." << "\n";
			return;
		}
		std::cout << "Mesh data exported successfully to " << outputFileName << "\n";
		if (lodIndex > 5)
			if (!ExportTimeVectorWithPointCountsInSeconds(apollonTimeTicks, dataOutPath + "IncrementalMeshBuilder_" + meshName + "/Apollon111M_Timings.txt"))
				throw std::logic_error("File not exported!");
		++lodIndex;
	};
	auto& meshBuilder = IMB::IncrementalMeshBuilder::GetInstance();
	meshBuilder.Init(
		inputFileName,
		nUpdates,
		//IMB::ReconstructionFunctionType::LagrangianShrinkWrapping,
		//IMB::ReconstructionFunctionType::BallPivoting, 
		//IMB::ReconstructionFunctionType::None,
		IMB::ReconstructionFunctionType::Poisson,
		IMB::VertexSelectionType::UniformRandom,
		//IMB::VertexSelectionType::Sequential,
		//exportPtsToPLY,
		//exportToOBJ,
		exportToVTK,
		40000
	);
	constexpr unsigned int seed = 4999;
	constexpr unsigned int nThreads = 4;
	meshBuilder.DispatchAndSyncWorkers(seed, nThreads);
}

void NanoflannDistanceTests()
{
	const std::vector<std::string> importedPtCloudNames{
		//"bunnyPts_3"
		"bunnyPts_Minimal"
	};

	for (const auto& ptCloudName : importedPtCloudNames)
	{
		std::cout << "----------------------------------------------------------------------\n";
		std::cout << "Procedure: " << ptCloudName << " Nanoflann inter-vertex evaluation ...\n";
		std::cout << "----------------------------------------------------------------------\n";
		const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataDirPath + ptCloudName + ".ply", true);
		if (!ptCloudOpt.has_value())
		{
			std::cerr << "ptCloudOpt == nullopt!\n";
			break;
		}
		const auto& ptCloud = ptCloudOpt.value();

		const auto minDist = Geometry::ComputeMinInterVertexDistance(ptCloud);
		std::cout << "minDist = " << minDist << "\n";
		const auto minDistBrute = Geometry::ComputeMinInterVertexDistanceBruteForce(ptCloud);
		std::cout << "minDistBrute = " << minDistBrute << "\n";
		std::cout << "..................................................................\n";

		const auto maxDistBrute = Geometry::ComputeMaxInterVertexDistanceBruteForce(ptCloud);
		std::cout << "maxDistBrute = " << maxDistBrute << "\n";
		std::cout << "..................................................................\n";

		const auto meanDist = Geometry::ComputeNearestNeighborMeanInterVertexDistance(ptCloud);
		std::cout << "meanDist = " << meanDist << "\n";
		const auto meanDistBrute = Geometry::ComputeMeanInterVertexDistanceBruteForce(ptCloud);
		std::cout << "meanDistBrute = " << meanDistBrute << "\n";
		std::cout << "..................................................................\n";
	}
}

void ApollonLSWSaliencyEval()
{
	const std::vector<std::string> importedMeshNames{
		"Apollon_ArtecEva_IMB"
	};

	constexpr size_t nLOD = 8;

	for (const auto& meshName : importedMeshNames)
	{
		for (size_t i = 0; i <= nLOD; ++i)
		{
			const auto meshNameFull = "/IncrementalMeshBuilder_Apollon_ArtecEva/" + meshName + "_LOD" + std::to_string(i);
			std::cout << "SaliencyEval: " << meshNameFull << "...\n";
			pmp::SurfaceMesh mesh;
			mesh.read(dataOutPath + meshNameFull + ".vtk");

			if (!Geometry::EvaluatePMPSurfaceMeshSaliency(mesh))
			{
				std::cout << "Error!\n";
				continue;
			}

			mesh.write(dataOutPath + meshNameFull + "_Saliency.vtk");
		}
	}

	const auto origMeshName = "Apollon_ArtecEva";
	std::cout << "SaliencyEval: " << origMeshName << ".ply ...\n";
	pmp::SurfaceMesh mesh;
	mesh.read(dataDirPath + origMeshName + ".ply");

	if (Geometry::EvaluatePMPSurfaceMeshSaliency(mesh, 10.5149))
	{
		mesh.write(dataOutPath + origMeshName + "_Saliency.vtk");
	}
}

void IncrementalMeshBuilderHausdorffEval()
{
	const std::vector<std::string> imbMeshNames{
		"armadillo",
		"bunny",
		"CaesarBust",
		//"maxPlanck",
		"nefertiti",
		//"Apollon_ArtecEva"
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;

	for (const auto& meshName : imbMeshNames)
	{
		// Load the original mesh
		pmp::SurfaceMesh origMesh;
		origMesh.read(dataDirPath + meshName + ".ply");

		// Create a mesh adapter for the original mesh
		Geometry::BaseMeshGeometryData meshData = Geometry::ConvertPMPSurfaceMeshToBaseMeshGeometryData(origMesh);
		Geometry::BaseMeshAdapter meshAdapter(std::make_shared<Geometry::BaseMeshGeometryData>(meshData));

		// Compute the distance field for the original mesh
		auto meshBBox = meshAdapter.GetBounds();
		auto meshBBoxSize = meshBBox.max() - meshBBox.min();
		pmp::Scalar meshMinSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		pmp::Scalar meshCellSize = meshMinSize / static_cast<pmp::Scalar>(nVoxelsPerMinDimension);

		SDF::DistanceFieldSettings meshDfSettings{
			meshCellSize,
			1.0,  // volExpansionFactor
			Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::KDTreeSplitType::Center,
			SDF::SignComputation::None,  // Unsigned distance field
			SDF::BlurPostprocessingType::None,
			SDF::PreprocessingType::Octree
		};

		auto meshDf = SDF::DistanceFieldGenerator::Generate(meshAdapter, meshDfSettings);

		// Process each LOD for this mesh
		int lodIndex = 0;
		while (true) {
			std::string filePath = dataOutPath + "IncrementalMeshBuilder_" + meshName + "/" + meshName + "_IMB_LOD" + std::to_string(lodIndex) + ".vtk";
			if (!std::filesystem::exists(filePath)) {
				std::cout << "No more LOD files found for " << meshName << " at LOD " << lodIndex << std::endl;
				break;
			}

			// Import LOD mesh
			try {
				auto lodMeshData = ImportVTK(filePath);
				Geometry::BaseMeshAdapter lodMeshAdapter(std::make_shared<Geometry::BaseMeshGeometryData>(lodMeshData));

				// Compute the Hausdorff distance between the original mesh and the LOD mesh
				auto hausdorffDistance = Geometry::ComputeMeshToMeshHausdorffDistance(lodMeshAdapter.GetBaseMesh(), meshAdapter.GetBaseMesh(), meshDf, nVoxelsPerMinDimension);
				std::cout << "Hausdorff Distance for " << meshName << " LOD " << lodIndex << ": " << (hausdorffDistance ? std::to_string(*hausdorffDistance) : "Failed") << std::endl;
			}
			catch (const std::exception& e) {
				std::cerr << "Error processing " << filePath << ": " << e.what() << '\n';
			}

			++lodIndex;  // Next LOD
		}
	}
}

void ApollonArtecEvaLSWHausdorffEval()
{
	const std::vector<std::string> imbMeshNames{
		"Apollon_ArtecEva"
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;

	for (const auto& meshName : imbMeshNames)
	{
		// Load the original mesh
		pmp::SurfaceMesh origMesh;
		origMesh.read(dataDirPath + meshName + ".ply");

		// Create a mesh adapter for the original mesh
		Geometry::BaseMeshGeometryData meshData = Geometry::ConvertPMPSurfaceMeshToBaseMeshGeometryData(origMesh);
		Geometry::BaseMeshAdapter meshAdapter(std::make_shared<Geometry::BaseMeshGeometryData>(meshData));

		// Compute the distance field for the original mesh
		auto meshBBox = meshAdapter.GetBounds();
		auto meshBBoxSize = meshBBox.max() - meshBBox.min();
		pmp::Scalar meshMinSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		pmp::Scalar meshCellSize = meshMinSize / static_cast<pmp::Scalar>(nVoxelsPerMinDimension);

		SDF::DistanceFieldSettings meshDfSettings{
			meshCellSize,
			1.0,  // volExpansionFactor
			Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::KDTreeSplitType::Center,
			SDF::SignComputation::None,  // Unsigned distance field
			SDF::BlurPostprocessingType::None,
			SDF::PreprocessingType::Octree
		};

		auto meshDf = SDF::DistanceFieldGenerator::Generate(meshAdapter, meshDfSettings);

		// Process each LOD for this mesh
		int lodIndex = 0;
		while (true) {
			std::string filePath = dataOutPath + "IncrementalMeshBuilder_" + meshName + "/" + meshName + "_IMB_LOD" + std::to_string(lodIndex) + ".vtk";
			if (!std::filesystem::exists(filePath)) {
				std::cout << "No more LOD files found for " << meshName << " at LOD " << lodIndex << std::endl;
				break;
			}

			// Import LOD mesh
			try {
				auto lodMeshData = ImportVTK(filePath);
				Geometry::BaseMeshAdapter lodMeshAdapter(std::make_shared<Geometry::BaseMeshGeometryData>(lodMeshData));

				// Compute the Hausdorff distance between the original mesh and the LOD mesh
				auto hausdorffDistance = Geometry::ComputeMeshToMeshHausdorffDistance(lodMeshAdapter.GetBaseMesh(), meshAdapter.GetBaseMesh(), meshDf, nVoxelsPerMinDimension);
				std::cout << "Hausdorff Distance for " << meshName << " LOD " << lodIndex << ": " << (hausdorffDistance ? std::to_string(*hausdorffDistance) : "Failed") << std::endl;
			}
			catch (const std::exception& e) {
				std::cerr << "Error processing " << filePath << ": " << e.what() << '\n';
			}

			++lodIndex;  // Next LOD
		}
	}
}

void VertexNormalSampling()
{
	const std::vector<std::string> importedMeshNames{
		"bunny_normals",
		//"CaesarBust"
	};

	constexpr size_t samplingLevel = 2;
	constexpr size_t nSamplings = 10;
	constexpr size_t minVerts = 9; // Minimum number of vertices to sample
	constexpr unsigned int seed = 4999;

	for (const auto& meshName : importedMeshNames)
	{
		std::cout << "==================================================================\n";
		std::cout << "Mesh To Pt Cloud: " << meshName << ".obj -> " << meshName << "Pts_" << samplingLevel << ".ply\n";
		std::cout << "------------------------------------------------------------------\n";
		const auto baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", false);
		if (!baseDataOpt.has_value())
		{
			std::cerr << "baseDataOpt == nullopt!\n";
			break;
		}
		std::cout << "meshName.obj" << " imported as BaseMeshGeometryData.\n";
		const auto& baseData = baseDataOpt.value();
		const size_t maxVerts = baseData.Vertices.size(); // Maximum number of vertices available
		size_t nVerts = minVerts + (maxVerts - minVerts) * samplingLevel / (nSamplings - 1);
		nVerts = std::max(minVerts, std::min(nVerts, maxVerts));

		std::cout << "Sampling " << nVerts << "/" << maxVerts << " vertices...\n";

		// Export sampled vertices with normals to VTK
		//std::string filename = dataOutPath + meshName + "Pts_" + std::to_string(samplingLevel) + ".vtk";
		//if (!ExportSampledVerticesWithNormalsToVTK(baseData, nVerts, filename, seed))
		//{
		//	std::cerr << "ExportSampledVerticesWithNormalsToVTK failed!\n";
		//	break;
		//}
		std::string filename = dataOutPath + meshName + "Pts_" + std::to_string(samplingLevel) + ".ply";
		if (!ExportSampledVerticesWithNormalsToPLY(baseData, nVerts, filename, seed))
		{
			std::cerr << "ExportSampledVerticesWithNormalsToPLY failed!\n";
			break;
		}
	}
}