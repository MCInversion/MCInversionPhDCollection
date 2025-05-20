// ---------------------------- ProjectOfDissertation ----------------------------
//     Spring of 2023
// ...............................................................................
/*
 * Some additional tests for verification and evaluation of some experimental ideas
 * written in my "project of dissertation" proposal. If you dive into these experiments,
 * you're about to see some results with voxelization and parallel mesh file sampling
 */
 // --------------------------------------------------------------------------------

#include "pmp/algorithms/Remeshing.h"

#include "geometry/GridUtil.h"
#include "geometry/MeshAnalysis.h"
#include "geometry/MobiusStripBuilder.h"

#include "utils/TimingUtils.h"

#include "core/ConversionUtils.h"

#include "IOEnvironment.h"
#include "../Experiments.h"

// NOTE: This one was used to override export during different stages of uniform remeshing for illustration purposes
void RemeshingTests()
{
	pmp::SurfaceMesh mesh;
	mesh.read(dataDirPath + "bunny_no_holes2.obj");

	pmp::Scalar meanEdgeLength = 0.0;
	for (const auto& e : mesh.edges())
	{
		meanEdgeLength += mesh.edge_length(e);
	}
	meanEdgeLength /= mesh.n_edges();

	/*constexpr int targetDecimPercentage = 10;
	constexpr int normalDeviation = 180;
	constexpr int aspectRatio = 10;
	pmp::Decimation decim(mesh);
	decim.initialize(aspectRatio, 0.0, 0.0, normalDeviation, 0.0);
	decim.decimate(mesh.n_vertices() * 0.01 * targetDecimPercentage);*/


	pmp::Remeshing remeshing(mesh);
	remeshing.uniform_remeshing(8.5, 1, true);
}

void MobiusStripVoxelization()
{
	constexpr Geometry::MobiusStripSettings mSettings{
		    1.0,
			1.0,
			40,
			10,
			false,
			true
	};
	Geometry::MobiusStripBuilder mb(mSettings);
	mb.BuildBaseData();
	mb.BuildPMPSurfaceMesh();
	auto mMesh = mb.GetPMPSurfaceMeshResult();

	//mMesh.write(dataOutPath + "mobius.vtk");
	mMesh.write(dataOutPath + "mobius.obj");
	auto bbox = mMesh.bounds();
	const auto bboxSize = bbox.max() - bbox.min();
	bbox.expand(0.1 * bboxSize[0], 0.1 * bboxSize[1], bboxSize[2]);
	Geometry::ScalarGrid grid(0.02, bbox);
	Geometry::ComputeInteriorExteriorSignFromMeshNormals(grid, mMesh);

	ExportToVTI(dataOutPath + "MobiusSignVals", grid);
}

void MetaballTest()
{
	constexpr double initVal = 0.0;
	// grid containing both balls
	//Geometry::ScalarGrid grid(0.05, pmp::BoundingBox{ pmp::vec3{}, pmp::vec3{10.0, 10.0, 10.0} }, initVal);

	// grid containing a clipped voxel field of the balls
	/**/Geometry::ScalarGrid grid(1.0, pmp::BoundingBox{
		pmp::vec3{2.1, 3.0, 1.6},
			pmp::vec3{7.3, 8.3, 6.2} }, initVal);

	const Geometry::ScalarGridBoolOpFunction opFnc = Geometry::SimpleUnion;

	// apply balls
	const Geometry::MetaBallParams mBall1Params = {
		pmp::vec3{3.0, 4.0, 4.0}, 4.0, opFnc
	};
	ApplyMetaBallToGrid(grid, mBall1Params);
	const Geometry::MetaBallParams mBall2Params = {
		pmp::vec3{4.0, 5.0, 4.0}, 5.0, opFnc
	};
	ApplyMetaBallToGrid(grid, mBall2Params);

	ExportToVTI(dataOutPath + "MetaBallVals", grid);

	/*constexpr double isoLevel = 0.1;
	const auto mcMesh = IlatsikMC::GetMarchingCubesMesh<double>(
		grid.Values().data(),
		grid.Dimensions().Nx, grid.Dimensions().Ny, grid.Dimensions().Nz,
		isoLevel);
	auto mcPMPMesh = Geometry::ConvertMCMeshToPMPSurfaceMesh(mcMesh);

	pmp::Remeshing remeshing(mcPMPMesh);
	remeshing.uniform_remeshing(1.5, 10, false);

	mcPMPMesh.write(dataOutPath + "MetaBallMC.vtk");*/
}

void ImportedObjMetricsEval()
{
	const std::vector<std::string> importedMeshNames{
		//"ArmadilloSWBlender_NearestSurfPt",
		//"ArmadilloSWBlender_ProjectNeg"
		"bunnyDanielLSW150"
	};

	for (const auto& meshName : importedMeshNames)
	{
		//try
		//{
		std::cout << "MetricsEval: " << meshName << "...\n";
		pmp::SurfaceMesh mesh;
		mesh.read(dataDirPath + meshName + ".obj");

		if (!Geometry::ComputeEquilateralTriangleJacobianConditionNumbers(mesh))
		{
			std::cout << "Error!\n";
			continue;
		}

		mesh.write(dataOutPath + meshName + "_Metric.vtk");
		//}
		//catch(...)
		//{
		//	std::cerr << "> > > > > > MetricsEval subroutine has thrown an exception! Continue... < < < < < \n";
		//}

	}
}

void MMapImportTest()
{
	const std::vector<std::string> importedMeshNames{
		"nefertiti"
	};

	for (const auto& meshName : importedMeshNames)
	{
		constexpr size_t nRuns = 10;

		// load parallel
		pmp::SurfaceMesh parImportedMesh;
		std::optional<Geometry::BaseMeshGeometryData> baseDataOpt;

		AVERAGE_TIMING(parImported, nRuns, {
			baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", true);
			if (!baseDataOpt.has_value())
			{
				std::cerr << "baseDataOpt == nullopt!\n";
				break;
			}
			}, true);

		// verify by export
		parImportedMesh = ConvertBufferGeomToPMPSurfaceMesh(baseDataOpt.value());
		parImportedMesh.write(dataOutPath + meshName + "_parallelImp.obj");

		// load single-threaded			
		pmp::SurfaceMesh stImportedMesh;
		std::optional<Geometry::BaseMeshGeometryData> stBaseDataOpt;

		AVERAGE_TIMING(singleThreadImported, nRuns, {
			stBaseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", false);
			if (!stBaseDataOpt.has_value())
			{
				std::cerr << "stBaseDataOpt == nullopt!\n";
				break;
			}
			}, true);

		// verify by export
		stImportedMesh = ConvertBufferGeomToPMPSurfaceMesh(stBaseDataOpt.value());
		stImportedMesh.write(dataOutPath + meshName + "_stImp.obj");

	}
}

void MMapOBJChunkMarkingTest()
{
	const std::vector<std::string> importedMeshNames{
		//"armadillo", /* ! non-manifold !? */
		//"BentChair",
		"blub",
		"bunny",
		"maxPlanck",
		"nefertiti",
		"ogre",
		"spot",
		"3holes", // messed up, also mutex issue?
		"fertility",
		//"happyBuddha", /* ! non-manifold !? */
		"rockerArm" // messed up, also mutex issue?
	};

	for (const auto& meshName : importedMeshNames)
	{
		std::cout << "Parallel loading mesh: " << meshName << ".obj ... ";

		// load parallel
		pmp::SurfaceMesh parImportedMesh;
		std::optional<Geometry::BaseMeshGeometryData> baseDataOpt;

		std::vector<pmp::Scalar> threadIds;
		baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", true, &threadIds);
		if (!baseDataOpt.has_value())
		{
			std::cerr << "baseDataOpt == nullopt!\n";
			break;
		}
		std::cout << "done.\nExporting ... ";

		// verify by export
		parImportedMesh = ConvertBufferGeomToPMPSurfaceMesh(baseDataOpt.value());
		if (parImportedMesh.n_vertices() != threadIds.size())
		{
			std::cerr << "parImportedMesh.n_vertices() != threadIds.size()!\n";
			break;
		}
		auto vThreadIdProp = parImportedMesh.add_vertex_property<pmp::Scalar>("v:threadId");
		for (const auto& v : parImportedMesh.vertices())
		{
			vThreadIdProp[v] = threadIds[v.idx()];
		}
		parImportedMesh.write(dataOutPath + meshName + "_parallelImp.vtk");
		std::cout << "done.\n";
	}
}

void SimpleBunnyOBJSamplingDemo()
{
	const auto baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + "bunny.obj", true);
	assert(baseDataOpt.has_value());
	const auto& baseData = baseDataOpt.value();

	constexpr size_t nSamplings = 10;
	constexpr size_t minVerts = 9; // Minimum number of vertices to sample
	const size_t maxVerts = baseData.Vertices.size(); // Maximum number of vertices available

	for (size_t i = 0; i < nSamplings; ++i) {
		// Determine the number of vertices to sample for this iteration
		size_t nVerts = minVerts + (maxVerts - minVerts) * i / (nSamplings - 1);

		// Ensure nVerts is within the valid range
		nVerts = std::max(minVerts, std::min(nVerts, maxVerts));
		std::cout << "Sampling " << nVerts << " vertices in iteration " << i << "\n";

		// Export sampled vertices to PLY
		std::string filename = dataOutPath + "bunnyPts_" + std::to_string(i) + ".ply";
		const auto bSuccess = ExportSampledVerticesToPLY(baseData, nVerts, filename);
		assert(bSuccess);
	}
}