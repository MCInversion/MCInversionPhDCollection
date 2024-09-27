// ------------------------------- SubdivisionSSGG --------------------------------
//     Fall of 2023, SSGG 2023
// ...............................................................................
/*
 * The origin of this output dates back to CESCG 2022 (Spring of 2022. See [Cavarga, 2022]).
 * The original evolution model had to be stabilized for general use by scaling evolving
 * meshes. The scale factor is computed from the ratio of the time step size and average
 * co-volume measure (area) which is estimated by assuming the distribution of mesh vertices
 * uniform on an icosphere. For this we need to know how many vertices there are. Of course,
 * we can initialize any mesh with N_V vertices and know directly, but my curiosity has
 * lead me to parametrize this for an icosphere constructed via 4-to-1 subdivision
 * using recurrence equations. After discovering the subdivision counting formula,
 * the local SSGG conference in September of 2023 was a suitable venue for this output
 * (see [Cavarga, 2023]).
 *
 * [Cavarga, 2022]
 * Cavarga, M. (2022). Advection-Driven Shrink-Wrapping of Triangulated Surfaces.
 * Proceedings of the 26th CESCG 2022. 95-104.
 *
 * [Cavarga, 2023]
 * Cavarga, M. (2023). 
 * Mesh primitive counting formulas for subdivision surfaces.
 * Proceedings of CSCG 2023. 9, 67-76.
 */
 // --------------------------------------------------------------------------------

#include "pmp/algorithms/Subdivision.h"
#include "pmp/algorithms/Decimation.h"
#include "pmp/algorithms/Remeshing.h"

#include "geometry/IcoSphereBuilder.h"
#include "geometry/MeshAnalysis.h"
#include "geometry/TorusBuilder.h"

#include "IOEnvironment.h"
#include "../Experiments.h"

void SubdivisionTests1()
{
	Geometry::IcoSphereBuilder ico({ 0 });
	ico.BuildBaseData();
	ico.BuildPMPSurfaceMesh();
	auto icoMesh = ico.GetPMPSurfaceMeshResult();
	for (unsigned int i = 0; i < 4; i++)
		icoMesh.delete_face(pmp::Face(i));
	icoMesh.garbage_collection();

	auto nEdges = icoMesh.n_edges();
	auto nVerts = icoMesh.n_vertices();
	auto nFaces = icoMesh.n_faces();
	int euler = nVerts - nEdges + nFaces;
	size_t nBdEdgesTheoretical = std::max(0, 2 - euler);
	auto nBdEdges0 = icoMesh.n_boundary_edges();
	std::cout << "s = 0, nBoundaryEdges = " << nBdEdges0 << ", nBdEdgesTheoretical = " << nBdEdgesTheoretical << "\n";

	pmp::Subdivision subdiv(icoMesh);

	for (unsigned int s = 1; s < 6; s++)
	{
		const auto theoreticalCount = (N_ICO_EDGES_0 * static_cast<unsigned int>(pow(4, s) - 1) + 3 * N_ICO_VERTS_0) / 3;
		subdiv.loop();
		const auto actualCount = icoMesh.n_vertices();
		std::cout << "s = " << s << ", theoreticalCount = " << theoreticalCount << ", actualCount = " << actualCount << "\n";

		nEdges = icoMesh.n_edges();
		nVerts = icoMesh.n_vertices();
		nFaces = icoMesh.n_faces();
		euler = nVerts - nEdges + nFaces;
		nBdEdgesTheoretical = std::max(0, 2 - euler);
		nBdEdges0 = icoMesh.n_boundary_edges();
		std::cout << "s = " << s << ", nBoundaryEdges = " << nBdEdges0 << ", nBdEdgesTheoretical = " << nBdEdgesTheoretical << "\n";

		icoMesh.write(dataOutPath + "ico_Loop" + std::to_string(s) + ".vtk"); /**/
	}
}

void SubdivisionTests2()
{
	constexpr int targetDecimPercentage = 50;
	constexpr int normalDeviation = 180;
	constexpr int aspectRatio = 10;

	Geometry::IcoSphereBuilder ico({ 3 });
	ico.BuildBaseData();
	ico.BuildPMPSurfaceMesh();
	auto icoMesh = ico.GetPMPSurfaceMeshResult();

	const pmp::mat4 matrixGeomScale{
		2.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
	};
	icoMesh *= matrixGeomScale;

	pmp::Decimation decim(icoMesh);
	decim.initialize(aspectRatio, 0, 0, normalDeviation, 0.0f);
	decim.decimate(icoMesh.n_vertices() * 0.01 * targetDecimPercentage);

	pmp::Remeshing remeshing(icoMesh);
	remeshing.uniform_remeshing(0.2f, 3);
	icoMesh.write(dataOutPath + "ico_Decimated0.vtk");

	// icoMesh is now an elongated decimated ellipsoid
	const size_t nVerts0 = icoMesh.n_vertices();
	const size_t nEdges0 = icoMesh.n_edges();

	pmp::Subdivision subdiv(icoMesh);

	for (unsigned int s = 1; s < 6; s++)
	{
		const auto theoreticalCount = (nEdges0 * static_cast<unsigned int>(pow(4, s) - 1) + 3 * nVerts0) / 3;
		subdiv.loop();
		const auto actualCount = icoMesh.n_vertices();
		std::cout << "s = " << s << ", theoreticalCount = " << theoreticalCount << ", actualCount = " << actualCount << "\n";

		icoMesh.write(dataOutPath + "ico_Decimated" + std::to_string(s) + ".vtk"); /**/
	}
}

void SubdivisionTests3()
{
	constexpr Geometry::TorusSettings tSettings{
		1.0f,
			0.4f,
			5,
			3,
			false
	};
	Geometry::TorusBuilder tb(tSettings);
	tb.BuildBaseData();
	tb.BuildPMPSurfaceMesh();
	auto tMesh = tb.GetPMPSurfaceMeshResult();

	tMesh.write(dataOutPath + "torus0.vtk");

	const size_t nVerts0 = tMesh.n_vertices();
	const size_t nEdges0 = tMesh.n_edges();

	pmp::Subdivision subdiv(tMesh);

	for (unsigned int s = 1; s < 6; s++)
	{
		const auto theoreticalCount = (nEdges0 * static_cast<unsigned int>(pow(4, s) - 1) + 3 * nVerts0) / 3;
		subdiv.loop();
		const auto actualCount = tMesh.n_vertices();
		std::cout << "s = " << s << ", theoreticalCount = " << theoreticalCount << ", actualCount = " << actualCount << "\n";

		tMesh.write(dataOutPath + "torus" + std::to_string(s) + ".vtk"); /**/
	}
}

void SubdivisionTest4()
{
	pmp::SurfaceMesh mesh;
	mesh.read(dataOutPath + "bunnyToSubdiv.obj");

	pmp::Subdivision subdiv(mesh);
	subdiv.loop();

	mesh.write(dataOutPath + "bunnySubdiv.vtk");
}

void SubdivTestsBoundary()
{
	std::cout << "performSubdivTestsBoundary...\n";
	Geometry::IcoSphereBuilder ico({ 1 });
	ico.BuildBaseData();
	ico.BuildPMPSurfaceMesh();
	auto icoMesh = ico.GetPMPSurfaceMeshResult();

	constexpr bool deleteSomeFaces = true;

	if (deleteSomeFaces)
	{
		std::vector<unsigned int> facesToDeleteIds{
			0, 1, 3, 10, 11
		};
		for (const auto i : facesToDeleteIds)
			icoMesh.delete_face(pmp::Face(i));
		icoMesh.garbage_collection();
	}

	icoMesh.write(dataOutPath + "icoMeshDeleteFaces0.obj");

	constexpr size_t maxSubdivLevel = 6;

	// estimate edge & vertex counts
	const auto [edgeCounts, vertCounts] = Geometry::GetEdgeVertCountsTheoreticalEstimate(icoMesh, maxSubdivLevel, true);

	pmp::Subdivision subdiv(icoMesh);

	for (size_t s = 1; s < maxSubdivLevel; s++)
	{
		subdiv.loop();
		const auto nEdges = icoMesh.n_edges();
		const auto nVerts = icoMesh.n_vertices();
		std::cout << "========= Edge Count (" << s << "): ==========\n";
		std::cout << "Actual: " << nEdges << ", Theoretical: " << edgeCounts[s] << ".\n";
		std::cout << "========= Vertex Count (" << s << "): ==========\n";
		std::cout << "Actual: " << nVerts << ", Theoretical: " << vertCounts[s] << ".\n";
		std::cout << "------------------------------------------------\n";

		icoMesh.write(dataOutPath + "icoMeshDeleteFaces" + std::to_string(s) + ".obj");
	}
}

void SubdivTestsMultiTorus()
{
	std::cout << "performSubdivTestsTorus...\n";

	for (size_t g = 1; g < 6; g++)
	{
		std::cout << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n";
		std::cout << "Genus : " << g << "\n";
		pmp::SurfaceMesh mesh;
		mesh.read(dataDirPath + std::to_string(g) + "Torus_Simple.obj");
		mesh.write(dataOutPath + std::to_string(g) + "Torus_Subdiv0.vtk");

		constexpr size_t maxSubdivLevel = 6;

		// estimate edge & vertex counts
		const auto [edgeCounts, vertCounts] = Geometry::GetEdgeVertCountsTheoreticalEstimate(mesh, maxSubdivLevel, true);

		pmp::Subdivision subdiv(mesh);

		for (size_t s = 1; s < maxSubdivLevel; s++)
		{
			subdiv.loop();
			const auto nEdges = mesh.n_edges();
			const auto nVerts = mesh.n_vertices();
			std::cout << "========= Edge Count (" << s << "): ==========\n";
			std::cout << "Actual: " << nEdges << ", Theoretical: " << edgeCounts[s] << ".\n";
			std::cout << "========= Vertex Count (" << s << "): ==========\n";
			std::cout << "Actual: " << nVerts << ", Theoretical: " << vertCounts[s] << ".\n";
			std::cout << "------------------------------------------------\n";

			mesh.write(dataOutPath + std::to_string(g) + "Torus_Subdiv" + std::to_string(s) + ".vtk");
		}
		std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
	}
}

void SubdivPreallocationTests()
{
	std::cout << " ... Preallocation Loop Subdivision Tests ..... \n";

	const std::vector<std::string> subdivMeshNames{
		/* 1 */ "armadillo_Simple",
		/* 2 */ "blub_Simple",
		/* 3 */ "bunny_Simple",
		/* 4 */ "maxPlanck_Simple",
		/* 5 */ "3holes",
		/* 6 */ "rockerArm_Simple"
	};

	constexpr size_t maxSubdivLevel = 6;

	for (const auto& meshName : subdivMeshNames)
	{
		std::cout << "meshName: " << meshName << "\n";
		// Load mesh
		pmp::SurfaceMesh mesh;
		mesh.read(dataDirPath + meshName + ".obj");

		double simpleTiming = 0.0;
		double preallocTiming = 0.0;
		constexpr size_t nTimings = 10;

		for (size_t i = 0; i < nTimings; i++)
		{
			std::cout << "timing " << i << "\n";
			// =================================================
			// ......... Plain Subdivision .....................

			auto meshForSubdiv0 = mesh;

			const auto startSimpleSubdiv = std::chrono::high_resolution_clock::now();
			pmp::Subdivision subdivSimple(meshForSubdiv0);

			for (size_t s = 1; s < maxSubdivLevel; s++)
			{
				subdivSimple.loop();
			}
			const auto endSimpleSubdiv = std::chrono::high_resolution_clock::now();
			const std::chrono::duration<double> timeDiffSimpleSubdiv = endSimpleSubdiv - startSimpleSubdiv;
			simpleTiming += timeDiffSimpleSubdiv.count();

			// export result for verification
			//meshForSubdiv0.write(dataOutPath + meshName + "_simpleSubdiv" + std::to_string(maxSubdivLevel - 1) + "timesResult.vtk");

			// =================================================
			// ......... Preallocated Subdivision .....................

			auto meshForSubdiv1 = mesh;

			const auto startPreallocSubdiv = std::chrono::high_resolution_clock::now();
			pmp::Subdivision subdivPrealloc(meshForSubdiv1);
			subdivPrealloc.loop_prealloc(maxSubdivLevel - 1);
			const auto endPreallocSubdiv = std::chrono::high_resolution_clock::now();
			const std::chrono::duration<double> timeDiffPreallocSubdiv = endPreallocSubdiv - startPreallocSubdiv;
			preallocTiming += timeDiffPreallocSubdiv.count();

			// export result for verification
			//meshForSubdiv1.write(dataOutPath + meshName + "_preallocSubdiv" + std::to_string(maxSubdivLevel - 1) + "timesResult.vtk");
		}

		simpleTiming /= nTimings;
		preallocTiming /= nTimings;

		// Report
		std::cout << "Simple Subdiv: " << simpleTiming << " s, Prealloc Subdiv: " << preallocTiming << " s\n";
	}
}

void NewIcosphereTests()
{
	std::cout << "performNewIcosphereTests...\n";
	Geometry::IcoSphereBuilder ico({ 5, 1.0f, true, false });
	ico.BuildBaseData();

	// test out BaseMeshGeometryData.
	const auto bSuccess = ExportBaseMeshGeometryDataToOBJ(ico.GetBaseResult(), dataOutPath + "icoPreallocatedBase.obj");
	assert(bSuccess);

	ico.BuildPMPSurfaceMesh();
	auto icoMesh = ico.GetPMPSurfaceMeshResult();

	icoMesh.write(dataOutPath + "icoPreallocated.obj");
}

void IcospherePerformanceTests()
{
	std::cout << "performIcosphereSpeedTests...\n";
	constexpr size_t maxSubdivLevel = 7;
	constexpr size_t nSphereRuns = 10;

	double simpleTiming = 0.0;
	double preallocTiming = 0.0;
	constexpr size_t nTimings = 10;

	for (size_t s = 1; s < maxSubdivLevel; s++)
	{
		std::cout << "s = " << s << ":\n";
		for (size_t j = 0; j < nTimings; j++)
		{
			//std::cout << "timing " << i << "\n";
			// =================================================
			// ......... Plain Subdivision .....................

			const unsigned int subdiv = s;

			const auto startSimpleSubdiv = std::chrono::high_resolution_clock::now();

			for (size_t i = 0; i < nSphereRuns; i++)
			{
				Geometry::IcoSphereBuilder ico0({ subdiv, 1.0f, true, true });
				ico0.BuildBaseData();
			}

			const auto endSimpleSubdiv = std::chrono::high_resolution_clock::now();
			const std::chrono::duration<double> timeDiffSimpleSubdiv = endSimpleSubdiv - startSimpleSubdiv;
			simpleTiming += timeDiffSimpleSubdiv.count();

			// =================================================
			// ......... Preallocated Subdivision .....................

			const auto startPreallocSubdiv = std::chrono::high_resolution_clock::now();

			for (size_t i = 0; i < nSphereRuns; i++)
			{
				Geometry::IcoSphereBuilder ico1({ subdiv, 1.0f, true, false });
				ico1.BuildBaseData();
			}

			const auto endPreallocSubdiv = std::chrono::high_resolution_clock::now();
			const std::chrono::duration<double> timeDiffPreallocSubdiv = endPreallocSubdiv - startPreallocSubdiv;
			preallocTiming += timeDiffPreallocSubdiv.count();
		}

		simpleTiming /= nTimings;
		preallocTiming /= nTimings;

		// Report
		std::cout << "Simple Icosphere Subdiv: " << simpleTiming << " s, Preallocated Icosphere Subdiv: " << preallocTiming << " s\n";
	}
}

void CatmullClarkCounting()
{
	// Load mesh
	pmp::SurfaceMesh mesh;
	mesh.read(dataDirPath + "CubeSphere.obj");

	constexpr size_t maxSubdivLevel = 6;
	pmp::Subdivision subdiv(mesh);

	for (size_t s = 1; s < maxSubdivLevel; s++)
	{
		subdiv.catmull_clark();
		const auto nEdges = mesh.n_edges();
		const auto nVerts = mesh.n_vertices();
		std::cout << "========= Edge Count (" << s << "): ==========\n";
		std::cout << "Actual: " << nEdges << ".\n";
		std::cout << "========= Vertex Count (" << s << "): ==========\n";
		std::cout << "Actual: " << nVerts << ".\n";
		std::cout << "------------------------------------------------\n";

		mesh.write(dataOutPath + "CubeSphereCC" + std::to_string(s) + ".obj");
	}
}