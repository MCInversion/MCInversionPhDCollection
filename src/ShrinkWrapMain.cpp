
#include "pmp/SurfaceMesh.h"

#include "sdf/SDF.h"
#include "geometry/GridUtil.h"

#include "SurfaceEvolver.h"
#include "ConversionUtils.h"

#include <filesystem>
#include <chrono>
#include <map>

#include "BrainSurfaceEvolver.h"
#include "IsosurfaceEvolver.h"
#include "SphereTest.h"
#include "geometry/IcoSphereBuilder.h"
#include "geometry/MarchingCubes.h"
#include "geometry/TorusBuilder.h"
#include "geometry/MobiusStripBuilder.h"
#include "pmp/algorithms/Decimation.h"
#include "pmp/algorithms/Normals.h"
#include "pmp/algorithms/Remeshing.h"
#include "pmp/algorithms/Subdivision.h"
//#include "geometry/MeshAnalysis.h"

// set up root directory
const std::filesystem::path fsRootPath = DROOT_DIR;
const auto fsDataDirPath = fsRootPath / "data\\";
const auto fsDataOutPath = fsRootPath / "output\\";
const std::string dataDirPath = fsDataDirPath.string();
const std::string dataOutPath = fsDataOutPath.string();

constexpr bool performSDFTests = false;
constexpr bool performSphereTest = false;
constexpr bool performEvolverTests = false;
constexpr bool performIsosurfaceEvolverTests = true;
// constexpr bool performNiftiTests = true; // TODO: nifti import not supported yet
constexpr bool performBrainEvolverTests = false;
constexpr bool performSubdivisionTests1 = false;
constexpr bool performSubdivisionTests2 = false;
constexpr bool performSubdivisionTests3 = false;
constexpr bool performSubdivisionTest4 = false;
constexpr bool performRemeshingTests = false;
constexpr bool performMobiusStripVoxelization = false;
constexpr bool performMetaballTest = false;

[[nodiscard]] size_t CountBoundaryEdges(const pmp::SurfaceMesh& mesh)
{
	size_t result = 0;
	for (const auto e : mesh.edges())
	{
		if (!mesh.is_boundary(e))
			continue;

		result++;
	}
	return result;
}

int main()
{
    // DISCLAIMER: the names need to match the models in "DROOT_DIR/data" except for the extension (which is always *.obj)
    const std::vector<std::string> meshNames{
        "armadillo",
        //"BentChair",
    	//"blub",
    	//"bunny",
        //"maxPlanck",
        //"nefertiti",
        //"ogre",
        //"spot"
    };

	if (performSDFTests)
	{
	    constexpr unsigned int nVoxelsPerMinDimension = 10;
		constexpr bool computeGradients = false;

	    for (const auto& name : meshNames)
	    {
		    pmp::SurfaceMesh mesh;
		    mesh.read(dataDirPath + name + ".obj");

			std::cout << "I break here!\n";

	        const auto meshBBox = mesh.bounds();
	        const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
	        const float minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
			const float cellSize = minSize / nVoxelsPerMinDimension;
		    const SDF::DistanceFieldSettings sdfSettings{
		        cellSize,
		        1.0f,
				DBL_MAX,
		        SDF::KDTreeSplitType::Center,
		        SDF::SignComputation::VoxelFloodFill,
		        SDF::BlurPostprocessingType::None,
		        SDF::PreprocessingType::Octree
		    };
			SDF::ReportInput(mesh, sdfSettings, std::cout);

		    const auto startSDF = std::chrono::high_resolution_clock::now();
		    const auto sdf = SDF::DistanceFieldGenerator::Generate(mesh, sdfSettings);
		    const auto endSDF = std::chrono::high_resolution_clock::now();

			SDF::ReportOutput(sdf, std::cout);
		    const std::chrono::duration<double> timeDiff = endSDF - startSDF;
		    std::cout << "SDF Time: " << timeDiff.count() << " s\n";
			std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
		    ExportToVTI(dataOutPath + name + "SDF", sdf);
			if (computeGradients)
			{
				std::cout << "Geometry::ComputeGradient(sdf) ...";
				const auto gradSdf = Geometry::ComputeGradient(sdf);
				ExportToVTK(dataOutPath + name + "gradSDF", gradSdf);
				std::cout << "... done\n";
				std::cout << "Geometry::ComputeNormalizedGradient(sdf) ...";
				const auto normGradSdf = Geometry::ComputeNormalizedGradient(sdf);
				ExportToVTK(dataOutPath + name + "normGradSDF", normGradSdf);
				std::cout << "... done\n";
				std::cout << "Geometry::ComputeNormalizedNegativeGradient(sdf) ...";
				const auto negNormGradSdf = Geometry::ComputeNormalizedNegativeGradient(sdf);
				ExportToVTK(dataOutPath + name + "negNormGradSDF", negNormGradSdf);
				std::cout << "... done\n";			
			}

			/*std::cout << "---------------------------------------------------\n";
			std::cout << "SDF - Angle Weighted Pseudonormal Approach:\n";
			std::cout << "---------------------------------------------------\n";

			pmp::BoundingBox sdfBox(meshBBox);
			const float expansion = 1.0f * minSize;
			sdfBox.expand(expansion, expansion, expansion);
			Geometry::ScalarGrid sdf2(cellSize, sdfBox);

			const auto startSDF2 = std::chrono::high_resolution_clock::now();
			pmp::Normals::compute_vertex_normals(mesh);
			ComputeMeshSignedDistanceFromNormals(sdf2, mesh);
			const auto endSDF2 = std::chrono::high_resolution_clock::now();
			SDF::ReportOutput(sdf2, std::cout);
			const std::chrono::duration<double> timeDiff2 = endSDF2 - startSDF2;
			std::cout << "SDF (Pseudonormal) Time: " << timeDiff2.count() << " s\n";
			std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
			ExportToVTI(dataOutPath + name + "SDF2", sdf2);*/
	    }
	} // endif performSDFTests

	if (performSphereTest)
	{
		/* { // Setup 1: No remeshing, No tangential redistribution
			const ST_MeshTopologySettings topoSettings{};

			const SphereTestEvolutionSettings stSettings{
				topoSettings,
				false, false, dataOutPath,
				ST_MeshLaplacian::Barycentric, 0.0f, false};

			SphereTest st(stSettings);
			st.PerformTest(4);			
		}

		std::cout << "=====================================================\n";

		{ // Setup 2: No remeshing, tangential redistribution with weight 0.25
			const ST_MeshTopologySettings topoSettings{};

			const SphereTestEvolutionSettings stSettings{
				topoSettings,
				false, false, dataOutPath,
				ST_MeshLaplacian::Barycentric, 0.25f, false };

			SphereTest st(stSettings);
			st.PerformTest(4);
		}*/

		std::cout << "=====================================================\n";

		{ // Setup 3: Remeshing, No tangential redistribution
			const ST_MeshTopologySettings topoSettings{};

			const SphereTestEvolutionSettings stSettings{
				topoSettings,
				false, false, dataOutPath,
				ST_MeshLaplacian::Barycentric, 0.0f, true };

			SphereTest st(stSettings);
			st.PerformTest(4);
		}

		std::cout << "=====================================================\n";

		{ // Setup 4: Remeshing, tangential redistribution with weight 0.25
			const ST_MeshTopologySettings topoSettings{};

			const SphereTestEvolutionSettings stSettings{
				topoSettings,
				false, false, dataOutPath,
				ST_MeshLaplacian::Barycentric, 0.25f, true };

			SphereTest st(stSettings);
			st.PerformTest(4);
		}
		
	} // endif performSphereTest

	if (performEvolverTests)
	{
		constexpr unsigned int nVoxelsPerMinDimension = 40;
        constexpr double defaultTimeStep = 0.05;
		const std::map<std::string, double> timeStepSizesForMeshes{
			{"armadillo", 0.05 },
			{"BentChair", 0.05 },
			{"blub", 0.05 },
			{"bunny", 0.0025 },
			{"maxPlanck", 0.05 },
			{"nefertiti", 0.05 },
			{"ogre", 0.05 },
			{"spot", 0.05 }
		};

		for (const auto& name : meshNames)
		{
			pmp::SurfaceMesh mesh;
			mesh.read(dataDirPath + name + ".obj");
			const auto meshBBox = mesh.bounds();
			const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
			const float minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
			const float maxSize = std::max({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
			const float cellSize = minSize / nVoxelsPerMinDimension;
			constexpr float volExpansionFactor = 1.0f;
			const SDF::DistanceFieldSettings sdfSettings{
				cellSize,
				volExpansionFactor,
				Geometry::DEFAULT_SCALAR_GRID_INIT_VAL, // 0.2, TODO: zero gradient values lead to slow MCF outside of the truncated SDF region
				SDF::KDTreeSplitType::Center,
				SDF::SignComputation::VoxelFloodFill,
				SDF::BlurPostprocessingType::None,
				SDF::PreprocessingType::Octree
			};
			SDF::ReportInput(mesh, sdfSettings, std::cout);

			const auto startSDF = std::chrono::high_resolution_clock::now();
			const auto sdf = SDF::DistanceFieldGenerator::Generate(mesh, sdfSettings);
			const auto endSDF = std::chrono::high_resolution_clock::now();

			SDF::ReportOutput(sdf, std::cout);
			const std::chrono::duration<double> timeDiff = endSDF - startSDF;
			std::cout << "SDF Time: " << timeDiff.count() << " s\n";
			std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
			ExportToVTI(dataOutPath + name + "SDF", sdf);

			const auto& sdfBox = sdf.Box();
			const auto sdfBoxSize = sdfBox.max() - sdfBox.min();
			const auto sdfBoxMaxDim = std::max<double>({ sdfBoxSize[0], sdfBoxSize[1], sdfBoxSize[2] });

			const double fieldIsoLevel = sqrt(3.0) / 2.0 * static_cast<double>(cellSize);
                        
			const double tau = (timeStepSizesForMeshes.contains(name) ? timeStepSizesForMeshes.at(name) : defaultTimeStep); // time step
			SurfaceEvolutionSettings seSettings{
				name,
				80,
				tau,
				fieldIsoLevel,
				3,
				PreComputeAdvectionDiffusionParams(0.5 * sdfBoxMaxDim, minSize),
				{},
				minSize, maxSize,
				meshBBox.center(),
				true, false,
				dataOutPath,
				MeshLaplacian::Voronoi,
				{"minAngle", "maxAngle", "jacobianConditionNumber", "equilateralJacobianCondition",/* "stiffnessMatrixConditioning" */},
				0.05f,
				true
			};
			ReportInput(seSettings, std::cout);
			SurfaceEvolver evolver(sdf, volExpansionFactor, seSettings);

			try
			{
				evolver.Evolve();
			}
			catch(...)
			{
				std::cerr << "> > > > > > > > > > > > > > SurfaceEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	} // endif performEvolverTests

	if (performIsosurfaceEvolverTests)
	{
		const std::vector<std::string> higherGenusMeshNames{
			"3holes",
			"fertility",
			"happyBuddha",
			"rockerArm"
		};

		constexpr unsigned int nVoxelsPerMinDimension = 40;
		constexpr double defaultTimeStep = 0.05;
		const std::map<std::string, double> timeStepSizesForMeshes{
			{"3holes", defaultTimeStep },
			{"fertility", defaultTimeStep },
			{"happyBuddha", defaultTimeStep },
			{"rockerArm", defaultTimeStep }
		};
		const std::map<std::string, double> effectiveIsolevelsForMeshes{
			{"3holes", 0.02 },
			{"fertility", 4.0 },
			{"happyBuddha", 1.5e-3 },
			{"rockerArm", 0.06 }
		};
		const std::map<std::string, float> resamplingFactors{
			{"3holes", 3.0f },
			{"fertility", 2.0f },
			{"happyBuddha", 1.0f },
			{"rockerArm", 2.0f }
		};

		for (const auto& name : higherGenusMeshNames)
		{
			pmp::SurfaceMesh mesh;
			mesh.read(dataDirPath + name + ".obj");
			const auto meshBBox = mesh.bounds();
			const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
			const float minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
			const float maxSize = std::max({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
			const float cellSize = minSize / nVoxelsPerMinDimension;
			constexpr float volExpansionFactor = 1.0f;
			const SDF::DistanceFieldSettings sdfSettings{
				cellSize,
				volExpansionFactor,
				//0.2, // TODO: will this truncation be OK?
				Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
				SDF::KDTreeSplitType::Center,
				SDF::SignComputation::VoxelFloodFill,
				SDF::BlurPostprocessingType::None,
				SDF::PreprocessingType::Octree
			};
			SDF::ReportInput(mesh, sdfSettings, std::cout);

			const auto startSDF = std::chrono::high_resolution_clock::now();
			const auto sdf = SDF::DistanceFieldGenerator::Generate(mesh, sdfSettings);
			const auto endSDF = std::chrono::high_resolution_clock::now();

			SDF::ReportOutput(sdf, std::cout);
			const std::chrono::duration<double> timeDiff = endSDF - startSDF;
			std::cout << "SDF Time: " << timeDiff.count() << " s\n";
			std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
			ExportToVTI(dataOutPath + name + "SDF", sdf);

			//constexpr double fieldIsoLevel = 0.0;
			const double fieldIsoLevel = sqrt(3.0) / 2.0 * static_cast<double>(cellSize);
			const double isoLevel = (effectiveIsolevelsForMeshes.contains(name) ? 
				(fieldIsoLevel < effectiveIsolevelsForMeshes.at(name) ? effectiveIsolevelsForMeshes.at(name) : 1.1 * fieldIsoLevel) : 5.0);

			const MeshTopologySettings topoParams{
				0.4f,
				0.0,
				1.0f,
				0.0,
				2,
				0.0,
				3,
				5,
				false,
				FeatureDetectionType::MeanCurvature,
				1.0 * M_PI_2 * 180.0, 2.0 * M_PI_2 * 180.0,
				2.0f,
				0.8f * static_cast<float>(M_PI_2),
				true
			};

			const double tau = (timeStepSizesForMeshes.contains(name) ? timeStepSizesForMeshes.at(name) : defaultTimeStep); // time step
			const float resamplingFactor = (resamplingFactors.contains(name) ? resamplingFactors.at(name) : 1.5f);
			IsoSurfaceEvolutionSettings seSettings{
				name,
				20,
				tau,
				fieldIsoLevel,
				isoLevel,
				cellSize * resamplingFactor,
				PreComputeAdvectionDiffusionParams(2.0, minSize),
				topoParams,
				true, false,
				dataOutPath,
				MeshLaplacian::Voronoi,
				{"minAngle", "maxAngle", "jacobianConditionNumber", "equilateralJacobianCondition",/* "stiffnessMatrixConditioning" */},
				0.05f,
				true
			};
			ReportInput(seSettings, std::cout);
			IsoSurfaceEvolver evolver(sdf, volExpansionFactor, seSettings);

			try
			{
				evolver.Evolve();
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > SurfaceEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	} // endif performIsosurfaceEvolverTests

	// DISCLAIMER: the names need to match the models in "DROOT_DIR/data" except for the extension (which is always *.vti)
	const std::vector<std::string> brainNames{
		"talairach",
		//"actual_brain" // TODO: use git lfs to upload ascii vti file of size > 50 MB or implement nifti import
	};

	/*if (performNiftiTests)
	{
		// TODO: nifti import not supported yet
	} // endif performNiftiTests
	*/


	if (performBrainEvolverTests)
	{
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// NOTE: these values are copy-pasted from bet2 under the same image data inputs.
		//       the true evaluation of threshold settings as well as radius and center are
		//       nearly impossible to reverse-engineer from bet2.
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// TODO: reverse-engineer bet2 threshold, and ico-sphere params evaluation

		constexpr BE_ThresholdSettings talairachThresholdSettings{
			7, 3, // min and max intensity depths for normal search
			0.0, // 2nd percentile threshold
			1024.33496, // 98th percentile threshold
			102.43349, // effective threshold
			348.0,// effective median threshold
		};
		constexpr BE_ThresholdSettings actualBrainThresholdSettings{
			7, 3, // min and max intensity depths for normal search
			0.0, // 2nd percentile threshold
			668.25, // 98th percentile threshold
			66.825, // effective threshold
			317.0, // effective median threshold
		};
		const std::map<std::string, BE_ThresholdSettings> bet2ThresholdSettings{
			{"talairach", talairachThresholdSettings},
			{"actual_brain", actualBrainThresholdSettings },
		};

		// ico-sphere params
		const pmp::vec3 talairachCenter{ 69.278477120258884f, 81.210907276033296f, 69.224956401243205f };
		const pmp::vec3 actualBrainCenter{ 96.536378893717142f, 126.13349816811417f, 116.99736547254018f };

		constexpr float talairachRadius = 66.061572538428337f;
		constexpr float actualBrainRadius = 102.09133074846271f;

		const BE_IcoSphereSettings talairachIcoSettings{ talairachCenter, talairachRadius };
		const BE_IcoSphereSettings actualBrainIcoSettings{ actualBrainCenter, actualBrainRadius };
		const std::map<std::string, BE_IcoSphereSettings> bet2IcoSphereSettings{
			{"talairach", talairachIcoSettings},
			{"actual_brain", actualBrainIcoSettings },
		};

		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		const BE_CurvatureSettings cSettings{
			// TODO: try values
		};

		for (const auto& name : brainNames)
		{
			std::string imagePathIn = dataDirPath + name + ".vti";
			const auto gridData = ImportVTI(imagePathIn);
			// ExportToVTI(dataOutPath + name + "_test", gridData);

			const auto& thresholdSettings = bet2ThresholdSettings.at(name);
			const auto& icoSphereSettings = bet2IcoSphereSettings.at(name);
			BrainExtractionSettings beSettings{
				name,
				80,
				0.01,
				3, // ico-sphere subdivision level: bet2 uses 5 by default
				cSettings,
				thresholdSettings,
				icoSphereSettings,
				{},
				true, false,
				dataOutPath,
				BE_MeshLaplacian::Voronoi,
				{"minAngle", "maxAngle", "jacobianConditionNumber",/* "stiffnessMatrixConditioning" */},
				true
			};

			BrainSurfaceEvolver evolver(gridData, beSettings);
			ReportInput(beSettings, std::cout);
			try
			{
				evolver.Evolve();
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > BrainSurfaceEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}


		
	} // endif performBrainEvolverTests

	if (performSubdivisionTests1)
	{
		Geometry::IcoSphereBuilder ico({0});
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
		auto nBdEdges0 = CountBoundaryEdges(icoMesh);
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
			nBdEdges0 = CountBoundaryEdges(icoMesh);
			std::cout << "s = " << s << ", nBoundaryEdges = " << nBdEdges0 << ", nBdEdgesTheoretical = " << nBdEdgesTheoretical << "\n";

			icoMesh.write(dataOutPath + "ico_Loop" + std::to_string(s) + ".vtk"); /**/
		}
	}

	if (performSubdivisionTests2)
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
		decim.initialize(aspectRatio, 0, 0, normalDeviation,0.0f);
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

	if (performSubdivisionTests3)
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

	if (performSubdivisionTest4)
	{
		pmp::SurfaceMesh mesh;
		mesh.read(dataOutPath + "bunnyToSubdiv.obj");

		pmp::Subdivision subdiv(mesh);
		subdiv.loop();

		mesh.write(dataOutPath + "bunnySubdiv.vtk");
	}

	if (performRemeshingTests)
	{
		pmp::SurfaceMesh mesh;
		mesh.read(dataDirPath + "bunny_no_holes2.obj");

		float meanEdgeLength = 0.0f;
		for (const auto& e : mesh.edges())
		{
			meanEdgeLength += mesh.edge_length(e);
		}
		meanEdgeLength /= mesh.n_edges();

		/*constexpr int targetDecimPercentage = 10;
		constexpr int normalDeviation = 180;
		constexpr int aspectRatio = 10;
		pmp::Decimation decim(mesh);
		decim.initialize(aspectRatio, 0.0, 0.0, normalDeviation, 0.0f);
		decim.decimate(mesh.n_vertices() * 0.01 * targetDecimPercentage);*/


		pmp::Remeshing remeshing(mesh);
		remeshing.uniform_remeshing(8.5, 1, true);
	}

	if (performMobiusStripVoxelization)
	{
		constexpr Geometry::MobiusStripSettings mSettings{
			1.0f,
			1.0f,
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
		bbox.expand(0.1f * bboxSize[0], 0.1f * bboxSize[1], bboxSize[2]);
		Geometry::ScalarGrid grid(0.02f, bbox);
		Geometry::ComputeInteriorExteriorSignFromMeshNormals(grid, mMesh);

		ExportToVTI(dataOutPath + "MobiusSignVals", grid);
	}

	if (performMetaballTest)
	{
		constexpr double initVal = 0.0;
		// grid containing both balls
		//Geometry::ScalarGrid grid(0.05f, pmp::BoundingBox{ pmp::vec3{}, pmp::vec3{10.0f, 10.0f, 10.0f} }, initVal);

		// grid containing a clipped voxel field of the balls
		/**/Geometry::ScalarGrid grid(1.0f, pmp::BoundingBox{
			pmp::vec3{2.1f, 3.0f, 1.6f},
			pmp::vec3{7.3f, 8.3f, 6.2f} }, initVal);

		const Geometry::ScalarGridBoolOpFunction opFnc = Geometry::SimpleUnion;

		// apply balls
		const Geometry::MetaBallParams mBall1Params = {
			pmp::vec3{3.0f, 4.0f, 4.0f}, 4.0f, opFnc
		};
		ApplyMetaBallToGrid(grid, mBall1Params);
		const Geometry::MetaBallParams mBall2Params = {
			pmp::vec3{4.0f, 5.0f, 4.0f}, 5.0f, opFnc
		};
		ApplyMetaBallToGrid(grid, mBall2Params);

		ExportToVTI(dataOutPath + "MetaBallVals", grid);

		/*constexpr double isoLevel = 0.1;
		const auto mcMesh = GetMarchingCubesMesh<double>(
			grid.Values().data(),
			grid.Dimensions().Nx, grid.Dimensions().Ny, grid.Dimensions().Nz,
			isoLevel);
		auto mcPMPMesh = Geometry::ConvertMCMeshToPMPSurfaceMesh(mcMesh);

		pmp::Remeshing remeshing(mcPMPMesh);
		remeshing.uniform_remeshing(1.5, 10, false);

		mcPMPMesh.write(dataOutPath + "MetaBallMC.vtk");*/
	}
}
