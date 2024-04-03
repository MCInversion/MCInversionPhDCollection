#include "BrainSurfaceEvolver.h"
#include "ConversionUtils.h"
#include "IsosurfaceEvolver.h"
#include "SurfaceEvolver.h"
#include "SheetMembraneEvolver.h"
#include "SphereTest.h"

#include "geometry/GridUtil.h"
#include "geometry/IcoSphereBuilder.h"
#include "geometry/MarchingCubes.h"
#include "geometry/MeshAnalysis.h"
#include "geometry/MobiusStripBuilder.h"
#include "geometry/PlaneBuilder.h"
#include "geometry/TorusBuilder.h"
#include "geometry/GeometryUtil.h"
#include "sdf/SDF.h"
#include "utils/TimingUtils.h"
#include "utils/StringUtils.h"

#include "pmp/SurfaceMesh.h"
#include "pmp/algorithms/Decimation.h"
#include "pmp/algorithms/Normals.h"
#include "pmp/algorithms/Remeshing.h"
#include "pmp/algorithms/Subdivision.h"

#include <filesystem>
#include <chrono>
#include <map>

// set up root directory
const std::filesystem::path fsRootPath = DROOT_DIR;
const auto fsDataDirPath = fsRootPath / "data\\";
const auto fsDataOutPath = fsRootPath / "output\\";
const std::string dataDirPath = fsDataDirPath.string();
const std::string dataOutPath = fsDataOutPath.string();

// test flags. The actual execution is just wrapped inside an if-statement.
constexpr bool performSDFTests = false;
constexpr bool performSphereTest = false;
constexpr bool performEvolverTests = false;
constexpr bool performOldArmadilloLSWTest = false;
constexpr bool performIsosurfaceEvolverTests = false;
constexpr bool performSheetEvolverTest = false;
// constexpr bool performNiftiTests = true; // TODO: nifti import not supported yet
constexpr bool performBrainEvolverTests = false;
constexpr bool performSubdivisionTests1 = false;
constexpr bool performSubdivisionTests2 = false;
constexpr bool performSubdivisionTests3 = false;
constexpr bool performSubdivisionTest4 = false;
constexpr bool performSubdivTestsBoundary = false;
constexpr bool performSubdivTestsMultiTorus = false;
constexpr bool performSubdivPreallocationTests = false;
constexpr bool performNewIcosphereTests = false;
constexpr bool performIcospherePerformanceTests = false;
constexpr bool pefrormCatmullClarkCounting = false;
constexpr bool performRemeshingTests = false;
constexpr bool performMobiusStripVoxelization = false;
constexpr bool performMetaballTest = false;
constexpr bool performImportedObjMetricsEval = false;
constexpr bool performMMapImportTest = false;
constexpr bool performMMapOBJChunkMarkingTest = false;
constexpr bool performSimpleBunnyOBJSamplingDemo = false;
constexpr bool performPDanielPtCloudPLYExport = false;
constexpr bool performPtCloudToDF = false;
constexpr bool performPDanielPtCloudComparisonTest = false;
constexpr bool performRepulsiveSurfResultEvaluation = false;
constexpr bool performHistogramResultEvaluation = false;
constexpr bool performOldResultJacobianMetricEval = false;
constexpr bool performHausdorffDistanceMeasurementsPerTimeStep = false;
constexpr bool performDirectHigherGenusPtCloudSampling = false;
constexpr bool performHigherGenusPtCloudLSW = false;
constexpr bool performTriTriIntersectionTests = false;
constexpr bool performMeshSelfIntersectionTests = false;
constexpr bool performHurtadoMeshesIsosurfaceEvolverTests = false;
constexpr bool performHurtadoTrexIcosphereLSW = false;
constexpr bool performImportVTIDebugTests = false;
constexpr bool performConvexHullTests = true;

int main()
{
	// DISCLAIMER: the names need to match the models in "DROOT_DIR/data" except for the extension (which is always *.obj)
	const std::vector<std::string> meshNames{
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
				3, // IcoSphereSubdivisionLevel
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
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > SurfaceEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	} // endif performEvolverTests

	if (performOldArmadilloLSWTest)
	{
		constexpr unsigned int nVoxelsPerMinDimension = 40;
		constexpr double defaultTimeStep = 0.05;
		const std::string name = "armadillo";
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

		const double tau = defaultTimeStep; // time step

		MeshTopologySettings topoParams;
		topoParams.EdgeLengthDecayFactor = 0.95f;
		topoParams.ExcludeEdgesWithoutBothFeaturePts = true;
		topoParams.NRemeshingIters = 5;
		topoParams.PrincipalCurvatureFactor = 2.0f;
		topoParams.FeatureType = FeatureDetectionType::MeanCurvature;
		topoParams.UseBackProjection = false;
		AdvectionDiffusionParameters adParams = PreComputeAdvectionDiffusionParams(0.5 * sdfBoxMaxDim, minSize);
		adParams.AdvectionSineMultiplier *= 1.5;
		SurfaceEvolutionSettings seSettings{
			name,
			80,
			tau,
			fieldIsoLevel,
			3,
			adParams,
			topoParams,
			minSize, maxSize,
			meshBBox.center(),
			true, false,
			dataOutPath,
			MeshLaplacian::Voronoi,
			{"minAngle", "maxAngle", "jacobianConditionNumber", "equilateralJacobianCondition",/* "stiffnessMatrixConditioning" */},
			0.05f,
			true,
			//false
		};
		ReportInput(seSettings, std::cout);
		SurfaceEvolver evolver(sdf, volExpansionFactor, seSettings);

		try
		{
			evolver.Evolve();
		}
		catch (...)
		{
			std::cerr << "> > > > > > > > > > > > > > SurfaceEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
		}
	} // endif performOldArmadilloLSWTest

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
				true,
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
				30,
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

	if (performSubdivTestsBoundary)
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

	if (performSubdivTestsMultiTorus)
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

	if (performSubdivPreallocationTests)
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

	if (performNewIcosphereTests)
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

	if (performIcospherePerformanceTests)
	{
		std::cout << "performIcospherePerformanceTests...\n";
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

	if (pefrormCatmullClarkCounting)
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

	if (performSheetEvolverTest)
	{
#define PERFORM_7PT_EXAMPLE false
		/**/
		constexpr float roiHalfDim = 5.0f;
		constexpr float roiDim = 2.0f * roiHalfDim;

#if PERFORM_7PT_EXAMPLE
		constexpr unsigned int nXSegments = 50;
		constexpr unsigned int nYSegments = 50;
#else
		constexpr unsigned int nXSegments = 40;
		constexpr unsigned int nYSegments = 40;
#endif

		constexpr Geometry::PlaneSettings mSettings{
			pmp::vec3{},
			roiDim,
			roiDim,
			nXSegments,
			nYSegments,
			true,
			true
		};
		Geometry::PlaneBuilder pb(mSettings);
		pb.BuildBaseData();
		pb.BuildPMPSurfaceMesh();
		auto pMesh = pb.GetPMPSurfaceMeshResult();

		pMesh.write(dataOutPath + "plane.vtk");
		//pMesh.write(dataOutPath + "plane.obj");

		constexpr float cellSize = 0.1f;
		constexpr float columnWeight = 0.5f;
#if PERFORM_7PT_EXAMPLE
		const auto gridBox = pmp::BoundingBox{ pmp::vec3{-5.0f, -5.0f, -roiHalfDim}, pmp::vec3{16.1f, 15.0f, roiHalfDim} };
		const auto grid = GetDistanceFieldWithSupportColumns(cellSize, gridBox, {
			{pmp::vec2{4.0f, 6.0f}, 0.5f * columnWeight},
			{pmp::vec2{0.0f, 0.0f}, 0.5f * columnWeight},
			{pmp::vec2{5.0f, 0.0f}, 0.5f * columnWeight},
			{pmp::vec2{11.1f, 0.1f}, 0.5f * columnWeight},
			{pmp::vec2{9.0f, 2.0f}, 0.5f * columnWeight},
			{pmp::vec2{7.0f, 2.0f}, 0.5f * columnWeight},
			{pmp::vec2{6.0f, 10.0f}, 0.5f * columnWeight}
			});
#else // 5-point example:
		const auto gridBox = pmp::BoundingBox{ pmp::vec3{0.0f, 0.0f, -roiHalfDim}, pmp::vec3{roiDim, roiDim, roiHalfDim} };
		const auto grid = GetDistanceFieldWithSupportColumns(cellSize, gridBox, {
			{pmp::vec2{2.5f, 2.5f}, 0.5f * columnWeight},
			{pmp::vec2{7.5f, 2.5f}, 0.5f * columnWeight},
			{pmp::vec2{7.5f, 7.5f}, 0.5f * columnWeight},
			{pmp::vec2{5.0f, 8.0f}, 0.5f * columnWeight},
			{pmp::vec2{2.5f, 7.5f}, 0.5f * columnWeight}
			});
#endif
		ExportToVTI(dataOutPath + "CapsuleVals", grid);

		const auto& sdfBox = grid.Box();
		const auto sdfBoxSize = sdfBox.max() - sdfBox.min();

		const double fieldIsoLevel = sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		const float startZHeight = sdfBox.min()[2] + 0.9f * sdfBoxSize[2];
		const float endZHeight = sdfBox.min()[2] + 0.5f * sdfBoxSize[2];

		constexpr double tau = 0.02;

		const MeshTopologySettings topoSettings{
			true,
			0.45f,
			0.0,
			1.0
		};

		const AdvectionDiffusionParameters adParams{
			1.0, 1.0,
			1.0, 0.5
		};

		SheetMembraneEvolutionSettings seSettings{
			"SheetMembrane",
			50,
			tau,
			fieldIsoLevel,
			startZHeight,
			endZHeight,
			nXSegments,
			nYSegments,
			adParams,
			topoSettings,
			true, false,
			dataOutPath,
			MeshLaplacian::Voronoi,
			{/*"minAngle", "maxAngle", "jacobianConditionNumber", "equilateralJacobianCondition", "stiffnessMatrixConditioning" */},
			0.05f,
			true,
			true
		};
		ReportInput(seSettings, std::cout);
		SheetMembraneEvolver evolver(grid, seSettings);

		try
		{
			evolver.Evolve();
		}
		catch (...)
		{
			std::cerr << "> > > > > > > > > > > > > > SurfaceEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
		}
	}

	if (performImportedObjMetricsEval)
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

	if (performMMapImportTest)
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

	if (performMMapOBJChunkMarkingTest)
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

			std::vector<float> threadIds;
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
			auto vThreadIdProp = parImportedMesh.add_vertex_property<float>("v:threadId");
			for (const auto& v : parImportedMesh.vertices())
			{
				vThreadIdProp[v] = threadIds[v.idx()];
			}
			parImportedMesh.write(dataOutPath + meshName + "_parallelImp.vtk");
			std::cout << "done.\n";
		}
	}

	if (performSimpleBunnyOBJSamplingDemo)
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

	if (performPDanielPtCloudPLYExport)
	{
		const std::vector<std::string> importedMeshNames{
			"bunny",
			"CaesarBust"
		};

		constexpr size_t samplingLevel = 3;
		constexpr size_t nSamplings = 10;
		constexpr size_t minVerts = 9; // Minimum number of vertices to sample

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

			// Export sampled vertices to PLY
			std::string filename = dataOutPath + meshName + "Pts_" + std::to_string(samplingLevel) + ".ply";
			if (!ExportSampledVerticesToPLY(baseData, nVerts, filename))
			{
				std::cerr << "ExportSampledVerticesToPLY failed!\n";
				break;
			}
		}
	}

	if (performPtCloudToDF)
	{
		const std::vector<std::string> importedPtCloudNames{
			"bunnyPts_3",
			"CaesarBustPts_3"
		};

		constexpr unsigned int nVoxelsPerMinDimension = 20;
		for (const auto& ptCloudName : importedPtCloudNames)
		{
			// const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataDirPath + ptCloudName + ".ply", true);
			const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataOutPath + ptCloudName + ".ply", true);
			//const auto ptCloudOpt = Geometry::ImportPLYPointCloudDataMainThread(dataOutPath + ptCloudName + ".ply");
			if (!ptCloudOpt.has_value())
			{
				std::cerr << "ptCloudOpt == nullopt!\n";
				break;
			}

			const auto& ptCloud = ptCloudOpt.value();
			const pmp::BoundingBox ptCloudBBox(ptCloud);
			const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
			const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
			const float cellSize = minSize / nVoxelsPerMinDimension;
			const SDF::PointCloudDistanceFieldSettings dfSettings{
						cellSize,
						1.0f,
						DBL_MAX,
						SDF::BlurPostprocessingType::None
			};

			std::cout << "==================================================================\n";
			std::cout << "Pt Cloud to DF: " << ptCloudName << ".ply -> " << ptCloudName << "_DF_" << nVoxelsPerMinDimension << "voxPerMinDim.vti\n";
			std::cout << "------------------------------------------------------------------\n";

			const auto startSDF = std::chrono::high_resolution_clock::now();
			const auto sdf = SDF::PointCloudDistanceFieldGenerator::Generate(ptCloud, dfSettings);
			const auto endSDF = std::chrono::high_resolution_clock::now();

			SDF::ReportOutput(sdf, std::cout);
			const std::chrono::duration<double> timeDiff = endSDF - startSDF;
			std::cout << "DF Time: " << timeDiff.count() << " s\n";
			std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
			ExportToVTI(dataOutPath + ptCloudName + "DF", sdf);
		}
	}

	if (performPDanielPtCloudComparisonTest)
	{
		const std::vector<std::string> importedPtCloudNames{
			"bunnyPts_3"//,
			//"CaesarBustPts_3"
		};
		const std::map<std::string, double> timeStepSizesForPtClouds{
			{"bunnyPts_3", 0.07 },
			{"CaesarBustPts_3", 0.05 }
		};
		const std::map<std::string, double> isoLevelOffsetFactors{
			{"bunnyPts_3", 2.0 },
			{"CaesarBustPts_3", 0.5 }
		};

		SetRemeshingAdjustmentTimeIndices({ 3, 10, 20, 50, 100, 120, 140, 145 });

		constexpr unsigned int nVoxelsPerMinDimension = 40;
		constexpr double defaultTimeStep = 0.05;
		constexpr double defaultOffsetFactor = 1.5;
		for (const auto& ptCloudName : importedPtCloudNames)
		{
			// const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataDirPath + ptCloudName + ".ply", true);
			const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataOutPath + ptCloudName + ".ply", true);
			// const auto ptCloudOpt = Geometry::ImportPLYPointCloudDataMainThread(dataOutPath + ptCloudName + ".ply");
			if (!ptCloudOpt.has_value())
			{
				std::cerr << "ptCloudOpt == nullopt!\n";
				break;
			}

			const auto& ptCloud = ptCloudOpt.value();
			const pmp::BoundingBox ptCloudBBox(ptCloud);
			const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
			const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
			const float maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
			const float cellSize = minSize / nVoxelsPerMinDimension;
			constexpr float volExpansionFactor = 1.0f;
			const SDF::PointCloudDistanceFieldSettings dfSettings{
						cellSize,
						volExpansionFactor,
						Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
						SDF::BlurPostprocessingType::None
			};

			std::cout << "==================================================================\n";
			std::cout << "Pt Cloud to DF: " << ptCloudName << ".ply -> " << ptCloudName << "_DF_" << nVoxelsPerMinDimension << "voxPerMinDim.vti\n";
			std::cout << "------------------------------------------------------------------\n";

			const auto startSDF = std::chrono::high_resolution_clock::now();
			auto sdf = SDF::PointCloudDistanceFieldGenerator::Generate(ptCloud, dfSettings);
			const auto endSDF = std::chrono::high_resolution_clock::now();

			//NormalizeScalarGridValues(sdf);

			SDF::ReportOutput(sdf, std::cout);
			const std::chrono::duration<double> timeDiff = endSDF - startSDF;
			std::cout << "DF Time: " << timeDiff.count() << " s\n";
			std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
			//ExportToVTI(dataOutPath + ptCloudName + "DF", sdf);

			//const auto& sdfBox = sdf.Box();
			//const auto sdfBoxSize = sdfBox.max() - sdfBox.min();
			//const auto sdfBoxMaxDim = std::max<double>({ sdfBoxSize[0], sdfBoxSize[1], sdfBoxSize[2] });

			const double isoLvlOffsetFactor = (timeStepSizesForPtClouds.contains(ptCloudName) ? isoLevelOffsetFactors.at(ptCloudName) : defaultOffsetFactor);
			const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

			const double tau = (timeStepSizesForPtClouds.contains(ptCloudName) ? timeStepSizesForPtClouds.at(ptCloudName) : defaultTimeStep); // time step

			MeshTopologySettings topoParams;
			topoParams.FixSelfIntersections = true;
			topoParams.MinEdgeMultiplier = 0.14f;
			topoParams.UseBackProjection = false;
			topoParams.PrincipalCurvatureFactor = 3.2f;
			topoParams.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);
			topoParams.EdgeLengthDecayFactor = 0.7f;
			topoParams.ExcludeEdgesWithoutBothFeaturePts = true;
			topoParams.FeatureType = FeatureDetectionType::MeanCurvature;

			AdvectionDiffusionParameters adParams{
				1.0, 1.0,
				2.0, 1.0
			};

			SurfaceEvolutionSettings seSettings{
				ptCloudName,
				146, // 150,
				tau,
				fieldIsoLevel,
				2, // IcoSphereSubdivisionLevel
				adParams,
				topoParams,
				minSize, maxSize,
				ptCloudBBox.center(),
				true, false,
				dataOutPath,
				MeshLaplacian::Voronoi,
				{"minAngle", "maxAngle", "jacobianConditionNumber", "equilateralJacobianCondition",/* "stiffnessMatrixConditioning" */},
				0.05f,
				true,
				false
			};
			ReportInput(seSettings, std::cout);
			SurfaceEvolver evolver(sdf, volExpansionFactor, seSettings);

			try
			{
				evolver.Evolve();
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > SurfaceEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	}

	if (performRepulsiveSurfResultEvaluation)
	{
		const std::vector<std::string> importedMeshNames{
			//"spot_RepulsiveResult220",
			"bunny_RepulsiveResult220",
			"bunnyLSW150_Obstacle",
			"bunnyLSW150_FullWrap",
			"bunnyLSW200_Daniel30k"
		};

		for (const auto& meshName : importedMeshNames)
		{
			// ===================================================================================
			// Triangle quality metrics eval for repulsive surfaces results [Yu, et al., 2021],
			// and for point cloud LSW [Daniel, et al., 2015]
			// -----------------------------------------------------------------------------------
			try
			{
				std::cout << "MetricsEval: " << meshName << "...\n";
				pmp::SurfaceMesh mesh;
				mesh.read(dataDirPath + meshName + ".obj");

				{
					const auto meshSurfArea = surface_area(mesh);
					assert(meshSurfArea > 0.0f);
					const size_t nVerts = mesh.n_vertices();
					const auto vertexDensity = static_cast<float>(nVerts) / meshSurfArea;
					std::cout << "Evaluated vertex density = (nVerts / meshSurfArea) = " << nVerts << "/" << meshSurfArea << " = " << vertexDensity << " verts / unit^2.\n";
					const auto bbox = mesh.bounds();
					const auto bboxVolume = bbox.volume();
					const auto vertexDensityPerUnitVolume = static_cast<float>(nVerts) / meshSurfArea * pow(bboxVolume, 2.0f/3.0f);
					std::cout << "Evaluated vertex density normalized per unit volume = nVerts / meshSurfArea * meshBBoxVolume^(2/3) = " <<
						nVerts << " / " << meshSurfArea << " * (" << bboxVolume << ")^(2/3) = " << vertexDensityPerUnitVolume << " verts.\n";
				}

				if (!Geometry::ComputeEquilateralTriangleJacobianConditionNumbers(mesh))
				{
					std::cout << "Error!\n";
					continue;
				}

				std::cout << "Writing to " << meshName << "_Metric.vtk ...";
				mesh.write(dataOutPath + meshName + "_Metric.vtk");
				std::cout << "done\n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > MetricsEval subroutine has thrown an exception! Continue... < < < < < \n";
			}
			std::cout << "---------------------------------------------\n";
		}
	}

	if (performHistogramResultEvaluation)
	{
		const std::vector<std::string> importedMeshNames{
			"bunny_RepulsiveResult220",
			"bunnyLSW150_Obstacle",
			"bunnyLSW150_FullWrap",
			"bunnyLSW200_Daniel30k"
		};

		const std::map<std::string, std::string> correspondingPtCloudNames{
			{"bunny_RepulsiveResult220", "bunny_PtCloudRep"},
			{"bunnyLSW150_Obstacle", "bunnyPts_3"},
			{"bunnyLSW150_FullWrap", "bunnyPts_3"},
			{"bunnyLSW200_Daniel30k", "bunnyPts_3_Daniel"}
		};

		const std::map<std::string, pmp::vec3> correspondingCorrectionTranslations{
			{"bunny_RepulsiveResult220", pmp::vec3{}},
			{"bunnyLSW150_Obstacle", pmp::vec3{0.001f, 0.001f, 0.001f}},
			{"bunnyLSW150_FullWrap", pmp::vec3{0.001f, 0.001f, 0.001f}},
			{"bunnyLSW200_Daniel30k", pmp::vec3{}}
		};

		constexpr unsigned int nVoxelsPerMinDimension = 40;

		for (const auto& meshName : importedMeshNames)
		{
			try
			{
				std::cout << "DistanceMetricsEval: " << meshName << "...\n";
				pmp::SurfaceMesh mesh;
				mesh.read(dataDirPath + meshName + ".obj");

				const auto& ptCloudName = correspondingPtCloudNames.at(meshName);
				const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataDirPath + ptCloudName + ".ply", true);
				if (!ptCloudOpt.has_value())
				{
					std::cerr << "ptCloudOpt == nullopt!\n";
					break;
				}
				const auto& ptCloudTranslationVec = correspondingCorrectionTranslations.at(meshName);
				const auto& ptCloud = ptCloudOpt.value();
				auto adjustedPtCloud = std::vector<pmp::Point>(ptCloud);
				for (auto& pt : adjustedPtCloud)
					pt += ptCloudTranslationVec;

				const auto distHistogramOpt = Geometry::ComputeMeshDistanceToPointCloudPerVertexHistogram(mesh, adjustedPtCloud, 30);

				if (!distHistogramOpt.has_value())
				{
					std::cerr << "distHistogramOpt evaluation error!\n";
					break;
				}

				const auto& distHistogram = distHistogramOpt.value();
				Geometry::PrintHistogramResultData(distHistogram, std::cout);
			}
			catch (...)
			{
				std::cerr << "> > > > > > MetricsEval subroutine has thrown an exception! Continue... < < < < < \n";
			}
			std::cout << "---------------------------------------------\n";
		}
	}

	if (performOldResultJacobianMetricEval)
	{
		const std::vector<std::string> importedMeshNames{
			"armadillo_ResultOLD",
			"rockerArm_ResultOLD"
		};

		for (const auto& meshName : importedMeshNames)
		{
			std::cout << "MetricsEval: " << meshName << "...\n";
			pmp::SurfaceMesh mesh;
			mesh.read(dataDirPath + meshName + ".obj");

			if (!Geometry::ComputeEquilateralTriangleJacobianConditionNumbers(mesh))
			{
				std::cout << "Error!\n";
				continue;
			}

			mesh.write(dataOutPath + meshName + "_Metric.vtk");
		}
	}

	if (performHausdorffDistanceMeasurementsPerTimeStep)
	{
		const std::vector<std::string> procedureNames{
			"bunnyRepulsive",
			"bunnyDanielLSW",
			"bunnyLSWObstacle",
			"bunnyLSWFullWrap"
		};

		const std::map<std::string, std::string> correspondingDataset{
			{"bunnyRepulsive", dataOutPath + "bunny_RepulsiveObstacleEvol/frame"},
			{"bunnyDanielLSW", dataOutPath + "bunnyPts_3_DanielEvol/ico02_bunnyPts_3_ply_df300__Evolution_time"},
			{"bunnyLSWObstacle", dataOutPath + "bunnyPts_3_40kVertZeroSineTwoAdvect/bunnyPts_3_Evol_"},
			{"bunnyLSWFullWrap", dataOutPath + "bunnyPts_3_40kVertFullWrap/bunnyPts_3_Evol_"}
		};

		const std::map<std::string, std::string> correspondingDatasetFormat{
			{"bunnyRepulsive", ".obj"},
			{"bunnyDanielLSW", ".vtk"},
			{"bunnyLSWObstacle", ".vtk"},
			{"bunnyLSWFullWrap", ".vtk"}
		};

		const std::map<std::string, std::string> correspondingPtCloudNames{
			{"bunnyRepulsive", "bunny_PtCloudRep"},
			{"bunnyDanielLSW", "bunnyPts_3_Daniel"},
			{"bunnyLSWObstacle", "bunnyPts_3"},
			{"bunnyLSWFullWrap", "bunnyPts_3"}
		};

		const std::map<std::string, pmp::vec3> correspondingCorrectionTranslations{
			{"bunnyRepulsive", pmp::vec3{}},
			{"bunnyDanielLSW", pmp::vec3{0.001f, 0.001f, 0.001f}},
			{"bunnyLSWObstacle", pmp::vec3{0.001f, 0.001f, 0.001f}},
			{"bunnyLSWFullWrap", pmp::vec3{}}
		};

		const std::map<std::string, unsigned int> correspondingNSteps{
			{"bunnyRepulsive", 220},
			{"bunnyDanielLSW", 200},
			{"bunnyLSWObstacle", 146},
			{"bunnyLSWFullWrap", 146}
		};
		//const std::map<std::string, unsigned int> correspondingNSteps{
		//	{"bunnyRepulsive", 5},
		//	{"bunnyDanielLSW", 5},
		//	{"bunnyLSWObstacle", 5},
		//	{"bunnyLSWFullWrap", 5}
		//};

		constexpr unsigned int nVoxelsPerMinDimension = 40;

		for (const auto& procedureName : procedureNames)
		{
			std::cout << "----------------------------------------------------------------------\n";
			std::cout << "Procedure: " << procedureName << " Hausdorff Distance per time step...\n";
			std::cout << "----------------------------------------------------------------------\n";
			const auto& ptCloudName = correspondingPtCloudNames.at(procedureName);
			const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataDirPath + ptCloudName + ".ply", true);
			if (!ptCloudOpt.has_value())
			{
				std::cerr << "ptCloudOpt == nullopt!\n";
				break;
			}
			const auto& ptCloudTranslationVec = correspondingCorrectionTranslations.at(procedureName);
			const auto& ptCloud = ptCloudOpt.value();
			auto adjustedPtCloud = std::vector(ptCloud);
			for (auto& pt : adjustedPtCloud)
				pt += ptCloudTranslationVec;

			const std::function timeIndexFormatFunction = (procedureName == "bunnyRepulsive") ?
				Utils::FormatIndex4DigitFill :
				Utils::FormatIndexSimple;

			// Compute distance field for the point cloud
			const pmp::BoundingBox ptCloudBBox(ptCloud);
			const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
			const float ptCloudMinSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
			const float ptCloudCellSize = ptCloudMinSize / static_cast<float>(nVoxelsPerMinDimension);
			const SDF::PointCloudDistanceFieldSettings ptCloudDfSettings{
				ptCloudCellSize,
				1.0f, // volExpansionFactor
				Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
				SDF::BlurPostprocessingType::None
			};
			const auto ptCloudDf = SDF::PointCloudDistanceFieldGenerator::Generate(ptCloud, ptCloudDfSettings);

			try
			{
				const unsigned int nSteps = correspondingNSteps.at(procedureName);
				const auto& format = correspondingDatasetFormat.at(procedureName);
				for (unsigned int i = 0; i < nSteps; ++i)
				{
					const auto timeStepIdString = timeIndexFormatFunction(i);
					const auto meshName = correspondingDataset.at(procedureName) + timeStepIdString;
					pmp::SurfaceMesh mesh;
					mesh.read(meshName + format);

					const auto hDistOpt = Geometry::ComputeMeshToPointCloudHausdorffDistance(mesh, ptCloud, ptCloudDf, nVoxelsPerMinDimension);
					if (!hDistOpt.has_value())
					{
						std::cerr << "hDistOpt == nullopt!\n";
						break;
					}

					std::cout << "step " << i << ": " << hDistOpt.value() << "\n";
				}				
			}
			catch (...)
			{
				std::cerr << "> Procedure \"" << procedureName << "\" has thrown an exception! Continue... < < < < < \n";
			}
		}
	}

	if (performDirectHigherGenusPtCloudSampling)
	{
		const std::vector<std::string> importedMeshNames{
			"1Torus",
			"2Torus",
			"3Torus",
			"4Torus",
			"5Torus"
		};

		for (const auto& meshName : importedMeshNames)
		{
			std::cout << "==================================================================\n";
			std::cout << "Mesh To Pt Cloud: " << meshName << ".obj -> " << meshName << "Pts.ply\n";
			std::cout << "------------------------------------------------------------------\n";
			const auto baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", false);
			if (!baseDataOpt.has_value())
			{
				std::cerr << "baseDataOpt == nullopt!\n";
				break;
			}
			std::cout << "meshName.obj" << " imported as BaseMeshGeometryData.\n";
			const auto& baseData = baseDataOpt.value();
			const auto nVerts = baseData.Vertices.size();
			std::cout << "Sampling " << nVerts << "/" << nVerts << " vertices...\n";

			// Export sampled vertices to PLY
			std::string filename = dataOutPath + meshName + "Pts.ply";
			if (!ExportSampledVerticesToPLY(baseData, nVerts, filename))
			{
				std::cerr << "ExportSampledVerticesToPLY failed!\n";
				break;
			}
		}
	}

	if (performHigherGenusPtCloudLSW)
	{
		const std::vector<std::string> importedPtCloudNames{
			//"1TorusPts",
			"2TorusPts",
			//"3TorusPts",
			//"4TorusPts",
			//"5TorusPts"
		};
		const std::map<std::string, double> timeStepSizesForPtClouds{
			{"1TorusPts", 0.05 },
			{"2TorusPts", 0.05 },
			{"3TorusPts", 0.05 },
			{"4TorusPts", 0.05 },
			{"5TorusPts", 0.05 }
		};
		const std::map<std::string, double> isoLevelOffsetFactors{
			{"1TorusPts", 1.5 },
			{"2TorusPts", 1.5 },
			{"3TorusPts", 1.5 },
			{"4TorusPts", 1.5 },
			{"5TorusPts", 1.5 }
		};

		SetRemeshingAdjustmentTimeIndices({ 3, 10, 20/*, 50 , 100, 120, 140, 145*/});

		constexpr unsigned int NTimeSteps = 65;
		constexpr unsigned int nVoxelsPerMinDimension = 40;
		constexpr double defaultTimeStep = 0.05;
		constexpr double defaultOffsetFactor = 1.5;
		for (const auto& ptCloudName : importedPtCloudNames)
		{
			const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataOutPath + ptCloudName + ".ply", true);
			if (!ptCloudOpt.has_value())
			{
				std::cerr << "ptCloudOpt == nullopt!\n";
				break;
			}

			const auto& ptCloud = ptCloudOpt.value();
			const pmp::BoundingBox ptCloudBBox(ptCloud);
			const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
			const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
			const float maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
			const float cellSize = minSize / nVoxelsPerMinDimension;
			constexpr float volExpansionFactor = 1.0f;
			const SDF::PointCloudDistanceFieldSettings dfSettings{
						cellSize,
						volExpansionFactor,
						Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
						SDF::BlurPostprocessingType::None
			};

			std::cout << "==================================================================\n";
			std::cout << "Pt Cloud to DF: " << ptCloudName << ".ply -> " << ptCloudName << "_DF_" << nVoxelsPerMinDimension << "voxPerMinDim.vti\n";
			std::cout << "------------------------------------------------------------------\n";

			const auto startSDF = std::chrono::high_resolution_clock::now();
			auto sdf = SDF::PointCloudDistanceFieldGenerator::Generate(ptCloud, dfSettings);
			const auto endSDF = std::chrono::high_resolution_clock::now();

			//NormalizeScalarGridValues(sdf);

			SDF::ReportOutput(sdf, std::cout);
			const std::chrono::duration<double> timeDiff = endSDF - startSDF;
			std::cout << "DF Time: " << timeDiff.count() << " s\n";
			std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
			//ExportToVTI(dataOutPath + ptCloudName + "DF", sdf);

			//const auto& sdfBox = sdf.Box();
			//const auto sdfBoxSize = sdfBox.max() - sdfBox.min();
			//const auto sdfBoxMaxDim = std::max<double>({ sdfBoxSize[0], sdfBoxSize[1], sdfBoxSize[2] });

			const double isoLvlOffsetFactor = (timeStepSizesForPtClouds.contains(ptCloudName) ? isoLevelOffsetFactors.at(ptCloudName) : defaultOffsetFactor);
			const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

			const double tau = (timeStepSizesForPtClouds.contains(ptCloudName) ? timeStepSizesForPtClouds.at(ptCloudName) : defaultTimeStep); // time step

			MeshTopologySettings topoParams;
			topoParams.FixSelfIntersections = true;
			topoParams.MinEdgeMultiplier = 0.14f;
			topoParams.UseBackProjection = false;
			topoParams.PrincipalCurvatureFactor = 3.2f;
			topoParams.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);
			topoParams.EdgeLengthDecayFactor = 0.7f;
			topoParams.ExcludeEdgesWithoutBothFeaturePts = true;
			topoParams.FeatureType = FeatureDetectionType::MeanCurvature;

			AdvectionDiffusionParameters adParams{
				1.0, 1.0,
				2.0, 2.0
			};

			SurfaceEvolutionSettings seSettings{
				ptCloudName,
				NTimeSteps,
				tau,
				fieldIsoLevel,
				2, // IcoSphereSubdivisionLevel
				adParams,
				topoParams,
				minSize, maxSize,
				ptCloudBBox.center(),
				true, false,
				dataOutPath,
				MeshLaplacian::Voronoi,
				{"minAngle", "maxAngle", "jacobianConditionNumber", "equilateralJacobianCondition",/* "stiffnessMatrixConditioning" */},
				0.05f,
				true,
				false
			};
			ReportInput(seSettings, std::cout);
			SurfaceEvolver evolver(sdf, volExpansionFactor, seSettings);

			try
			{
				evolver.Evolve();
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > SurfaceEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	}

	if (performTriTriIntersectionTests)
	{
		const std::vector tri0Vertices{
			pmp::vec3{0.0572398156f, - 0.0717622489f, - 1.03659534f},
			pmp::vec3{0.198035538f, -0.0351596251f, -1.18061316f},
			pmp::vec3{0.250069439f, -0.122829974f, -1.05019867f}
		};
		const std::vector tri1Vertices{
			pmp::vec3{0.0789409578f, -0.0649886578f, -1.11610115f},
			pmp::vec3{0.211945787f, 0.000253308594f, -1.07397616f},
			pmp::vec3{0.233549535f, -0.0453912392f, -1.20364273f}
		};

		if (Geometry::TriangleIntersectsTriangle(tri0Vertices, tri1Vertices))
			std::cout << "Geometry::TriangleIntersectsTriangle: tri0 intersects tri1.\n";
		else
			std::cout << "Geometry::TriangleIntersectsTriangle: tri0 does not intersect tri1.\n";

		const auto iLineOpt = Geometry::ComputeTriangleTriangleIntersectionLine(tri0Vertices, tri1Vertices);
		if (iLineOpt.has_value())
		{
			std::cout << "Geometry::ComputeTriangleTriangleIntersectionLine: tri0 intersects tri1 at (" << iLineOpt.value().first << ")->(" << iLineOpt.value().second << ").\n";
			const std::vector<std::vector<pmp::vec3>> cutPolylines = { {iLineOpt.value().first, iLineOpt.value().second} };
			if (!Geometry::ExportPolylinesToOBJ(cutPolylines, dataOutPath + "dummy_cutPolylines.obj"))
			{
				std::cerr << "Geometry::ExportPolylinesToOBJ failed!\n";
			}
		}
		else
			std::cout << "Geometry::ComputeTriangleTriangleIntersectionLine: tri0 does not intersect tri1.\n";
	}

	if (performMeshSelfIntersectionTests)
	{
		const std::vector<std::string> importedMeshNames{
			//"3holes",
			"SelfIntersection2TorusTest_1",
			"SelfIntersection2TorusTest_2",
			"SelfIntersection2TorusTest_3"
		};

		for (const auto& meshName : importedMeshNames)
		{
			std::cout << "==================================================================\n";
			std::cout << "Mesh Self-Intersection Test: " << meshName << ".obj\n";
			std::cout << "------------------------------------------------------------------\n";
			pmp::SurfaceMesh mesh;
			mesh.read(dataDirPath + meshName + ".obj");

			const auto hasSelfIntersections = Geometry::PMPSurfaceMeshHasSelfIntersections(mesh);
			std::cout << (hasSelfIntersections ? "Self-intersections detected!\n" : "No self-intersections were detected.\n");
			const auto nSelfIntFaces = Geometry::CountPMPSurfaceMeshSelfIntersectingFaces(mesh, true);
			std::cout << "Detected " << nSelfIntFaces << " faces intersecting another face.\n";
			//const auto fIntMMap = Geometry::ExtractPMPSurfaceMeshFaceIntersectionMultimap(mesh);
			//Utils::PrintFaceIntersectionsMultimap(fIntMMap);
			const auto cutPolylines = Geometry::ComputeSurfaceMeshSelfIntersectionPolylines(mesh);
			if (!cutPolylines.empty())
			{
				std::cout << "Writing " << cutPolylines.size() << " cutting polylines to " << meshName << "_cutPolylines.obj ...";
				if (!Geometry::ExportPolylinesToOBJ(cutPolylines, dataOutPath + meshName + "_cutPolylines.obj"))
				{
					std::cerr << "Geometry::ExportPolylinesToOBJ failed!\n";
					break;
				}
				std::cout << "done\n";
			}

			Geometry::ConvertPMPSurfaceMeshBoolFacePropertyToScalarVertexProperty(mesh, "f:isSelfIntersecting");
			std::cout << "Writing to " << meshName << "_SelfIntersections.vtk ...";
			mesh.write(dataOutPath + meshName + "_SelfIntersections.vtk");
			std::cout << "done\n";

		}
	}

	if (performHurtadoMeshesIsosurfaceEvolverTests)
	{
		const std::vector<std::string> meshNames{
			"drone",
			//"engine",
			//"trex"
		};

		//constexpr unsigned int nVoxelsPerMinDimension = 40;
		constexpr double defaultTimeStep = 0.02;
		const std::map<std::string, double> timeStepSizesForMeshes{
			{"drone", 0.02 },
			{"engine", 0.02 },
			{"trex", defaultTimeStep }
		};
		const std::map<std::string, double> effectiveIsolevelsForMeshes{
			{"drone", 18.0 },
			{"engine", 20.0 },
			{"trex", 200.0 }
		};
		const std::map<std::string, float> resamplingFactors{
			{"drone", 2.0f },
			{"engine", 4.0f },
			{"trex", 0.3f }
		};

		for (const auto& name : meshNames)
		{
			//pmp::SurfaceMesh mesh;
			//mesh.read(dataDirPath + name + ".obj");
			//const auto meshBBox = mesh.bounds();
			//const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
			//const float minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
			//const float maxSize = std::max({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
			//const float cellSize = minSize / nVoxelsPerMinDimension;
			//constexpr float volExpansionFactor = 1.0f;
			//const SDF::DistanceFieldSettings sdfSettings{
			//	cellSize,
			//	volExpansionFactor,
			//	//0.2, // TODO: will this truncation be OK?
			//	Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
			//	SDF::KDTreeSplitType::Center,
			//	SDF::SignComputation::VoxelFloodFill,
			//	SDF::BlurPostprocessingType::None,
			//	SDF::PreprocessingType::Octree
			//};
			//SDF::ReportInput(mesh, sdfSettings, std::cout);

			//const auto startSDF = std::chrono::high_resolution_clock::now();
			//const auto sdf = SDF::DistanceFieldGenerator::Generate(mesh, sdfSettings);
			//const auto endSDF = std::chrono::high_resolution_clock::now();

			//SDF::ReportOutput(sdf, std::cout);
			//const std::chrono::duration<double> timeDiff = endSDF - startSDF;
			//std::cout << "SDF Time: " << timeDiff.count() << " s\n";
			//std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
			//ExportToVTI(dataOutPath + name + "SDF", sdf);

			// DistanceFieldGenerator::Generate only supports pmp::Surface mesh which cannot have non-manifoldness
			// TODO: be able to process BaseMeshGeometryData with DistanceFieldGenerator::Generate

			std::cout << "Opening " << name << "_voxFieldSDF90_FSM_90.vti ...";
			constexpr float volExpansionFactor = 1.0f;
			const auto sdf = ImportVTI(dataOutPath + name + "_voxFieldSDF90_FSM_90.vti");
			std::cout << " done.\n";
			const auto cellSize = sdf.CellSize();

			//ExportToVTI(dataOutPath + name + "_reexport", sdf);
			//continue;

			//constexpr double fieldIsoLevel = 0.0;
			const double fieldIsoLevel = 0.01 * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);
			const double isoLevel = (effectiveIsolevelsForMeshes.contains(name) ?
				(fieldIsoLevel < effectiveIsolevelsForMeshes.at(name) ? effectiveIsolevelsForMeshes.at(name) : 1.1 * fieldIsoLevel) : 5.0);

			const MeshTopologySettings topoParams{
				true,
				0.24f, //0.1f, //0.4f,
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
				1.0f * static_cast<float>(M_PI_2),
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
				{1.0f, 1.0, 2.0, 1.0},
				topoParams,
				true, false,
				dataOutPath,
				MeshLaplacian::Voronoi,
				{"minAngle", "maxAngle", "jacobianConditionNumber", "equilateralJacobianCondition",/* "stiffnessMatrixConditioning" */},
				0.05f,
				true,
				false
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
	} // endif performHurtadoMeshesIsosurfaceEvolverTests

	if (performHurtadoTrexIcosphereLSW)
	{
		constexpr double defaultTimeStep = 0.05;
		SetRemeshingAdjustmentTimeIndices({ 3, 8, 10, 20, 50, 100, /* 120, 140, 145*/ });

		constexpr float minSize = 3996.9329f;
		constexpr float maxSize = 11613.2236f;
		const pmp::vec3 center{ 5806.6118f, 2353.8142f, 2005.4388f };

		std::cout << "Opening " << "trex_voxFieldSDF90_FSM_90.vti ...";
		const auto sdf = ImportVTI(dataOutPath + "trex_voxFieldSDF90_FSM_90.vti");
		std::cout << " done.\n";
		const auto cellSize = sdf.CellSize();

		const double isoLvlOffsetFactor = 1.5;
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		const double tau = defaultTimeStep; // time step

		MeshTopologySettings topoParams;
		topoParams.FixSelfIntersections = true;
		topoParams.MinEdgeMultiplier = 0.14f;
		topoParams.UseBackProjection = false;
		topoParams.PrincipalCurvatureFactor = 3.2f;
		topoParams.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);
		topoParams.EdgeLengthDecayFactor = 0.7f;
		topoParams.ExcludeEdgesWithoutBothFeaturePts = true;
		topoParams.FeatureType = FeatureDetectionType::MeanCurvature;

		AdvectionDiffusionParameters adParams{
			1.0, 1.0,
			1.0, 1.5
		};

		SurfaceEvolutionSettings seSettings{
			"TRex",
			146, // 150,
			tau,
			fieldIsoLevel,
			3, // IcoSphereSubdivisionLevel
			adParams,
			topoParams,
			minSize, maxSize,
			center,
			true, false,
			dataOutPath,
			MeshLaplacian::Voronoi,
			{"minAngle", "maxAngle", "jacobianConditionNumber", "equilateralJacobianCondition",/* "stiffnessMatrixConditioning" */},
			0.05f,
			true,
			false
		};
		ReportInput(seSettings, std::cout);
		SurfaceEvolver evolver(sdf, 1.0f, seSettings);

		try
		{
			evolver.Evolve();
		}
		catch (...)
		{
			std::cerr << "> > > > > > > > > > > > > > SurfaceEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
		}


	} // endif performHurtadoTrexIcosphereLSW

	if (performImportVTIDebugTests)
	{
		std::cout << "Importing MetaBallVals.vti ...";
		const auto grid = ImportVTI(dataOutPath + "MetaBallVals.vti");
		std::cout << " done.\n";
		std::cout << "Exporting MetaBallVals_reexport.vti ...";
		ExportToVTI(dataOutPath + "MetaBallVals_reexport", grid);
		std::cout << " done.\n";
	} // endif performImportVTIDebugTests

	if (performConvexHullTests)
	{
		const std::vector<std::string> importedPtCloudNames{
			"bunnyPts_3"//,
			//"CaesarBustPts_3"
		};

		for (const auto& ptCloudName : importedPtCloudNames)
		{
			const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataOutPath + ptCloudName + ".ply", true);
			if (!ptCloudOpt.has_value())
			{
				std::cerr << "ptCloudOpt == nullopt!\n";
				break;
			}

			const auto& ptCloud = ptCloudOpt.value();

			const auto convexHullMeshOpt = Geometry::ComputePointCloudConvexHull(ptCloud);
			if (!convexHullMeshOpt.has_value())
			{
				std::cerr << "convexHullMeshOpt == nullopt!\n";
				break;
			}
			const auto& convexHull = convexHullMeshOpt.value();

			convexHull.write(dataOutPath + ptCloudName + "_convexHull.obj");
		}
	} // performConvexHullTests
}
