// -------------------------- RevisitingOldEvolverWSCG ----------------------------
//  Spring of 2024, WSCG 2024
// ................................................................................
/*
 * The resubmission of the GRAPP 2023 paper to GRAPP 2024 was unsuccessful. The article
 * was reviewed far more strictly, and more work needed to be done. Most importantly, it
 * lacked application and well-tuned experiments as well a comparison with existing work.
 * By comparing the ability of point cloud surface reconstruction (approximation) with the
 * earlier work of [Daniel et al. 2015] and [Yu, et al. 2021], and the CAD model simplification options tested by
 * a model of [Hurtado et al. 2022] using only the isosurface evolver, this work got published
 * and presented at WSCG 2024 in Pilsen (see [Cavarga, 2024]).
 *
 * [Daniel et al. 2015]
 * Daniel, P., Medla, M., Mikula, K., & Remesikova, M. (2015, April).
 * Reconstruction of surfaces from point clouds using a lagrangian surface evolution model.
 * SSVM 2015, 589-600.
 *
 * [Yu, et al. 2021]
 * Yu, C., Brakensiek, C., Schumacher, H., & Crane, K. (2021).
 * Repulsive surfaces. ACM Transactions on Graphics, 40(6). ACM.
 *
 * [Hurtado et al. 2022]
 * Hurtado, J., Montenegro, A., Gattass, M., Carvalho, F., & Raposo, A. (2022).
 * Enveloping CAD models for visualization and interaction in XR applications. 
 * Engineering with Computers, 1-19.
 *
 * [Cavarga, 2024]
 * Cavarga, M. (2024).
 * Automated Surface Extraction: Adaptive Remeshing Meets Lagrangian Shrink-Wrapping.
 * Journal of WSCG 2024. 32, 21-30.
 */
 // --------------------------------------------------------------------------------

#include "pmp/algorithms/Remeshing.h"

#include "geometry/GeometryConversionUtils.h"
#include "geometry/GeometryUtil.h"
#include "geometry/GridUtil.h"

#include "sdf/SDF.h"

#include "utils/StringUtils.h"

#include "core/ConversionUtils.h"
#include "core/ConvexHullEvolver.h"
#include "core/EvolverUtilsCommon.h"
#include "core/IcoSphereEvolver.h"
#include "core/IsosurfaceEvolver.h"
#include "core/MeasurementUtils.h"
#include "core/SurfaceEvolver.h"

#include "IOEnvironment.h"
#include "../Experiments.h"

void PDanielPtCloudPLYExport()
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

void PtCloudToDF()
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
		const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
		const SDF::PointCloudDistanceFieldSettings dfSettings{
			cellSize,
			1.0,
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

void PDanielPtCloudComparisonTest()
{
	const std::vector<std::string> importedPtCloudNames{
		"bunnyPts_3"//,
		//"CaesarBustPts_3"
	};
	const std::map<std::string, double> timeStepSizesForPtClouds{
		{"bunnyPts_3", 0.07 },
		{ "CaesarBustPts_3", 0.05 }
	};
	const std::map<std::string, double> isoLevelOffsetFactors{
		{"bunnyPts_3", 2.0 },
		{ "CaesarBustPts_3", 0.5 }
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
		const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
		constexpr pmp::Scalar volExpansionFactor = 1.0;
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

		//ormalizeScalarGridValues<Geometry::ScalarGrid>(sdf);

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
			0.05,
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

void RepulsiveSurfResultEvaluation()
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
				assert(meshSurfArea > 0.0);
				const size_t nVerts = mesh.n_vertices();
				const auto vertexDensity = static_cast<pmp::Scalar>(nVerts) / meshSurfArea;
				std::cout << "Evaluated vertex density = (nVerts / meshSurfArea) = " << nVerts << "/" << meshSurfArea << " = " << vertexDensity << " verts / unit^2.\n";
				const auto bbox = mesh.bounds();
				const auto bboxVolume = bbox.volume();
				const auto vertexDensityPerUnitVolume = static_cast<pmp::Scalar>(nVerts) / meshSurfArea * pow(bboxVolume, 2.0 / 3.0);
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

void HistogramResultEvaluation()
{
	const std::vector<std::string> importedMeshNames{
		"bunny_RepulsiveResult220",
		"bunnyLSW150_Obstacle",
		"bunnyLSW150_FullWrap",
		"bunnyLSW200_Daniel30k"
	};

	const std::map<std::string, std::string> correspondingPtCloudNames{
		{"bunny_RepulsiveResult220", "bunny_PtCloudRep"},
		{ "bunnyLSW150_Obstacle", "bunnyPts_3" },
		{ "bunnyLSW150_FullWrap", "bunnyPts_3" },
		{ "bunnyLSW200_Daniel30k", "bunnyPts_3_Daniel" }
	};

	const std::map<std::string, pmp::vec3> correspondingCorrectionTranslations{
		{"bunny_RepulsiveResult220", pmp::vec3{}},
		{ "bunnyLSW150_Obstacle", pmp::vec3{0.001, 0.001, 0.001} },
		{ "bunnyLSW150_FullWrap", pmp::vec3{0.001, 0.001, 0.001} },
		{ "bunnyLSW200_Daniel30k", pmp::vec3{} }
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

void OldResultJacobianMetricEval()
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

void HausdorffDistanceMeasurementsPerTimeStep()
{
	const std::vector<std::string> procedureNames{
		"bunnyRepulsive",
		"bunnyDanielLSW",
		"bunnyLSWObstacle",
		"bunnyLSWFullWrap"
	};

	const std::map<std::string, std::string> correspondingDataset{
		{"bunnyRepulsive", dataOutPath + "bunny_RepulsiveObstacleEvol/frame"},
		{ "bunnyDanielLSW", dataOutPath + "bunnyPts_3_DanielEvol/ico02_bunnyPts_3_ply_df300__Evolution_time" },
		{ "bunnyLSWObstacle", dataOutPath + "bunnyPts_3_40kVertZeroSineTwoAdvect/bunnyPts_3_Evol_" },
		{ "bunnyLSWFullWrap", dataOutPath + "bunnyPts_3_40kVertFullWrap/bunnyPts_3_Evol_" }
	};

	const std::map<std::string, std::string> correspondingDatasetFormat{
		{"bunnyRepulsive", ".obj"},
		{ "bunnyDanielLSW", ".vtk" },
		{ "bunnyLSWObstacle", ".vtk" },
		{ "bunnyLSWFullWrap", ".vtk" }
	};

	const std::map<std::string, std::string> correspondingPtCloudNames{
		{"bunnyRepulsive", "bunny_PtCloudRep"},
		{ "bunnyDanielLSW", "bunnyPts_3_Daniel" },
		{ "bunnyLSWObstacle", "bunnyPts_3" },
		{ "bunnyLSWFullWrap", "bunnyPts_3" }
	};

	const std::map<std::string, pmp::vec3> correspondingCorrectionTranslations{
		{"bunnyRepulsive", pmp::vec3{}},
		{ "bunnyDanielLSW", pmp::vec3{0.001, 0.001, 0.001} },
		{ "bunnyLSWObstacle", pmp::vec3{0.001, 0.001, 0.001} },
		{ "bunnyLSWFullWrap", pmp::vec3{} }
	};

	const std::map<std::string, unsigned int> correspondingNSteps{
		{"bunnyRepulsive", 220},
		{ "bunnyDanielLSW", 200 },
		{ "bunnyLSWObstacle", 146 },
		{ "bunnyLSWFullWrap", 146 }
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
		const pmp::Scalar ptCloudMinSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar ptCloudCellSize = ptCloudMinSize / static_cast<pmp::Scalar>(nVoxelsPerMinDimension);
		const SDF::PointCloudDistanceFieldSettings ptCloudDfSettings{
			ptCloudCellSize,
				1.0, // volExpansionFactor
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

void DirectHigherGenusPtCloudSampling()
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

void HigherGenusPtCloudLSW()
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
		{ "2TorusPts", 0.05 },
		{ "3TorusPts", 0.05 },
		{ "4TorusPts", 0.05 },
		{ "5TorusPts", 0.05 }
	};
	const std::map<std::string, double> isoLevelOffsetFactors{
		{"1TorusPts", 1.5 },
		{ "2TorusPts", 1.5 },
		{ "3TorusPts", 1.5 },
		{ "4TorusPts", 1.5 },
		{ "5TorusPts", 1.5 }
	};

	SetRemeshingAdjustmentTimeIndices({ 3, 10, 20/*, 50 , 100, 120, 140, 145*/ });

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
		const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
		constexpr pmp::Scalar volExpansionFactor = 1.0;
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
		topoParams.MinEdgeMultiplier = 0.14;
		topoParams.UseBackProjection = false;
		topoParams.PrincipalCurvatureFactor = 3.2;
		topoParams.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);
		topoParams.EdgeLengthDecayFactor = 0.7;
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
			0.05,
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

void TriTriIntersectionTests()
{
	const std::vector tri0Vertices{
		pmp::vec3{0.0572398156, -0.0717622489, -1.03659534},
			pmp::vec3{0.198035538, -0.0351596251, -1.18061316},
			pmp::vec3{0.250069439, -0.122829974, -1.05019867}
	};
	const std::vector tri1Vertices{
		pmp::vec3{0.0789409578, -0.0649886578, -1.11610115},
			pmp::vec3{0.211945787, 0.000253308594, -1.07397616},
			pmp::vec3{0.233549535, -0.0453912392, -1.20364273}
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

void MeshSelfIntersectionTests()
{
	const std::vector<std::string> importedMeshNames{
		//"3holes",
		"SelfIntersection2TorusTest_1",
		//"SelfIntersection2TorusTest_2",
		//"SelfIntersection2TorusTest_3"
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
		const auto fIntMMap = Geometry::ExtractPMPSurfaceMeshFaceIntersectionMultimap(mesh);
		Utils::PrintFaceIntersectionsMultimap(fIntMMap);
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

void HurtadoMeshesIsosurfaceEvolverTests()
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
		{ "engine", 0.02 },
		{ "trex", defaultTimeStep }
	};
	const std::map<std::string, double> effectiveIsolevelsForMeshes{
		{"drone", 18.0 },
		{ "engine", 20.0 },
		{ "trex", 200.0 }
	};
	const std::map<std::string, pmp::Scalar> resamplingFactors{
		{"drone", 2.0 },
		{ "engine", 4.0 },
		{ "trex", 0.3 }
	};

	for (const auto& name : meshNames)
	{
		//pmp::SurfaceMesh mesh;
		//mesh.read(dataDirPath + name + ".obj");
		//const auto meshBBox = mesh.bounds();
		//const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
		//const pmp::Scalar minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		//const pmp::Scalar maxSize = std::max({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		//const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
		//constexpr pmp::Scalar volExpansionFactor = 1.0;
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

		std::cout << "Opening " << name << "_voxFieldSDF90_FSM_90.vti ...";
		constexpr pmp::Scalar volExpansionFactor = 1.0;
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
			0.24, //0.1, //0.4,
			0.0,
			1.0,
			0.0,
			2,
			0.0,
			3,
			5,
			false,
			FeatureDetectionType::MeanCurvature,
			1.0 * M_PI_2 * 180.0, 2.0 * M_PI_2 * 180.0,
			2.0,
			1.0 * static_cast<pmp::Scalar>(M_PI_2),
			true
		};

		const double tau = (timeStepSizesForMeshes.contains(name) ? timeStepSizesForMeshes.at(name) : defaultTimeStep); // time step
		const pmp::Scalar resamplingFactor = (resamplingFactors.contains(name) ? resamplingFactors.at(name) : 1.5);
		IsoSurfaceEvolutionSettings seSettings{
			name,
			20,
			tau,
			fieldIsoLevel,
			isoLevel,
			cellSize * resamplingFactor,
			{1.0, 1.0, 2.0, 1.0},
			topoParams,
			true, false,
			dataOutPath,
			MeshLaplacian::Voronoi,
			{"minAngle", "maxAngle", "jacobianConditionNumber", "equilateralJacobianCondition",/* "stiffnessMatrixConditioning" */},
			0.05,
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
}

void HurtadoTrexIcosphereLSW()
{
	constexpr double defaultTimeStep = 0.05;
	SetRemeshingAdjustmentTimeIndices({ 3, 8, 10, 20, 50, 100, /* 120, 140, 145*/ });

	constexpr pmp::Scalar minSize = 3996.9329;
	constexpr pmp::Scalar maxSize = 11613.2236;
	const pmp::vec3 center{ 5806.6118, 2353.8142, 2005.4388 };

	std::cout << "Opening " << "trex_voxFieldSDF90_FSM_90.vti ...";
	const auto sdf = ImportVTI(dataOutPath + "trex_voxFieldSDF90_FSM_90.vti");
	std::cout << " done.\n";
	const auto cellSize = sdf.CellSize();

	const double isoLvlOffsetFactor = 1.5;
	const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

	const double tau = defaultTimeStep; // time step

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
		0.05,
		true,
		false
	};
	ReportInput(seSettings, std::cout);
	SurfaceEvolver evolver(sdf, 1.0, seSettings);

	try
	{
		evolver.Evolve();
	}
	catch (...)
	{
		std::cerr << "> > > > > > > > > > > > > > SurfaceEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
	}
}

void ImportVTIDebugTests()
{
	std::cout << "Importing MetaBallVals.vti ...";
	const auto grid = ImportVTI(dataOutPath + "MetaBallVals.vti");
	std::cout << " done.\n";
	std::cout << "Exporting MetaBallVals_reexport.vti ...";
	ExportToVTI(dataOutPath + "MetaBallVals_reexport", grid);
	std::cout << " done.\n";
}


// More focused LSW experiments were carried out to see whether they'll be applicable
// for the IncrementalMeshBuilder. The result was that they were kind of hard to work with.

void ConvexHullTests()
{
	const std::vector<std::string> importedPtCloudNames{
		"bunnyPts_3"//,
		//"CaesarBustPts_3"
	};

	const std::map<std::string, pmp::Point> slicingPlaneRefPts{
		{"bunnyPts_3", pmp::Point{-0.01684039831161499, 0.11015420407056808, 0.0012007840834242693} },
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"bunnyPts_3", pmp::vec3{0.0, 0.0, 1.0} },
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

		const auto convexHullMeshOpt = Geometry::ComputePMPConvexHullFromPoints(ptCloud);
		if (!convexHullMeshOpt.has_value())
		{
			std::cerr << "convexHullMeshOpt == nullopt!\n";
			break;
		}
		const auto& convexHull = convexHullMeshOpt.value();

		convexHull.write(dataOutPath + ptCloudName + "_convexHull.obj");

		const pmp::BoundingBox ptCloudBBox(ptCloud);
		const auto center = ptCloudBBox.center();
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		//const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar distTolerance = 0.01 * minSize;

		const auto planeRefPt = (slicingPlaneRefPts.contains(ptCloudName) ? slicingPlaneRefPts.at(ptCloudName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(ptCloudName) ? slicingPlaneNormals.at(ptCloudName) : pmp::vec3{ -1.0, 0.0, 0.0 });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for point cloud " << ptCloudName << "!\n";
			continue;
		}

		const auto delaunayMeshOpt = Geometry::ComputeDelaunayMeshFrom2DPoints(pts2D);
		if (!delaunayMeshOpt.has_value())
		{
			std::cerr << "delaunayMeshOpt == nullopt!\n";
			break;
		}

		if (!ExportBaseMeshGeometryDataToOBJ(*delaunayMeshOpt, dataOutPath + ptCloudName + "_2DSliceDelaunay.obj"))
		{
			std::cerr << "ExportBaseMeshGeometryDataToOBJ failed!\n";
			break;
		}
	}
}

void ConvexHullRemeshingTests()
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

		auto convexHullMeshOpt = Geometry::ComputePMPConvexHullFromPoints(ptCloud);
		if (!convexHullMeshOpt.has_value())
		{
			std::cerr << "convexHullMeshOpt == nullopt!\n";
			break;
		}
		auto& convexHull = convexHullMeshOpt.value();
		const auto [lengthMin, lengthMean, lengthMax] = Geometry::ComputeEdgeLengthMinAverageAndMax(convexHull);
		pmp::Remeshing remeshing(convexHull);
		remeshing.convex_hull_adaptive_remeshing({
			static_cast<pmp::Scalar>(2.0) * lengthMin, 
			static_cast<pmp::Scalar>(8.0) * lengthMin, 
			static_cast<pmp::Scalar>(0.5) * lengthMin,
			3, 10, true
			});
		const auto [newLengthMin, newLengthMean, newLengthMax] = Geometry::ComputeEdgeLengthMinAverageAndMax(convexHull);
		std::cout << "Edge lengths stats: [" << lengthMin << ", " << lengthMean << ", " << lengthMax << "] -> [" << newLengthMin << ", " << newLengthMean << ", " << newLengthMax << "]\n";
		convexHull.write(dataOutPath + ptCloudName + "_convexHullRemeshed.obj");
	}
}

void ConvexHullEvolverTests()
{
	const std::vector<std::string> meshForPtCloudNames{
		//"armadillo",
		//"blub",
		"bunny",
		//"maxPlanck",
		//"nefertiti",
		//"ogre",
		//"spot"
	};
	const std::map<std::string, double> timeStepSizesForPtClouds{
		{"armadillo", 0.05 },
		{ "blub", 0.05 },
		{ "bunny", 0.05 },
		{ "maxPlanck", 0.05 },
		{ "nefertiti", 0.05 },
		{ "ogre", 0.05 },
		{ "spot", 0.05 }
	};
	const std::map<std::string, double> isoLevelOffsetFactors{
		{"armadillo", 0.5 },
		{ "blub", 0.5 },
		{ "bunny", 0.5 },
		{ "maxPlanck", 0.5 },
		{ "nefertiti", 0.5 },
		{ "ogre", 0.5 },
		{ "spot", 0.5 }
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;

	constexpr size_t samplingLevel = 3;
	constexpr size_t nSamplings = 10;
	constexpr size_t minVerts = 9; // Minimum number of vertices to sample

	constexpr unsigned int seed = 5000; // seed for the pt cloud sampling RNG

	SetRemeshingAdjustmentTimeIndices({}); // no remeshing adjustment
	//SetRemeshingAdjustmentTimeIndices({ 3, 10, 20, 50, 100, 120, 140, 145 });

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
		std::cout << "meshName.obj" << " imported as BaseMeshGeometryData.\n";
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
		const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
		const double isoLvlOffsetFactor = (timeStepSizesForPtClouds.contains(ptCloudName) ? isoLevelOffsetFactors.at(ptCloudName) : defaultOffsetFactor);
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);
		const double tau = (timeStepSizesForPtClouds.contains(ptCloudName) ? timeStepSizesForPtClouds.at(ptCloudName) : defaultTimeStep); // time step

		MeshTopologySettings topoParams;
		topoParams.EdgeLengthDecayFactor = 0.7;
		topoParams.ExcludeEdgesWithoutBothFeaturePts = true;

		AdvectionDiffusionParameters adParams{
			1.0, 1.0,
			2.0, 1.0
		};

		ConvexHullSurfaceEvolutionSettings seSettings{
			meshName,
			80,
			tau,
			fieldIsoLevel,
			nVoxelsPerMinDimension,
			adParams,
			topoParams,
			ptCloudBBox.center(),
			maxSize,
			true, true,
			dataOutPath,
			MeshLaplacian::Voronoi,
			{"minAngle", "maxAngle", "jacobianConditionNumber", "equilateralJacobianCondition",/* "stiffnessMatrixConditioning" */},
			0.05,
			true,
			false,
			false,
			true
		};
		ReportCHEvolverInput(seSettings, std::cout);
		ConvexHullEvolver evolver(ptCloud, seSettings);

		try
		{
			evolver.Evolve();
		}
		catch (...)
		{
			std::cerr << "> > > > > > > > > > > > > > ConvexHullEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
		}
	}
}

void IcoSphereEvolverTests()
{
	const std::vector<std::string> meshForPtCloudNames{
		//"armadillo",
		//"blub",
		"bunny",
		//"maxPlanck",
		//"nefertiti",
		//"ogre",
		//"spot"
	};
	const std::map<std::string, double> timeStepSizesForPtClouds{
		{"armadillo", 0.05 },
		{ "blub", 0.05 },
		{ "bunny", 0.05 },
		{ "maxPlanck", 0.05 },
		{ "nefertiti", 0.05 },
		{ "ogre", 0.05 },
		{ "spot", 0.05 }
	};
	const std::map<std::string, double> isoLevelOffsetFactors{
		{"armadillo", 0.5 },
		{ "blub", 0.5 },
		{ "bunny", 1.5 },
		{ "maxPlanck", 0.5 },
		{ "nefertiti", 0.5 },
		{ "ogre", 0.5 },
		{ "spot", 0.5 }
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;

	constexpr size_t samplingLevel = 3;
	constexpr size_t nSamplings = 10;
	constexpr size_t minVerts = 9; // Minimum number of vertices to sample

	constexpr unsigned int seed = 5000; // seed for the pt cloud sampling RNG

	SetRemeshingAdjustmentTimeIndices({}); // no remeshing adjustment
	//SetRemeshingAdjustmentTimeIndices({ 3, 7, 11, 15, 20, 50, 100 });

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
		std::cout << "meshName.obj" << " imported as BaseMeshGeometryData.\n";
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
		const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
		const double isoLvlOffsetFactor = (timeStepSizesForPtClouds.contains(ptCloudName) ? isoLevelOffsetFactors.at(ptCloudName) : defaultOffsetFactor);
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);
		const double tau = (timeStepSizesForPtClouds.contains(ptCloudName) ? timeStepSizesForPtClouds.at(ptCloudName) : defaultTimeStep); // time step

		MeshTopologySettings topoParams;
		topoParams.MinEdgeMultiplier = 0.14;
		topoParams.UseBackProjection = false;
		topoParams.EdgeLengthDecayFactor = 0.7;
		topoParams.ExcludeEdgesWithoutBothFeaturePts = true;

		AdvectionDiffusionParameters adParams{
			1.0, 1.0,
			2.0, 1.0
		};

		IcoSphereEvolutionSettings seSettings{
			meshName,
			80,
			tau,
			fieldIsoLevel,
			nVoxelsPerMinDimension,
			2,
			adParams,
			topoParams,
			ptCloudBBox.center(),
			maxSize,
			true, true,
			dataOutPath,
			MeshLaplacian::Voronoi,
			{"minAngle", "maxAngle", "jacobianConditionNumber", "equilateralJacobianCondition",/* "stiffnessMatrixConditioning" */},
			0.05,
			true,
			false,
			false,
			true
		};
		ReportIcoEvolverInput(seSettings, std::cout);
		IcoSphereEvolver evolver(ptCloud, seSettings);

		try
		{
			evolver.Evolve();
		}
		catch (...)
		{
			std::cerr << "> > > > > > > > > > > > > > IcoSphereEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
		}
	}
}
