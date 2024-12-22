// ------------------------------ OldEvolverInPMP --------------------------------
//     Late 2022, GRAPP 2023.
// ...............................................................................
/*
 * This work was done from the summer of 2022 through fall,
 * while initially abandoning old custom project "Symplektis" (And yes, if you read
 * this comment, I am imposing a right on this name for future use (see the license)).
 * Symplektis was a reinvention of the wheel and I switched to the PMP library
 * as the foundation of this project. This first set of experiments was described
 * in a preprint of a submission to GRAPP 2023. The paper submission was abandoned
 * due to personal reasons.
 */
 // --------------------------------------------------------------------------------

#include "pmp/SurfaceMesh.h"

#include "geometry/GridUtil.h"
#include "geometry/PlaneBuilder.h"

#include "sdf/SDF.h"

#include "core/BrainSurfaceEvolver.h"
#include "core/ConversionUtils.h"
#include "core/IsosurfaceEvolver.h"
#include "core/SheetMembraneEvolver.h"
#include "core/SphereTest.h"
#include "core/SurfaceEvolver.h"

#include "IOEnvironment.h"
#include "../Experiments.h"

// DISCLAIMER: the names need to match the models in "DROOT_DIR/data" except for the extension (which is always *.obj)
const std::vector<std::string> meshNames{
	    //"BentChair",
		//"blub",
		"bunny",
		//"maxPlanck",
		//"nefertiti",
		//"ogre",
		//"spot"
};

void SDFTests()
{
	constexpr unsigned int nVoxelsPerMinDimension = 10;
	constexpr bool computeGradients = false;

	for (const auto& name : meshNames)
	{
		pmp::SurfaceMesh mesh;
		mesh.read(dataDirPath + name + ".obj");

		const auto meshBBox = mesh.bounds();
		const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
		const pmp::Scalar minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
		const SDF::DistanceFieldSettings sdfSettings{
			cellSize,
				1.0,
				DBL_MAX,
				SDF::KDTreeSplitType::Center,
				SDF::SignComputation::VoxelFloodFill,
				SDF::BlurPostprocessingType::None,
				SDF::PreprocessingType::Octree
		};
		SDF::ReportInput(mesh, sdfSettings, std::cout);

		const auto startSDF = std::chrono::high_resolution_clock::now();
		const Geometry::PMPSurfaceMeshAdapter meshAdapter(std::make_shared<pmp::SurfaceMesh>(mesh));
		const auto sdf = SDF::DistanceFieldGenerator::Generate(meshAdapter, sdfSettings);
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
		const pmp::Scalar expansion = 1.0 * minSize;
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
}

void SphereTestOLD()
{
	{ // Setup 1: No remeshing, No tangential redistribution
		const ST_MeshTopologySettings topoSettings{};

		const SphereTestEvolutionSettings stSettings{
			topoSettings,
			false, false, dataOutPath,
			ST_MeshLaplacian::Barycentric, 0.0, false };

		SphereTest st(stSettings);
		st.PerformTest(4);
	}

	std::cout << "=====================================================\n";

	{ // Setup 2: No remeshing, tangential redistribution with weight 0.25
		const ST_MeshTopologySettings topoSettings{};

		const SphereTestEvolutionSettings stSettings{
			topoSettings,
			false, false, dataOutPath,
			ST_MeshLaplacian::Barycentric, 0.25, false };

		SphereTest st(stSettings);
		st.PerformTest(4);
	}

	/*std::cout << "=====================================================\n";

	{ // Setup 3: Remeshing, No tangential redistribution
		const ST_MeshTopologySettings topoSettings{};

		const SphereTestEvolutionSettings stSettings{
			topoSettings,
			false, false, dataOutPath,
			ST_MeshLaplacian::Barycentric, 0.0, true };

		SphereTest st(stSettings);
		st.PerformTest(4);
	}

	std::cout << "=====================================================\n";

	{ // Setup 4: Remeshing, tangential redistribution with weight 0.25
		const ST_MeshTopologySettings topoSettings{};

		const SphereTestEvolutionSettings stSettings{
			topoSettings,
			false, false, dataOutPath,
			ST_MeshLaplacian::Barycentric, 0.25, true };

		SphereTest st(stSettings);
		st.PerformTest(4);
	}*/
}

void EvolverTests()
{
	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	const std::map<std::string, double> timeStepSizesForMeshes{
		{"armadillo", 0.05 },
		{ "BentChair", 0.05 },
		{ "blub", 0.05 },
		{ "bunny", 0.0025 },
		{ "maxPlanck", 0.05 },
		{ "nefertiti", 0.05 },
		{ "ogre", 0.05 },
		{ "spot", 0.05 }
	};

	for (const auto& name : meshNames)
	{
		pmp::SurfaceMesh mesh;
		mesh.read(dataDirPath + name + ".obj");
		const auto meshBBox = mesh.bounds();
		const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
		const pmp::Scalar minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		const pmp::Scalar maxSize = std::max({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
		constexpr pmp::Scalar volExpansionFactor = 1.0;
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
		const Geometry::PMPSurfaceMeshAdapter meshAdapter(std::make_shared<pmp::SurfaceMesh>(mesh));
		const auto sdf = SDF::DistanceFieldGenerator::Generate(meshAdapter, sdfSettings);
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
			0.05,
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
}

void OldArmadilloLSWTest()
{
	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	const std::string name = "armadillo";
	pmp::SurfaceMesh mesh;
	mesh.read(dataDirPath + name + ".obj");
	const auto meshBBox = mesh.bounds();
	const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
	const pmp::Scalar minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
	const pmp::Scalar maxSize = std::max({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
	const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
	constexpr pmp::Scalar volExpansionFactor = 1.0;
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
	const Geometry::PMPSurfaceMeshAdapter meshAdapter(std::make_shared<pmp::SurfaceMesh>(mesh));
	const auto sdf = SDF::DistanceFieldGenerator::Generate(meshAdapter, sdfSettings);
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
	topoParams.EdgeLengthDecayFactor = 0.95;
	topoParams.ExcludeEdgesWithoutBothFeaturePts = true;
	topoParams.NRemeshingIters = 5;
	topoParams.PrincipalCurvatureFactor = 2.0;
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
		0.05,
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
}

void IsosurfaceEvolverTests()
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
		{ "fertility", defaultTimeStep },
		{ "happyBuddha", defaultTimeStep },
		{ "rockerArm", defaultTimeStep }
	};
	const std::map<std::string, double> effectiveIsolevelsForMeshes{
		{"3holes", 0.02 },
		{ "fertility", 4.0 },
		{ "happyBuddha", 1.5e-3 },
		{ "rockerArm", 0.06 }
	};
	const std::map<std::string, pmp::Scalar> resamplingFactors{
		{"3holes", 3.0 },
		{ "fertility", 2.0 },
		{ "happyBuddha", 1.0 },
		{ "rockerArm", 2.0 }
	};

	for (const auto& name : higherGenusMeshNames)
	{
		pmp::SurfaceMesh mesh;
		mesh.read(dataDirPath + name + ".obj");
		const auto meshBBox = mesh.bounds();
		const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
		const pmp::Scalar minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		const pmp::Scalar maxSize = std::max({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
		constexpr pmp::Scalar volExpansionFactor = 1.0;
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
		const Geometry::PMPSurfaceMeshAdapter meshAdapter(std::make_shared<pmp::SurfaceMesh>(mesh));
		const auto sdf = SDF::DistanceFieldGenerator::Generate(meshAdapter, sdfSettings);
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
			0.4,
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
			0.8 * static_cast<pmp::Scalar>(M_PI_2),
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
			PreComputeAdvectionDiffusionParams(2.0, minSize),
			topoParams,
			true, false,
			dataOutPath,
			MeshLaplacian::Voronoi,
			{"minAngle", "maxAngle", "jacobianConditionNumber", "equilateralJacobianCondition",/* "stiffnessMatrixConditioning" */},
			0.05,
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
}

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

void BrainEvolverTests()
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
		{ "actual_brain", actualBrainThresholdSettings },
	};

	// ico-sphere params
	const pmp::vec3 talairachCenter{ 69.278477120258884, 81.210907276033296, 69.224956401243205 };
	const pmp::vec3 actualBrainCenter{ 96.536378893717142, 126.13349816811417, 116.99736547254018 };

	constexpr pmp::Scalar talairachRadius = 66.061572538428337;
	constexpr pmp::Scalar actualBrainRadius = 102.09133074846271;

	const BE_IcoSphereSettings talairachIcoSettings{ talairachCenter, talairachRadius };
	const BE_IcoSphereSettings actualBrainIcoSettings{ actualBrainCenter, actualBrainRadius };
	const std::map<std::string, BE_IcoSphereSettings> bet2IcoSphereSettings{
		{"talairach", talairachIcoSettings},
		{ "actual_brain", actualBrainIcoSettings },
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
}


void SheetEvolverTest()
{
#define PERFORM_7PT_EXAMPLE false
	/**/
	constexpr pmp::Scalar roiHalfDim = 5.0;
	constexpr pmp::Scalar roiDim = 2.0 * roiHalfDim;

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

	constexpr pmp::Scalar cellSize = 0.1;
	constexpr pmp::Scalar columnWeight = 0.5;
#if PERFORM_7PT_EXAMPLE
	const auto gridBox = pmp::BoundingBox{ pmp::vec3{-5.0, -5.0, -roiHalfDim}, pmp::vec3{16.1, 15.0, roiHalfDim} };
	const auto grid = GetDistanceFieldWithSupportColumns(cellSize, gridBox, {
		{pmp::vec2{4.0, 6.0}, 0.5 * columnWeight},
		{pmp::vec2{0.0, 0.0}, 0.5 * columnWeight},
		{pmp::vec2{5.0, 0.0}, 0.5 * columnWeight},
		{pmp::vec2{11.1, 0.1}, 0.5 * columnWeight},
		{pmp::vec2{9.0, 2.0}, 0.5 * columnWeight},
		{pmp::vec2{7.0, 2.0}, 0.5 * columnWeight},
		{pmp::vec2{6.0, 10.0}, 0.5 * columnWeight}
		});
#else // 5-point example:
	const auto gridBox = pmp::BoundingBox{ pmp::vec3{0.0, 0.0, -roiHalfDim}, pmp::vec3{roiDim, roiDim, roiHalfDim} };
	const auto grid = GetDistanceFieldWithSupportColumns(cellSize, gridBox, {
		{pmp::vec2{2.5, 2.5}, 0.5 * columnWeight},
		{pmp::vec2{7.5, 2.5}, 0.5 * columnWeight},
		{pmp::vec2{7.5, 7.5}, 0.5 * columnWeight},
		{pmp::vec2{5.0, 8.0}, 0.5 * columnWeight},
		{pmp::vec2{2.5, 7.5}, 0.5 * columnWeight}
		});
#endif
	ExportToVTI(dataOutPath + "CapsuleVals", grid);

	const auto& sdfBox = grid.Box();
	const auto sdfBoxSize = sdfBox.max() - sdfBox.min();

	const double fieldIsoLevel = sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

	const pmp::Scalar startZHeight = sdfBox.min()[2] + 0.9 * sdfBoxSize[2];
	const pmp::Scalar endZHeight = sdfBox.min()[2] + 0.5 * sdfBoxSize[2];

	constexpr double tau = 0.02;

	const MeshTopologySettings topoSettings{
		true,
		0.45,
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
		0.05,
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