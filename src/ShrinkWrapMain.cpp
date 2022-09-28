
#include "pmp/SurfaceMesh.h"

#include "sdf/SDF.h"
#include "geometry/GridUtil.h"

#include "SurfaceEvolver.h"
#include "ConversionUtils.h"

#include <filesystem>
#include <chrono>
#include <map>

// set up root directory
const std::filesystem::path fsRootPath = DROOT_DIR;
const auto fsDataDirPath = fsRootPath / "data\\";
const auto fsDataOutPath = fsRootPath / "output\\";
const std::string dataDirPath = fsDataDirPath.string();
const std::string dataOutPath = fsDataOutPath.string();

constexpr bool performSDFTests = false;
constexpr bool performEvolverTests = true;

int main()
{
    // DISCLAIMER: the names need to match the models in "DROOT_DIR/data" except for the extension (which is always *.obj)
    const std::vector<std::string> meshNames{
        //"armadillo",
        //"BentChair",
        //"blub",
        "bunny",
        //"maxPlanck",
        //"nefertiti", /* <<<< unstable >>>> */
        //"ogre", /* <<<< unstable >>>> */
        //"spot"
    };

	if (performSDFTests)
	{
	    constexpr unsigned int nVoxelsPerMinDimension = 40;
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
		        0.2,
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
	    }		
	} // endif performSDFTests

	if (performEvolverTests)
	{
		constexpr unsigned int nVoxelsPerMinDimension = 40;

		const std::map<std::string, double> timeStepSizesForMeshes{
			{"armadillo", 0.05 },
			{"BentChair", 0.05 },
			{"blub", 0.01 },
			{"bunny", 0.002 },
			{"maxPlanck", 0.05 },
			{"nefertiti", 0.015 },
			{"ogre", 0.01 },
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

			const auto& sdfBox = sdf.Box();
			const auto sdfBoxSize = sdfBox.max() - sdfBox.min();
			const auto sdfBoxMaxDim = std::max<double>({ sdfBoxSize[0], sdfBoxSize[1], sdfBoxSize[2] });

			const double tau = timeStepSizesForMeshes.at(name); // time step
			SurfaceEvolutionSettings seSettings{
				name,
				80,
				tau,
				0.00 * minSize,
				3,
				PreComputeAdvectionDiffusionParams(0.5 * sdfBoxMaxDim, minSize),
				{},
				minSize, maxSize,
				meshBBox.center(),
				true, false,
				dataOutPath,
				true
			};
			ReportInput(seSettings, std::cout);
			SurfaceEvolver evolver(sdf, seSettings);

			try
			{
				evolver.Evolve();				
			}
			catch(...)
			{
				std::cerr << "> > > > > > > > > > > > > > SurfaceEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	} // endif performSDFTests
}