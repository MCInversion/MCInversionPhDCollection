
#include "pmp/SurfaceMesh.h"

#include "sdf/SDF.h"

#include "ConversionUtils.h"

#include <filesystem>
#include <chrono>

// set up root directory
const std::filesystem::path fsRootPath = DROOT_DIR;
const auto fsDataDirPath = fsRootPath / "data\\";
const auto fsDataOutPath = fsRootPath / "output\\";
const std::string dataDirPath = fsDataDirPath.string();
const std::string dataOutPath = fsDataOutPath.string();

int main()
{
    constexpr float cellSize = 0.1207f / 40.0f;
    constexpr SDF::DistanceFieldSettings sdfSettings{
        cellSize,
        SDF::KDTreeSplitType::Center,
        1.0f,
        0.2,
        SDF::SignComputation::None,
        SDF::BlurPostprocessingType::None,
        SDF::PreprocessingType::Octree
    };

    pmp::SurfaceMesh mesh;
    mesh.read(dataDirPath + "bunny.obj");

    const auto startSDF = std::chrono::high_resolution_clock::now();
    const auto sdf = SDF::DistanceFieldGenerator::Generate(mesh, sdfSettings);

    const auto endSDF = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> timeDiff = endSDF - startSDF;
    std::cout << "SDF Time: " << timeDiff.count() << "s\n";

    ExportToVTI(dataOutPath + "bunnySDF", sdf);
}