
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
    pmp::SurfaceMesh mesh;
    mesh.read(dataDirPath + "bunny.obj");

    const auto startSDF = std::chrono::high_resolution_clock::now();
    constexpr float cellSize = 0.1207f / 40.0f;
    const auto sdf = SDF::ComputeDistanceField(mesh, { cellSize });

    const auto endSDF = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> timeDiff = endSDF - startSDF;
    std::cout << "SDF Time: " << timeDiff.count() << "s\n";

    ExportToVTI(dataOutPath + "bunnySDF", sdf);
}