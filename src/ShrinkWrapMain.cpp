
#include "pmp/SurfaceMesh.h"
#include "sdf/SDF.h"
#include "ConversionUtils.h"

#include <chrono>


const std::string dataDirPath =
    "D:\\mcava\\Repos\\ImplicitSurfaceWrap\\data\\"; //"..\\..\\data\\";

// D:\\mcava\\Repos\\ImplicitSurfaceWrap\\data
// C:\\Users\\Martin\\source\\repos\\ImplicitSurfaceWrap

const std::string dataOutPath =
	"D:\\mcava\\Repos\\ImplicitSurfaceWrap\\output\\";

int main()
{
    pmp::SurfaceMesh mesh;
    mesh.read(dataDirPath + "bunny.obj");

    const auto startSDF = std::chrono::high_resolution_clock::now();
    constexpr float cellSize = 0.1207f / 40.0f;
    const auto sdf = SDF::ComputeDistanceField(mesh, cellSize);

    const auto endSDF = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> timeDiff = endSDF - startSDF;
    std::cout << "SDF Time: " << timeDiff.count() << "s\n";

    ExportToVTI(dataOutPath + "bunnySDF", sdf);
}