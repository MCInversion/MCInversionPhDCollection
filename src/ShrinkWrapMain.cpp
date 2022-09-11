
#include "pmp/SurfaceMesh.h"
#include "pmp/SurfaceMeshIO.h"

#include "sdfgen/makelevelset3.h"

#include "ConversionUtils.h"

#include <chrono>

const std::string dataDirPath =
    "C:\\Users\\Martin\\source\\repos\\ImplicitSurfaceWrap\\data\\"; //"..\\..\\data\\";

const std::string dataOutPath =
	"C:\\Users\\Martin\\source\\repos\\ImplicitSurfaceWrap\\output\\";

int main()
{
    pmp::SurfaceMesh mesh;
    mesh.read(dataDirPath + "bunny.obj");

    const auto startSDF = std::chrono::high_resolution_clock::now();

    const float cellSize = 0.1207f / 50.0f;
    auto sdfBBox = mesh.bounds();
    sdfBBox.expand(0.1207f, 0.1207f, 0.1207f);
    const auto newBoxSize = sdfBBox.max() - sdfBBox.min();

    const auto triMesh = GetTriangulatedSurfaceMesh(mesh);
    const auto triIds = GetTriangleIndices(triMesh);
    const auto vertPos = GetVertexPositions(triMesh);

    const int nx = std::floor(newBoxSize[0] / cellSize);
    const int ny = std::floor(newBoxSize[1] / cellSize);
    const int nz = std::floor(newBoxSize[2] / cellSize);

    const auto bboxMin = sdfBBox.min();
    const Vec3f origin(bboxMin[0], bboxMin[1], bboxMin[2]);

    Array3f distVals;
    make_level_set3(triIds, vertPos, origin, cellSize, nx, ny, nz, distVals);

    // TODO: Don't use make_level_set3. Recycle the one from SurfaceEvolver or utilize pmp::TriangleKdTree

    const auto endSDF = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> timeDiff = endSDF - startSDF;
    std::cout << "SDF Time: " << timeDiff.count() << "s\n";

    ExportToVTI(dataOutPath + "bunnySDF", origin, cellSize, nx, ny, nz, distVals);
}