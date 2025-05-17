// ------------------------- IncrementalMeshBuilderTUWien ------------------------------
// Late Spring / Summer of 2025. TU Wien June 6th Talk.
// ..................................................................................
/*
 * This talk builds on the previously introduced Incremental Mesh Builder (IMB) framework for progressive 
 * Level-of-Detail (LOD) mesh reconstruction from large 3D datasets. The current focus is on evaluating 
 * the viability of a geometry-aware stochastic vertex sampling process as a foundation for adaptive mesh 
 * reconstruction — balancing simplification speed with feature preservation.
 *
 * The latest developments test this sampling approach on both structured and noisy 3D datasets, using 
 * simplified artificial heightmaps and real-world scan-derived meshes. Performance results include 
 * comparisons between uniform and softmax-based sampling, along with LOD convergence times and visual 
 * fidelity metrics.
 *
 * The talk will also outline the ongoing effort to shape this framework into a lightweight, open-source 
 * application tailored to experts in 3D scanning, offering an efficient and scalable way to browse 
 * massive raw mesh files without full preprocessing.
 *
 * DISCLAIMER: This is work-in-progress!
 */
 // --------------------------------------------------------------------------------

#include "geometry/GeometryIOUtils.h"

#include "core/PointCloudMeshingStrategies.h"

#include "IOEnvironment.h"

#include "../Experiments.h"

void TestPoissonMeshingStrategy()
{
	const auto ptCloudName = "bunnyPts_3";
	const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataOutPath + ptCloudName + ".ply", true);
	if (!ptCloudOpt.has_value())
	{
		std::cerr << "ptCloudOpt == nullopt!\n";
		return;
	}
	const auto& pts = *ptCloudOpt;

	const auto meshingStrategy = IMB::GetReconstructionStrategy(IMB::ReconstructionFunctionType::Poisson);
	Geometry::BaseMeshGeometryData mesh;
	mesh.Vertices = pts;
	meshingStrategy->Process(mesh.Vertices, mesh.PolyIndices);

	const std::string outputFileName = dataOutPath + ptCloudName + "_PoissonRecon.vtk";
	if (!Geometry::ExportBaseMeshGeometryDataToVTK(mesh, outputFileName))
	{
		std::cout << "Failed to export mesh data." << "\n";
		return;
	}
}

void TestIMBShrinkWrapperNormalEstimation()
{
	const auto ptCloudName = "bunnyPts_3";
	const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataOutPath + ptCloudName + ".ply", true);
	if (!ptCloudOpt.has_value())
	{
		std::cerr << "ptCloudOpt == nullopt!\n";
		return;
	}
	const auto& pts = *ptCloudOpt;

	const auto meshingStrategy = IMB::GetReconstructionStrategy(IMB::ReconstructionFunctionType::LagrangianShrinkWrapping);
	Geometry::BaseMeshGeometryData mesh;
	mesh.Vertices = pts;
	meshingStrategy->Process(mesh.Vertices, mesh.PolyIndices);

	const std::string outputFileName = dataOutPath + ptCloudName + "_LSWRecon.vtk";
	if (!Geometry::ExportBaseMeshGeometryDataToVTK(mesh, outputFileName))
	{
		std::cout << "Failed to export mesh data." << "\n";
		return;
	}
}