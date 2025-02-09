# pragma once

# include <map>
# include <functional>
# include <string>

// ======================================================
// Declarations for all experiment functions
// ------------------------------------------------------

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
void SDFTests();
void SphereTestOLD();
void EvolverTests();
void OldArmadilloLSWTest();
void IsosurfaceEvolverTests();
void SheetEvolverTest();
void BrainEvolverTests();

// ------------------------------- SubdivisionSSGG --------------------------------
//     Fall of 2023, SSGG 2023
// ...............................................................................
/*
 * The origin of this output dates back to CESCG 2022 (Spring of 2022. See [Cavarga, 2022]).
 * The original evolution model had to be stabilized for general use by scaling evolving
 * meshes. The scale factor is computed from the ratio of the time step size and average
 * co-volume measure (area) which is estimated by assuming the distribution of mesh vertices
 * uniform on an icosphere. For this we need to know how many vertices there are. Of course,
 * we can initialize any mesh with N_V vertices and know directly, but my curiosity has
 * lead me to parametrize this for an icosphere constructed via 4-to-1 subdivision
 * using recurrence equations. After discovering the subdivision counting formula,
 * the local SSGG conference in September of 2023 was a suitable venue for this output
 * (see [Cavarga, 2023]).
 *
 * [Cavarga, 2022]
 * Cavarga, M. (2022). Advection-Driven Shrink-Wrapping of Triangulated Surfaces.
 * Proceedings of the 26th CESCG 2022. 95-104.
 *
 * [Cavarga, 2023]
 * Cavarga, M. (2023). 
 * Mesh primitive counting formulas for subdivision surfaces.
 * Proceedings of CSCG 2023. 9, 67-76.
 */
// --------------------------------------------------------------------------------
void SubdivisionTests1();
void SubdivisionTests2();
void SubdivisionTests3();
void SubdivisionTest4();
void SubdivTestsBoundary();
void SubdivTestsMultiTorus();
void SubdivPreallocationTests();
void NewIcosphereTests();
void IcospherePerformanceTests();
void CatmullClarkCounting();

// ---------------------------- ProjectOfDissertation ----------------------------
//     Spring of 2023
// ...............................................................................
/*
 * Some additional tests for verification and evaluation of some experimental ideas
 * written in my "project of dissertation" proposal. If you dive into these experiments,
 * you're about to see some results with voxelization and parallel mesh file sampling
 */
// --------------------------------------------------------------------------------
void RemeshingTests();
void MobiusStripVoxelization();
void MetaballTest();
void ImportedObjMetricsEval();
void MMapImportTest();
void MMapOBJChunkMarkingTest();
void SimpleBunnyOBJSamplingDemo();

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
void PDanielPtCloudPLYExport();
void PtCloudToDF();
void PDanielPtCloudComparisonTest();
void RepulsiveSurfResultEvaluation();
void HistogramResultEvaluation();
void OldResultJacobianMetricEval();
void HausdorffDistanceMeasurementsPerTimeStep();
void DirectHigherGenusPtCloudSampling();
void HigherGenusPtCloudLSW();
void TriTriIntersectionTests();
void MeshSelfIntersectionTests();
void HurtadoMeshesIsosurfaceEvolverTests();
void HurtadoTrexIcosphereLSW();
void ImportVTIDebugTests();

// More focused LSW experiments were carried out to see whether they'll be applicable
// for the IncrementalMeshBuilder. The result was that they were kind of hard to work with.
void ConvexHullTests();
void ConvexHullRemeshingTests();
void ConvexHullEvolverTests();
void IcoSphereEvolverTests();

// ------------------------- IncrementalMeshBuilderCESCG ----------------------------
// Spring of 2024. CESCG 2024 doctoral colloquium.
// ..................................................................................
/*
 * This was a chance to revisit an idea from the project of dissertation proposal about
 * progressively sampling very large mesh files. This lead to the development of the
 * IMB (Incremental Mesh Builder) framework. The results of these experiments were
 * presented on CESCG 2024 doctoral colloquium.
 */
 // --------------------------------------------------------------------------------
void BPATest();
void IncrementalMeshBuilderTests();
void TwoGBApollonMeshBuilderTest();
void NanoflannDistanceTests();
void ApollonLSWSaliencyEval();
void IncrementalMeshBuilderHausdorffEval();
void ApollonArtecEvaLSWHausdorffEval();
void VertexNormalSampling();

// ----------------------------- TerrainGeneration ----------------------------------
//   Summer of 2024
// ..................................................................................
/*
 * Experiments for constructing a triangulated terrain generated by Perlin noise.
 * Multiple libraries are compared, including the free version of commercial Fade library.
 */
 // --------------------------------------------------------------------------------
void TerrainPtGenerationTest();
void TerrainTriangulationTest();

// ------------------------- NewEvolverInnerOuterLSW ------------------------------
// Summer/Fall of 2024. MDPI Mathematics Journal.
// ..................................................................................
/*
 * This work focuses on extending the framework used while publishing [Cavarga, 2024].
 * Its core idea is to solve the stopping problem, i.e.: what is the right time to stop
 * Lagrangian shrink-wrapping (LSW)? Evolving meshes flowing into gaps in the target
 * point cloud prevent the use of a straight-forward heuristic of: if it doesn't move,
 * it's done! So we counterbalance LSW by introducing an inner manifold.
 *
 * DISCLAIMER: This is work-in-progress!
 *
 * [Cavarga, 2024]
 * Cavarga, M. (2024).
 * Automated Surface Extraction: Adaptive Remeshing Meets Lagrangian Shrink-Wrapping.
 * Journal of WSCG 2024. 32, 21-30.
 */
 // --------------------------------------------------------------------------------
void DistanceFieldHashTest();
void DistanceField2DHashTest();
void ManifoldCurve2DTests();
void PointCloud2DSliceTests();
void CurveReorientTests();
void MeshReorientTests();
void PropertyPerformanceTests();
void OldVsNewLSWTests();
void PairedLSWTests();
void PairedLSWRepulsionTests();
void OutwardEvolvingInnerCircleTest();
void AdditionalCurveShapesTests();
void AdvectionDrivenInnerCircleTests();
void ConcentricCirclesTests();
void NonConcentricCirclesTest();
void EquilibriumPairedManifoldTests();
void EquilibriumPairedConcaveManifoldTests();
void VisualizeMCF();
void VisualizeMultipleInnerCurves();
void ExportSlicingPlanes();
void AdvectionDrivenInnerOuterCircleTests();
void OuterOnlySimpleShapeTests();
void SimpleMeshesIOLSWTests();
void StandardMeshesIOLSWTests();
void MedialAxisTests();
void HyperellipseEllipsoidEquilibriumTests();
void JunkCan2DTests();
void InscribedCircleCalculatorVisualization();
void StandardMeshesExportWithNormals();
void ExportBallPivotingBoundaryEdges();
void TestProblematicMedialAxisPtClouds();
void TestDFDivergence2D();
void TestArcLengthCalculation();
void TestCurve2DRotation();
void TestSmoothingAdvectionEquilibrium();
void TestImageSegmentation();

//
// =================================================================================
//
// ------------------------ experiment identification ------------------------------
//

/// \brief global registered experiments map
inline std::map<std::string, std::function<void()>, std::less<>> REGISTERED_EXPERIMENTS{
// ------------------------------ OldEvolverInPMP --------------------------------
//     Late 2022, GRAPP 2023.
// ...............................................................................
    {"SDFTests", SDFTests},
    {"SphereTest", SphereTestOLD},
    {"EvolverTests", EvolverTests},
    {"OldArmadilloLSWTest", OldArmadilloLSWTest},
    {"IsosurfaceEvolverTests", IsosurfaceEvolverTests},
    {"SheetEvolverTest", SheetEvolverTest},
    {"BrainEvolverTests", BrainEvolverTests},

// ------------------------------- SubdivisionSSGG --------------------------------
//     Fall of 2023, SSGG 2023
// ...............................................................................
    {"SubdivisionTests1", SubdivisionTests1},
    {"SubdivisionTests2", SubdivisionTests2},
    {"SubdivisionTests3", SubdivisionTests3},
    {"SubdivisionTest4", SubdivisionTest4},
    {"SubdivTestsBoundary", SubdivTestsBoundary},
    {"SubdivTestsMultiTorus", SubdivTestsMultiTorus},
    {"SubdivPreallocationTests", SubdivPreallocationTests},
    {"NewIcosphereTests", NewIcosphereTests},
    {"IcospherePerformanceTests", IcospherePerformanceTests},
    {"CatmullClarkCounting", CatmullClarkCounting},

// ---------------------------- ProjectOfDissertation ----------------------------
//     Spring of 2023
// ...............................................................................
    {"RemeshingTests", RemeshingTests},
    {"MobiusStripVoxelization", MobiusStripVoxelization},
    {"MetaballTest", MetaballTest},
    {"ImportedObjMetricsEval", ImportedObjMetricsEval},
    {"MMapImportTest", MMapImportTest},
    {"MMapOBJChunkMarkingTest", MMapOBJChunkMarkingTest},
    {"SimpleBunnyOBJSamplingDemo", SimpleBunnyOBJSamplingDemo},

// -------------------------- RevisitingOldEvolverWSCG ----------------------------
//  Spring of 2024, WSCG 2024
// ................................................................................
    {"PDanielPtCloudPLYExport", PDanielPtCloudPLYExport},
    {"PtCloudToDF", PtCloudToDF},
    {"PDanielPtCloudComparisonTest", PDanielPtCloudComparisonTest},
    {"RepulsiveSurfResultEvaluation", RepulsiveSurfResultEvaluation},
    {"HistogramResultEvaluation", HistogramResultEvaluation},
    {"OldResultJacobianMetricEval", OldResultJacobianMetricEval},
    {"HausdorffDistanceMeasurementsPerTimeStep", HausdorffDistanceMeasurementsPerTimeStep},
    {"DirectHigherGenusPtCloudSampling", DirectHigherGenusPtCloudSampling},
    {"HigherGenusPtCloudLSW", HigherGenusPtCloudLSW},
    {"TriTriIntersectionTests", TriTriIntersectionTests},
    {"MeshSelfIntersectionTests", MeshSelfIntersectionTests},
    {"HurtadoMeshesIsosurfaceEvolverTests", HurtadoMeshesIsosurfaceEvolverTests},
    {"HurtadoTrexIcosphereLSW", HurtadoTrexIcosphereLSW},
    {"ImportVTIDebugTests", ImportVTIDebugTests},

    {"ConvexHullTests", ConvexHullTests},
    {"ConvexHullRemeshingTests", ConvexHullRemeshingTests},
    {"ConvexHullEvolverTests", ConvexHullEvolverTests},
    {"IcoSphereEvolverTests", IcoSphereEvolverTests},

// ------------------------- IncrementalMeshBuilderCESCG ----------------------------
// Spring of 2024. CESCG 2024 doctoral colloquium.
// ..................................................................................
    {"BPATest", BPATest},
    {"IncrementalMeshBuilderTests", IncrementalMeshBuilderTests},
    {"TwoGBApollonMeshBuilderTest", TwoGBApollonMeshBuilderTest},
    {"NanoflannDistanceTests", NanoflannDistanceTests},
    {"ApollonLSWSaliencyEval", ApollonLSWSaliencyEval},
    {"IncrementalMeshBuilderHausdorffEval", IncrementalMeshBuilderHausdorffEval},
    {"ApollonArtecEvaLSWHausdorffEval", ApollonArtecEvaLSWHausdorffEval},
    {"VertexNormalSampling", VertexNormalSampling},

// ----------------------------- TerrainGeneration ----------------------------------
//   Summer of 2024
// ..................................................................................
    {"TerrainPtGenerationTest", TerrainPtGenerationTest},
    {"TerrainTriangulationTest", TerrainTriangulationTest},

// ------------------------- NewEvolverInnerOuterLSW ------------------------------
// Summer/Fall of 2024. MDPI Mathematics Journal.
// ..................................................................................
    {"DistanceFieldHashTest", DistanceFieldHashTest},
    {"DistanceField2DHashTest", DistanceField2DHashTest},
    {"ManifoldCurve2DTests", ManifoldCurve2DTests},
    {"PointCloud2DSliceTests", PointCloud2DSliceTests},
    {"CurveReorientTests", CurveReorientTests},
    {"MeshReorientTests", MeshReorientTests},
    {"PropertyPerformanceTests", PropertyPerformanceTests},
    {"OldVsNewLSWTests", OldVsNewLSWTests},
    {"PairedLSWTests", PairedLSWTests},
    {"PairedLSWRepulsionTests", PairedLSWRepulsionTests},
    {"OutwardEvolvingInnerCircleTest", OutwardEvolvingInnerCircleTest},
    {"AdditionalCurveShapesTests", AdditionalCurveShapesTests},
    {"AdvectionDrivenInnerCircleTests", AdvectionDrivenInnerCircleTests},
    {"ConcentricCirclesTests", ConcentricCirclesTests},
    {"NonConcentricCirclesTest", NonConcentricCirclesTest},
    {"EquilibriumPairedManifoldTests", EquilibriumPairedManifoldTests},
    {"EquilibriumPairedConcaveManifoldTests", EquilibriumPairedConcaveManifoldTests},
    {"VisualizeMCF", VisualizeMCF},
    {"VisualizeMultipleInnerCurves", VisualizeMultipleInnerCurves},
    {"ExportSlicingPlanes", ExportSlicingPlanes},
    {"AdvectionDrivenInnerOuterCircleTests", AdvectionDrivenInnerOuterCircleTests},
    {"OuterOnlySimpleShapeTests", OuterOnlySimpleShapeTests },
    {"SimpleMeshesIOLSWTests", SimpleMeshesIOLSWTests },
    {"StandardMeshesIOLSWTests", StandardMeshesIOLSWTests },
    {"MedialAxisTests", MedialAxisTests },
    {"HyperellipseEllipsoidEquilibriumTests", HyperellipseEllipsoidEquilibriumTests },
    {"JunkCan2DTests", JunkCan2DTests },
    {"InscribedCircleCalculatorVisualization", InscribedCircleCalculatorVisualization },
    {"StandardMeshesExportWithNormals", StandardMeshesExportWithNormals },
    {"ExportBallPivotingBoundaryEdges", ExportBallPivotingBoundaryEdges },
    {"TestProblematicMedialAxisPtClouds", TestProblematicMedialAxisPtClouds },
    {"TestDFDivergence2D", TestDFDivergence2D },
    {"TestArcLengthCalculation", TestArcLengthCalculation },
    {"TestCurve2DRotation", TestCurve2DRotation },
    {"TestSmoothingAdvectionEquilibrium", TestSmoothingAdvectionEquilibrium },
    {"TestImageSegmentation", TestImageSegmentation },

// ------------------------- NEW_EXPERIMENTS ------------------------------
// Register new experiments here
// ........................................................................
};