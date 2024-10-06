# include <iostream>
# include <map>
# include <functional>
# include <string>
# include <vector>

# include "Experiments.h"  // Include the global registered experiments map

int main()
{
    // Available experiment identifiers (comment/uncomment as needed)
    const std::vector<std::string> enabledExperiments = {
// ------------------------------ OldEvolverInPMP --------------------------------
//     Late 2022, GRAPP 2023.
// ...............................................................................

        // "SDFTests",
        // "SphereTestOLD",
        // "EvolverTests",
        // "OldArmadilloLSWTest",
        // "IsosurfaceEvolverTests",
        // "SheetEvolverTest",
        // "BrainEvolverTests",

// ------------------------------- SubdivisionSSGG --------------------------------
//     Fall of 2023, SSGG 2023
// ...............................................................................
        // "SubdivisionTests1",
        // "SubdivisionTests2",
        // "SubdivisionTests3",
        // "SubdivisionTest4",
        // "SubdivTestsBoundary",
        // "SubdivTestsMultiTorus",
        // "SubdivPreallocationTests",
        // "NewIcosphereTests",
        // "IcospherePerformanceTests",
        // "CatmullClarkCounting",

// ---------------------------- ProjectOfDissertation ----------------------------
//     Spring of 2023
// ...............................................................................
        // "RemeshingTests",
        // "MobiusStripVoxelization",
        // "MetaballTest",
        // "ImportedObjMetricsEval",
        // "MMapImportTest",
        // "MMapOBJChunkMarkingTest",
        // "SimpleBunnyOBJSamplingDemo",

// -------------------------- RevisitingOldEvolverWSCG ----------------------------
//  Spring of 2024, WSCG 2024
// ................................................................................
        // "PDanielPtCloudPLYExport",
        // "PtCloudToDF",
        // "PDanielPtCloudComparisonTest",
        // "RepulsiveSurfResultEvaluation",
        // "HistogramResultEvaluation",
        // "OldResultJacobianMetricEval",
        // "HausdorffDistanceMeasurementsPerTimeStep",
        // "DirectHigherGenusPtCloudSampling",
        // "HigherGenusPtCloudLSW",
        // "TriTriIntersectionTests",
        // "MeshSelfIntersectionTests",
        // "HurtadoMeshesIsosurfaceEvolverTests",
        // "HurtadoTrexIcosphereLSW",
        // "ImportVTIDebugTests",

        // More focused LSW experiments were carried out to see whether they'll be applicable
		// for the IncrementalMeshBuilder. The result was that they were kind of hard to work with.

        // "ConvexHullTests",
        // "ConvexHullRemeshingTests",
        // "ConvexHullEvolverTests",
        // "IcoSphereEvolverTests",

// ------------------------- IncrementalMeshBuilderCESCG ----------------------------
// Spring of 2024. CESCG 2024 doctoral colloquium.
// ..................................................................................
        // "BPATest",
        // "IncrementalMeshBuilderTests",
        // "2GBApollonMeshBuilderTest",
        // "NanoflannDistanceTests",
        // "ApollonLSWSaliencyEval",
        // "IncrementalMeshBuilderHausdorffEval",
        // "ApollonArtecEvaLSWHausdorffEval",
        // "VertexNormalSampling",

// ----------------------------- TerrainGeneration ----------------------------------
//   Summer of 2024
// ..................................................................................
        // "TerrainPtGenerationTest",
        // "TerrainTriangulationTest",

// ------------------------- NewEvolverInnerOuterLSW ------------------------------
// Summer/Fall of 2024. MDPI Mathematics Journal.
// ..................................................................................
        // "DistanceFieldHashTest",
        // "DistanceField2DHashTest",
        // "ManifoldCurve2DTests",
        // "PointCloud2DSliceTests",
        // "CurveReorientTests",
        // "MeshReorientTests",
        // "PropertyPerformanceTests",
        // "OldVsNewLSWTests",
        // "PairedLSWTests",
        // "PairedLSWRepulsionTests",
        // "OutwardEvolvingInnerCircleTest",
        // "AdditionalCurveShapesTests",
        // "AdvectionDrivenInnerCircleTests",
        // "ConcentricCirclesTests",
        // "NonConcentricCirclesTest",
        // "EquilibriumPairedManifoldTests",
        // "VisualizeMCF",
        // "VisualizeMultipleInnerCurves",
        // "ExportSlicingPlanes",
        // "AdvectionDrivenInnerOuterCircleTests",
        // "OuterOnlySimpleShapeTests",
        "StandardMeshesIOLSWTests",
		// "MedialAxisTests",
		// "HyperellipseEllipsoidEquilibriumTests",
        // "JunkCan2DTests"

// ------------------------- NEW_EXPERIMENTS ------------------------------
// Allow new experiments here
// ........................................................................
    };

    if (enabledExperiments.empty())
    {
        std::cout << "main: No experiments have been chosen. Uncomment experiments in enabledExperiments vector.\n";
    }

    for (const auto& experimentId : enabledExperiments)
    {
        if (!REGISTERED_EXPERIMENTS.contains(experimentId))
        {
            std::cout << "main: Experiment \"" << experimentId << "\" not found in REGISTERED_EXPERIMENTS!" << "\n";
            continue;
        }
        std::cout << "\nvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n";
        std::cout << "main: Executing experiment \"" << experimentId << "\" ... \n";
        std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
        REGISTERED_EXPERIMENTS[experimentId]();  // Call the corresponding function
    }

    return 0;
}
