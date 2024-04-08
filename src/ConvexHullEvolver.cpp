#include "ConvexHullEvolver.h"

#include "geometry/GeometryConversionUtils.h"
#include "sdf/SDF.h"

/// \brief if true individual steps of surface evolution will be printed out into a given stream.
#define REPORT_EVOL_STEPS true // Note: may affect performance

/// \brief if true, upon computing trilinear system solution, new vertices are verified for belonging in the field bounds.
#define VERIFY_SOLUTION_WITHIN_BOUNDS false // Note: useful for detecting numerical explosions of the solution.

ConvexHullEvolver::ConvexHullEvolver(const std::vector<pmp::Point>& pointCloud, const ConvexHullSurfaceEvolutionSettings& settings)
    : m_PointCloud(pointCloud),
	  m_EvolSettings(settings)
{
	m_ImplicitLaplacianFunction =
		(m_EvolSettings.LaplacianType == MeshLaplacian::Barycentric ?
			pmp::laplace_implicit_barycentric : pmp::laplace_implicit_voronoi);
	m_LaplacianAreaFunction =
		(m_EvolSettings.LaplacianType == MeshLaplacian::Barycentric ?
			pmp::voronoi_area_barycentric : pmp::voronoi_area);
}

void ConvexHullEvolver::Evolve()
{
    Preprocess();
    // Evolution logic goes here
}

void ConvexHullEvolver::Preprocess()
{
    ConstructConvexHull();
    ComputeDistanceField();
    // Other preprocessing tasks
}

void ConvexHullEvolver::ConstructConvexHull()
{
    auto convexHullMeshOpt = Geometry::ComputePMPConvexHullFromPoints(m_PointCloud);
    if (!convexHullMeshOpt.has_value())
        throw std::logic_error("ConvexHullEvolver::ConstructConvexHull: m_PointCloud ComputePMPConvexHullFromPoints error! Terminating!\n");

    m_EvolvingSurface = std::make_shared<pmp::SurfaceMesh>(convexHullMeshOpt.value());
}

void ConvexHullEvolver::ComputeDistanceField()
{
	const pmp::BoundingBox ptCloudBBox(m_PointCloud);
	const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
	const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
	const float cellSize = minSize / m_EvolSettings.NVoxelsPerMinDimension;
	constexpr float volExpansionFactor = 1.0f;
	const SDF::PointCloudDistanceFieldSettings dfSettings{
				cellSize,
				volExpansionFactor,
				Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
				SDF::BlurPostprocessingType::None
	};

	m_Field = std::make_shared<Geometry::ScalarGrid>(
		SDF::PointCloudDistanceFieldGenerator::Generate(m_PointCloud, dfSettings));
}

void ReportInput(const ConvexHullSurfaceEvolutionSettings& evolSettings, std::ostream& os)
{
	os << "======================================================================\n";
	os << "> > > > > > > > > > Initiating ConvexHullEvolver: < < < < < < < < < < < <\n";
	os << "Target Name: " << evolSettings.ProcedureName << ",\n";
	os << "NSteps: " << evolSettings.NSteps << ",\n";
	os << "TimeStep: " << evolSettings.TimeStep << ",\n";
	os << "FieldIsoLevel: " << evolSettings.FieldIsoLevel << ",\n";
	os << "......................................................................\n";
	os << "NVoxelsPerMinDimension: " << evolSettings.NVoxelsPerMinDimension << ",\n";
	os << "......................................................................\n";
	const auto& c1 = evolSettings.ADParams.MCFMultiplier;
	const auto& c2 = evolSettings.ADParams.MCFVariance;
	os << "Curvature diffusion weight: " << c1 << " * (1 - exp(d^2 / " << c2 << ")),\n";
	const auto& d1 = evolSettings.ADParams.AdvectionMultiplier;
	const auto& d2 = evolSettings.ADParams.AdvectionSineMultiplier;
	os << "Advection weight: " << d1 << " * d * ((-grad(d) . N) - " << d2 << " * sqrt(1 - (grad(d) . N)^2)),\n";
	os << "......................................................................\n";
	os << "Min. Target Size: " << evolSettings.MinTargetSize << ",\n";
	os << "Max. Target Size: " << evolSettings.MaxTargetSize << ",\n";
	os << "Export Surface per Time Step: " << (evolSettings.ExportSurfacePerTimeStep ? "true" : "false") << ",\n";
	os << "Output Path: " << evolSettings.OutputPath << ",\n";
	os << "Do Remeshing: " << (evolSettings.DoRemeshing ? "true" : "false") << ",\n";
	os << "Do Feature Detection: " << (evolSettings.DoFeatureDetection ? "true" : "false") << ",\n";
	os << "----------------------------------------------------------------------\n";
}

float GetStabilizationScalingFactor(const double& timeStep, const pmp::SurfaceMesh& convexHullMesh, const float& stabilizationFactor)
{
	const auto surfaceArea = surface_area(convexHullMesh);
	const auto nVertices = convexHullMesh.n_vertices();
	return 1.0f;
}
