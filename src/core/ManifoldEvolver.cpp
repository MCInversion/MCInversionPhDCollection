#include "ManifoldEvolver.h"

#include "InscribedManifold.h"
#include "geometry/GridUtil.h"
#include "geometry/IcoSphereBuilder.h"
#include "pmp/algorithms/CurveFactory.h"
#include "pmp/algorithms/SurfaceFactory.h"
#include "sdf/SDF.h"

/// \brief A factor by which the radius of any constructed outer/inner sphere is shrunken.
constexpr float SPHERE_RADIUS_FACTOR = 0.8f;

//
// ======================================================================================
//                    The strategy for 1D Curves in 2D space
// ---------------------------------------------------------------------------------------
//

void ManifoldCurveEvolutionStrategy::Preprocess(double timeStep)
{
	const auto [minTargetSize, maxTargetSize] = ComputeAmbientFields();

	ConstructInitialManifolds(minTargetSize, maxTargetSize);
}

void CustomManifoldCurveEvolutionStrategy::Preprocess(double timeStep)
{
	if (!GetOuterCurve() && GetInnerCurves().empty())
		throw std::invalid_argument("CustomManifoldCurveEvolutionStrategy::Preprocess: There's nothing to evolve!\n");

	std::tie(std::ignore, std::ignore) = ComputeAmbientFields();
}

void ManifoldCurveEvolutionStrategy::PerformEvolutionStep(unsigned int step)
{
}

bool ManifoldCurveEvolutionStrategy::ShouldRemesh()
{
	return false;
}

void ManifoldCurveEvolutionStrategy::Remesh()
{
}

void ManifoldCurveEvolutionStrategy::ExportCurrentState(unsigned int step)
{
}

void ManifoldCurveEvolutionStrategy::ExportFinalResult()
{
}

// -------------------------------------------------------------------------------------

std::pair<float, float> ManifoldCurveEvolutionStrategy::ComputeAmbientFields()
{
	if (!m_TargetPointCloud)
	{
		std::cerr << "ManifoldCurveEvolutionStrategy::ComputeAmbientFields: No m_TargetPointCloud found! Initializing empty fields: m_DistanceField and m_DFNegNormalizedGradient.\n";
		m_DistanceField = std::make_shared<Geometry::ScalarGrid2D>(1.0f, pmp::BoundingBox2{});
		m_DFNegNormalizedGradient = std::make_shared<Geometry::VectorGrid2D>(1.0f, pmp::BoundingBox2{});
		return {};
	}

	const pmp::BoundingBox2 ptCloudBBox(*m_TargetPointCloud);
	const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
	const float minSize = std::min(ptCloudBBoxSize[0], ptCloudBBoxSize[1]);
	const float maxSize = std::max(ptCloudBBoxSize[0], ptCloudBBoxSize[1]);
	const float cellSize = minSize / static_cast<float>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
	const SDF::PointCloudDistanceField2DSettings dfSettings{
		cellSize,
		GetSettings().FieldSettings.FieldExpansionFactor,
		Geometry::DEFAULT_SCALAR_GRID_INIT_VAL
	};
	m_DistanceField = std::make_shared<Geometry::ScalarGrid2D>(
		SDF::PlanarPointCloudDistanceFieldGenerator::Generate(*m_TargetPointCloud, dfSettings));
	m_DFNegNormalizedGradient = std::make_shared<Geometry::VectorGrid2D>(ComputeNormalizedNegativeGradient(*m_DistanceField));
	return { minSize, maxSize };
}

void ManifoldCurveEvolutionStrategy::ConstructInitialManifolds(float minTargetSize, float maxTargetSize)
{
	const float circleRadius = 0.5f * SPHERE_RADIUS_FACTOR *
		(minTargetSize + (0.5f + GetSettings().FieldSettings.FieldExpansionFactor) * maxTargetSize);

	const auto nSegments = 5 * (1 + 2 * GetSettings().LevelOfDetail);
	m_OuterCurve = std::make_shared<pmp::ManifoldCurve2D>(pmp::CurveFactory::circle(pmp::Point2(0, 0), circleRadius, nSegments));

	if (!GetSettings().UseInnerManifolds || !m_TargetPointCloud || !m_DistanceField)
		return;

	const InscribedCircleInputData calcData{
		*m_TargetPointCloud,
		std::make_shared<Geometry::ScalarGrid2D>(*m_DistanceField) // clone
	};
	ParticleSwarmDistanceFieldInscribedCircleCalculator inscribedCircleCalculator;
	const auto circles = inscribedCircleCalculator.Calculate(calcData);

	for (const auto& circle : circles)
	{
		m_InnerCurves.emplace_back(std::make_shared<pmp::ManifoldCurve2D>(pmp::CurveFactory::circle(
			circle.Center,
			circle.Radius * SPHERE_RADIUS_FACTOR, 
			nSegments
		)));
	}
}

//
// ======================================================================================
//                    The strategy for 2D Surfaces in 3D space
// ---------------------------------------------------------------------------------------
//

void ManifoldSurfaceEvolutionStrategy::Preprocess(double timeStep)
{
	const auto [minTargetSize, maxTargetSize] = ComputeAmbientFields();

	ConstructInitialManifolds(minTargetSize, maxTargetSize);
}

void CustomManifoldSurfaceEvolutionStrategy::Preprocess(double timeStep)
{
	if (!GetOuterSurface() && GetInnerSurfaces().empty())
		throw std::invalid_argument("CustomManifoldSurfaceEvolutionStrategy::Preprocess: There's nothing to evolve!\n");

	std::tie(std::ignore, std::ignore) = ComputeAmbientFields();
}

void ManifoldSurfaceEvolutionStrategy::PerformEvolutionStep(unsigned int step)
{
}

bool ManifoldSurfaceEvolutionStrategy::ShouldRemesh()
{
	return false;
}

void ManifoldSurfaceEvolutionStrategy::Remesh()
{
}

void ManifoldSurfaceEvolutionStrategy::ExportCurrentState(unsigned int step)
{
}

void ManifoldSurfaceEvolutionStrategy::ExportFinalResult()
{
}

// ------------------------------------------------

std::pair<float, float> ManifoldSurfaceEvolutionStrategy::ComputeAmbientFields()
{
	if (!m_TargetPointCloud)
	{
		std::cerr << "ManifoldSurfaceEvolutionStrategy::ComputeAmbientFields: No m_TargetPointCloud found! Initializing empty fields: m_DistanceField and m_DFNegNormalizedGradient.\n";
		m_DistanceField = std::make_shared<Geometry::ScalarGrid>(1.0f, pmp::BoundingBox{});
		m_DFNegNormalizedGradient = std::make_shared<Geometry::VectorGrid>(1.0f, pmp::BoundingBox{});
		return {};
	}

	const pmp::BoundingBox ptCloudBBox(*m_TargetPointCloud);
	const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
	const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
	const float maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
	const float cellSize = minSize / static_cast<float>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
	const SDF::PointCloudDistanceFieldSettings dfSettings{
		cellSize,
		GetSettings().FieldSettings.FieldExpansionFactor,
		Geometry::DEFAULT_SCALAR_GRID_INIT_VAL
	};
	m_DistanceField = std::make_shared<Geometry::ScalarGrid>(
		SDF::PointCloudDistanceFieldGenerator::Generate(*m_TargetPointCloud, dfSettings));
	m_DFNegNormalizedGradient = std::make_shared<Geometry::VectorGrid>(ComputeNormalizedNegativeGradient(*m_DistanceField));
	return { minSize, maxSize };
}

void ManifoldSurfaceEvolutionStrategy::ConstructInitialManifolds(float minTargetSize, float maxTargetSize)
{
	const float icoSphereRadius = 0.5f * SPHERE_RADIUS_FACTOR *
		(minTargetSize + (0.5f + GetSettings().FieldSettings.FieldExpansionFactor) * maxTargetSize);

	Geometry::IcoSphereBuilder icoBuilder({ GetSettings().LevelOfDetail, icoSphereRadius });
	icoBuilder.BuildBaseData();
	icoBuilder.BuildPMPSurfaceMesh();
	m_OuterSurface = std::make_shared<pmp::SurfaceMesh>(icoBuilder.GetPMPSurfaceMeshResult());

	if (!GetSettings().UseInnerManifolds || !m_TargetPointCloud || !m_DistanceField)
		return;

	// TODO: implement 3D version
}

// ================================================
