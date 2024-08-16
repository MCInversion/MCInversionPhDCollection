#include "ManifoldEvolver.h"

#include "sdf/SDF.h"

//
// ======================================================================================
//                    The strategy for 1D Curves in 2D space
// ---------------------------------------------------------------------------------------
//

void ManifoldCurveEvolutionStrategy::Preprocess()
{

}

void CustomManifoldCurveEvolutionStrategy::Preprocess()
{
	if (!GetOuterCurve() && GetInnerCurves().empty())
		throw std::invalid_argument("CustomManifoldCurveEvolutionStrategy::Preprocess: There's nothing to evolve!\n");

	ComputeAmbientFields();
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

void ManifoldCurveEvolutionStrategy::ComputeAmbientFields()
{
	if (!m_TargetPointCloud)
	{
		std::cerr << "ManifoldCurveEvolutionStrategy::ComputeAmbientFields: No m_TargetPointCloud found! Initializing empty fields: m_DistanceField and m_DFNegNormalizedGradient.\n";
		m_DistanceField = std::make_shared<Geometry::ScalarGrid2D>(1.0f, pmp::BoundingBox2{});
		m_DFNegNormalizedGradient = std::make_shared<Geometry::VectorGrid2D>(1.0f, pmp::BoundingBox2{});
		return;
	}

	const pmp::BoundingBox2 ptCloudBBox(*m_TargetPointCloud);
	const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
	const float minSize = std::min(ptCloudBBoxSize[0], ptCloudBBoxSize[1]);
	const float cellSize = minSize / static_cast<float>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
	const SDF::PointCloudDistanceField2DSettings dfSettings{
		cellSize,
		GetSettings().FieldSettings.FieldExpansionFactor,
		Geometry::DEFAULT_SCALAR_GRID_INIT_VAL
	};
	m_DistanceField = std::make_shared<Geometry::ScalarGrid2D>(
		SDF::PlanarPointCloudDistanceFieldGenerator::Generate(*m_TargetPointCloud, dfSettings));
}

//
// ======================================================================================
//                    The strategy for 2D Surfaces in 3D space
// ---------------------------------------------------------------------------------------
//

void ManifoldSurfaceEvolutionStrategy::Preprocess()
{
}

void CustomManifoldSurfaceEvolutionStrategy::Preprocess()
{
	if (!GetOuterSurface() && GetInnerSurfaces().empty())
		throw std::invalid_argument("CustomManifoldSurfaceEvolutionStrategy::Preprocess: There's nothing to evolve!\n");
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

void ManifoldSurfaceEvolutionStrategy::ComputeAmbientFields()
{
	if (!m_TargetPointCloud)
	{
		std::cerr << "ManifoldSurfaceEvolutionStrategy::ComputeAmbientFields: No m_TargetPointCloud found! Initializing empty fields: m_DistanceField and m_DFNegNormalizedGradient.\n";
		m_DistanceField = std::make_shared<Geometry::ScalarGrid>(1.0f, pmp::BoundingBox{});
		m_DFNegNormalizedGradient = std::make_shared<Geometry::VectorGrid>(1.0f, pmp::BoundingBox{});
		return;
	}

	const pmp::BoundingBox ptCloudBBox(*m_TargetPointCloud);
	const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
	const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
	const float cellSize = minSize / static_cast<float>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
	const SDF::PointCloudDistanceFieldSettings dfSettings{
		cellSize,
		GetSettings().FieldSettings.FieldExpansionFactor,
		Geometry::DEFAULT_SCALAR_GRID_INIT_VAL
	};
	m_DistanceField = std::make_shared<Geometry::ScalarGrid>(
		SDF::PointCloudDistanceFieldGenerator::Generate(*m_TargetPointCloud, dfSettings));
}

// ================================================
