#include "ManifoldEvolver.h"

#include "pmp/algorithms/CurveFactory.h"

#include "geometry/GridUtil.h"
#include "geometry/IcoSphereBuilder.h"
#include "geometry/MeshAnalysis.h"

#include "sdf/SDF.h"

#include "InscribedManifold.h"

/// \brief A factor by which the radius of any constructed outer/inner sphere is shrunken.
constexpr float SPHERE_RADIUS_FACTOR = 0.8f;

//
// ======================================================================================
//                    The strategy for 1D Curves in 2D space
// ---------------------------------------------------------------------------------------
//

void ManifoldCurveEvolutionStrategy::Preprocess(double timeStep)
{
	const auto [minTargetSize, maxTargetSize, targetCenter] = ComputeAmbientFields();
	const auto outerRadius = ConstructInitialManifolds(minTargetSize, maxTargetSize, targetCenter);
	StabilizeGeometries(timeStep, outerRadius);
}

void CustomManifoldCurveEvolutionStrategy::Preprocess(double timeStep)
{
	if (!GetOuterCurve() && GetInnerCurves().empty())
		throw std::invalid_argument("CustomManifoldCurveEvolutionStrategy::Preprocess: There's nothing to evolve!\n");
	
	if (!HasValidInnerOuterManifolds())
		throw std::invalid_argument("CustomManifoldCurveEvolutionStrategy::Preprocess: Invalid inner /outer manifold geometry! all custom inner curves are contained within the custom outer curve.\n");

	std::tie(std::ignore, std::ignore, std::ignore) = ComputeAmbientFields();

	const auto [minLength, maxLength] = CalculateCoVolumeRange();
	StabilizeCustomGeometries(timeStep, minLength, maxLength);
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

void ManifoldCurveEvolutionStrategy::ExportCurrentState(unsigned int step, const std::string& baseOutputFilename)
{
	const std::string connectingName = "_Evol_" + std::to_string(step);
	auto exportedOuterCurve = *m_OuterCurve;
	exportedOuterCurve *= m_TransformToOriginal;
	if (!write_to_ply(exportedOuterCurve, baseOutputFilename + "_Outer" + connectingName + ".ply"))
		std::cerr << "ManifoldCurveEvolutionStrategy::ExportCurrentState: error writing " << (baseOutputFilename + "_Outer" + connectingName + ".ply") << "!\n";

	for (size_t i = 0; const auto & innerCurve : m_InnerCurves)
	{
		auto exportedInnerCurve = *innerCurve;
		exportedInnerCurve *= m_TransformToOriginal;
		if (!write_to_ply(exportedInnerCurve, baseOutputFilename + "_Inner" + std::to_string(i++) + connectingName + ".ply"))
			std::cerr << "ManifoldCurveEvolutionStrategy::ExportCurrentState: error writing " << (baseOutputFilename + "_Inner" + std::to_string(i++) + connectingName + ".ply") << "!\n";
	}
}

void ManifoldCurveEvolutionStrategy::ExportFinalResult(const std::string& baseOutputFilename)
{
	const std::string connectingName = "_Evol_Result";
	auto exportedOuterCurve = *m_OuterCurve;
	exportedOuterCurve *= m_TransformToOriginal;
	if (!write_to_ply(exportedOuterCurve, baseOutputFilename + "_Outer" + connectingName + ".ply"))
		std::cerr << "ManifoldCurveEvolutionStrategy::ExportCurrentState: error writing " << (baseOutputFilename + "_Outer" + connectingName + ".ply") << "!\n";

	for (size_t i = 0; const auto & innerCurve : m_InnerCurves)
	{
		auto exportedInnerCurve = *innerCurve;
		exportedInnerCurve *= m_TransformToOriginal;
		if (!write_to_ply(exportedInnerCurve, baseOutputFilename + "_Inner" + std::to_string(i++) + connectingName + ".ply"))
			std::cerr << "ManifoldCurveEvolutionStrategy::ExportCurrentState: error writing " << (baseOutputFilename + "_Inner" + std::to_string(i++) + connectingName + ".ply") << "!\n";
	}
}

std::shared_ptr<pmp::ManifoldCurve2D> ManifoldCurveEvolutionStrategy::GetOuterCurveInOrigScale() const
{
	if (!m_OuterCurve)
		return nullptr;

	auto outerCurveTransformed = std::make_shared<pmp::ManifoldCurve2D>(*m_OuterCurve);
	*outerCurveTransformed *= m_TransformToOriginal;
	return outerCurveTransformed;
}

std::vector<std::shared_ptr<pmp::ManifoldCurve2D>> ManifoldCurveEvolutionStrategy::GetInnerCurvesInOrigScale() const
{
	std::vector<std::shared_ptr<pmp::ManifoldCurve2D>> innerCurvesTransformed{};
	innerCurvesTransformed.reserve(m_InnerCurves.size());
	for (const auto& innerCurve : m_InnerCurves)
	{
		if (!innerCurve)
			continue;

		auto innerCurveTransformed = std::make_shared<pmp::ManifoldCurve2D>(*innerCurve);
		*innerCurveTransformed *= m_TransformToOriginal;
		innerCurvesTransformed.push_back(innerCurveTransformed);
	}
	innerCurvesTransformed.shrink_to_fit();
	return innerCurvesTransformed;
}

// -------------------------------------------------------------------------------------

void ManifoldCurveEvolutionStrategy::SemiImplicitIntegrationStep(unsigned int step)
{
}

void ManifoldCurveEvolutionStrategy::ExplicitIntegrationStep(unsigned int step)
{
}

std::tuple<float, float, pmp::Point2> ManifoldCurveEvolutionStrategy::ComputeAmbientFields()
{
	if (!m_TargetPointCloud)
	{
		std::cerr << "ManifoldCurveEvolutionStrategy::ComputeAmbientFields: No m_TargetPointCloud found! Initializing empty fields: m_DistanceField and m_DFNegNormalizedGradient.\n";
		m_DistanceField = std::make_shared<Geometry::ScalarGrid2D>(1.0f, pmp::BoundingBox2{});
		m_DFNegNormalizedGradient = std::make_shared<Geometry::VectorGrid2D>(1.0f, pmp::BoundingBox2{});
		return { FLT_MAX, FLT_MAX, pmp::Point2(0, 0)};
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
	return { minSize, maxSize, ptCloudBBox.center() };
}

/// \brief the smallest allowed number of vertices in a manifold curve.
constexpr unsigned int N_CIRCLE_VERTS_0{ 5 };

float ManifoldCurveEvolutionStrategy::ConstructInitialManifolds(float minTargetSize, float maxTargetSize, const pmp::Point2& targetBoundsCenter)
{
	const float outerCircleRadius = 0.5f * SPHERE_RADIUS_FACTOR *
		(minTargetSize + (0.5f + GetSettings().FieldSettings.FieldExpansionFactor) * maxTargetSize);

	const auto nSegments = static_cast<unsigned int>(pow(2, GetSettings().LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
	m_OuterCurve = std::make_shared<pmp::ManifoldCurve2D>(pmp::CurveFactory::circle(targetBoundsCenter, outerCircleRadius, nSegments));

	if (!GetSettings().UseInnerManifolds || !m_TargetPointCloud || !m_DistanceField)
		return outerCircleRadius;

	const InscribedCircleInputData calcData{
		*m_TargetPointCloud,
		std::make_shared<Geometry::ScalarGrid2D>(*m_DistanceField) // clone
	};
	ParticleSwarmDistanceFieldInscribedCircleCalculator inscribedCircleCalculator;
	const auto circles = inscribedCircleCalculator.Calculate(calcData);
	m_InnerCurves.reserve(circles.size());

	for (const auto& circle : circles)
	{
		// keep the same vertex density for inner circles
		const auto nInnerSegments = static_cast<unsigned int>(static_cast<pmp::Scalar>(nSegments) * circle.Radius / outerCircleRadius);
		m_InnerCurves.emplace_back(std::make_shared<pmp::ManifoldCurve2D>(pmp::CurveFactory::circle(
			circle.Center,
			circle.Radius * SPHERE_RADIUS_FACTOR, 
			nInnerSegments
		)));
	}

	return outerCircleRadius;
}

/// \brief The power of the stabilizing scale factor.
constexpr float SCALE_FACTOR_POWER_1D = 1.0f;
/// \brief the reciprocal value of how many times the surface area element shrinks during evolution.
constexpr float INV_SHRINK_FACTOR_1D = 5.0f;

void ManifoldCurveEvolutionStrategy::StabilizeGeometries(double timeStep, float outerRadius, float stabilizationFactor)
{
	const auto expectedVertexCount = static_cast<unsigned int>(pow(2, GetSettings().LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
	const auto expectedMeanCoVolLength = stabilizationFactor * (2.0f * static_cast<float>(M_PI) * outerRadius / static_cast<float>(expectedVertexCount));
	const auto scalingFactor = pow(static_cast<float>(timeStep) / expectedMeanCoVolLength * INV_SHRINK_FACTOR_1D, SCALE_FACTOR_POWER_1D);
	GetScalingFactor() = scalingFactor;

	const pmp::mat3 transfMatrixGeomScale{
		scalingFactor, 0.0f, 0.0f,
		0.0f, scalingFactor, 0.0f,
		0.0f, 0.0f, 1.0f
	};
	m_TransformToOriginal = inverse(transfMatrixGeomScale);

	// transform geometries
	(*m_OuterCurve) *= transfMatrixGeomScale;
	(*m_DistanceField) *= transfMatrixGeomScale;
	(*m_DFNegNormalizedGradient) *= transfMatrixGeomScale;
	(*m_DistanceField) *= static_cast<double>(scalingFactor); // scale also the distance values.
}

bool CustomManifoldCurveEvolutionStrategy::HasValidInnerOuterManifolds() const
{
	// check self-intersections
	if (Geometry::PMPManifoldCurve2DHasSelfIntersections(*GetOuterCurve()))
		return false;

	// check inner curves
	for (const auto& innerCurve : GetInnerCurves())
	{
		// check self-intersections
		if (Geometry::PMPManifoldCurve2DHasSelfIntersections(*innerCurve))
			return false;

		// check whether the inner curve is contained inside the outer curve.
		const auto& outerCurve = *GetOuterCurve();
		for (const auto& p : innerCurve->positions())
		{
			if (!Geometry::IsPointInsidePMPManifoldCurve(p, outerCurve))
				return false;
		}
	}

	return true;
}

std::pair<float, float> CustomManifoldCurveEvolutionStrategy::CalculateCoVolumeRange() const
{
	float minCoVolLength = FLT_MAX;
	float maxCoVolLength = -FLT_MAX;

	if (GetOuterCurve())
	{
		const auto& outerCurve = *GetOuterCurve();
		for (const auto v : outerCurve.vertices())
		{
			const auto [eTo, eFrom] = outerCurve.edges(v);
			const auto currentCoVolLength = 0.5f * (outerCurve.edge_length(eTo) + outerCurve.edge_length(eFrom));

			if (currentCoVolLength > maxCoVolLength) maxCoVolLength = currentCoVolLength;
			if (currentCoVolLength < minCoVolLength) minCoVolLength = currentCoVolLength;
		}
	}

	for (const auto& curve : GetInnerCurves())
	{
		const auto& innerCurve = *curve;
		for (const auto v : innerCurve.vertices())
		{
			const auto [eTo, eFrom] = innerCurve.edges(v);
			const auto currentCoVolLength = 0.5f * (innerCurve.edge_length(eTo) + innerCurve.edge_length(eFrom));

			if (currentCoVolLength > maxCoVolLength) maxCoVolLength = currentCoVolLength;
			if (currentCoVolLength < minCoVolLength) minCoVolLength = currentCoVolLength;
		}
	}

	return { minCoVolLength, maxCoVolLength };
}

void CustomManifoldCurveEvolutionStrategy::StabilizeCustomGeometries(double timeStep, float minLength, float maxLength, float stabilizationFactor)
{
	const float expectedMeanCoVolLength = (1.0f - stabilizationFactor) * minLength + stabilizationFactor * maxLength;
	const float scalingFactor = pow(static_cast<float>(timeStep) / expectedMeanCoVolLength * INV_SHRINK_FACTOR_1D, SCALE_FACTOR_POWER_1D);
	GetScalingFactor() = scalingFactor;
	const pmp::mat3 transfMatrixGeomScale{
		scalingFactor, 0.0f, 0.0f,
		0.0f, scalingFactor, 0.0f,
		0.0f, 0.0f, 1.0f
	};
	GetTransformToOriginal() = inverse(transfMatrixGeomScale);

	// Transform the geometries
	(*GetOuterCurve()) *= transfMatrixGeomScale;
	for (auto& innerCurve : GetInnerCurves())
	{
		(*innerCurve) *= transfMatrixGeomScale;
	}
	(*GetDistanceField()) *= transfMatrixGeomScale;
	(*GetDFNegNormalizedGradient()) *= transfMatrixGeomScale;
	(*GetDistanceField()) *= static_cast<double>(scalingFactor); // Scale the distance values
}

//
// ======================================================================================
//                    The strategy for 2D Surfaces in 3D space
// ---------------------------------------------------------------------------------------
//

void ManifoldSurfaceEvolutionStrategy::Preprocess(double timeStep)
{
	const auto [minTargetSize, maxTargetSize, targetCenter] = ComputeAmbientFields();
	const auto outerRadius = ConstructInitialManifolds(minTargetSize, maxTargetSize, targetCenter);
	StabilizeGeometries(timeStep, outerRadius);
}

void CustomManifoldSurfaceEvolutionStrategy::Preprocess(double timeStep)
{
	if (!GetOuterSurface() && GetInnerSurfaces().empty())
		throw std::invalid_argument("CustomManifoldSurfaceEvolutionStrategy::Preprocess: There's nothing to evolve!\n");

	if (!HasValidInnerOuterManifolds())
		throw std::invalid_argument("CustomManifoldSurfaceEvolutionStrategy::Preprocess: Invalid inner /outer manifold geometry! all custom inner curves are contained within the custom outer curve.\n");

	std::tie(std::ignore, std::ignore, std::ignore) = ComputeAmbientFields();

	const auto [minArea, maxArea] = CalculateCoVolumeRange();
	StabilizeCustomGeometries(timeStep, minArea, maxArea);
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

void ManifoldSurfaceEvolutionStrategy::ExportCurrentState(unsigned int step, const std::string& baseOutputFilename)
{
	const std::string connectingName = "_Evol_" + std::to_string(step);
	auto exportedOuterSurface = *m_OuterSurface;
	exportedOuterSurface *= m_TransformToOriginal;
	exportedOuterSurface.write(baseOutputFilename + "_Outer" + connectingName + ".vtk");

	for (size_t i = 0; const auto& innerSurface : m_InnerSurfaces)
	{
		auto exportedInnerSurface = *innerSurface;
		exportedInnerSurface *= m_TransformToOriginal;
		exportedInnerSurface.write(baseOutputFilename + "_Inner" + std::to_string(i++) + connectingName + ".vtk");
	}
}

void ManifoldSurfaceEvolutionStrategy::ExportFinalResult(const std::string& baseOutputFilename)
{
	const std::string connectingName = "_Evol_Result";
	auto exportedOuterSurface = *m_OuterSurface;
	exportedOuterSurface *= m_TransformToOriginal;
	exportedOuterSurface.write(baseOutputFilename + "_Outer" + connectingName + ".vtk");

	for (size_t i = 0; const auto & innerSurface : m_InnerSurfaces)
	{
		auto exportedInnerSurface = *innerSurface;
		exportedInnerSurface *= m_TransformToOriginal;
		exportedInnerSurface.write(baseOutputFilename + "_Inner" + std::to_string(i++) + connectingName + ".vtk");
	}
}

std::shared_ptr<pmp::SurfaceMesh> ManifoldSurfaceEvolutionStrategy::GetOuterSurfaceInOrigScale() const
{
	if (!m_OuterSurface)
		return nullptr;

	auto outerSurfaceTransformed = std::make_shared<pmp::SurfaceMesh>(*m_OuterSurface);
	*outerSurfaceTransformed *= m_TransformToOriginal;
	return outerSurfaceTransformed;
}

std::vector<std::shared_ptr<pmp::SurfaceMesh>> ManifoldSurfaceEvolutionStrategy::GetInnerSurfacesInOrigScale() const
{
	std::vector<std::shared_ptr<pmp::SurfaceMesh>> innerSurfacesTransformed{};
	innerSurfacesTransformed.reserve(m_InnerSurfaces.size());
	for (const auto& innerSurface : m_InnerSurfaces)
	{
		if (!innerSurface)
			continue;

		auto innerCurveTransformed = std::make_shared<pmp::SurfaceMesh>(*innerSurface);
		*innerCurveTransformed *= m_TransformToOriginal;
		innerSurfacesTransformed.push_back(innerCurveTransformed);
	}
	innerSurfacesTransformed.shrink_to_fit();
	return innerSurfacesTransformed;
}

// ------------------------------------------------

void ManifoldSurfaceEvolutionStrategy::SemiImplicitIntegrationStep(unsigned int step)
{
}

void ManifoldSurfaceEvolutionStrategy::ExplicitIntegrationStep(unsigned int step)
{
}

std::tuple<float, float, pmp::Point> ManifoldSurfaceEvolutionStrategy::ComputeAmbientFields()
{
	if (!m_TargetPointCloud)
	{
		std::cerr << "ManifoldSurfaceEvolutionStrategy::ComputeAmbientFields: No m_TargetPointCloud found! Initializing empty fields: m_DistanceField and m_DFNegNormalizedGradient.\n";
		m_DistanceField = std::make_shared<Geometry::ScalarGrid>(1.0f, pmp::BoundingBox{});
		m_DFNegNormalizedGradient = std::make_shared<Geometry::VectorGrid>(1.0f, pmp::BoundingBox{});
		return { FLT_MAX, FLT_MAX, pmp::Point(0, 0, 0)};
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
	return { minSize, maxSize, ptCloudBBox.center() };
}

float ManifoldSurfaceEvolutionStrategy::ConstructInitialManifolds(float minTargetSize, float maxTargetSize, const pmp::Point& targetBoundsCenter)
{
	const float outerSphereRadius = 0.5f * SPHERE_RADIUS_FACTOR *
		(minTargetSize + (0.5f + GetSettings().FieldSettings.FieldExpansionFactor) * maxTargetSize);

	Geometry::IcoSphereBuilder icoBuilder({ GetSettings().LevelOfDetail, outerSphereRadius });
	icoBuilder.BuildBaseData();
	icoBuilder.BuildPMPSurfaceMesh();
	m_OuterSurface = std::make_shared<pmp::SurfaceMesh>(icoBuilder.GetPMPSurfaceMeshResult());
	const pmp::mat4 transfMatrixGeomMove{
		1.0f, 0.0f, 0.0f, -targetBoundsCenter[0],
		0.0f, 1.0f, 0.0f, -targetBoundsCenter[1],
		0.0f, 0.0f, 1.0f, -targetBoundsCenter[2],
		0.0f, 0.0f, 0.0f, 1.0f
	};
	(*m_OuterSurface) *= transfMatrixGeomMove; // center to target bounds

	if (!GetSettings().UseInnerManifolds || !m_TargetPointCloud || !m_DistanceField)
		return outerSphereRadius;

	// TODO: implement 3D version

	return outerSphereRadius;
}

/// \brief The power of the stabilizing scale factor.
constexpr float SCALE_FACTOR_POWER_2D = 1.0f / 2.0f;
/// \brief the reciprocal value of how many times the surface area element shrinks during evolution.
constexpr float INV_SHRINK_FACTOR_2D = 5.0f;

void ManifoldSurfaceEvolutionStrategy::StabilizeGeometries(double timeStep, float outerRadius, float stabilizationFactor)
{
	const auto expectedVertexCount = static_cast<unsigned int>(pow(2, GetSettings().LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
	const auto expectedMeanCoVolLength = stabilizationFactor * (2.0f * static_cast<float>(M_PI) * outerRadius / static_cast<float>(expectedVertexCount));
	const auto scalingFactor = pow(static_cast<float>(timeStep) / expectedMeanCoVolLength * INV_SHRINK_FACTOR_2D, SCALE_FACTOR_POWER_2D);
	GetScalingFactor() = scalingFactor;

	const pmp::mat4 transfMatrixGeomScale{
		scalingFactor, 0.0f, 0.0f, 0.0f,
		0.0f, scalingFactor, 0.0f, 0.0f,
		0.0f, 0.0f, scalingFactor, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	m_TransformToOriginal = inverse(transfMatrixGeomScale);

	// transform geometries
	(*m_OuterSurface) *= transfMatrixGeomScale;
	(*m_DistanceField) *= transfMatrixGeomScale;
	(*m_DFNegNormalizedGradient) *= transfMatrixGeomScale;
	(*m_DistanceField) *= static_cast<double>(scalingFactor); // scale also the distance values.
}

bool CustomManifoldSurfaceEvolutionStrategy::HasValidInnerOuterManifolds() const
{
	// check self-intersections
	if (Geometry::PMPSurfaceMeshHasSelfIntersections(*GetOuterSurface()))
		return false;

	// check inner surfaces
	for (const auto& innerSurface : GetInnerSurfaces())
	{
		// check self-intersections
		if (Geometry::PMPSurfaceMeshHasSelfIntersections(*innerSurface))
			return false;

		// check whether the inner surface is contained inside the outer surface.
		const auto& outerSurface = *GetOuterSurface();
		const Geometry::PMPSurfaceMeshAdapter outerSurfaceAdapter(std::make_shared<pmp::SurfaceMesh>(outerSurface));
		const auto outerSurfaceKdTree = std::make_shared<Geometry::CollisionKdTree>(outerSurfaceAdapter, Geometry::CenterSplitFunction);
		for (const auto& p : innerSurface->positions())
		{
			if (!IsPointInsidePMPSurfaceMesh(p, outerSurfaceKdTree))
				return false;
		}
	}

	return true;
}

std::pair<float, float> CustomManifoldSurfaceEvolutionStrategy::CalculateCoVolumeRange() const
{
	float minCoVolArea = FLT_MAX;
	float maxCoVolArea = -FLT_MAX;

	if (GetOuterSurface())
	{
		const auto& outerSurface = *GetOuterSurface();
		for (const auto v : outerSurface.vertices())
		{
			const auto currentCoVolArea = GetLaplacianAreaFunction()(outerSurface, v);

			if (currentCoVolArea > maxCoVolArea) maxCoVolArea = currentCoVolArea;
			if (currentCoVolArea < minCoVolArea) minCoVolArea = currentCoVolArea;
		}
	}

	for (const auto& surface : GetInnerSurfaces())
	{
		const auto& innerSurface = *surface;
		for (const auto v : innerSurface.vertices())
		{
			const auto currentCoVolArea = GetLaplacianAreaFunction()(innerSurface, v);

			if (currentCoVolArea > maxCoVolArea) maxCoVolArea = currentCoVolArea;
			if (currentCoVolArea < minCoVolArea) minCoVolArea = currentCoVolArea;
		}
	}

	return { minCoVolArea, maxCoVolArea };
}

void CustomManifoldSurfaceEvolutionStrategy::StabilizeCustomGeometries(double timeStep, float minLength, float maxLength, float stabilizationFactor)
{
	const float expectedMeanCoVolLength = (1.0f - stabilizationFactor) * minLength + stabilizationFactor * maxLength;
	const float scalingFactor = pow(static_cast<float>(timeStep) / expectedMeanCoVolLength * INV_SHRINK_FACTOR_1D, SCALE_FACTOR_POWER_1D);
	GetScalingFactor() = scalingFactor;
	const pmp::mat4 transfMatrixGeomScale{
		scalingFactor, 0.0f, 0.0f, 0.0f,
		0.0f, scalingFactor, 0.0f, 0.0f,
		0.0f, 0.0f, scalingFactor, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	GetTransformToOriginal() = inverse(transfMatrixGeomScale);

	// Transform the geometries
	(*GetOuterSurface()) *= transfMatrixGeomScale;
	for (auto& innerSurface : GetInnerSurfaces())
	{
		(*innerSurface) *= transfMatrixGeomScale;
	}
	(*GetDistanceField()) *= transfMatrixGeomScale;
	(*GetDFNegNormalizedGradient()) *= transfMatrixGeomScale;
	(*GetDistanceField()) *= static_cast<double>(scalingFactor); // Scale the distance values
}

// ================================================
