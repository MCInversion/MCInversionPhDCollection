#include "BrainSurfaceEvolver.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "pmp/algorithms/Remeshing.h"
#include "pmp/algorithms/Normals.h"
#include "pmp/algorithms/Decimation.h"
#include "pmp/algorithms/Features.h"

#include "geometry/GridUtil.h"
#include "geometry/IcoSphereBuilder.h"
#include "geometry/MeshAnalysis.h"

// ================================================================================================

/// \brief a magic multiplier computing the radius of an ico-sphere that fits into the field's box.
constexpr float ICO_SPHERE_RADIUS_FACTOR = 0.4f;

/// \brief if true individual steps of surface evolution will be printed out into a given stream.
#define REPORT_EVOL_STEPS true // Note: may affect performance

/// \brief if true, upon computing trilinear system solution, new vertices are verified for belonging in the field bounds.
#define VERIFY_SOLUTION_WITHIN_BOUNDS false // Note: useful for detecting numerical explosions of the solution.

/// \brief repairs infinities and nan values of the input grid.
#define REPAIR_INPUT_GRID true //Note: very useful! You never know what field we use! Infinities and nans lead to the Eigen::NoConvergence error code whose cause is difficult to find!


BrainSurfaceEvolver::BrainSurfaceEvolver(const Geometry::ScalarGrid& field, const BrainExtractionSettings& settings)
	: m_EvolSettings(settings), m_Field(std::make_shared<Geometry::ScalarGrid>(field))
{
#if REPAIR_INPUT_GRID
	Geometry::RepairScalarGrid(*m_Field); // repair needed in case of invalid cell values.
#endif
	m_ImplicitLaplacianFunction =
		(m_EvolSettings.LaplacianType == BE_MeshLaplacian::Barycentric ?
			pmp::laplace_implicit_barycentric : pmp::laplace_implicit_voronoi);
	m_LaplacianAreaFunction =
		(m_EvolSettings.LaplacianType == BE_MeshLaplacian::Barycentric ?
			pmp::voronoi_area_barycentric : pmp::voronoi_area);
}

/*
 * \brief Computes scaling factor for stabilizing the finite volume method on assumed spherical surface meshes based on time step.
 * \param timeStep               time step size.
 * \param icoRadius              radius of an evolving geodesic icosahedron.
 * \param icoSubdiv              subdivision level of an evolving geodesic icosahedron.
 * \param stabilizationFactor    a multiplier for stabilizing mean co-volume area.
 * \return scaling factor for mesh and scalar grid.
 */
[[nodiscard]] float GetStabilizationScalingFactor(const double& timeStep, const float& icoRadius, const unsigned int& icoSubdiv, const float& stabilizationFactor = 1.0f)
{
	const unsigned int expectedVertexCount = (N_ICO_EDGES_0 * static_cast<unsigned int>(pow(4, icoSubdiv) - 1) + 3 * N_ICO_VERTS_0) / 3;
	const float weighedIcoRadius = icoRadius;
	const float expectedMeanCoVolArea = stabilizationFactor * (4.0f * static_cast<float>(M_PI) * weighedIcoRadius * weighedIcoRadius / static_cast<float>(expectedVertexCount));
	return pow(static_cast<float>(timeStep) / expectedMeanCoVolArea, 1.0f / 2.0f);
}


void BrainSurfaceEvolver::Preprocess()
{
	auto& field = *m_Field;
	const auto& fieldBox = field.Box();
	const auto fieldBoxSize = fieldBox.max() - fieldBox.min();
	const float minDim = std::min({ fieldBoxSize[0], fieldBoxSize[1], fieldBoxSize[2] });
	const float maxDim = std::max({ fieldBoxSize[0], fieldBoxSize[1], fieldBoxSize[2] });

	// build ico-sphere
	const float icoSphereRadius = ICO_SPHERE_RADIUS_FACTOR * (minDim + maxDim);
	m_StartingSurfaceRadius = icoSphereRadius;
#if REPORT_EVOL_STEPS
	std::cout << "Ico-Sphere Radius: " << icoSphereRadius << ",\n";
#endif
	const unsigned int icoSphereSubdiv = m_EvolSettings.IcoSphereSubdivisionLevel;
	Geometry::IcoSphereBuilder icoBuilder({ m_EvolSettings.IcoSphereSubdivisionLevel, icoSphereRadius });
	icoBuilder.BuildBaseData();
	icoBuilder.BuildPMPSurfaceMesh();
	m_EvolvingSurface = std::make_shared<pmp::SurfaceMesh>(icoBuilder.GetPMPSurfaceMeshResult());

	// transform mesh and grid
	// >>> uniform scale to ensure numerical method's stability.
	// >>> translation to origin for fields not centered at (0,0,0).
	// >>> scaling factor value is intended for stabilization of the numerical method.
	const float scalingFactor = GetStabilizationScalingFactor(m_EvolSettings.TimeStep, icoSphereRadius, icoSphereSubdiv);
	m_ScalingFactor = scalingFactor;

	// origin needs to be computed from bet2
	const auto origin = m_EvolSettings.TargetOrigin;


#if REPORT_EVOL_STEPS
	std::cout << "Stabilization Scaling Factor: " << scalingFactor << ",\n";
	std::cout << "Target Origin: {" << origin[0] << ", " << origin[1] << ", " << origin[2] << "},\n";
#endif
	const pmp::mat4 transfMatrixGeomScale{
		scalingFactor, 0.0f, 0.0f, 0.0f,
		0.0f, scalingFactor, 0.0f, 0.0f,
		0.0f, 0.0f, scalingFactor, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	const pmp::mat4 transfMatrixGeomMove{
		1.0f, 0.0f, 0.0f, -origin[0],
		0.0f, 1.0f, 0.0f, -origin[1],
		0.0f, 0.0f, 1.0f, -origin[2],
		0.0f, 0.0f, 0.0f, 1.0f
	};
	const auto transfMatrixFull = transfMatrixGeomScale * transfMatrixGeomMove;
	m_TransformToOriginal = inverse(transfMatrixFull);

	(*m_EvolvingSurface) *= transfMatrixGeomScale; // ico sphere is already centered at (0,0,0).
	field *= transfMatrixFull; // field needs to be moved to (0,0,0) and also scaled.
}

void BrainSurfaceEvolver::Evolve()
{
	if (!m_Field)
		throw std::invalid_argument("SurfaceEvolver::Evolve: m_Field not set! Terminating!\n");
	if (!m_Field->IsValid())
		throw std::invalid_argument("SurfaceEvolver::Evolve: m_Field is invalid! Terminating!\n");

	const auto& field = *m_Field;

	Preprocess();

	// do stuff



	if (!m_EvolvingSurface)
		throw std::invalid_argument("SurfaceEvolver::Evolve: m_EvolvingSurface not set! Terminating!\n");
}

void BrainSurfaceEvolver::ExportSurface(const unsigned int& tId, const bool& isResult, const bool& transformToOriginal) const
{
	const std::string connectingName = (isResult ? "_BE_Result" : "_BE_Evol_" + std::to_string(tId));
	if (!transformToOriginal)
	{
		m_EvolvingSurface->write(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + connectingName + m_OutputMeshExtension);
		return;
	}
	auto exportedSurface = *m_EvolvingSurface;
	exportedSurface *= m_TransformToOriginal;
	exportedSurface.write(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + connectingName + m_OutputMeshExtension);
}

void BrainSurfaceEvolver::ComputeTriangleMetrics() const
{
	for (const auto& metricName : m_EvolSettings.TriMetrics)
	{
		if (!Geometry::IsMetricRegistered(metricName))
			continue;

		const auto metricFunction = Geometry::IdentifyMetricFunction(metricName);

		if (!metricFunction(*m_EvolvingSurface))
		{
			std::cerr << "SurfaceEvolver::ComputeTriangleMetrics: [WARNING] Computation of metric " << metricName << " finished with errors!\n";
		}
	}
}
