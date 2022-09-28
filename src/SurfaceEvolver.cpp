#include "SurfaceEvolver.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "ConversionUtils.h"
#include "pmp/algorithms/Remeshing.h"
#include "pmp/algorithms/Normals.h"
#include "pmp/algorithms/DifferentialGeometry.h"

#include "geometry/GridUtil.h"
#include "geometry/IcoSphereBuilder.h"
#include "pmp/algorithms/Decimation.h"
#include "pmp/algorithms/Features.h"

/// \brief a magic multiplier computing the radius of an ico-sphere that fits into the field's box.
constexpr float ICO_SPHERE_RADIUS_FACTOR = 0.4f;

constexpr unsigned int N_ICO_VERTS_0 = 12; // number of vertices in an icosahedron.
constexpr unsigned int N_ICO_EDGES_0 = 30; // number of edges in an icosahedron.

#define REPORT_EVOL_STEPS true // Note: may affect performance

/// \brief needs no explanation because it does not have one.
constexpr float UNEXPLAINABLE_MAGIC_CONSTANT = 1.0f;
// constexpr float UNEXPLAINABLE_MAGIC_CONSTANT = 0.93f;
// constexpr float UNEXPLAINABLE_MAGIC_CONSTANT = 8.0f;

/**
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
	const float weighedIcoRadius = icoRadius * UNEXPLAINABLE_MAGIC_CONSTANT;
	const float expectedMeanCoVolArea = stabilizationFactor * (4.0f * static_cast<float>(M_PI) * weighedIcoRadius * weighedIcoRadius / static_cast<float>(expectedVertexCount));
	return pow(static_cast<float>(timeStep) / expectedMeanCoVolArea, 1.0f / 3.0f);
}

void SurfaceEvolver::Preprocess()
{
	// build ico-sphere
	const float icoSphereRadius = ICO_SPHERE_RADIUS_FACTOR * (m_EvolSettings.MinTargetSize + 1.5f * m_EvolSettings.MaxTargetSize);
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
	(*m_Field) *= transfMatrixFull; // field needs to be moved to (0,0,0) and also scaled.
	(*m_Field) *= (static_cast<double>(scalingFactor)); // scale also distance values.

	// >>>>> Scaled geometry & field (use when debugging) <<<<<<
	//ExportToVTI(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + "_scaledField", *m_Field);
	//pmp::SurfaceMesh scaledTargetMesh = *m_EvolSettings.DO_NOT_KEEP_ptrTargetSurface;
	//scaledTargetMesh *= transfMatrixFull;
	//scaledTargetMesh.write(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + "_stableScale.obj");
}

// ================================================================================================

/// \brief identifier for sparse matrix.
using SparseMatrix = Eigen::SparseMatrix<double>;

/// \brief A utility for converting Eigen::ComputationInfo to a string message.
[[nodiscard]] std::string InterpretSolverErrorCode(const Eigen::ComputationInfo& cInfo)
{
	if (cInfo == Eigen::Success)
		return "Eigen::Success";

	if (cInfo == Eigen::NumericalIssue)
		return "Eigen::NumericalIssue";

	if (cInfo == Eigen::NoConvergence)
		return "Eigen::NoConvergence";

	return "Eigen::InvalidInput";
}

// ================================================================================================

/// \brief a stats wrapper for co-volume measures affecting the stability of the finite volume method.
struct CoVolumeStats
{
	double Mean{ 0.0 };
	double Max{ -DBL_MAX };
	double Min{ DBL_MAX };
};

const std::string coVolMeasureVertexPropertyName{ "v:coVolumeMeasure" };

/**
 * \brief Analyzes the stats of co-volumes around each mesh vertex, and creates a vertex property for the measure values.
 * \param mesh         input mesh.
 * \return co-volume stats.
 */
[[nodiscard]] CoVolumeStats AnalyzeMeshCoVolumes(pmp::SurfaceMesh& mesh)
{
	// vertex property for co-volume measures.
	if (!mesh.has_vertex_property(coVolMeasureVertexPropertyName))
		mesh.add_vertex_property<pmp::Scalar>(coVolMeasureVertexPropertyName);
	auto vCoVols = mesh.get_vertex_property<pmp::Scalar>(coVolMeasureVertexPropertyName);

	CoVolumeStats stats{};
	for (const auto v : mesh.vertices())
	{
		const double measure = pmp::voronoi_area(mesh, v);
		vCoVols[v] = static_cast<pmp::Scalar>(measure);

		if (stats.Max < measure) stats.Max = measure;
		if (stats.Min > measure) stats.Min = measure;
		stats.Mean += measure;
	}
	stats.Mean /= static_cast<double>(mesh.n_vertices());
	return stats;
}

// ================================================================================================

double SurfaceEvolver::LaplacianDistanceWeightFunction(const double& distanceAtVertex) const
{
	if (distanceAtVertex < 0.0 && m_EvolSettings.ADParams.MCFSupportPositive)
		return 0.0;
	const auto& c1 = m_EvolSettings.ADParams.MCFMultiplier;
	const auto& c2 = m_EvolSettings.ADParams.MCFVariance;
	return c1 * (1.0 - exp(-(distanceAtVertex * distanceAtVertex) / c2));
}

double SurfaceEvolver::AdvectionDistanceWeightFunction(const double& distanceAtVertex,
	const pmp::dvec3& negDistanceGradient, const pmp::Point& vertexNormal) const
{
	if (distanceAtVertex < 0.0 && m_EvolSettings.ADParams.AdvectionSupportPositive)
		return 0.0;
	const auto& d1 = m_EvolSettings.ADParams.AdvectionMultiplier;
	const auto& d2 = m_EvolSettings.ADParams.AdvectionSineMultiplier;
	const auto negGradDotNormal = pmp::ddot(negDistanceGradient, vertexNormal);
	return d1 * distanceAtVertex * (negGradDotNormal - d2 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
}

void SurfaceEvolver::ExportSurface(const unsigned int& tId, const bool& isResult, const bool& transformToOriginal) const
{
	const std::string connectingName = (isResult ? "_Result" : "_Evol_" + std::to_string(tId));
	if (!transformToOriginal)
	{
		m_EvolvingSurface->write(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + connectingName + m_OutputMeshExtension);
	    return;		
	}
	auto exportedSurface = *m_EvolvingSurface;
	exportedSurface *= m_TransformToOriginal;
	exportedSurface.write(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + connectingName + m_OutputMeshExtension);
}

// ================================================================================================


void SurfaceEvolver::Evolve()
{
	if (!m_Field)
		throw std::invalid_argument("SurfaceEvolver::Evolve: m_Field not set! Terminating!\n");
	if (!m_Field->IsValid())
		throw std::invalid_argument("SurfaceEvolver::Evolve: m_Field is invalid! Terminating!\n");

	Preprocess();

	if (!m_EvolvingSurface)
		throw std::invalid_argument("SurfaceEvolver::Evolve: m_EvolvingSurface not set! Terminating!\n");

	const auto& fieldBox = m_Field->Box();
	const auto fieldNegGradient = Geometry::ComputeNormalizedNegativeGradient(*m_Field);

	const auto& NSteps = m_EvolSettings.NSteps;
	const auto& tStep = m_EvolSettings.TimeStep;

	// ........ evaluate edge lengths for remeshing ....................
	const float phi = (1.0f + sqrt(5.0f)) / 2.0f; /// golden ratio.
	const auto subdiv = static_cast<float>(m_EvolSettings.IcoSphereSubdivisionLevel);
	const float r = m_StartingSurfaceRadius * m_ScalingFactor;
	const float minEdgeMultiplier = m_EvolSettings.TopoParams.MinEdgeMultiplier;
	auto minEdgeLength = minEdgeMultiplier * (2.0f * r / (sqrt(phi * sqrt(5.0f)) * subdiv)); // from icosahedron edge length
	auto maxEdgeLength = 4.0f * minEdgeLength;
#if REPORT_EVOL_STEPS
	std::cout << "minEdgeLength for remeshing: " << minEdgeLength << "\n";
#endif
	// .................................................................

	// DISCLAIMER: dhe dimensionality of the system depends on the number of mesh vertices which can change if remeshing is used.
	auto NVertices = static_cast<unsigned int>(m_EvolvingSurface->n_vertices());
	SparseMatrix sysMat(NVertices, NVertices);
	Eigen::MatrixXd sysRhs(NVertices, 3);
	auto vDistance = m_EvolvingSurface->add_vertex_property<pmp::Scalar>("v:distance"); // vertex property for distance field values.

	// property container for surface vertex normals
	pmp::VertexProperty<pmp::Point> vNormalsProp{};

	// ----------- System fill function --------------------------------
	const auto fillMatrixAndRHSTriplesFromMesh = [&]()
	{
		for (const auto v : m_EvolvingSurface->vertices())
		{
			const auto vPosToUpdate = m_EvolvingSurface->position(v);

			if (m_EvolvingSurface->is_boundary(v))
			{
				// freeze boundary/feature vertices
				const Eigen::Vector3d vertexRhs = vPosToUpdate;
				sysRhs.row(v.idx()) = vertexRhs;
				sysMat.coeffRef(v.idx(), v.idx()) = 1.0;
				continue;
			}

			const auto vNegGradDistanceToTarget = Geometry::TrilinearInterpolateVectorValue(vPosToUpdate, fieldNegGradient);

			const auto vNormal = vNormalsProp[v]; // vertex unit normal

			const double epsilonCtrlWeight = LaplacianDistanceWeightFunction(static_cast<double>(vDistance[v]) - m_EvolSettings.FieldIsoLevel);
			const double etaCtrlWeight = AdvectionDistanceWeightFunction(static_cast<double>(vDistance[v]) - m_EvolSettings.FieldIsoLevel, vNegGradDistanceToTarget, vNormal);


			const Eigen::Vector3d vertexRhs = vPosToUpdate + tStep * etaCtrlWeight * vNormal /* + tStep * tanRedistWeight * vTanVelocity */;
			sysRhs.row(v.idx()) = vertexRhs;

			const auto laplaceWeightInfo = pmp::laplace_implicit(*m_EvolvingSurface, v); // Laplacian weights
			sysMat.coeffRef(v.idx(), v.idx()) = 1.0 + tStep * epsilonCtrlWeight * static_cast<double>(laplaceWeightInfo.weightSum);

			for (const auto& [w, weight] : laplaceWeightInfo.vertexWeights)
			{
				sysMat.coeffRef(v.idx(), w.idx()) = -1.0 * tStep * epsilonCtrlWeight * static_cast<double>(weight);
			}
		}
	};
	// -----------------------------------------------------------------

	// write initial surface
	auto coVolStats = AnalyzeMeshCoVolumes(*m_EvolvingSurface);
#if REPORT_EVOL_STEPS
	std::cout << "Co-Volume Measure Stats: { Mean: " << coVolStats.Mean << ", Min: " << coVolStats.Min << ", Max: " << coVolStats.Max << "},\n";
#endif
	// set initial surface vertex properties
	for (const auto v : m_EvolvingSurface->vertices())
	{
		const auto vPos = m_EvolvingSurface->position(v);
		const double vDistanceToTarget = Geometry::TrilinearInterpolateScalarValue(vPos, *m_Field);
		vDistance[v] = static_cast<pmp::Scalar>(vDistanceToTarget);
	}
	if (m_EvolSettings.ExportSurfacePerTimeStep)
		ExportSurface(0);

	// main loop
	for (unsigned int ti = 1; ti <= NSteps; ti++)
	{
#if REPORT_EVOL_STEPS
		std::cout << "time step id: " << ti << "/" << NSteps << ", time: " << tStep * ti << "/" << tStep * NSteps << "\n";
		std::cout << "pmp::Normals::compute_vertex_normals ... ";
#endif
		pmp::Normals::compute_vertex_normals(*m_EvolvingSurface);
		vNormalsProp = m_EvolvingSurface->vertex_property<pmp::Point>("v:normal");
#if REPORT_EVOL_STEPS
		std::cout << "done\n";
		std::cout << "fillMatrixAndRHSTriplesFromMesh for " << NVertices << " vertices ... ";
#endif

		// prepare matrix & rhs
		fillMatrixAndRHSTriplesFromMesh();

#if REPORT_EVOL_STEPS
		std::cout << "done\n";
		std::cout << "Solving linear system ... ";
#endif
		// solve
		Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double>> solver(sysMat);
		Eigen::MatrixXd x = solver.solve(sysRhs);
		if (solver.info() != Eigen::Success)
		{
			const std::string msg = "\nSurfaceEvolver::Evolve: solver.info() != Eigen::Success for time step id: "
				+ std::to_string(ti) + ", Error code: " + InterpretSolverErrorCode(solver.info()) + "\n";
			std::cerr << msg;
			throw std::runtime_error(msg);
		}
#if REPORT_EVOL_STEPS
		std::cout << "done\n";
		std::cout << "Updating vertex positions ... ";
#endif

		// update vertex positions & verify mesh within bounds
		size_t nVertsOutOfBounds = 0;
		for (unsigned int i = 0; i < NVertices; i++)
		{
			const auto newPos = x.row(i);
			if (!fieldBox.Contains(newPos))
			{
				if (nVertsOutOfBounds == 0) std::cerr << "\n";
				std::cerr << "SurfaceEvolver::Evolve: vertex " << i << " out of field bounds!\n";
				nVertsOutOfBounds++;
			}
			m_EvolvingSurface->position(pmp::Vertex(i)) = x.row(i);
		}
		if (nVertsOutOfBounds > 0)
		{
			std::cerr << "SurfaceEvolver::Evolve: found " << nVertsOutOfBounds << " vertices out of bounds! Terminating!\n";
			break;
		}

		if (m_EvolSettings.DoRemeshing && ti > NSteps * m_EvolSettings.TopoParams.RemeshingStartTimeFactor)
		{
			// remeshing
#if REPORT_EVOL_STEPS
			std::cout << "done\n";
			std::cout << "pmp::Remeshing::adaptive_remeshing(minEdgeLength: " << minEdgeLength << ", maxEdgeLength: " << maxEdgeLength << ") ... ";
			//std::cout << "pmp::Remeshing::uniform_remeshing(targetEdgeLength: " << targetEdgeLength << ") ... ";
#endif
			if (m_EvolSettings.DoFeatureDetection && ti > NSteps * m_EvolSettings.TopoParams.FeatureDetectionStartTimeFactor)
			{
				// detect features
				pmp::Features feat(*m_EvolvingSurface);
				const auto minDihedralAngle = static_cast<pmp::Scalar>(m_EvolSettings.TopoParams.MinDihedralAngle);
				const auto maxDihedralAngle = static_cast<pmp::Scalar>(m_EvolSettings.TopoParams.MaxDihedralAngle);
				feat.detect_angle_within_bounds(minDihedralAngle, maxDihedralAngle);
			}
			pmp::Remeshing remeshing(*m_EvolvingSurface);
			remeshing.adaptive_remeshing(
				minEdgeLength, maxEdgeLength, 2.0f * minEdgeLength,
				m_EvolSettings.TopoParams.NRemeshingIters,
				m_EvolSettings.TopoParams.UseBackProjection
			);
			//remeshing.uniform_remeshing(targetEdgeLength);
#if REPORT_EVOL_STEPS
			std::cout << "done\n";
#endif
			if (ti % m_EvolSettings.TopoParams.StepStrideForEdgeDecay == 0 &&
				ti > NSteps * m_EvolSettings.TopoParams.RemeshingSizeDecayStartTimeFactor)
			{
				minEdgeLength *= 0.97f;
				maxEdgeLength *= 0.97f;
			}
		}

		coVolStats = AnalyzeMeshCoVolumes(*m_EvolvingSurface);
#if REPORT_EVOL_STEPS
		std::cout << "Co-Volume Measure Stats: { Mean: " << coVolStats.Mean << ", Min: " << coVolStats.Min << ", Max: " << coVolStats.Max << "},\n";
#endif
		// set surface vertex properties
		for (const auto v : m_EvolvingSurface->vertices())
		{
			const auto vPos = m_EvolvingSurface->position(v);
			const double vDistanceToTarget = Geometry::TrilinearInterpolateScalarValue(vPos, *m_Field);
			vDistance[v] = static_cast<pmp::Scalar>(vDistanceToTarget);
		}

		if (m_EvolSettings.ExportSurfacePerTimeStep)
			ExportSurface(ti);

		// update linear system dims for next time step:
		if (ti < NSteps && NVertices != m_EvolvingSurface->n_vertices())
		{
			NVertices = static_cast<unsigned int>(m_EvolvingSurface->n_vertices());
			sysMat = SparseMatrix(NVertices, NVertices);
			sysRhs = Eigen::MatrixXd(NVertices, 3);
		}

#if REPORT_EVOL_STEPS
		std::cout << ">>> Time step " << ti << " finished.\n";
		std::cout << "----------------------------------------------------------------------\n";
#endif

	} // end main loop

	if (m_EvolSettings.ExportResultSurface)
		ExportSurface(NSteps, true);
}

void ReportInput(const SurfaceEvolutionSettings& evolSettings, std::ostream& os)
{
	os << "======================================================================\n";
	os << "> > > > > > > > > > Initiating SurfaceEvolver: < < < < < < < < < < < <\n";
	os << "Target Name: " << evolSettings.ProcedureName << ",\n";
	os << "NSteps: " << evolSettings.NSteps << ",\n";
	os << "TimeStep: " << evolSettings.TimeStep << ",\n";
	os << "FieldIsoLevel: " << evolSettings.FieldIsoLevel << ",\n";
	os << "IcoSphereSubdivisionLevel: " << evolSettings.IcoSphereSubdivisionLevel << ",\n";
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
	os << "----------------------------------------------------------------------\n";
}

/// \brief a unit speed of a shrink-wrapping sphere sufficiently far away from target.
constexpr double BASE_DISTANCE_MULTIPLIER = 1.0;

/// \brief a factor by which distance field variance is multiplied.
constexpr double BASE_DISTANCE_VARIANCE = 1.0;

AdvectionDiffusionParameters PreComputeAdvectionDiffusionParams(const double& distanceMax, const double& targetMinDimension)
{
	// the radius within which target object starts slowing down shrink wrapping.
	const double distanceVariance = BASE_DISTANCE_VARIANCE; //* 0.125 * targetMinDimension * targetMinDimension;
	const double distanceMultiplier = BASE_DISTANCE_MULTIPLIER / (1.0 - exp(-distanceMax * distanceMax / distanceVariance));
	return { distanceMultiplier, distanceVariance, distanceMultiplier, 1.0 };
}