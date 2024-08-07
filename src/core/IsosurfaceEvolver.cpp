#include "IsosurfaceEvolver.h"

#include "pmp/algorithms/Remeshing.h"
#include "pmp/algorithms/Normals.h"
#include "pmp/algorithms/Decimation.h"
#include "pmp/algorithms/Features.h"

#include "geometry/GridUtil.h"
#include "geometry/MeshAnalysis.h"
#include "geometry/MarchingCubes.h"

#include <fstream>

#include "ConversionUtils.h"
#include "geometry/GeometryConversionUtils.h"

// ================================================================================================

/// \brief if true individual steps of surface evolution will be printed out into a given stream.
#define REPORT_EVOL_STEPS true // Note: may affect performance

/// \brief if true, upon computing trilinear system solution, new vertices are verified for belonging in the field bounds.
#define VERIFY_SOLUTION_WITHIN_BOUNDS false // Note: useful for detecting numerical explosions of the solution.

/// \brief repairs infinities and nan values of the input grid.
#define REPAIR_INPUT_GRID true //Note: very useful! You never know what field we use! Infinities and nans lead to the Eigen::NoConvergence error code whose cause is difficult to find!

// ================================================================================================

size_t IsoSurfaceEvolver::DetectFeatures(const FeatureDetectionType& type) const
{
	pmp::Features featuresDetector(*m_EvolvingSurface);
	if (type == FeatureDetectionType::Angle)
	{
		const auto maxDihedralAngle = static_cast<pmp::Scalar>(m_EvolSettings.TopoParams.MaxDihedralAngle);
		return featuresDetector.detect_angle(maxDihedralAngle);
	}

	if (type == FeatureDetectionType::AngleWithinBounds)
	{
		const auto minDihedralAngle = static_cast<pmp::Scalar>(m_EvolSettings.TopoParams.MinDihedralAngle);
		const auto maxDihedralAngle = static_cast<pmp::Scalar>(m_EvolSettings.TopoParams.MaxDihedralAngle);
		return featuresDetector.detect_angle_within_bounds(minDihedralAngle, maxDihedralAngle);
	}

	if (type == FeatureDetectionType::PrincipalCurvatures)
	{
		const auto curvatureFactor = m_EvolSettings.TopoParams.PrincipalCurvatureFactor;
		const bool exclude = m_EvolSettings.TopoParams.ExcludeEdgesWithoutBothFeaturePts;
		return featuresDetector.detect_vertices_with_curvatures_imbalance(curvatureFactor, exclude);
	}

	const auto curvatureAngle = m_EvolSettings.TopoParams.CriticalMeanCurvatureAngle;
	const auto curvatureFactor = m_EvolSettings.TopoParams.PrincipalCurvatureFactor;
	const bool exclude = m_EvolSettings.TopoParams.ExcludeEdgesWithoutBothFeaturePts;
	return featuresDetector.detect_vertices_with_high_curvature(curvatureAngle, curvatureFactor, exclude);
}

// ================================================================================================

IsoSurfaceEvolver::IsoSurfaceEvolver(const Geometry::ScalarGrid& field, const float& fieldExpansionFactor, const IsoSurfaceEvolutionSettings& settings)
	: m_EvolSettings(settings), m_Field(std::make_shared<Geometry::ScalarGrid>(field)), m_ExpansionFactor(fieldExpansionFactor)
{
#if REPAIR_INPUT_GRID
	Geometry::RepairScalarGrid(*m_Field); // repair needed in case of invalid cell values.
#endif
	m_ImplicitLaplacianFunction =
		(m_EvolSettings.LaplacianType == MeshLaplacian::Barycentric ?
			pmp::laplace_implicit_barycentric : pmp::laplace_implicit_voronoi);
	m_LaplacianAreaFunction =
		(m_EvolSettings.LaplacianType == MeshLaplacian::Barycentric ?
			pmp::voronoi_area_barycentric : pmp::voronoi_area);
}

// ================================================================================================

void IsoSurfaceEvolver::Preprocess()
{
	auto& field = *m_Field;

	// build isosurface
	const double isoLevel = m_EvolSettings.FieldIsoLevelOffset;
	const auto cellSize = m_EvolSettings.ReSampledGridCellSize;
	const auto reSampledField = ExtractReSampledGrid(cellSize, field);

	const auto& dim = reSampledField.Dimensions();
	const auto Nx = static_cast<unsigned int>(dim.Nx);
	const auto Ny = static_cast<unsigned int>(dim.Ny);
	const auto Nz = static_cast<unsigned int>(dim.Nz);
	const auto mcMesh = MarchingCubes::GetMarchingCubesMesh<double>(reSampledField.Values().data(), Nx, Ny, Nz, isoLevel);
	m_EvolvingSurface = std::make_shared<pmp::SurfaceMesh>(Geometry::ConvertMCMeshToPMPSurfaceMesh(mcMesh));

	// Ilatsik's Marching cubes implementation outputs a mesh from a grid with unit cell size, and zero origin (0,0,0)
	const auto& orig = reSampledField.Box().min();
	const pmp::mat4 voxelTransformMat{
		cellSize, 0.0f, 0.0f, orig[0],
		0.0f, cellSize, 0.0f, orig[1],
		0.0f, 0.0f, cellSize, orig[2],
		0.0f, 0.0f, 0.0f, 1.0f
	};
	(*m_EvolvingSurface) *= voxelTransformMat;

	// basic 1-iter remesh for bad quality mesh from marching cubes
	const float minEdgeLength =static_cast<float>(M_SQRT2) * cellSize * m_EvolSettings.TopoParams.MinEdgeMultiplier;
	const float maxEdgeLength = 4.0f * minEdgeLength;
	const float approxError = 0.25f * (minEdgeLength + maxEdgeLength);
	pmp::Remeshing remeshing(*m_EvolvingSurface);
	remeshing.adaptive_remeshing({
		minEdgeLength, maxEdgeLength, approxError,
		1,
		m_EvolSettings.TopoParams.NTanSmoothingIters,
		true });

	// transform mesh and grid
	// >>> uniform scale to ensure numerical method's stability.
	// >>> translation to origin for fields not centered at (0,0,0).
	// >>> scaling factor value is intended for stabilization of the numerical method.
	const float scalingFactor = GetStabilizationScalingFactor(m_EvolSettings.TimeStep, cellSize);
	m_ScalingFactor = scalingFactor;
	m_EvolSettings.FieldIsoLevel *= static_cast<double>(scalingFactor);
#if REPORT_EVOL_STEPS
	std::cout << "Stabilization Scaling Factor: " << scalingFactor << ",\n";
#endif
	const pmp::mat4 transfMatrixGeomScale{
		scalingFactor, 0.0f, 0.0f, 0.0f,
		0.0f, scalingFactor, 0.0f, 0.0f,
		0.0f, 0.0f, scalingFactor, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	m_TransformToOriginal = inverse(transfMatrixGeomScale);

	(*m_EvolvingSurface) *= transfMatrixGeomScale; // ico sphere is already centered at (0,0,0).
	field *= transfMatrixGeomScale; // field needs to be moved to (0,0,0) and also scaled.
	field *= static_cast<double>(scalingFactor); // scale also distance values.

	// >>>>> Scaled geometry & field (use when debugging) <<<<<<
	//ExportToVTI(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + "_resampledField", reSampledField);
	//pmp::SurfaceMesh scaledTargetMesh = *m_EvolvingSurface;
	//scaledTargetMesh *= transfMatrixFull;
	//scaledTargetMesh.write(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + "_stableScale.obj");

}

// ================================================================================================

double IsoSurfaceEvolver::LaplacianDistanceWeightFunction(const double& distanceAtVertex) const
{
	if (distanceAtVertex < 0.0 && m_EvolSettings.ADParams.MCFSupportPositive)
		return 0.0;
	const auto& c1 = m_EvolSettings.ADParams.MCFMultiplier;
	const auto& c2 = m_EvolSettings.ADParams.MCFVariance;
	return c1 * (1.0 - exp(-(distanceAtVertex * distanceAtVertex) / c2));
}

double IsoSurfaceEvolver::AdvectionDistanceWeightFunction(const double& distanceAtVertex,
	const pmp::dvec3& negDistanceGradient, const pmp::Point& vertexNormal) const
{
	if (distanceAtVertex < 0.0 && m_EvolSettings.ADParams.AdvectionSupportPositive)
		return 0.0;
	const auto& d1 = m_EvolSettings.ADParams.AdvectionMultiplier;
	const auto& d2 = m_EvolSettings.ADParams.AdvectionSineMultiplier;
	const auto negGradDotNormal = pmp::ddot(negDistanceGradient, vertexNormal);
	return d1 * distanceAtVertex * (negGradDotNormal - d2 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
}

void IsoSurfaceEvolver::ExportSurface(const unsigned int& tId, const bool& isResult, const bool& transformToOriginal) const
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

void IsoSurfaceEvolver::ComputeTriangleMetrics() const
{
	for (const auto& metricName : m_EvolSettings.TriMetrics)
	{
		if (!Geometry::IsMetricRegistered(metricName))
			continue;

		const auto metricFunction = Geometry::IdentifyMetricFunction(metricName);

		if (!metricFunction(*m_EvolvingSurface))
		{
			std::cerr << "IsoSurfaceEvolver::ComputeTriangleMetrics: [WARNING] Computation of metric " << metricName << " finished with errors!\n";
		}
	}
}

// ================================================================================================

void IsoSurfaceEvolver::Evolve()
{
	if (!m_Field)
		throw std::invalid_argument("IsoSurfaceEvolver::Evolve: m_Field not set! Terminating!\n");
	if (!m_Field->IsValid())
		throw std::invalid_argument("IsoSurfaceEvolver::Evolve: m_Field is invalid! Terminating!\n");

	const auto& field = *m_Field;

	Preprocess();

	if (!m_EvolvingSurface)
		throw std::invalid_argument("IsoSurfaceEvolver::Evolve: m_EvolvingSurface not set! Terminating!\n");

#if VERIFY_SOLUTION_WITHIN_BOUNDS
	const auto& fieldBox = field.Box();
#endif
	const auto fieldNegGradient = Geometry::ComputeNormalizedNegativeGradient(field);

	const auto& NSteps = m_EvolSettings.NSteps;
	const auto& tStep = m_EvolSettings.TimeStep;

	// ........ evaluate edge lengths for remeshing ....................
	const float cellSize = m_EvolSettings.ReSampledGridCellSize * m_ScalingFactor;
	float minEdgeLength = static_cast<float>(M_SQRT2) * cellSize * m_EvolSettings.TopoParams.MinEdgeMultiplier;
	float maxEdgeLength = 4.0f * minEdgeLength;
	float approxError = 0.25f * (minEdgeLength + maxEdgeLength);
	//auto approxError = 2.0f * minEdgeLength;
#if REPORT_EVOL_STEPS
	std::cout << "minEdgeLength for remeshing: " << minEdgeLength << "\n";
#endif
	// .................................................................

	// DISCLAIMER: the dimensionality of the system depends on the number of mesh vertices which can change if remeshing is used.
	auto NVertices = static_cast<unsigned int>(m_EvolvingSurface->n_vertices());
	SparseMatrix sysMat(NVertices, NVertices);
	Eigen::MatrixXd sysRhs(NVertices, 3);
	auto vDistance = m_EvolvingSurface->add_vertex_property<pmp::Scalar>("v:distance"); // vertex property for distance field values.
	auto vFeature = m_EvolvingSurface->vertex_property<bool>("v:feature", false);

	// property container for surface vertex normals
	pmp::VertexProperty<pmp::Point> vNormalsProp{};

	// ----------- System fill function --------------------------------
	const auto fillMatrixAndRHSTriplesFromMesh = [&]()
	{
		for (const auto v : m_EvolvingSurface->vertices())
		{
			const auto vPosToUpdate = m_EvolvingSurface->position(v);

			if ((m_EvolSettings.IdentityForBoundaryVertices && m_EvolvingSurface->is_boundary(v)) ||
				(m_EvolSettings.IdentityForFeatureVertices && vFeature[v]))
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

			const Eigen::Vector3d vertexRhs = vPosToUpdate + tStep * etaCtrlWeight * vNormal;
			sysRhs.row(v.idx()) = vertexRhs;
			const float tanRedistWeight = m_EvolSettings.TangentialVelocityWeight * epsilonCtrlWeight;
			if (tanRedistWeight > 0.0f)
			{
				// compute tangential velocity
				const auto vTanVelocity = ComputeTangentialUpdateVelocityAtVertex(*m_EvolvingSurface, v, vNormal, tanRedistWeight);
				sysRhs.row(v.idx()) += tStep * Eigen::Vector3d(vTanVelocity);
			}

			const auto laplaceWeightInfo = m_ImplicitLaplacianFunction(*m_EvolvingSurface, v); // Laplacian weights
			sysMat.coeffRef(v.idx(), v.idx()) = 1.0 + tStep * epsilonCtrlWeight * static_cast<double>(laplaceWeightInfo.weightSum);

			for (const auto& [w, weight] : laplaceWeightInfo.vertexWeights)
			{
				sysMat.coeffRef(v.idx(), w.idx()) = -1.0 * tStep * epsilonCtrlWeight * static_cast<double>(weight);
			}
		}
	};
	// -----------------------------------------------------------------

	// write initial surface
	auto coVolStats = AnalyzeMeshCoVolumes(*m_EvolvingSurface, m_LaplacianAreaFunction);
#if REPORT_EVOL_STEPS
	std::ofstream fileOStreamMins(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + "_CoVolMins.txt");
	std::ofstream fileOStreamMeans(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + "_CoVolMeans.txt");
	std::ofstream fileOStreamMaxes(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + "_CoVolMaxes.txt");
	std::cout << "Co-Volume Measure Stats: { Mean: " << coVolStats.Mean << ", Min: " << coVolStats.Min << ", Max: " << coVolStats.Max << "},\n";
	fileOStreamMins << coVolStats.Min << ", ";
	fileOStreamMeans << coVolStats.Mean << ", ";
	fileOStreamMaxes << coVolStats.Max << ", ";
#endif
	// set initial surface vertex properties
	for (const auto v : m_EvolvingSurface->vertices())
	{
		const auto vPos = m_EvolvingSurface->position(v);
		const double vDistanceToTarget = Geometry::TrilinearInterpolateScalarValue(vPos, field);
		vDistance[v] = static_cast<pmp::Scalar>(vDistanceToTarget);
	}
	Geometry::ComputeEdgeDihedralAngles(*m_EvolvingSurface);
	Geometry::ComputeVertexCurvaturesAndRelatedProperties(*m_EvolvingSurface);
	ComputeTriangleMetrics();
	if (m_EvolSettings.ExportSurfacePerTimeStep)
		ExportSurface(0);

	// -------------------------------------------------------------------------------------------------------------
	// ........................................ main loop ..........................................................
	// -------------------------------------------------------------------------------------------------------------
	for (unsigned int ti = 1; ti <= NSteps; ti++)
	{
#if REPORT_EVOL_STEPS
		std::cout << "time step id: " << ti << "/" << NSteps << ", time: " << tStep * ti << "/" << tStep * NSteps
			<< ", Procedure Name: " << m_EvolSettings.ProcedureName << "\n";
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
			const std::string msg = "\nIsoSurfaceEvolver::Evolve: solver.info() != Eigen::Success for time step id: "
				+ std::to_string(ti) + ", Error code: " + InterpretSolverErrorCode(solver.info()) + "\n";
			std::cerr << msg;
			throw std::runtime_error(msg);
		}
#if REPORT_EVOL_STEPS
		std::cout << "done\n";
		std::cout << "Updating vertex positions ... ";
#endif

		// update vertex positions & verify mesh within bounds
#if VERIFY_SOLUTION_WITHIN_BOUNDS
		size_t nVertsOutOfBounds = 0;
#endif
		for (unsigned int i = 0; i < NVertices; i++)
		{
			const auto newPos = x.row(i);
#if VERIFY_SOLUTION_WITHIN_BOUNDS
			if (!fieldBox.Contains(newPos))
			{
				if (nVertsOutOfBounds == 0) std::cerr << "\n";
				std::cerr << "IsoSurfaceEvolver::Evolve: vertex " << i << " out of field bounds!\n";
				nVertsOutOfBounds++;
			}
#endif
			m_EvolvingSurface->position(pmp::Vertex(i)) = x.row(i);
		}
#if VERIFY_SOLUTION_WITHIN_BOUNDS
		if ((m_EvolSettings.DoRemeshing && nVertsOutOfBounds > static_cast<double>(NVertices) * m_EvolSettings.MaxFractionOfVerticesOutOfBounds) ||
			(!m_EvolSettings.DoRemeshing && nVertsOutOfBounds > 0))
		{
			std::cerr << "IsoSurfaceEvolver::Evolve: found " << nVertsOutOfBounds << " vertices out of bounds! Terminating!\n";
			break;
		}
#endif
#if REPORT_EVOL_STEPS
		std::cout << "done\n";
#endif

		if (m_EvolSettings.DoRemeshing /* && ti > NSteps * m_EvolSettings.TopoParams.RemeshingStartTimeFactor*/)
		{
			// remeshing
#if REPORT_EVOL_STEPS
			std::cout << "Detecting Features ...";
#endif
			if (m_EvolSettings.DoFeatureDetection && ti > NSteps * m_EvolSettings.TopoParams.FeatureDetectionStartTimeFactor)
			{
				const auto nEdges = DetectFeatures(m_EvolSettings.TopoParams.FeatureType);
#if REPORT_EVOL_STEPS
				std::cout << "done. " << nEdges << " feature edges detected.\n";
#endif
			}
#if REPORT_EVOL_STEPS
			std::cout << "pmp::Remeshing::adaptive_remeshing(minEdgeLength: " << minEdgeLength << ", maxEdgeLength: " << maxEdgeLength << ") ... ";
#endif
			//std::cout << "pmp::Remeshing::uniform_remeshing(targetEdgeLength: " << targetEdgeLength << ") ... ";
			pmp::Remeshing remeshing(*m_EvolvingSurface);
			remeshing.adaptive_remeshing({
				minEdgeLength, maxEdgeLength, approxError,
				m_EvolSettings.TopoParams.NRemeshingIters,
				m_EvolSettings.TopoParams.NTanSmoothingIters,
				m_EvolSettings.TopoParams.UseBackProjection });
			//remeshing.uniform_remeshing(targetEdgeLength);
#if REPORT_EVOL_STEPS
			std::cout << "done\n";
#endif
			if (ti % m_EvolSettings.TopoParams.StepStrideForEdgeDecay == 0 &&
				ti > NSteps * m_EvolSettings.TopoParams.RemeshingSizeDecayStartTimeFactor)
			{
				// shorter edges are needed for features close to the target.
				minEdgeLength *= m_EvolSettings.TopoParams.EdgeLengthDecayFactor;
				maxEdgeLength *= m_EvolSettings.TopoParams.EdgeLengthDecayFactor;
				//approxError *= m_EvolSettings.TopoParams.EdgeLengthDecayFactor;
			}
		}

		coVolStats = AnalyzeMeshCoVolumes(*m_EvolvingSurface, m_LaplacianAreaFunction);
#if REPORT_EVOL_STEPS
		std::cout << "Co-Volume Measure Stats: { Mean: " << coVolStats.Mean << ", Min: " << coVolStats.Min << ", Max: " << coVolStats.Max << "},\n";
		fileOStreamMins << coVolStats.Min << (ti < NSteps ? ", " : "");
		fileOStreamMeans << coVolStats.Mean << (ti < NSteps ? ", " : "");
		fileOStreamMaxes << coVolStats.Max << (ti < NSteps ? ", " : "");
#endif
		// set surface vertex properties
		for (const auto v : m_EvolvingSurface->vertices())
		{
			const auto vPos = m_EvolvingSurface->position(v);
			const double vDistanceToTarget = Geometry::TrilinearInterpolateScalarValue(vPos, field);
			vDistance[v] = static_cast<pmp::Scalar>(vDistanceToTarget);
		}
		Geometry::ComputeEdgeDihedralAngles(*m_EvolvingSurface);
		Geometry::ComputeVertexCurvaturesAndRelatedProperties(*m_EvolvingSurface);
		ComputeTriangleMetrics();

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
	// -------------------------------------------------------------------------------------------------------------

	if (m_EvolSettings.ExportResultSurface)
		ExportSurface(NSteps, true);

#if REPORT_EVOL_STEPS
	fileOStreamMins.close();
	fileOStreamMaxes.close();
	fileOStreamMeans.close();
#endif
}

void ReportInput(const IsoSurfaceEvolutionSettings& evolSettings, std::ostream& os)
{
	os << "======================================================================\n";
	os << "> > > > > > > > > > Initiating IsoSurfaceEvolver: < < < < < < < < < < < <\n";
	os << "Target Name: " << evolSettings.ProcedureName << ",\n";
	os << "NSteps: " << evolSettings.NSteps << ",\n";
	os << "TimeStep: " << evolSettings.TimeStep << ",\n";
	os << "FieldIsoLevel: " << evolSettings.FieldIsoLevel << ",\n";
	os << "ReSampledGridCellSize: " << evolSettings.ReSampledGridCellSize << ",\n";
	os << "......................................................................\n";
	const auto& c1 = evolSettings.ADParams.MCFMultiplier;
	const auto& c2 = evolSettings.ADParams.MCFVariance;
	os << "Curvature diffusion weight: " << c1 << " * (1 - exp(d^2 / " << c2 << ")),\n";
	const auto& d1 = evolSettings.ADParams.AdvectionMultiplier;
	const auto& d2 = evolSettings.ADParams.AdvectionSineMultiplier;
	os << "Advection weight: " << d1 << " * d * ((-grad(d) . N) - " << d2 << " * sqrt(1 - (grad(d) . N)^2)),\n";
	os << "......................................................................\n";
	os << "Export Surface per Time Step: " << (evolSettings.ExportSurfacePerTimeStep ? "true" : "false") << ",\n";
	os << "Output Path: " << evolSettings.OutputPath << ",\n";
	os << "Do Remeshing: " << (evolSettings.DoRemeshing ? "true" : "false") << ",\n";
	os << "Do Feature Detection: " << (evolSettings.DoFeatureDetection ? "true" : "false") << ",\n";
	os << "----------------------------------------------------------------------\n";
}

/// \brief The power of the stabilizing scale factor.
constexpr float SCALE_FACTOR_POWER = 1.0f / 2.0f;
/// \brief the reciprocal value of how many times the surface area element shrinks during evolution.
constexpr float INV_SHRINK_FACTOR = 5.0f;

float GetStabilizationScalingFactor(const double& timeStep, const float& cellSize, const float& stabilizationFactor)
{
	// const float expectedMeanCoVolArea = stabilizationFactor * 2.0f * sqrt(3.0f) * (cellSize * cellSize) / 3.0f;
	const float expectedMeanCoVolArea = stabilizationFactor * 4.0f * sqrt(2.0f) * (cellSize * cellSize) / 3.0f;
	return pow(static_cast<float>(timeStep) / expectedMeanCoVolArea * INV_SHRINK_FACTOR, SCALE_FACTOR_POWER);
}

