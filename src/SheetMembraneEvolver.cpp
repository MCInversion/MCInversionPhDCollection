#include "SheetMembraneEvolver.h"

#include "pmp/algorithms/Remeshing.h"
#include "pmp/algorithms/Normals.h"
#include "pmp/algorithms/Decimation.h"
#include "pmp/algorithms/Features.h"

#include "geometry/GridUtil.h"
#include "geometry/MeshAnalysis.h"

#include <fstream>

#include "ConversionUtils.h"
#include "geometry/GeometryConversionUtils.h"
#include "geometry/PlaneBuilder.h"

// ================================================================================================

/// \brief if true individual steps of surface evolution will be printed out into a given stream.
#define REPORT_EVOL_STEPS true // Note: may affect performance

/// \brief if true, upon computing trilinear system solution, new vertices are verified for belonging in the field bounds.
#define VERIFY_SOLUTION_WITHIN_BOUNDS false // Note: useful for detecting numerical explosions of the solution.

/// \brief repairs infinities and nan values of the input grid.
#define REPAIR_INPUT_GRID true //Note: very useful! You never know what field we use! Infinities and nans lead to the Eigen::NoConvergence error code whose cause is difficult to find!

// ================================================================================================

size_t SheetMembraneEvolver::DetectFeatures(const FeatureDetectionType& type) const
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

SheetMembraneEvolver::SheetMembraneEvolver(const Geometry::ScalarGrid& field, const SheetMembraneEvolutionSettings& settings)
	: m_EvolSettings(settings), m_Field(std::make_shared<Geometry::ScalarGrid>(field))
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

void SheetMembraneEvolver::Preprocess()
{
	auto& field = *m_Field;
	const auto& fieldBox = field.Box();

	const float startZHeight = m_EvolSettings.StartZHeight;
	const float endZHeight = m_EvolSettings.EndZHeight;
	m_SheetSurfaceVelocity = static_cast<double>(startZHeight) - static_cast<double>(endZHeight);

	// build plane surface
	const pmp::vec3 fieldOrig = fieldBox.min();
	const pmp::vec3 planeOrig = pmp::vec3{ fieldOrig[0], fieldOrig[1], startZHeight };
	const pmp::vec3 fieldSize = fieldBox.max() - fieldBox.min();
	const Geometry::PlaneSettings mSettings{
		planeOrig,
		fieldSize[0],
		fieldSize[1],
		m_EvolSettings.nXSegments,
		m_EvolSettings.nYSegments,
		false,
		true
	};
	Geometry::PlaneBuilder pb(mSettings);
	pb.BuildBaseData();
	pb.BuildPMPSurfaceMesh();
	m_EvolvingSurface = std::make_shared<pmp::SurfaceMesh>(pb.GetPMPSurfaceMeshResult());

	// transform mesh and grid
	// >>> uniform scale to ensure numerical method's stability.
	// >>> translation to origin for fields not centered at (0,0,0).
	// >>> scaling factor value is intended for stabilization of the numerical method.

	const float cellSizeX = fieldSize[0] / static_cast<float>(m_EvolSettings.nXSegments);
	const float cellSizeY = fieldSize[1] / static_cast<float>(m_EvolSettings.nYSegments);
	m_MeanEdgeLength = (cellSizeX + cellSizeY + sqrt(cellSizeX * cellSizeX + cellSizeY * cellSizeY)) / 3.0f;

	const float scalingFactor = GetStabilizationScalingFactor(m_EvolSettings.TimeStep, cellSizeX, cellSizeY);
	m_ScalingFactor = scalingFactor;
#if REPORT_EVOL_STEPS
	std::cout << "Stabilization Scaling Factor: " << scalingFactor << ",\n";
#endif

	m_EvolSettings.FieldIsoLevel *= static_cast<double>(scalingFactor);
	m_SheetSurfaceVelocity *= static_cast<double>(scalingFactor);
	m_MeanEdgeLength *= scalingFactor;

	const pmp::mat4 transfMatrixGeomScale{
		scalingFactor, 0.0f, 0.0f, 0.0f,
		0.0f, scalingFactor, 0.0f, 0.0f,
		0.0f, 0.0f, scalingFactor, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	const auto planeCenter = pmp::vec3{ fieldOrig[0] + 0.5f * fieldSize[0], fieldOrig[1] + 0.5f * fieldSize[1], planeOrig[2] };
	const pmp::mat4 transfMatrixGeomMove{
		1.0f, 0.0f, 0.0f, -planeCenter[0],
		0.0f, 1.0f, 0.0f, -planeCenter[1],
		0.0f, 0.0f, 1.0f, -planeCenter[2],
		0.0f, 0.0f, 0.0f, 1.0f
	};
	const auto transfMatrixFull = transfMatrixGeomScale * transfMatrixGeomMove;
	m_TransformToOriginal = inverse(transfMatrixFull);

	(*m_EvolvingSurface) *= transfMatrixFull; // ico sphere is already centered at (0,0,0).
	field *= transfMatrixFull; // field needs to be moved to (0,0,0) and also scaled.
	field *= static_cast<double>(scalingFactor); // scale also distance values.

	// >>>>> Scaled geometry & field (use when debugging) <<<<<<
	//ExportToVTI(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + "_resampledField", reSampledField);
	//pmp::SurfaceMesh scaledTargetMesh = *m_EvolvingSurface;
	//scaledTargetMesh *= transfMatrixFull;
	//scaledTargetMesh.write(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + "_stableScale.obj");

}

// ================================================================================================

double SheetMembraneEvolver::LaplacianDistanceWeightFunction(const double& distanceAtVertex) const
{
	if (distanceAtVertex < 0.0 && m_EvolSettings.ADParams.MCFSupportPositive)
		return 0.0;
	const auto& c1 = m_EvolSettings.ADParams.MCFMultiplier;
	const auto& c2 = m_EvolSettings.ADParams.MCFVariance;
	return c1 * (1.0 - exp(-(distanceAtVertex * distanceAtVertex) / c2));
}

double SheetMembraneEvolver::AdvectionDistanceWeightFunction(const double& distanceAtVertex,
	const pmp::dvec3& negDistanceGradient, const pmp::Point& vertexNormal) const
{
	if (distanceAtVertex < 0.0 && m_EvolSettings.ADParams.AdvectionSupportPositive)
		return 0.0;
	const auto& d1 = m_EvolSettings.ADParams.AdvectionMultiplier;
	const auto& d2 = m_EvolSettings.ADParams.AdvectionSineMultiplier;
	const auto negGradDotNormal = pmp::ddot(negDistanceGradient, vertexNormal);
	return d1 * distanceAtVertex * (negGradDotNormal - d2 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
}

void SheetMembraneEvolver::ExportSurface(const unsigned int& tId, const bool& isResult, const bool& transformToOriginal) const
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

void SheetMembraneEvolver::ComputeTriangleMetrics() const
{
	for (const auto& metricName : m_EvolSettings.TriMetrics)
	{
		if (!Geometry::IsMetricRegistered(metricName))
			continue;

		const auto metricFunction = Geometry::IdentifyMetricFunction(metricName);

		if (!metricFunction(*m_EvolvingSurface))
		{
			std::cerr << "SheetMembraneEvolver::ComputeTriangleMetrics: [WARNING] Computation of metric " << metricName << " finished with errors!\n";
		}
	}
}

// ================================================================================================

void SheetMembraneEvolver::Evolve()
{
	if (!m_Field)
		throw std::invalid_argument("SheetMembraneEvolver::Evolve: m_Field not set! Terminating!\n");
	if (!m_Field->IsValid())
		throw std::invalid_argument("SheetMembraneEvolver::Evolve: m_Field is invalid! Terminating!\n");

	const auto& field = *m_Field;

	Preprocess();

	if (!m_EvolvingSurface)
		throw std::invalid_argument("SheetMembraneEvolver::Evolve: m_EvolvingSurface not set! Terminating!\n");

#if VERIFY_SOLUTION_WITHIN_BOUNDS
	const auto& fieldBox = field.Box();
#endif
	const auto fieldNegGradient = Geometry::ComputeNormalizedNegativeGradient(field);

	const auto& NSteps = m_EvolSettings.NSteps;
	const auto& tStep = m_EvolSettings.TimeStep;

	// ........ evaluate edge lengths for remeshing ....................
	float minEdgeLength = m_MeanEdgeLength * m_EvolSettings.TopoParams.MinEdgeMultiplier;
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
		const Eigen::Vector3d downVec{ 0.0, 0.0, -1.0 };
		for (const auto v : m_EvolvingSurface->vertices())
		{
			const auto vPosToUpdate = m_EvolvingSurface->position(v);

			if ((m_EvolSettings.IdentityForBoundaryVertices && m_EvolvingSurface->is_boundary(v)) ||
				(m_EvolSettings.IdentityForFeatureVertices && vFeature[v]))
			{
				// move boundary/feature vertices along downVec
				const Eigen::Vector3d vertexRhs = vPosToUpdate;
				sysRhs.row(v.idx()) = vertexRhs + tStep * m_SheetSurfaceVelocity * downVec;
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
	// Geometry::ComputeEdgeDihedralAngles(*m_EvolvingSurface);
	Geometry::ComputeVertexCurvatures(*m_EvolvingSurface);
	Geometry::ComputeZLevelElevations(*m_EvolvingSurface);
	// ComputeTriangleMetrics();
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
			const std::string msg = "\nSheetMembraneEvolver::Evolve: solver.info() != Eigen::Success for time step id: "
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
				std::cerr << "SheetMembraneEvolver::Evolve: vertex " << i << " out of field bounds!\n";
				nVertsOutOfBounds++;
			}
#endif
			m_EvolvingSurface->position(pmp::Vertex(i)) = x.row(i);
		}
#if VERIFY_SOLUTION_WITHIN_BOUNDS
		if ((m_EvolSettings.DoRemeshing && nVertsOutOfBounds > static_cast<double>(NVertices) * m_EvolSettings.MaxFractionOfVerticesOutOfBounds) ||
			(!m_EvolSettings.DoRemeshing && nVertsOutOfBounds > 0))
		{
			std::cerr << "SheetMembraneEvolver::Evolve: found " << nVertsOutOfBounds << " vertices out of bounds! Terminating!\n";
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
		// Geometry::ComputeEdgeDihedralAngles(*m_EvolvingSurface);
		Geometry::ComputeVertexCurvatures(*m_EvolvingSurface);
		Geometry::ComputeZLevelElevations(*m_EvolvingSurface);
		// ComputeTriangleMetrics();

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

void ReportInput(const SheetMembraneEvolutionSettings& evolSettings, std::ostream& os)
{
	os << "======================================================================\n";
	os << "> > > > > > > > > > Initiating SheetMembraneEvolver: < < < < < < < < < < < <\n";
	os << "Target Name: " << evolSettings.ProcedureName << ",\n";
	os << "NSteps: " << evolSettings.NSteps << ",\n";
	os << "TimeStep: " << evolSettings.TimeStep << ",\n";
	os << "StartZHeight: " << evolSettings.StartZHeight << ",\n";
	os << "EndZHeight: " << evolSettings.EndZHeight << ",\n";
	os << "nXSegments: " << evolSettings.nXSegments << ",\n";
	os << "nYSegments: " << evolSettings.nYSegments << ",\n";
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
constexpr float INV_SHRINK_FACTOR = 2.0f;

float GetStabilizationScalingFactor(const double& timeStep, const float& cellSizeX, const float& cellSizeY, const float& stabilizationFactor)
{
	const float expectedMeanCoVolArea = stabilizationFactor * cellSizeX * cellSizeY;
	return pow(static_cast<float>(timeStep) / expectedMeanCoVolArea * INV_SHRINK_FACTOR, SCALE_FACTOR_POWER);
}

//
// ============================================================================================================
//

/**
 * \brief A verification function for the column 2D position within field box.
 * \param fieldBox      bounding box of the scalar field
 * \param pos           2D position in the x,y-plane
 * \return true if pos is within box's x,y-range
 */
static [[nodiscard]] bool IsColumnInField(const pmp::BoundingBox& fieldBox, const pmp::vec2& pos)
{
	if (pos[0] < fieldBox.min()[0])
		return false;

	if (pos[0] > fieldBox.max()[0])
		return false;

	if (pos[1] < fieldBox.min()[1])
		return false;

	return pos[1] <= fieldBox.max()[1];
}

Geometry::ScalarGrid GetDistanceFieldWithSupportColumns(
	const float& cellSize, const pmp::BoundingBox& box, 
	const std::vector<WeightedColumnPosition>& weightedColumnPositions, const float& supportZLevel)
{
	if (box.is_empty())
	{
		std::cerr << "GetDistanceFieldWithSupportColumns: box.is_empty()!\n";
		throw std::logic_error("GetDistanceFieldWithSupportColumns: box.is_empty()!\n");
	}

	const auto boxSize = box.max() - box.min();

	// input verification
	if (cellSize >= boxSize[0] || cellSize >= boxSize[1] || cellSize >= boxSize[2])
	{
		std::cerr << "GetDistanceFieldWithSupportColumns: cellSize too large!\n";
		throw std::logic_error("GetDistanceFieldWithSupportColumns: cellSize too large!\n");
	}
	if (supportZLevel < 0.0f || supportZLevel > 1.0f)
	{
		std::cerr << "GetDistanceFieldWithSupportColumns: supportZLevel must be a value between 0 and 1!\n";
		throw std::logic_error("GetDistanceFieldWithSupportColumns: supportZLevel must be a value between 0 and 1!\n");
	}

	const float maxColumnRadius = 0.5f * std::fmaxf(boxSize[0], boxSize[1]);
	const float preferredColumnRadius = 0.05f * maxColumnRadius;
	const float columnZPosition = box.min()[2] + (supportZLevel - 0.1f) * boxSize[2];
	const float maxColumnHeight = box.max()[2] - columnZPosition;

	constexpr double initVal = Geometry::DEFAULT_SCALAR_GRID_INIT_VAL;
	Geometry::ScalarGrid result(cellSize, box, initVal);
	Geometry::CapsuleParams cp{};
	cp.Radius = preferredColumnRadius;
	cp.BoolOpFunction = Geometry::DistanceUnion;

	for (const auto& wPos : weightedColumnPositions)
	{
		const auto& pos = wPos.Position;
		if (!IsColumnInField(box, pos))
			continue;

		const auto& weight = wPos.Weight;
		if (weight < 0.0f || weight > 1.0f)
		{
			std::cerr << "GetDistanceFieldWithSupportColums: Weight must be a value between 0 and 1!\n";
			continue;
		}

		cp.Position = pmp::vec3{ pos[0], pos[1], columnZPosition };
		const auto capsuleHeight = maxColumnHeight * weight;
		cp.Height = capsuleHeight;

		ApplyCapsuleDistanceFieldToGrid(result, cp);
	}

	return result;
}
