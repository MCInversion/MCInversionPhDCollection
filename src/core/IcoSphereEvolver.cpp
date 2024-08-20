#include "IcoSphereEvolver.h"

#include "pmp/algorithms/Normals.h"
#include "geometry/GeometryConversionUtils.h"
#include "geometry/GridUtil.h"
#include "geometry/MeshAnalysis.h"
#include "sdf/SDF.h"
#include "ConversionUtils.h"

#include <fstream>

#include "geometry/IcoSphereBuilder.h"


/// \brief if true individual steps of surface evolution will be printed out into a given stream.
#define REPORT_EVOL_STEPS true // Note: may affect performance

/// \brief if true, upon computing trilinear system solution, new vertices are verified for belonging in the field bounds.
#define VERIFY_SOLUTION_WITHIN_BOUNDS false // Note: useful for detecting numerical explosions of the solution.

IcoSphereEvolver::IcoSphereEvolver(const std::vector<pmp::Point>& pointCloud, const IcoSphereEvolutionSettings& settings)
	: m_PointCloud(pointCloud),
	m_EvolSettings(settings)
{
	m_ImplicitLaplacianFunction =
		(m_EvolSettings.LaplacianType == MeshLaplacian::Barycentric ?
			pmp::laplace_implicit_barycentric : pmp::laplace_implicit_voronoi);
	m_LaplacianAreaFunction =
		(m_EvolSettings.LaplacianType == MeshLaplacian::Barycentric ?
			pmp::barycentric_area : pmp::voronoi_area);
}

//
// ================================================================================================
//

void IcoSphereEvolver::ExportSurface(const unsigned int& tId, const bool& isResult, const bool& transformToOriginal) const
{
	if (!m_EvolvingSurface)
	{
		std::cerr << "IcoSphereEvolver::ExportSurface: m_EvolvingSurface == nullptr!\n";
		return;
	}
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

void IcoSphereEvolver::ExportField(const bool& transformToOriginal) const
{
	if (!m_Field)
	{
		std::cerr << "IcoSphereEvolver::ExportField: m_Field == nullptr!\n";
		return;
	}
	if (!transformToOriginal)
	{
		ExportToVTI(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + "_SDF", *m_Field);
		return;
	}
	auto exportedField = *m_Field;
	exportedField *= m_TransformToOriginal;
	ExportToVTI(m_EvolSettings.OutputPath + m_EvolSettings.ProcedureName + "_SDF", exportedField);
}

void IcoSphereEvolver::ComputeTriangleMetrics() const
{
	for (const auto& metricName : m_EvolSettings.TriMetrics)
	{
		if (!Geometry::IsMetricRegistered(metricName))
			continue;

		const auto metricFunction = Geometry::IdentifyMetricFunction(metricName);

		if (!metricFunction(*m_EvolvingSurface))
		{
			std::cerr << "IcoSphereEvolver::ComputeTriangleMetrics: [WARNING] Computation of metric " << metricName << " finished with errors!\n";
		}
	}
}

//
// ================================================================================================
//

void IcoSphereEvolver::ComputeDistanceField()
{
#if REPORT_EVOL_STEPS
	std::cout << "IcoSphereEvolver::ComputeDistanceField: with " << m_EvolSettings.NVoxelsPerMinDimension << " voxels per min dimension ... ";
#endif
	const pmp::BoundingBox ptCloudBBox(m_PointCloud);
	const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
	const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
	const float cellSize = minSize / static_cast<float>(m_EvolSettings.NVoxelsPerMinDimension);
	constexpr float volExpansionFactor = 1.0f;
	const SDF::PointCloudDistanceFieldSettings dfSettings{
				cellSize,
				volExpansionFactor,
				Geometry::DEFAULT_SCALAR_GRID_INIT_VAL, // m_EvolSettings.MaxDim * 1.0
				SDF::BlurPostprocessingType::None
	};
	m_Field = std::make_shared<Geometry::ScalarGrid>(
		SDF::PointCloudDistanceFieldGenerator::Generate(m_PointCloud, dfSettings));
#if REPORT_EVOL_STEPS
	std::cout << "done.\n";
	const auto& dims = m_Field->Dimensions();
	std::cout << "IcoSphereEvolver::ComputeDistanceField: field resolution: " << dims.Nx << " x " << dims.Ny << " x " << dims.Nz << " voxels.\n";
#endif
}

/// \brief a multiplier allowing some deviation from the bounding sphere.
constexpr float ICO_SPHERE_RADIUS_FACTOR = 1.001f;

void IcoSphereEvolver::ConstructIcoSphere()
{
	const auto [boundingSphereCenter, boundingSphereRadius] = Geometry::ComputePointCloudBoundingSphere(m_PointCloud);
	m_StartingSurfaceRadius = boundingSphereRadius * ICO_SPHERE_RADIUS_FACTOR;
	m_StartingSurfaceCenter = boundingSphereCenter;
#if REPORT_EVOL_STEPS
	std::cout << "Ico-Sphere Radius: " << m_StartingSurfaceRadius << ",\n";
#endif
	const unsigned int icoSphereSubdiv = m_EvolSettings.IcoSphereSubdivisionLevel;
	Geometry::IcoSphereBuilder icoBuilder({ icoSphereSubdiv, m_StartingSurfaceRadius });
	icoBuilder.BuildBaseData();
	icoBuilder.BuildPMPSurfaceMesh();
	m_EvolvingSurface = std::make_shared<pmp::SurfaceMesh>(icoBuilder.GetPMPSurfaceMeshResult());
}

void IcoSphereEvolver::DetectFeatureVerticesFromDistanceField() const
{
	if (!m_EvolvingSurface)
	{
		throw std::runtime_error("IcoSphereEvolver::DetectFeatureVerticesFromDistanceField: m_EvolvingSurface not initialized!\n");
	}

	if (!m_Field)
	{
		throw std::runtime_error("IcoSphereEvolver::DetectFeatureVerticesFromDistanceField: m_Field not initialized!\n");
	}

	const auto& field = *m_Field;
	auto vIsFeature = m_EvolvingSurface->add_vertex_property<bool>("v:feature", false);
#if REPORT_EVOL_STEPS
	unsigned int nFeatureVerts = 0;
#endif
	for (const auto v : m_EvolvingSurface->vertices())
	{
		const auto vPos = m_EvolvingSurface->position(v);
		const auto vDistanceToTarget = Geometry::TrilinearInterpolateScalarValue(vPos, field);

		if (vDistanceToTarget > 2 * m_EvolSettings.FieldIsoLevel)
			continue;

		vIsFeature[v] = true;
#if REPORT_EVOL_STEPS
		nFeatureVerts++;
#endif
	}
#if REPORT_EVOL_STEPS
	if (nFeatureVerts > 0)
	{
		std::cout << "IcoSphereEvolver::DetectFeatureVerticesFromDistanceField: detected " << nFeatureVerts << " vertices within " << m_EvolSettings.FieldIsoLevel << " units to the target point cloud.\n";
	}
#endif
}

void IcoSphereEvolver::StabilizeEvolutionData()
{
	if (!m_EvolvingSurface)
	{
		throw std::runtime_error("IcoSphereEvolver::StabilizeEvolutionData: m_EvolvingSurface not initialized!\n");
	}

	// transform mesh and grid
	const unsigned int& icoSphereSubdiv = m_EvolSettings.IcoSphereSubdivisionLevel;
	const float scalingFactor = GetStabilizationScalingFactor(m_EvolSettings.TimeStep, m_StartingSurfaceRadius, icoSphereSubdiv);
	m_ScalingFactor = scalingFactor;
	m_EvolSettings.FieldIsoLevel *= static_cast<double>(scalingFactor);
	const auto origin = m_StartingSurfaceCenter;
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
	(*m_Field) *= static_cast<double>(scalingFactor); // scale also distance values.

	// Export field for debugging purposes
	//ExportField();
}

// ================================================================================================

double IcoSphereEvolver::LaplacianDistanceWeightFunction(const double& distanceAtVertex) const
{
	if (distanceAtVertex < 0.0 && m_EvolSettings.ADParams.MCFSupportPositive)
		return 0.0;
	const auto& c1 = m_EvolSettings.ADParams.MCFMultiplier;
	const auto& c2 = m_EvolSettings.ADParams.MCFVariance;
	return c1 * (1.0 - exp(-(distanceAtVertex * distanceAtVertex) / c2));
}

double IcoSphereEvolver::AdvectionDistanceWeightFunction(const double& distanceAtVertex,
	const pmp::dvec3& negDistanceGradient, const pmp::Point& vertexNormal) const
{
	if (distanceAtVertex < 0.0 && m_EvolSettings.ADParams.AdvectionSupportPositive)
		return 0.0;
	const auto& d1 = m_EvolSettings.ADParams.AdvectionMultiplier;
	const auto& d2 = m_EvolSettings.ADParams.AdvectionSineMultiplier;
	const auto negGradDotNormal = pmp::ddot(negDistanceGradient, vertexNormal);
	return d1 * distanceAtVertex * (negGradDotNormal - d2 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
}

//
// ================================================================================================
//

void IcoSphereEvolver::Preprocess()
{
	ComputeDistanceField();
	ConstructIcoSphere();
	DetectFeatureVerticesFromDistanceField();
	StabilizeEvolutionData();
}

//
// ================================================================================================
//

void IcoSphereEvolver::Evolve()
{
	if (m_PointCloud.empty())
		throw std::invalid_argument("IcoSphereEvolver::Evolve: m_PointCloud.empty()!\n");

	Preprocess();

	if (!m_Field)
		throw std::invalid_argument("ConvexHullEvolver::Evolve: m_Field not set! Terminating!\n");
	if (!m_Field->IsValid())
		throw std::invalid_argument("ConvexHullEvolver::Evolve: m_Field is invalid! Terminating!\n");
	if (!m_EvolvingSurface)
		throw std::invalid_argument("ConvexHullEvolver::Evolve: m_EvolvingSurface not set! Terminating!\n");

	const auto& field = *m_Field;

#if VERIFY_SOLUTION_WITHIN_BOUNDS
	const auto& fieldBox = field.Box();
#endif
	const auto fieldNegGradient = Geometry::ComputeNormalizedNegativeGradient(field);

	const auto& NSteps = m_EvolSettings.NSteps;
	auto tStep = m_EvolSettings.TimeStep;

	// ........ evaluate edge lengths for remeshing ....................
	const auto subdiv = static_cast<float>(m_EvolSettings.IcoSphereSubdivisionLevel);
	const float r = m_StartingSurfaceRadius * m_ScalingFactor;
	constexpr float baseIcoHalfAngle = 2.0f * M_PI / 10.0f;
	const float minEdgeMultiplier = m_EvolSettings.TopoParams.MinEdgeMultiplier;

	auto minEdgeLength = minEdgeMultiplier * 2.0f * r * sin(baseIcoHalfAngle * pow(2.0f, -subdiv)); // from icosahedron edge length
	auto maxEdgeLength = 4.0f * minEdgeLength;
	auto approxError = 0.25f * (minEdgeLength + maxEdgeLength);
	m_Remesher = std::make_shared<pmp::Remeshing>(*m_EvolvingSurface);

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
	auto vIsFeatureVal = m_EvolvingSurface->vertex_property<pmp::Scalar>("v:isFeature", -1.0f);

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
		vIsFeatureVal[v] = (vFeature[v] ? 1.0f : -1.0f);
	}
	Geometry::ComputeEdgeDihedralAngles(*m_EvolvingSurface);
	Geometry::ComputeVertexCurvaturesAndRelatedProperties(*m_EvolvingSurface, m_EvolSettings.TopoParams.PrincipalCurvatureFactor);
	ComputeTriangleMetrics();
	if (m_EvolSettings.ExportSurfacePerTimeStep)
		ExportSurface(0);

	// -------------------------------------------------------------------------------------------------------------
	// ........................................ main loop ..........................................................
	// -------------------------------------------------------------------------------------------------------------
	for (unsigned int ti = 1; ti <= NSteps; ti++)
	{
#if REPORT_EVOL_STEPS
		std::cout << "time step id: " << ti << "/" << NSteps << ", time: " << tStep * ti << " "
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
				std::cerr << "SurfaceEvolver::Evolve: vertex " << i << " out of field bounds!\n";
				nVertsOutOfBounds++;
			}
#endif
			m_EvolvingSurface->position(pmp::Vertex(i)) = x.row(i);
		}
#if VERIFY_SOLUTION_WITHIN_BOUNDS
		if ((m_EvolSettings.DoRemeshing && nVertsOutOfBounds > static_cast<double>(NVertices) * m_EvolSettings.MaxFractionOfVerticesOutOfBounds) ||
			(!m_EvolSettings.DoRemeshing && nVertsOutOfBounds > 0))
		{
			std::cerr << "SurfaceEvolver::Evolve: found " << nVertsOutOfBounds << " vertices out of bounds! Terminating!\n";
			break;
		}
#endif
#if REPORT_EVOL_STEPS
		std::cout << "done\n";
#endif

		// --------------------------------------------------------------------

		const auto meshQualityProp = m_EvolvingSurface->get_vertex_property<float>("v:equilateralJacobianCondition");
		if (m_EvolSettings.DoRemeshing && IsRemeshingNecessary(meshQualityProp.vector()))
		{
			// remeshing
#if REPORT_EVOL_STEPS
			std::cout << "Remeshing ...";
#endif
#if REPORT_EVOL_STEPS
			std::cout << "pmp::Remeshing::adaptive_remeshing(minEdgeLength: " << minEdgeLength << ", maxEdgeLength: " << maxEdgeLength << ", approxError: " << approxError << ") ... ";
#endif
			pmp::Remeshing remeshing(*m_EvolvingSurface);
			remeshing.adaptive_remeshing({
				minEdgeLength, maxEdgeLength, approxError,
				m_EvolSettings.TopoParams.NRemeshingIters,
				m_EvolSettings.TopoParams.NTanSmoothingIters,
				m_EvolSettings.TopoParams.UseBackProjection });
#if REPORT_EVOL_STEPS
			std::cout << "done\n";
#endif
		}

		// --------------------------------------------------------------------

		if (ShouldAdjustRemeshingLengths(ti))
		{
			// shorter edges are needed for features close to the target.
			AdjustRemeshingLengths(m_EvolSettings.TopoParams.EdgeLengthDecayFactor, minEdgeLength, maxEdgeLength, approxError);
#if REPORT_EVOL_STEPS
			std::cout << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n";
			std::cout << "Lengths for adaptive remeshing adjusted to:\n";
			std::cout << "min: " << minEdgeLength << ", max: " << maxEdgeLength << ", error: " << approxError << "\n";
			std::cout << "by a factor of " << m_EvolSettings.TopoParams.EdgeLengthDecayFactor << ".\n";
			std::cout << "time step adjustment: " << tStep << " -> ";
#endif
			tStep *= pow(m_EvolSettings.TopoParams.EdgeLengthDecayFactor, 2);
#if REPORT_EVOL_STEPS
			std::cout << tStep << ".\n";
			std::cout << "AdvectionMultiplier adjustment: " << m_EvolSettings.ADParams.AdvectionMultiplier << " -> ";
#endif
			m_EvolSettings.ADParams.AdvectionMultiplier /= m_EvolSettings.TopoParams.EdgeLengthDecayFactor;
			//m_EvolSettings.ADParams.AdvectionSineMultiplier /= m_EvolSettings.TopoParams.EdgeLengthDecayFactor;
#if REPORT_EVOL_STEPS
			std::cout << m_EvolSettings.ADParams.AdvectionMultiplier << "\n";
			std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
#endif
		}

		// --------------------------------------------------------------------

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
			vIsFeatureVal[v] = (vFeature[v] ? 1.0f : -1.0f);
		}
		Geometry::ComputeEdgeDihedralAngles(*m_EvolvingSurface);
		Geometry::ComputeVertexCurvaturesAndRelatedProperties(*m_EvolvingSurface, m_EvolSettings.TopoParams.PrincipalCurvatureFactor);
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

//
// ================================================================================================
//

void ReportIcoEvolverInput(const IcoSphereEvolutionSettings& evolSettings, std::ostream& os)
{
	os << "======================================================================\n";
	os << "> > > > > > > > > > Initiating IcoSphereEvolver: < < < < < < < < < < < <\n";
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
	os << "Target Origin: " << evolSettings.TargetOrigin << ",\n";
	os << "Export Surface per Time Step: " << (evolSettings.ExportSurfacePerTimeStep ? "true" : "false") << ",\n";
	os << "Output Path: " << evolSettings.OutputPath << ",\n";
	os << "Do Remeshing: " << (evolSettings.DoRemeshing ? "true" : "false") << ",\n";
	os << "Do Feature Detection: " << (evolSettings.DoFeatureDetection ? "true" : "false") << ",\n";
	os << "----------------------------------------------------------------------\n";
}
