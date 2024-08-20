#include "SphereTest.h"

#include "geometry/IcoSphereBuilder.h"
#include "geometry/MeshAnalysis.h"

#include "EvolverUtilsCommon.h"
#include "pmp/algorithms/Normals.h"
#include "pmp/algorithms/Remeshing.h"

/// \brief if true individual steps of surface evolution will be printed out into a given stream.
#define REPORT_EVOL_STEPS false // Note: may affect performance

/// \brief if true, upon computing trilinear system solution, new vertices are verified for belonging in the field bounds.
#define VERIFY_SOLUTION_WITHIN_BOUNDS true // Note: useful for detecting numerical explosions of the solution.

/// \brief fundamental sphere test settings
constexpr double STARTING_SPHERE_RADIUS = 1.0;
constexpr unsigned int STARTING_SPHERE_SUBDIVISION = 1;
constexpr unsigned int N_STARTING_SPHERE_TEST_TIME_STEPS = 6;
constexpr double STARTING_TIME_STEP_SIZE = 0.01;

// ================================================================================================

SphereTest::SphereTest(const SphereTestEvolutionSettings& settings)
	: m_EvolSettings(settings), m_TotalTime(N_STARTING_SPHERE_TEST_TIME_STEPS * STARTING_TIME_STEP_SIZE)
{
	m_ImplicitLaplacianFunction =
		(m_EvolSettings.LaplacianType == ST_MeshLaplacian::Barycentric ?
			pmp::laplace_implicit_barycentric : pmp::laplace_implicit_voronoi);
	m_LaplacianAreaFunction =
		(m_EvolSettings.LaplacianType == ST_MeshLaplacian::Barycentric ?
			pmp::barycentric_area : pmp::voronoi_area);
}

// ================================================================================================

void SphereTest::Preprocess(const unsigned int& iter)
{
	// build ico-sphere
	const unsigned int icoSphereSubdiv = STARTING_SPHERE_SUBDIVISION + iter;
	Geometry::IcoSphereBuilder icoBuilder({ icoSphereSubdiv, STARTING_SPHERE_RADIUS });
	icoBuilder.BuildBaseData();
	icoBuilder.BuildPMPSurfaceMesh();
	m_EvolvingSurface = std::make_shared<pmp::SurfaceMesh>(icoBuilder.GetPMPSurfaceMeshResult());
}

// ================================================================================================

void SphereTest::Evolve(const unsigned int& iter, const double& tStep)
{
	Preprocess(iter);

	if (!m_EvolvingSurface)
		throw std::invalid_argument("SphereTest::Evolve: m_EvolvingSurface not set! Terminating!\n");

#if VERIFY_SOLUTION_WITHIN_BOUNDS
	auto bbox = m_EvolvingSurface->bounds();
	bbox.expand(0.5f * STARTING_SPHERE_RADIUS, 0.5f * STARTING_SPHERE_RADIUS, 0.5f * STARTING_SPHERE_RADIUS);
#endif

	// ........ evaluate edge lengths for remeshing ....................
	//const float phi = (1.0f + sqrt(5.0f)) / 2.0f; /// golden ratio.
	const auto subdiv = static_cast<float>(STARTING_SPHERE_SUBDIVISION + iter);
	constexpr float r = STARTING_SPHERE_RADIUS;
	constexpr float baseIcoHalfAngle = 2.0f * M_PI / 10.0f;
	const float minEdgeMultiplier = m_EvolSettings.TopoParams.MinEdgeMultiplier;
	//auto minEdgeLength = minEdgeMultiplier * (r / (pow(2.0f, subdiv - 1) * sqrt(phi * sqrt(5.0f)))); // from icosahedron edge length
	auto minEdgeLength = minEdgeMultiplier * 2.0f * r * sin(baseIcoHalfAngle * pow(2.0f, -subdiv)); // from icosahedron edge length
	auto maxEdgeLength = 1.5f * minEdgeLength;
	auto approxError = 0.25f * (minEdgeLength + maxEdgeLength);
#if REPORT_EVOL_STEPS
	std::cout << "minEdgeLength for remeshing: " << minEdgeLength << "\n";
#endif
	// .................................................................

	const auto NSteps = static_cast<unsigned int>(m_TotalTime / tStep);
	// DISCLAIMER: the dimensionality of the system depends on the number of mesh vertices which can change if remeshing is used.
	auto NVertices = static_cast<unsigned int>(m_EvolvingSurface->n_vertices());
	SparseMatrix sysMat(NVertices, NVertices);
	Eigen::MatrixXd sysRhs(NVertices, 3);

	// property container for surface vertex normals
	pmp::VertexProperty<pmp::Point> vNormalsProp{};

	// ----------- System fill function --------------------------------
	const auto fillMatrixAndRHSTriplesFromMesh = [&]()
	{
		for (const auto v : m_EvolvingSurface->vertices())
		{
			const auto vPosToUpdate = m_EvolvingSurface->position(v);
			const auto vNormal = vNormalsProp[v]; // vertex unit normal

			const Eigen::Vector3d vertexRhs = vPosToUpdate;
			sysRhs.row(v.idx()) = vertexRhs;
			const float tanRedistWeight = m_EvolSettings.TangentialVelocityWeight;
			if (tanRedistWeight > 0.0f)
			{
				// compute tangential velocity
				const auto vTanVelocity = ComputeTangentialUpdateVelocityAtVertex(*m_EvolvingSurface, v, vNormal, tanRedistWeight);
				sysRhs.row(v.idx()) += tStep * Eigen::Vector3d(vTanVelocity);
			}

			const auto laplaceWeightInfo = m_ImplicitLaplacianFunction(*m_EvolvingSurface, v); // Laplacian weights
			sysMat.coeffRef(v.idx(), v.idx()) = 1.0 + tStep * static_cast<double>(laplaceWeightInfo.weightSum);

			for (const auto& [w, weight] : laplaceWeightInfo.vertexWeights)
			{
				sysMat.coeffRef(v.idx(), w.idx()) = -1.0 * tStep * static_cast<double>(weight);
			}
		}
	};
	// -----------------------------------------------------------------
	// starting surface error 
	double totalSurfaceSqError = 0.0;

	if (m_EvolSettings.DoRemeshing)
	{
		// init mean edge length
		m_MeanEdgeLength = ComputeMeanEdgeLength();
	}

	// -------------------------------------------------------------------------------------------------------------
	// ........................................ main loop ..........................................................
	// -------------------------------------------------------------------------------------------------------------
	for (unsigned int ti = 1; ti <= NSteps; ti++)
	{
		std::cout << "\riter" << iter << "[" << ti << "/" << NSteps << "]:";
#if REPORT_EVOL_STEPS
		std::cout << "SphereTest: time step id: " << ti << "/" << NSteps << ", time: " << tStep * ti << "/" << tStep * NSteps << "\n";
		std::cout << "pmp::Normals::compute_vertex_normals ... ";
#endif
		pmp::Normals::compute_vertex_normals(*m_EvolvingSurface);
#if REPORT_EVOL_STEPS
		std::cout << "done\n";
		std::cout << "fillMatrixAndRHSTriplesFromMesh for " << NVertices << " vertices ... ";
#endif
		vNormalsProp = m_EvolvingSurface->vertex_property<pmp::Point>("v:normal");

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
			const std::string msg = "\nSphereTest::Evolve: solver.info() != Eigen::Success for time step id: "
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
			if (!bbox.Contains(newPos))
			{
				if (nVertsOutOfBounds == 0) std::cerr << "\n";
				std::cerr << "SphereTest::Evolve: vertex " << i << " out of field bounds!\n";
				nVertsOutOfBounds++;
			}
#endif
			m_EvolvingSurface->position(pmp::Vertex(i)) = x.row(i);
		}
#if VERIFY_SOLUTION_WITHIN_BOUNDS
		if ((m_EvolSettings.DoRemeshing && nVertsOutOfBounds > static_cast<double>(NVertices) * m_EvolSettings.MaxFractionOfVerticesOutOfBounds) ||
			(!m_EvolSettings.DoRemeshing && nVertsOutOfBounds > 0))
		{
			std::cerr << "SphereTest::Evolve: found " << nVertsOutOfBounds << " vertices out of bounds! Terminating!\n";
			break;
		}
#endif
#if REPORT_EVOL_STEPS
		std::cout << "done\n";
#endif

		if (m_EvolSettings.DoRemeshing && ti > NSteps * m_EvolSettings.TopoParams.RemeshingStartTimeFactor)
		{
			// remeshing
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

			m_MeanEdgeLength += ComputeMeanEdgeLength();
			
			// update linear system dims for next time step:
			if (ti < NSteps && NVertices != m_EvolvingSurface->n_vertices())
			{
				NVertices = static_cast<unsigned int>(m_EvolvingSurface->n_vertices());
				sysMat = SparseMatrix(NVertices, NVertices);
				sysRhs = Eigen::MatrixXd(NVertices, 3);
			}
		}

		// evaluate error
		totalSurfaceSqError += EvaluateSquaredSurfaceError(ti * tStep);

#if REPORT_EVOL_STEPS
		std::cout << ">>> Time step " << ti << " finished.\n";
		std::cout << "----------------------------------------------------------------------\n";
#endif

	} // end main loop
	// -------------------------------------------------------------------------------------------------------------

	std::cout << "...done.";
	m_MeanErrors.push_back(sqrt(tStep * totalSurfaceSqError));

	if (m_EvolSettings.DoRemeshing)
	{
		m_MeanEdgeLength /= NSteps;
		m_EvolvingSurface->write(m_EvolSettings.OutputPath + "SphereTest" + (m_EvolSettings.TangentialVelocityWeight > 0.0 ? "_redist_" : "") + std::to_string(iter) + ".vtk");
	}
}

/// \brief exact ground truth solution.
[[nodiscard]] double ExactShrinkingSphereRadius(const double& time)
{
	return sqrt(STARTING_SPHERE_RADIUS * STARTING_SPHERE_RADIUS - 4.0 * time);
}

double SphereTest::EvaluateSquaredSurfaceError(const double& time) const
{
	const double exactRad = ExactShrinkingSphereRadius(time);
	double result = 0.0;
	for (const auto v : m_EvolvingSurface->vertices())
	{
		const auto vPos = m_EvolvingSurface->position(v);
		const auto numericalRad = static_cast<double>(pmp::norm(vPos));
		const double diff = numericalRad - exactRad;
		const double coVolMeasure = m_LaplacianAreaFunction(*m_EvolvingSurface, v);
		result += diff * diff * coVolMeasure;
	}
	return result;
}

// ================================================================================================

double SphereTest::ComputeMeanEdgeLength() const
{
	if (!m_EvolvingSurface)
		throw std::invalid_argument("SphereTest::Evolve: m_EvolvingSurface not set! Terminating!\n");

	double result = 0.0;
	for (const auto e : m_EvolvingSurface->edges())
		result += static_cast<double>(m_EvolvingSurface->edge_length(e));

	return result / static_cast<double>(m_EvolvingSurface->n_edges());
}

void SphereTest::PerformTest(const unsigned int& nIterations)
{
	std::cout << "......................................................................\n";
	std::cout << " < < < < < < < < < < < Sphere Test Results: > > > > > > > > > > > > > \n";
	std::cout << "......................................................................\n";
	double tStep = STARTING_TIME_STEP_SIZE;
	// -------------------------------
	/*
	const double phi = (1.0 + sqrt(5.0)) / 2.0; /// golden ratio.
	const auto subdiv = static_cast<double>(STARTING_SPHERE_SUBDIVISION);
	const auto meanR = sqrt(19.0) / 10.0 * STARTING_SPHERE_RADIUS;
	double prevMeanEdgeLength = (m_EvolSettings.DoRemeshing ? (meanR / (sqrt(phi * sqrt(5.0)) * pow(2.0, subdiv - 1))) : 2.0);*/
	double prevMeanEdgeLength{ 1.0 };
	for (unsigned int i = 0; i < nIterations; i++)
	{
		Evolve(i, tStep);
		if (m_EvolSettings.DoRemeshing && i == 0)
			prevMeanEdgeLength = m_MeanEdgeLength;

		const auto nVertices = m_EvolvingSurface->n_vertices();
		std::cout << ":: SurfaceError: " << m_MeanErrors[i] << ", nVertices : " << nVertices;
		if (i > 0)
		{
			// Evaluate the (experimental) order of convergence: EOC = log2( previousIterMeanError / currentIterMeanError );
			const double errorRatio = m_MeanErrors[i - 1] / m_MeanErrors[i];
			if (m_EvolSettings.DoRemeshing)
			{
				const double meanLengthRatio = prevMeanEdgeLength / m_MeanEdgeLength;
				const double baseChangeLog = log2(meanLengthRatio);
				m_EOCs.push_back(log2(errorRatio) / baseChangeLog);
				std::cout << ", EOC(" << meanLengthRatio << "): " << m_EOCs[i - 1] << "\n";
				prevMeanEdgeLength = m_MeanEdgeLength;
				tStep /= 4.0; // reduce time step by 2^2
				continue;
			}
			m_EOCs.push_back(log2(errorRatio));
			std::cout << ", EOC: " << m_EOCs[i - 1];
		}
		std::cout << "\n";
		tStep /= 4.0; // reduce time step by 2^2
	}

}
// ================================================================================================


void ReportInput(const SphereTestEvolutionSettings& settings, std::ostream& os)
{
	os << "======================================================================\n";
	os << "> > > > > > > > > > > Initiating SphereTest: < < < < < < < < < < < < <\n";
	os << "......................................................................\n";
	os << "Export Surface per Time Step: " << (settings.ExportSurfacePerTimeStep ? "true" : "false") << ",\n";
	os << "Output Path: " << settings.OutputPath << ",\n";
	os << "Do Remeshing: " << (settings.DoRemeshing ? "true" : "false") << ",\n";
	os << "----------------------------------------------------------------------\n";
}
