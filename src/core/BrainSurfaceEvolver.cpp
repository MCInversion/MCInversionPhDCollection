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

#include "EvolverUtilsCommon.h"
//#include "ConversionUtils.h"

// ================================================================================================

/// \brief a magic multiplier computing the radius of an ico-sphere that fits into the field's box.
constexpr float ICO_SPHERE_RADIUS_FACTOR = 0.2f;

/// \brief this value is verified for [Smith, 2002]
//constexpr double BET_NORMAL_INTENSITY_FACTOR = 0.05;

// Custom values
//constexpr double BET_NORMAL_INTENSITY_FACTOR = -15.0;
constexpr double BET_NORMAL_INTENSITY_FACTOR = 1.0;
//constexpr double BET_NORMAL_INTENSITY_FACTOR = 100.0;

constexpr double CURVATURE_INTENSITY_FACTOR = 0.5;

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

// ================================================================================================

void BrainSurfaceEvolver::Preprocess()
{
	auto& field = *m_Field;
	const auto& fieldBox = field.Box();
	const auto fieldBoxSize = fieldBox.max() - fieldBox.min();
	const float minDim = std::min({ fieldBoxSize[0], fieldBoxSize[1], fieldBoxSize[2] });
	const float maxDim = std::max({ fieldBoxSize[0], fieldBoxSize[1], fieldBoxSize[2] });

	// build ico-sphere
	const float icoSphereRadius = ICO_SPHERE_RADIUS_FACTOR * (minDim + maxDim);
	//const float icoSphereRadius = m_EvolSettings.IcoSphereSettings.Radius;
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
	m_EvolvingSurfaceRadiusEstimate = static_cast<double>(icoSphereRadius * scalingFactor);
	//m_EvolSettings.IcoSphereSettings.Radius *= scalingFactor;

	m_UnitNormalToGridScaleFactor = static_cast<float>(m_EvolvingSurfaceRadiusEstimate) / (maxDim * m_ScalingFactor);

	// origin needs to be computed from bet2
	const auto origin = m_EvolSettings.IcoSphereSettings.Center;

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

// ================================================================================================

void BrainSurfaceEvolver::UpdateRadiusEstimate()
{
	if (!m_EvolvingSurface)
		throw std::invalid_argument("BrainSurfaceEvolver::ComputeRadiusEstimate: m_EvolvingSurface is null! Terminating!\n");

	const auto surfaceBBox = m_EvolvingSurface->bounds();
	const auto surfaceBBoxSize = surfaceBBox.max() - surfaceBBox.min();
	const double minDim = std::min({ surfaceBBoxSize[0], surfaceBBoxSize[1], surfaceBBoxSize[2] });
	const double maxDim = std::max({ surfaceBBoxSize[0], surfaceBBoxSize[1], surfaceBBoxSize[2] });

	m_EvolvingSurfaceRadiusEstimate = 0.5 * (minDim + maxDim);
}

void BrainSurfaceEvolver::ComputeNormalIntensities()
{
}

// ================================================================================================

constexpr double MAGIC_RADIUS_CONSTANT = 0.0;

double BrainSurfaceEvolver::LaplacianWeightFunction() const
{
	return 1.0;

	//const double bet2Radius = m_EvolSettings.IcoSphereSettings.Radius;
	//return m_EvolvingSurfaceRadiusEstimate - MAGIC_RADIUS_CONSTANT * bet2Radius;
}

double BrainSurfaceEvolver::NormalIntensityWeightFunction(const pmp::vec3& vPos, const pmp::vec3& vNormal) const
{
	const auto cellSize = m_Field->CellSize();
	const auto& fieldBox = m_Field->Box();
	const auto fieldOrigin = fieldBox.min();
	const auto& fieldValues = m_Field->Values();
	const auto [Nx, Ny, Nz] = m_Field->Dimensions();
	const auto gridExtent = Nx * Ny * Nz;

	double Imin = m_EvolSettings.ThresholdSettings.ThresholdEffectiveMedian; // tm
	double Imax = m_EvolSettings.ThresholdSettings.ThresholdEffective; // t

	// need to use scaling factor to access the correct grid voxels from normals
	const auto vShiftedByNormal = vPos - m_ScalingFactor * vNormal;

	// transform vertex from real space to grid index space
	const auto ix0 = static_cast<int>(std::floor((vShiftedByNormal[0] - fieldOrigin[0]) / cellSize));
	const auto iy0 = static_cast<int>(std::floor((vShiftedByNormal[1] - fieldOrigin[1]) / cellSize));
	const auto iz0 = static_cast<int>(std::floor((vShiftedByNormal[2] - fieldOrigin[2]) / cellSize));
	if (ix0 < 0 || iy0 < 0 || iz0 < 0 || ix0 >= Nx || iy0 >= Ny || iz0 >= Nz)
		return 0.0; // Leave if indices are out of bounds.
	const unsigned int gridPos0 = Nx * Ny * iz0 + Nx * iy0 + ix0;
	const double im0 = fieldValues[gridPos0];

	if (im0 < Imin) Imin = im0;
	if (im0 > Imax) Imax = im0;

	const auto d1 = m_EvolSettings.ThresholdSettings.MinIntensitySearchDepth;
	const auto d2 = m_EvolSettings.ThresholdSettings.MaxIntensitySearchDepth;

	const float normalStepFactor = m_UnitNormalToGridScaleFactor / static_cast<float>(d2 - d1);
	const auto iNormalStepX = static_cast<int>(std::round(vNormal[0] * normalStepFactor));
	const auto iNormalStepY = static_cast<int>(std::round(vNormal[1] * normalStepFactor));
	const auto iNormalStepZ = static_cast<int>(std::round(vNormal[2] * normalStepFactor));

	const auto ix1 = ix0 - (d1 - 1) * iNormalStepX;
	const auto iy1 = iy0 - (d1 - 1) * iNormalStepY;
	const auto iz1 = iz0 - (d1 - 1) * iNormalStepZ;
	if (ix1 < 0 || iy1 < 0 || iz1 < 0 || ix1 >= Nx || iy1 >= Ny || iz1 >= Nz)
		return 0.0; // Leave if indices are out of bounds.
	const unsigned int gridPos1 = Nx * Ny * iz1 + Nx * iy1 + ix1;
	const double im1 = fieldValues[gridPos1];

	if (im1 < Imin) Imin = im1;

	// search further in normal direction
	auto ix = ix0;
	auto iy = iy0;
	auto iz = iz0;
	//const auto cellSizeCubed = pow(cellSize, 3);
	for (double gi = 2.0; gi < d1; gi += 1.0)
	{
		ix -= iNormalStepX;
		iy -= iNormalStepY;
		iz -= iNormalStepZ;
		// Check bounds to ensure indices are within the grid dimensions.
		if (ix < 0 || iy < 0 || iz < 0 || ix >= Nx || iy >= Ny || iz >= Nz)
			break; // Break out of the loop if indices are out of bounds.

		const unsigned int gridPos = Nx * Ny * iz + Nx * iy + ix;
		if (gridPos >= gridExtent)
			continue;
		const double im = fieldValues[gridPos];
		// Update Imin and Imax based on the intensity value at the current position.
		Imin = std::min(Imin, im);
		// Only update Imax within the range specified by d2.
		if (gi < d2)
		{
			Imax = std::max(Imax, im);
		}
	}

	const auto t2 = m_EvolSettings.ThresholdSettings.Threshold2ndPercentile;
	const auto tm = m_EvolSettings.ThresholdSettings.ThresholdEffectiveMedian;
	if (t2 < Imin) Imin = t2;
	if (tm > Imax) Imax = tm;

	const auto b = m_EvolSettings.ThresholdSettings.BetMainParam;
	const double tl = (Imax - t2) * b + t2; // local threshold value

	double f3 = 2.0 * (Imin - tl);
	if (Imax - t2 > 0.0) 
		f3 /= (Imax - t2);

	return f3;
}

// ================================================================================================

/// \brief a utility for computing the mean value of mean edge lengths for each vertex neighborhood.
[[nodiscard]] float ComputeMeanInterVertexDistance(const pmp::SurfaceMesh& mesh)
{
	float result{ 0.0f };
	for (const auto v : mesh.vertices())
	{
		const auto vPos = mesh.position(v);
		float neighborhoodResult{ 0.0 };
		for (const auto w : mesh.vertices(v))
		{
			const auto wPos = mesh.position(w);
			neighborhoodResult += pmp::norm(vPos - wPos);
		}
		const auto vValence = mesh.valence(v);
		neighborhoodResult /= static_cast<float>(vValence);
		result += neighborhoodResult;
	}
	result /= static_cast<float>(mesh.n_vertices());

	return result;
}

// ================================================================================================

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

// ================================================================================================

void BrainSurfaceEvolver::Evolve()
{
	if (!m_Field)
		throw std::invalid_argument("SurfaceEvolver::Evolve: m_Field not set! Terminating!\n");
	if (!m_Field->IsValid())
		throw std::invalid_argument("SurfaceEvolver::Evolve: m_Field is invalid! Terminating!\n");

	const auto& field = *m_Field;

	Preprocess();

	if (!m_EvolvingSurface)
		throw std::invalid_argument("SurfaceEvolver::Evolve: m_EvolvingSurface not set! Terminating!\n");

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

	// DISCLAIMER: the dimensionality of the system depends on the number of mesh vertices which can change if remeshing is used.
	auto NVertices = static_cast<unsigned int>(m_EvolvingSurface->n_vertices());
	SparseMatrix sysMat(NVertices, NVertices);
	Eigen::MatrixXd sysRhs(NVertices, 3);
	auto vIntensity = m_EvolvingSurface->vertex_property<pmp::Scalar>("v:normalIntensity", 0.0f);
	auto vFeature = m_EvolvingSurface->vertex_property<bool>("v:feature", false);

	// property container for surface vertex normals
	pmp::VertexProperty<pmp::Point> vNormalsProp{};

	// ----------- System fill function --------------------------------
	const auto fillMatrixAndRHSTriplesFromMesh = [&]()
	{
		const float meanInterVertexDistance = ComputeMeanInterVertexDistance(*m_EvolvingSurface);

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
			
			const auto vNormal = vNormalsProp[v]; // vertex unit normal

			//const double absIntensity = std::abs(vIntensity[v]);
			const double epsilonCtrlWeight = CURVATURE_INTENSITY_FACTOR; // (absIntensity > 0.0 ? absIntensity : 1.0);
			const double etaCtrlWeight = BET_NORMAL_INTENSITY_FACTOR * meanInterVertexDistance * vIntensity[v];

			const Eigen::Vector3d vertexRhs = vPosToUpdate + tStep * etaCtrlWeight * vNormal /* + tStep * tanRedistWeight * vTanVelocity */;
			sysRhs.row(v.idx()) = vertexRhs;

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
	std::cout << "Co-Volume Measure Stats: { Mean: " << coVolStats.Mean << ", Min: " << coVolStats.Min << ", Max: " << coVolStats.Max << "},\n";
#endif
	// set initial surface vertex properties
	{
		pmp::Normals::compute_vertex_normals(*m_EvolvingSurface);
		vNormalsProp = m_EvolvingSurface->vertex_property<pmp::Point>("v:normal");
		for (const auto v : m_EvolvingSurface->vertices())
		{
			const auto vPos = m_EvolvingSurface->position(v);
			vIntensity[v] = static_cast<pmp::Scalar>(NormalIntensityWeightFunction(vPos, vNormalsProp[v]));
		}
	}
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
			remeshing.adaptive_remeshing({
				minEdgeLength, maxEdgeLength, 2.0f * minEdgeLength,
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
			}
		}

		coVolStats = AnalyzeMeshCoVolumes(*m_EvolvingSurface, m_LaplacianAreaFunction);
#if REPORT_EVOL_STEPS
		std::cout << "Co-Volume Measure Stats: { Mean: " << coVolStats.Mean << ", Min: " << coVolStats.Min << ", Max: " << coVolStats.Max << "},\n";
#endif
		// set surface vertex properties
		for (const auto v : m_EvolvingSurface->vertices())
		{
			const auto vPos = m_EvolvingSurface->position(v);
			vIntensity[v] = static_cast<pmp::Scalar>(NormalIntensityWeightFunction(vPos, vNormalsProp[v]));
		}
		ComputeTriangleMetrics();
		UpdateRadiusEstimate();

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
}

// ================================================================================================

void ReportInput(const BrainExtractionSettings& evolSettings, std::ostream& os)
{
	os << "======================================================================\n";
	os << "> > > > > > > > Initiating BrainSurfaceEvolver: < < < < < < < < < < < \n";
	os << "Brain Name: " << evolSettings.ProcedureName << ",\n";
	os << "NSteps: " << evolSettings.NSteps << ",\n";
	os << "TimeStep: " << evolSettings.TimeStep << ",\n";
	os << "IcoSphereSubdivisionLevel: " << evolSettings.IcoSphereSubdivisionLevel << ",\n";
	os << "......................................................................\n";
	os << "WIP: not printing threshold and ico-sphere BET values ... \n";
	os << "......................................................................\n";
	os << "Export Surface per Time Step: " << (evolSettings.ExportSurfacePerTimeStep ? "true" : "false") << ",\n";
	os << "Output Path: " << evolSettings.OutputPath << ",\n";
	os << "Do Remeshing: " << (evolSettings.DoRemeshing ? "true" : "false") << ",\n";
	os << "----------------------------------------------------------------------\n";
}
