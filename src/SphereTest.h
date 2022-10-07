#pragma once

#include "pmp/SurfaceMesh.h"
#include "pmp/algorithms/DifferentialGeometry.h"

/**
 * \brief A wrapper for parameters related to mesh topology adjustments (remeshing etc.).
 * \struct ST_MeshTopologySettings
 */
struct ST_MeshTopologySettings
{
	float MinEdgeMultiplier{ 0.07f }; //>! multiplier for minimum edge length in adaptive remeshing.
	double RemeshingStartTimeFactor{ 0.1 }; //>! the fraction of total time steps after which remeshing should take place.
	float EdgeLengthDecayFactor{ 0.98f }; //>! decay factor for minimum (and consequently maximum) edge length.
	double RemeshingSizeDecayStartTimeFactor{ 0.2 }; //>! decay of edge length bounds should take place after (this value) * NSteps of evolution.
	unsigned int StepStrideForEdgeDecay{ 5 }; //>! the number of steps after which edge length bound decay takes place.
	double FeatureDetectionStartTimeFactor{ 0.4 }; //>! feature detection becomes relevant after (this value) * NSteps.
	unsigned int NRemeshingIters{ 2 }; //>! the number of iterations for pmp::Remeshing.
	unsigned int NTanSmoothingIters{ 5 }; //>! the number of tangential smoothing iterations for pmp::Remeshing.
	bool UseBackProjection{ false }; //>! if true surface kd-tree back-projection will be used for pmp::Remeshing.
};

/// \brief An enumerator for the choice of mesh Laplacian scheme [Meyer, Desbrun, Schroder, Barr, 2003].
enum class [[nodiscard]] ST_MeshLaplacian
{
	Voronoi = 0, //>! the finite volume for mesh Laplacian is generated from a true Voronoi neighborhood of each vertex.
	Barycentric = 1 //>! the finite volume for mesh Laplacian is generated from face barycenters.
};

/// \brief a list of triangle metrics to be computed.
using TriangleMetrics = std::vector<std::string>;

/**
 * \brief A wrapper for sphere test settings.
 * \struct SphereTestEvolutionSettings
 */
struct SphereTestEvolutionSettings
{
	ST_MeshTopologySettings TopoParams{}; //>! parameters for mesh topology adjustments.
	bool ExportSurfacePerTimeStep{ false }; //>! whether to export evolving surface for each time step.
	bool ExportResultSurface{ false }; //>! whether to export resulting evolving surface.
	std::string OutputPath{}; //>! path where output surfaces are to be exported.

	ST_MeshLaplacian LaplacianType{}; //>! type of mesh Laplacian.

	float TangentialVelocityWeight{ 0.0f }; //>! the weight of tangential velocity update vector for each time step.
	bool DoRemeshing{ true }; //>! if true, adaptive remeshing will be performed after the first 10-th of time steps.

	double MaxFractionOfVerticesOutOfBounds{ 0.02 }; //>! fraction of vertices allowed to be out of bounds (because it will be decimated).
};

/**
 * \brief A variant of surface evolution testing the the convergence to compact ground truth exact solution - shrinking sphere with radius r(t) = sqrt(r0 - 4t).
 * \class SphereTest
 */
class SphereTest
{
public:
	/**
	 * \brief Constructor. Initialize from settings.
	 * \param settings         sphere test settings.
	 */
	explicit SphereTest(const SphereTestEvolutionSettings& settings);

	/**
	 * \brief Main functionality.
	 * \param nIterations       the number of time step subdivisions to be tested.
	 */
	void PerformTest(const unsigned int& nIterations = 4);


private:
	/**
	 * \brief Preprocess for evolution, i.e.: generate m_EvolvingSurface, and transform both m_Field and m_EvolvingSurface for stabilization.
	 * \param iter    index of the iteration for which to preprocess.
	 */
	void Preprocess(const unsigned int& iter);

	/**
	 * \brief Perform a particular evolution iteration.
	 * \param iter    if 0, time step will be equal to base time step, and ico-sphere subdivision level will also correspond to base settings.
	 * \param tStep   time step size.
	 */
	void Evolve(const unsigned int& iter, const double& tStep);

	// ----------------------------------------------------------------

	/**
	 * \brief Evaluate the sum of squared differences between mesh surface and exact solution (for a particular time step) weighed by co-volume measures (m_LaplacianAreaFunction).
	 * \param time       input time (not time index, but ti * tStep).
	 * \return error for the particular time.
	 */
	[[nodiscard]] double EvaluateSquaredSurfaceError(const double& time) const;

	// ----------------------------------------------------------------

	SphereTestEvolutionSettings m_EvolSettings{}; //>! settings.
	const double m_TotalTime; //>! total time into which all time steps need to fit for the sphere test to be numerically valid.
	std::shared_ptr<pmp::SurfaceMesh> m_EvolvingSurface{ nullptr }; //>! (stabilized) evolving surface.

	std::function<pmp::ImplicitLaplaceInfo(const pmp::SurfaceMesh&, pmp::Vertex)> m_ImplicitLaplacianFunction{}; //>! a Laplacian function chosen from parameter MeshLaplacian.
	std::function<double(const pmp::SurfaceMesh&, pmp::Vertex)> m_LaplacianAreaFunction{}; //>! a Laplacian area function chosen from parameter MeshLaplacian.

	// export
	std::string m_OutputMeshExtension = ".vtk"; //>! extension of the exported mesh geometry.

	// error analysis
	std::vector<double> m_MeanErrors{}; //>! the mean absolute difference between mesh surface and exact solution for all time steps, for all iterations.
	std::vector<double> m_EOCs{}; //>! Experimental order of convergence, for all iterations.
};

/**
 * \brief Reports SphereTest's input to a given stream.
 * \param settings        settings for SphereTest.
 * \param os              output stream.
 */
void ReportInput(const SphereTestEvolutionSettings& settings, std::ostream& os);