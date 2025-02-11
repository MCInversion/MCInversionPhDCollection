#pragma once

#include "pmp/ManifoldCurve2D.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "geometry/MeshAnalysis.h"

#include "utils/nlohmann/json.hpp"

#include <ranges>
#include <fstream>
#include <iomanip>
#include <functional>
#include <unordered_set>

namespace pmp
{
	class SurfaceMesh;
	class ManifoldCurve2D;
	class Vertex;
	struct AdaptiveRemeshingSettings;
}

/// \brief a stats wrapper for co-volume measures affecting the stability of the finite volume method.
struct CoVolumeStats
{
	double Mean{ 0.0 };
	double Max{ -DBL_MAX };
	double Min{ DBL_MAX };
};

const std::string coVolMeasureVertexPropertyName{ "v:coVolumeMeasure" };

/// \brief a co-volume area evaluator.
using AreaFunction = std::function<double(const pmp::SurfaceMesh&, pmp::Vertex)>;

/**
 * \brief Analyzes the stats of co-volumes around each mesh vertex, and creates a vertex property for the measure values.
 * \param mesh             input mesh.
 * \param areaFunction     function to evaluate vertex co-volume area.
 * \return co-volume stats.
 */
[[nodiscard]] CoVolumeStats AnalyzeMeshCoVolumes(pmp::SurfaceMesh& mesh, const AreaFunction& areaFunction);

/*
 * \brief Computes scaling factor for stabilizing the finite volume method on assumed spherical surface meshes based on time step.
 * \param timeStep               time step size.
 * \param icoRadius              radius of an evolving geodesic icosahedron.
 * \param icoSubdiv              subdivision level of an evolving geodesic icosahedron.
 * \param stabilizationFactor    a multiplier for stabilizing mean co-volume area.
 * \return scaling factor for mesh and scalar grid.
 */
[[nodiscard]] pmp::Scalar GetStabilizationScalingFactor(const double& timeStep, const pmp::Scalar& icoRadius, const unsigned int& icoSubdiv, const pmp::Scalar& stabilizationFactor = 1.0);

/// \brief identifier for sparse matrix.
using SparseMatrix = Eigen::SparseMatrix<double>;

/**
 * @brief Utility to print the entire SparseMatrix, including zeros, in a matrix-like format.
 *
 * @param matrix The sparse matrix to print.
 */
inline void PrintSparseMatrix(const SparseMatrix& matrix)
{
	std::cout << "SparseMatrix (" << matrix.rows() << " x " << matrix.cols() << "):" << std::endl;

	// Set output formatting for better readability
	std::cout << std::fixed << std::setprecision(6);

	for (int i = 0; i < matrix.rows(); ++i)
	{
		for (int j = 0; j < matrix.cols(); ++j)
		{
			// Get the value at (i, j)
			double value = matrix.coeff(i, j);
			std::cout << std::setw(10) << value << " ";
		}
		std::cout << "\n";
	}
}

/**
 * @brief Utility to print the entire SparseMatrix, including zeros, in a matrix-like format,
 *        alongside the right-hand side vector/matrix (rhs).
 *
 * @param matrix The sparse matrix to print.
 * @param rhs The right-hand side matrix to print alongside (should have 2 columns).
 * @param matrixName    an optional identifier for matrix printout
 */
inline void PrintSparseMatrixAndRHS(const SparseMatrix& matrix, const Eigen::MatrixXd& rhs, const std::string& matrixName = "")
{
	// Validate that the RHS matrix has the correct dimensions
	if (rhs.rows() != matrix.rows() || rhs.cols() != 2)
	{
		throw std::invalid_argument("The RHS matrix must have the same number of rows as the sparse matrix and exactly 2 columns.");
	}

	std::cout << "\n-----------------------------------------------------------------\n";
	std::cout << (matrixName.empty() ? "Matrix" : matrixName) << " (" << matrix.rows() << " x " << matrix.cols() << "):\n";

	// Set output formatting for better readability
	std::cout << std::fixed << std::setprecision(6);

	for (int i = 0; i < matrix.rows(); ++i)
	{
		for (int j = 0; j < matrix.cols(); ++j)
		{
			// Get the value at (i, j)
			double value = matrix.coeff(i, j);
			std::cout << std::setw(10) << value << " ";
		}

		// Print the RHS values for the current row
		std::cout << "| " << std::setw(10) << rhs(i, 0) << " " << std::setw(10) << rhs(i, 1) << "\n";
	}
	std::cout << "\n-----------------------------------------------------------------\n";
}

/**
 * @brief Check if a sparse matrix is diagonally dominant.
 *
 * @param matrix The sparse matrix to check.
 * @return true if the matrix is diagonally dominant, false otherwise.
 */
inline bool IsDiagonallyDominant(const SparseMatrix& matrix)
{
	if (matrix.rows() != matrix.cols())
	{
		std::cerr << "Matrix must be square to check diagonal dominance!" << std::endl;
		return false;
	}

	for (int i = 0; i < matrix.rows(); ++i)
	{
		double diagonalValue = std::abs(matrix.coeff(i, i));
		double offDiagonalSum = 0.0;

		for (SparseMatrix::InnerIterator it(matrix, i); it; ++it)
		{
			if (it.row() != it.col()) // Skip diagonal element
			{
				offDiagonalSum += std::abs(it.value());
			}
		}

		if (diagonalValue < offDiagonalSum)
		{
			return false; // Not diagonally dominant
		}
	}

	return true; // All rows satisfy the condition
}

/**
 * @brief Regularize a sparse matrix in place by adding lambda times the identity matrix.
 *
 * @param matrix The sparse matrix to regularize (modified in place).
 * @param lambda The regularization parameter.
 */
inline void RegularizeMatrixInPlace(SparseMatrix& matrix, double lambda)
{
	// Ensure the matrix is square
	if (matrix.rows() != matrix.cols())
	{
		throw std::invalid_argument("Matrix must be square for regularization!");
	}

	// Add lambda * I to the diagonal in place
	for (int i = 0; i < matrix.rows(); ++i)
	{
		matrix.coeffRef(i, i) += lambda;
	}
}

/// \brief A utility for converting Eigen::ComputationInfo to a string message.
[[nodiscard]] std::string InterpretSolverErrorCode(const Eigen::ComputationInfo& cInfo);

/**
 * \brief Computes angle-based tangential update velocity for a mesh vertex with a given weight.
 * \param mesh       a surface mesh to whom the vertex belongs.
 * \param v          vertex handle of a vertex where the tangential velocity is to be computed.
 * \param vNormal    normal to vertex with handle v.
 * \param weight     weight of the velocity vector.
 * \return tangential velocity vector.
 */
[[nodiscard]] pmp::vec3 ComputeTangentialUpdateVelocityAtVertex(const pmp::SurfaceMesh& mesh, const pmp::Vertex& v, const pmp::vec3& vNormal, const pmp::Scalar& weight = 1.0);

/**
 * \brief Computes tangential update velocity for a mesh vertex with a given weight.
 * \param mesh       a manifold curve to whom the vertex belongs.
 * \param v          vertex handle of a vertex where the tangential velocity is to be computed.
 * \param vNormal    normal to vertex with handle v.
 * \param weight     weight of the velocity vector.
 * \return tangential velocity vector.
 */
[[nodiscard]] pmp::vec2 ComputeTangentialUpdateVelocityAtVertex(const pmp::ManifoldCurve2D& curve, const pmp::Vertex& v, const pmp::vec2& vNormal, const pmp::Scalar& weight = 1.0);

// ======================================================================================================================

/// \brief An enumerator for the choice of mesh Laplacian scheme [Meyer, Desbrun, Schroder, Barr, 2003].
enum class [[nodiscard]] MeshLaplacian
{
	Voronoi = 0, //>! the finite volume for mesh Laplacian is generated from a true Voronoi neighborhood of each vertex.
	Barycentric = 1 //>! the finite volume for mesh Laplacian is generated from face barycenters.
};

/// \brief a list of triangle metrics to be computed.
using TriangleMetrics = std::vector<std::string>;

/**
 * \brief A wrapper for surface evolution advection-diffusion parameters.
 * \struct AdvectionDiffusionParameters
 */
struct AdvectionDiffusionParameters
{
	// constants for weight function:
	// C1 * (1 - exp(d^2 / C2)), where d is distance
	double MCFMultiplier{ 1.0 }; // C1
	double MCFVariance{ 1.0 }; // C2

	// constants for weight function:
	// D1 * d * ((-grad(d) . N) - D2 * sqrt(1 - (grad(d) . N)^2)), where d is distance, and N unit normal.
	double AdvectionMultiplier{ 1.0 }; // D1
	double AdvectionSineMultiplier{ 1.0 }; // D2

	bool MCFSupportPositive{ true }; //>! if true, the diffusion weight function has non-zero support for positive values only.
	bool AdvectionSupportPositive{ true }; //>! if true, the advection weight function has non-zero support for positive values only.
};

/**
 * \brief An enumerator for the type of feature detection function used during evolution.
 * \enum FeatureDetectionType
 */
enum class [[nodiscard]] FeatureDetectionType
{
	Angle = 0, //>! edges with dihedral angle larger than a given threshold value are marked as features.
	AngleWithinBounds = 1, //>! edges with dihedral angle within a given range are marked as features.
	PrincipalCurvatures = 2, //>! vertices with too much imbalance in principal curvatures are marked as features.
	MeanCurvature = 3 //>! features are preferred for vertices with high positive mean curvature.
};

/**
 * \brief A wrapper for parameters related to mesh topology adjustments (remeshing etc.).
 * \struct MeshTopologySettings
 */
struct MeshTopologySettings
{
	bool FixSelfIntersections{ true }; //>! if true, self-intersecting faces within the evolving surface will be removed, and the holes will be patched by pmp::HoleFilling.
	pmp::Scalar MinEdgeMultiplier{ 0.14 }; //>! multiplier for minimum edge length in adaptive remeshing.
	double RemeshingStartTimeFactor{ 0.1 }; //>! the fraction of total time steps after which remeshing should take place.
	pmp::Scalar EdgeLengthDecayFactor{ 0.98 }; //>! decay factor for minimum (and consequently maximum) edge length.
	double RemeshingSizeDecayStartTimeFactor{ 0.2 }; //>! decay of edge length bounds should take place after (this value) * NSteps of evolution.
	unsigned int StepStrideForEdgeDecay{ 5 }; //>! the number of steps after which edge length bound decay takes place.
	double FeatureDetectionStartTimeFactor{ 0.4 }; //>! feature detection becomes relevant after (this value) * NSteps.
	unsigned int NRemeshingIters{ 2 }; //>! the number of iterations for pmp::Remeshing.
	unsigned int NTanSmoothingIters{ 5 }; //>! the number of tangential smoothing iterations for pmp::Remeshing.
	bool UseBackProjection{ true }; //>! if true surface kd-tree back-projection will be used for pmp::Remeshing.

	FeatureDetectionType FeatureType{ FeatureDetectionType::MeanCurvature }; //>! type of feature detection function.
	double MinDihedralAngle{ 1.0 * M_PI_2 * 180.0 }; //>! critical dihedral angle for feature detection
	double MaxDihedralAngle{ 2.0 * M_PI_2 * 180.0 }; //>! critical dihedral angle for feature detection
	pmp::Scalar PrincipalCurvatureFactor{ 2.0 }; //>! vertices with |Kmax| > \p principalCurvatureFactor * |Kmin| are marked as feature.
	pmp::Scalar CriticalMeanCurvatureAngle{ 1.0 * static_cast<pmp::Scalar>(M_PI_2) }; //>! vertices with curvature angles smaller than this value are feature vertices. 
	bool ExcludeEdgesWithoutBothFeaturePts{ false }; //>! if true, edges with only one vertex detected as feature will not be marked as feature.
};

/**
 * \brief Precomputes parameters for advection-diffusion model within (Iso)SurfaceEvolver.
 * \param distanceMax             maximum effective distance from target (affects diffusion term weight).
 * \param targetMinDimension      minimum target size (affects diffusion term variance).
 * \return desired advection-diffusion parameters.
 */
[[nodiscard]] AdvectionDiffusionParameters PreComputeAdvectionDiffusionParams(const double& distanceMax, const double& targetMinDimension);

/// \brief Evaluates whether remeshing is necessary from the co-volume stats and time step.
[[nodiscard]] bool IsRemeshingNecessary(const CoVolumeStats& stats, const double& tStep);

/// \brief Evaluates whether remeshing is necessary from the condition number metric for equilateral triangles.
[[nodiscard]] bool IsRemeshingNecessary(const std::vector<pmp::Scalar>& equilateralJacobianConditionNumbers);

/// \brief Evaluates whether remeshing is necessary from the condition number metric for equilateral triangles that do not have a feature vertex.
[[nodiscard]] bool IsNonFeatureRemeshingNecessary(const pmp::SurfaceMesh& mesh);

/// \brief Evaluates whether remeshing is necessary from the system matrix spectral radius
[[nodiscard]] bool IsRemeshingNecessary(const SparseMatrix& lswMatrix);

/// \brief Evaluates whether remeshing is necessary for a manifold curve.
///	Uses midpoint co-volumes around each vertex.
[[nodiscard]] bool IsRemeshingNecessary(const pmp::ManifoldCurve2D& curve, const pmp::AdaptiveRemeshingSettings& remeshingSettings);

/// \brief Evaluates whether remeshing is necessary for a surface mesh based on vertex density.
///	\param[in] mesh                  input mesh.
///	\param[in] remeshingSettings     settings to be used.
///	\param[in] areaFunc              control volume measure evaluation.
///	\return true if remeshing is necessary.
[[nodiscard]] bool IsRemeshingNecessary(const pmp::SurfaceMesh& mesh, const pmp::AdaptiveRemeshingSettings& remeshingSettings, const AreaFunction& areaFunc);

/// \brief Evaluates whether remeshing is necessary for a surface mesh.
[[nodiscard]] bool IsRemeshingNecessary(const pmp::SurfaceMesh& mesh, const Geometry::FaceQualityFunction& qualityFunc, const Geometry::FaceQualityRange& qualityRange);

/// \brief A (one-time) evaluation whether the distance to target reaches a lower bound.
///	\param distancePerVertexValues    a vector of distance values on the evolving surface.
///	\return true if the conditions for feature detection are satisfied.
[[nodiscard]] bool ShouldDetectFeatures(const std::vector<pmp::Scalar>& distancePerVertexValues);

/// \brief Sets a static container for time indices for a particular evolver setup.
void SetRemeshingAdjustmentTimeIndices(const std::unordered_set<unsigned int>& valuesSet);

/// \brief A getter for the global container for time indices during which remeshing sizing should be changed.
std::unordered_set<unsigned int>& GetRemeshingAdjustmentTimeIndices();

/// 
/// \brief Evaluates whether the target edge lengths for adaptive remeshing should be decreased.
/// \param ti               time index from 1 to NSteps.
/// \return true if the adjustment should take place.
///
[[nodiscard]] bool ShouldAdjustRemeshingLengths(const unsigned int& ti /*, const unsigned int& NSteps*/);

///
/// \brief Adjusts edge lengths for adaptive remeshing.
///	\param decayFactor      decay factor from [0, 1] for edge length.
///	\param minEdgeLength    the minimum edge length to be adjusted.
///	\param maxEdgeLength    the maximum edge length to be adjusted.
///	\param approxError      approximation error to be adjusted.
///
void AdjustRemeshingLengths(const pmp::Scalar& decayFactor, pmp::Scalar& minEdgeLength, pmp::Scalar& maxEdgeLength, pmp::Scalar& approxError);

/**
 * \brief Logs manifolds which need remeshing.
 * \tparam ManifoldType  either pmp::ManifoldCurve2D or pmp::SurfaceMesh.
 */
template<typename ManifoldType>
class ManifoldsToRemeshTracker
{
public:
	/// \brief Add a manifold that requires remeshing
	void AddManifold(ManifoldType* manifoldPtr)
	{
		if (!manifoldPtr)
		{
			std::cerr << "ManifoldsToRemeshTracker::AddManifold: attempting to log a null manifold!\n";
			return;
		}
		m_Manifolds.push_back(manifoldPtr);
	}

	/// \brief Get all manifolds that need remeshing
	[[nodiscard]] std::vector<ManifoldType*> GetManifoldsToRemesh() const
	{
		return m_Manifolds;
	}

	/// \brief Clear all stored manifolds (reset for a new time step)
	void Reset()
	{
		m_Manifolds.clear();
	}

private:
	// TODO: use something like: bool m_IsReset{ true };
	std::vector<ManifoldType*> m_Manifolds;
};

/**
 * \brief Input parameters for adaptive remeshing within manifold evolution strategy.
 * \struct ManifoldAdaptiveRemeshingParams
 */
struct ManifoldAdaptiveRemeshingParams
{
	pmp::Scalar MinEdgeMultiplier{ 0.14 }; //>! multiplier for minimum edge length in adaptive remeshing.
	unsigned int NRemeshingIters{ 3 }; //>! the number of iterations for pmp::Remeshing.
	unsigned int NTanSmoothingIters{ 5 }; //>! the number of tangential smoothing iterations for pmp::Remeshing.
	bool UseBackProjection{ true }; //>! if true surface kd-tree back-projection will be used for pmp::Remeshing.
};

/**
 * \brief A utility for estimating the edge sizing, error bound and other parameters for adaptive remeshing.
 * \param subdiv                 subdivision level of the initial icosphere (affects initial sizing).
 * \param radius                 radius of the initial icosphere (affects initial sizing).
 * \param remeshingParams        input parameters for manifold adaptive remeshing.
 * \return the result AdaptiveRemeshingSettings.
 */
[[nodiscard]] pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromIcoSphere_OLD(unsigned int subdiv, pmp::Scalar radius, const ManifoldAdaptiveRemeshingParams& remeshingParams);

/**
 * \brief A utility for estimating the edge sizing, error bound and other parameters for adaptive remeshing.
 * \param manifold               initial manifold (a pmp::ManifoldCurve2D or pmp::SurfaceMesh).
 * \param remeshingParams        input parameters for manifold adaptive remeshing.
 * \return the result AdaptiveRemeshingSettings.
 */
template <typename ManifoldType>
[[nodiscard]] pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsForManifold_EXPERIMENTAL(const std::shared_ptr<ManifoldType>& manifold, const ManifoldAdaptiveRemeshingParams& remeshingParams)
{
	if (!manifold)
	{
		throw std::invalid_argument("CollectRemeshingSettingsForManifold_EXPERIMENTAL: manifold == nullptr!\n");
	}

	pmp::AdaptiveRemeshingSettings settings;

	auto minEdgeLength = std::numeric_limits<pmp::Scalar>::max();
	for (const auto e : manifold->edges())
	{
		auto edgeLength = manifold->edge_length(e);
		minEdgeLength = std::min(minEdgeLength, edgeLength);
	}

	settings.MinEdgeLength = remeshingParams.MinEdgeMultiplier * minEdgeLength / 2.0;
	settings.MaxEdgeLength = 4.0 * settings.MinEdgeLength;
	settings.ApproxError = 0.25 * (settings.MinEdgeLength + settings.MaxEdgeLength);

	settings.NRemeshingIterations = remeshingParams.NRemeshingIters;
	settings.NTangentialSmoothingIters = remeshingParams.NTanSmoothingIters;
	settings.UseProjection = remeshingParams.UseBackProjection;

	return settings;
}

/**
 * \brief A utility for computing the edge sizing and error limits from an icosphere mesh.
 * \param icosphere      the input icosphere surface.
 * \param radius         the radius of the smooth sphere which the icosphere's supposed to approximate.
 * \param center         the center of the smooth sphere which the icosphere's supposed to approximate.
 * \return the result AdaptiveRemeshingSettings.
 */
[[nodiscard]] pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromIcoSphere(
	const std::shared_ptr<pmp::SurfaceMesh>& icosphere, pmp::Scalar radius, const pmp::Point& center);

/// \brief  A utility for computing the edge sizing and error limits from an arbitrary mesh.
[[nodiscard]] pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromMesh(const std::shared_ptr<pmp::SurfaceMesh>& mesh);

/**
 * \brief A utility for computing the edge sizing and error limits from an icosphere mesh.
 * \param circlePolyline    the input circle polyline curve.
 * \param radius            the radius of the smooth circle which the circlePolyline's supposed to approximate.
 * \param center            the center of the smooth circle which the circlePolyline's supposed to approximate.
 * \return the result AdaptiveRemeshingSettings.
 */
[[nodiscard]] pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromCircleCurve(
	const std::shared_ptr<pmp::ManifoldCurve2D>& circlePolyline, pmp::Scalar radius, const pmp::Point2& center);

/// \brief  A utility for computing the edge sizing and error limits from an arbitrary manifold curve.
[[nodiscard]] pmp::AdaptiveRemeshingSettings CollectRemeshingSettingsFromCurve(const std::shared_ptr<pmp::ManifoldCurve2D>& curve);


/**
 * \brief Logs settings for manifold remeshing.
 * \class ManifoldRemeshingSettingsWrapper
 * \tparam ManifoldType  either pmp::ManifoldCurve2D or pmp::SurfaceMesh.
 */
template<typename ManifoldType>
class ManifoldRemeshingSettingsWrapper
{
public:
	/// \brief Overload the [] operator for assigning and retrieving map values
	pmp::AdaptiveRemeshingSettings& operator[](ManifoldType* manifold)
	{
		return m_ManifoldSettings[manifold];
	}

	/// \brief Overload the [] operator for constant access (read-only)
	const pmp::AdaptiveRemeshingSettings& operator[](const ManifoldType* manifold) const
	{
		return m_ManifoldSettings.at(manifold);
	}

	/// \brief Adjust all remeshing settings using a provided resize factor
	void AdjustAllRemeshingLengths(pmp::Scalar resizeFactor)
	{
		for (auto& remeshingSettings : m_ManifoldSettings | std::views::values)
		{
			AdjustRemeshingLengths(resizeFactor,
				remeshingSettings.MinEdgeLength,
				remeshingSettings.MaxEdgeLength,
				remeshingSettings.ApproxError);
		}
	}

private:
	std::map<ManifoldType*, pmp::AdaptiveRemeshingSettings> m_ManifoldSettings{};
};

/**
 * \brief Logs initial sphere manifold settings.
 * \class InitialSphereSettingsWrapper
 * \tparam ManifoldType  either pmp::ManifoldCurve2D or pmp::SurfaceMesh.
 * \tparam SphereType    either Sphere2D or Sphere3D.
 */
template<typename ManifoldType, typename SphereType>
class InitialSphereSettingsWrapper
{
public:
	/// \brief Overload the [] operator for assigning and retrieving map values
	SphereType& operator[](ManifoldType* manifold)
	{
		return m_ManifoldSettings[manifold];
	}

	/// \brief Overload the [] operator for constant access (read-only)
	SphereType& operator[](const ManifoldType* manifold) const
	{
		return m_ManifoldSettings.at(manifold);
	}

	/// \brief Find the maximum radius among all stored SphereType objects
	[[nodiscard]] pmp::Scalar MaxRadius() const
	{
		if (m_ManifoldSettings.empty())
		{
			return -1.0; // Return a default value if the map is empty
		}

		// Use std::max_element to find the SphereType with the maximum radius
		auto maxIt = std::max_element(m_ManifoldSettings.begin(), m_ManifoldSettings.end(),
			[](const auto& lhs, const auto& rhs) {
				return lhs.second.Radius < rhs.second.Radius;
			});

		return maxIt->second.Radius;
	}

	/// \brief Find the minimum radius among all stored SphereType objects
	[[nodiscard]] pmp::Scalar MinRadius() const
	{
		if (m_ManifoldSettings.empty())
		{
			return -1.0; // Return a default value if the map is empty
		}

		// Use std::max_element to find the SphereType with the maximum radius
		auto maxIt = std::min_element(m_ManifoldSettings.begin(), m_ManifoldSettings.end(),
			[](const auto& lhs, const auto& rhs) {
				return lhs.second.Radius < rhs.second.Radius;
			});

		return maxIt->second.Radius;
	}

private:
	std::map<ManifoldType*, SphereType> m_ManifoldSettings{};
};

/**
 * \brief An enumerator for the type of selection used when the difference in distances is below critical interaction distance.
 * \enum FeatureDetectionType
 */
enum class [[nodiscard]] DistanceSelectionType
{
	PlainMinimum = 0, //>! plain minimum, produces a C^1 discontinuity in the interaction distance and a C^0 jummp in gradient.
	QuadricBlend = 1, //>! within a given critical radius, produces a C^2 discontinuity in the interaction distance and a C^1 jump in gradient.
};

/**
 * \brief Gathers interaction distance info between interacting evolving manifolds and target point cloud.
 * \class InteractionDistanceInfo
 * \tparam VectorType   either pmp::dvec2 or pmp::dvec3.
 */
template<typename VectorType>
class InteractionDistanceInfo
{
public:
	double Distance{ DBL_MAX };   //>! interaction distance
	VectorType NegGradient{}; //>! interaction vector

	/// \brief Stream operator for updating the interaction info
	InteractionDistanceInfo& operator<<(const InteractionDistanceInfo& newInfo)
	{
		if (newInfo.Distance < this->Distance)
		{
			// update both distance and gradient
			this->Distance = newInfo.Distance;
			this->NegGradient = newInfo.NegGradient;
		}
		return *this;
	}
};

///**
// * \brief Gathers vertex values from multiple manifolds and dumps them into a specified log file.
// * \class VertexValueLogger
// * \tparam ManifoldType   The type of the manifold (pmp::ManifoldCurve2D or pmp::SurfaceMesh).
// */
template<typename ManifoldType>
class VertexValueLogger
{
public:
	/// \brief Initialize the logger with a file name
	void Init(const std::string& fileName)
	{
		m_FileName = fileName;
	}

	/// \brief Reset the values for a new time step
	void StartNewTimeStep(unsigned int timeStep)
	{
		//if (!EXPERIMENTAL_CheckValues())
		//{
		//	std::cout << "VertexValueLogger::StartNewTimeStep: some values are incorrect!\n";
		//}

		m_CurrentTimeStep = timeStep;
		m_TimeStepData[timeStep] = nlohmann::json::object();
	}

	/// \brief Add a new manifold to the logger
	void AddManifold(const ManifoldType* manifold)
	{
		m_ManifoldKeys[manifold] = "Manifold_" + std::to_string(m_ManifoldCounter++);
	}

	///// \brief Reserve buffers for all value types for a manifold
	//void ReserveBuffers(const ManifoldType* manifold)
	//{
	//	size_t nVertices = manifold->n_vertices();
	//	std::string manifoldKey = m_ManifoldKeys[manifold];

	//	for (auto& [timeStep, data] : m_TimeStepData)
	//	{
	//		auto& timeStepJson = data;
	//		if (!timeStepJson.contains(manifoldKey))
	//			continue;

	//		auto& manifoldJson = timeStepJson[manifoldKey];
	//		for (auto& [valueId, vertexValues] : manifoldJson.items())
	//		{
	//			// Resize JSON arrays to match the updated vertex count
	//			while (vertexValues.size() < nVertices)
	//			{
	//				vertexValues[std::to_string(vertexValues.size())] = 0.0; // Default value
	//			}
	//			while (vertexValues.size() > nVertices)
	//			{
	//				vertexValues.erase(std::to_string(vertexValues.size() - 1)); // Remove excess
	//			}
	//		}
	//	}
	//}

	//double DBG_EPSILON{ 1e-6 };

	/// \brief Log a value for a specific manifold and vertex
	void LogValue(const ManifoldType* manifold, const std::string& valueId, unsigned int vertexIndex, double value)
	{
		if (!m_ManifoldKeys.contains(manifold))
			return;

		auto& timeStepJson = m_TimeStepData[m_CurrentTimeStep];
		std::string manifoldKey = m_ManifoldKeys[manifold];

		// Initialize structure for manifold if it doesn't exist
		if (!timeStepJson.contains(manifoldKey))
			timeStepJson[manifoldKey] = nlohmann::json::object();

		auto& manifoldJson = timeStepJson[manifoldKey];

		// Initialize structure for value type if it doesn't exist
		if (!manifoldJson.contains(valueId))
			manifoldJson[valueId] = nlohmann::json::object();

		//if (std::abs(value) < DBG_EPSILON)
		//{
		//	std::cout << "|" << value << "| < " << DBG_EPSILON << "\n";
		//}

		// Add value for vertex
		manifoldJson[valueId][std::to_string(vertexIndex)] = value;
	}

	/// \brief Save all logged data to JSON format
	// \param[in] omitLastTimeStep        a flag useful when an exception has been thrown during some time step.
	void Save(bool omitLastTimeStep = false)
	{
		nlohmann::json outputJson;
		const auto nTimeSteps = m_TimeStepData.size() - (omitLastTimeStep ? 1 : 0);
		for (size_t timeStep = 1; timeStep <= nTimeSteps; ++timeStep)
		{
			auto& timeStepJson = outputJson["TimeStep_" + std::to_string(timeStep)];

			for (const auto& [manifoldKey, manifoldData] : m_TimeStepData[timeStep].items())
			{
				auto& manifoldJson = timeStepJson[manifoldKey];

				for (const auto& [valueId, vertexValues] : manifoldData.items())
				{
					// Extract vertex indices and values, then sort them
					std::vector<std::pair<unsigned int, double>> sortedVertexValues;
					for (const auto& [vertexIndex, value] : vertexValues.items())
					{
						sortedVertexValues.emplace_back(std::stoi(vertexIndex), value.template get<double>());
					}
					std::ranges::sort(sortedVertexValues, [](const auto& a, const auto& b) { return a.first < b.first; });

					// Store sorted values in an array
					nlohmann::json sortedArray = nlohmann::json::array();
					for (const auto& [vertexIndex, value] : sortedVertexValues)
					{
						//if (std::abs(value) < DBG_EPSILON)
						//{
						//	std::cout << "|" << value << "| < " << DBG_EPSILON << "\n";
						//}

						sortedArray.push_back({ {"vertexIndex", vertexIndex}, {"value", value} });
					}

					manifoldJson[valueId] = sortedArray;
				}
			}
		}

		// Write the JSON to the output file
		std::ofstream file(m_FileName);
		if (!file.is_open())
		{
			throw std::runtime_error("VertexValueLogger::Save: Failed to open file: " + m_FileName);
		}

		file << outputJson.dump(4); // Pretty-print with 4 spaces
		file.close();
	}


private:

	/// Experimental: Check if any logged values in the previous time step are below a given epsilon
	//bool EXPERIMENTAL_CheckValues(double epsilon = 1e-5) const
	//{
	//	if (m_CurrentTimeStep == 0)
	//	{
	//		return true; // Nothing to check for the first time step
	//	}

	//	auto prevTimeStep = m_TimeStepData.find(m_CurrentTimeStep - 1);
	//	if (prevTimeStep == m_TimeStepData.end())
	//	{
	//		return true; // No data for the previous time step
	//	}

	//	const auto& timeStepJson = prevTimeStep->second;

	//	for (const auto& [manifoldKey, manifoldJson] : timeStepJson.items())
	//	{
	//		for (const auto& [valueId, valueJson] : manifoldJson.items())
	//		{
	//			for (const auto& [vertexIndex, value] : valueJson.items())
	//			{
	//				if (value.is_number())
	//				{
	//					double val = value.get<double>();
	//					if (std::abs(val) < epsilon)
	//					{
	//						std::cout << "CheckValues: Found value below epsilon (" << epsilon
	//							<< ") in TimeStep: " << m_CurrentTimeStep - 1
	//							<< ", Manifold: " << manifoldKey
	//							<< ", ValueId: " << valueId
	//							<< ", VertexIndex: " << vertexIndex
	//							<< ", Value: " << val << "\n";
	//						return false;
	//					}
	//				}
	//			}
	//		}
	//	}
	//	return true;
	//}


	//[[nodiscard]] std::vector<unsigned int> GetSortedVertexIds(const std::string& manifoldKey, )


	std::string m_FileName; //!< Log file name
	unsigned int m_CurrentTimeStep = 0; //!< Current time step being logged
	unsigned int m_ManifoldCounter = 0; //!< Counter for assigning manifold keys
	std::unordered_map<const ManifoldType*, std::string> m_ManifoldKeys; //!< Map from manifold pointer to unique key
	std::unordered_map<unsigned int, nlohmann::json> m_TimeStepData; //!< Map from time step to JSON data
	std::string m_SortByValueId{}; //>! if empty the vertex values will be pre-sorted by vertex index, otherwise by a given value id.
};

//
// ==========================================================================
//

/// \brief Extended control function wrapper
template<typename ControlFunctionSignature>
class ExtendedFunction;

/// \brief Extended control function wrapper definition
template<typename... Args>
class ExtendedFunction<std::function<double(double /*dist*/, Args...)>>
{
public:

	ExtendedFunction() = default;

	explicit ExtendedFunction(std::function<double(double /*dist*/, Args...)> func)
	{
		Bind(std::move(func));
	}

	void Bind(double percentage, std::function<double(double /*dist*/, Args...)> func)
	{
		m_InputPercentage = percentage;
		Bind(std::move(func));
	}

	void ComputeLimit(const double& maxDistance = 0.0)
	{
		if (maxDistance < 0.0 || m_InputPercentage < 0.0)
		{
			// invalid input
			m_LowerDistanceLimit = 0.0;
			return;
		}

		m_LowerDistanceLimit = maxDistance * m_InputPercentage;
	}

	double operator()(double dist, Args... args) const
	{
		if (!m_Lambda)
		{
			// invocation not available
			return 0.0;
		}

		if (dist <= m_LowerDistanceLimit)
			return 0.0;

		return m_Lambda(dist, std::forward<Args>(args)...);
	}

	ExtendedFunction& operator=(std::function<double(double /*dist*/, Args...)> func)
	{
		Bind(std::move(func));
		return *this;
	}

private:

	void Bind(std::function<double(double /*dist*/, Args...)> func)
	{
		m_Lambda = std::move(func);
	}

	double m_InputPercentage{ 0.0 };
	std::function<double(double /*dist*/, Args...)> m_Lambda;
	double m_LowerDistanceLimit{ 0.0 };
};