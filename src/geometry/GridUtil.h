# pragma once

#include "pmp/Types.h"

#include "Grid.h"

#include <functional>
#include <optional>
#include <random>

namespace pmp
{
	class SurfaceMesh;
}

namespace Geometry
{
	/**
	 * \brief Negate the values of a scalar grid.
	 * \param grid    2D/3D input scalar grid.
	 */
	template <typename Grid>
	void NegateGrid(Grid& grid)
	{
		if (!grid.IsValid())
			return;
		grid *= -1.0;
	}

	/// \brief Compute a hash from the provided 2D/3D scalar grid
	template <typename Grid>
	[[nodiscard]] size_t HashScalarGrid(const Grid& grid)
	{
		std::size_t hash = 0;
		for (const auto& value : grid.Values())
		{
			hash ^= std::hash<double>{}(value) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
		}
		return hash;
	}

	/**
	 * \brief Negates the relevant sub-volume of the scalar grid.
	 * \param grid        input scalar grid.
	 * \param subBox      box which should be the sub-box of the grid's box.
	 */
	void NegateGridSubVolume(ScalarGrid& grid, const pmp::BoundingBox& subBox);

	/**
	 * \brief Apply a 3x3x3 averaging kernel onto a given grid
	 * \param grid    input grid.
	 */
	void ApplyNarrowAveragingBlur(ScalarGrid& grid);

	/**
	 * \brief Apply a 5x5x5 averaging kernel onto a given grid
	 * \param grid    input grid.
	 */
	void ApplyWideAveragingBlur(ScalarGrid& grid);

	/**
	 * \brief Apply a 3x3x3 Gaussian kernel onto a given grid
	 * \param grid    input grid.
	 */
	void ApplyNarrowGaussianBlur(ScalarGrid& grid);

	/**
	 * \brief Apply a 5x5x5 Gaussian kernel onto a given grid
	 * \param grid    input grid.
	 */
	void ApplyWideGaussianBlur(ScalarGrid& grid);

	/**
	 * \brief Apply a 3x3 averaging kernel onto a given grid
	 * \param grid    input grid.
	 */
	void ApplyNarrowAveragingBlur2D(ScalarGrid2D& grid);

	/**
	 * \brief Apply a 5x5 averaging kernel onto a given grid
	 * \param grid    input grid.
	 */
	void ApplyWideAveragingBlur2D(ScalarGrid2D& grid);

	/**
	 * \brief Apply a 3x3 Gaussian kernel onto a given grid
	 * \param grid    input grid.
	 */
	void ApplyNarrowGaussianBlur2D(ScalarGrid2D& grid);

	/**
	 * \brief Apply a 5x5 Gaussian kernel onto a given grid
	 * \param grid    input grid.
	 */
	void ApplyWideGaussianBlur2D(ScalarGrid2D& grid);

	/**
	 * \brief Looks for nans and infinities in the grid, if the cell neighbors have valid values, averaged value is written for an invalid cell. Otherwise the cell value is set to default init value.
	 * \param grid    input grid.
	 */
	void RepairScalarGrid(ScalarGrid& grid, bool verbose = false);

	/**
	 * \brief Looks for nans and infinities in the grid, if the cell neighbors have valid values, averaged value is written for an invalid cell. Otherwise the cell value is set to default init value.
	 * \param grid    input grid.
	 */
	void RepairScalarGrid2D(ScalarGrid2D& grid, bool verbose = false);

	/**
	 * \brief Normalizes the values of the scalar grid to interval [-1, 1]
	 * \param grid    input grid.
	 */
	void NormalizeScalarGridValues(ScalarGrid& grid);

	/**
	 * \brief Computes a gradient from a given scalar grid.
	 * \param scalarGrid     input grid.
	 * \return gradient field.
	 *
	 * DISCLAIMER: This function uses central difference for approximating partial derivatives of scalarGrid values. Boundary voxels are skipped and contain default vector values.
	 */
	[[nodiscard]] VectorGrid ComputeGradient(const ScalarGrid& scalarGrid);

	/**
	 * \brief Computes a gradient from a given 2D scalar grid.
	 * \param scalarGrid     input grid.
	 * \return gradient field.
	 *
	 * DISCLAIMER: This function uses central difference for approximating partial derivatives of scalarGrid values. Boundary pixels are skipped and contain default vector values.
	 */
	[[nodiscard]] VectorGrid2D ComputeGradient(const ScalarGrid2D& scalarGrid);

	/**
	 * \brief Computes a normalized gradient from a given scalar grid.
	 * \param scalarGrid     input grid.
	 * \return normalized gradient field.
	 *
	 * DISCLAIMER: This function uses central difference for approximating partial derivatives of scalarGrid values. Boundary pixels are skipped and contain default vector values.
	 */
	[[nodiscard]] VectorGrid2D ComputeNormalizedGradient(const ScalarGrid2D& scalarGrid);

	/**
	 * \brief Computes a normalized gradient from a given scalar grid.
	 * \param scalarGrid     input grid.
	 * \return normalized gradient field.
	 *
	 * DISCLAIMER: This function uses central difference for approximating partial derivatives of scalarGrid values. Boundary voxels are skipped and contain default vector values.
	 */
	[[nodiscard]] VectorGrid ComputeNormalizedGradient(const ScalarGrid& scalarGrid);

	/**
	 * \brief Computes a normalized negative gradient from a given scalar grid.
	 * \param scalarGrid     input grid.
	 * \return normalized negative gradient field.
	 *
	 * DISCLAIMER: This function uses central difference for approximating partial derivatives of scalarGrid values. Boundary voxels are skipped and contain default vector values.
	 */
	[[nodiscard]] VectorGrid ComputeNormalizedNegativeGradient(const ScalarGrid& scalarGrid);

	/**
	 * \brief Computes a normalized negative gradient from a given scalar grid.
	 * \param scalarGrid     input grid.
	 * \return normalized negative gradient field.
	 *
	 * DISCLAIMER: This function uses central difference for approximating partial derivatives of scalarGrid values. Boundary voxels are skipped and contain default vector values.
	 */
	[[nodiscard]] VectorGrid2D ComputeNormalizedNegativeGradient(const ScalarGrid2D& scalarGrid);

	/**
	 * \brief Computes a divergence field from the given vector field.
	 * \param vectorGrid     input grid.
	 * \return divergence field.
	 *
	 * DISCLAIMER: This function uses central difference for approximating partial derivatives of scalarGrid values. Boundary voxels are skipped and contain default vector values.
	 */
	[[nodiscard]] ScalarGrid2D ComputeDivergenceField(const VectorGrid2D& vectorGrid);

	/**
	 * \brief A parameter container for the generation of streamline paths on a vector grid
	 * \struct StreamLineSettings
	 */
	struct StreamLineSettings
	{
		unsigned int NSamplePts{ 10 }; //>! determines the sampling density for Poisson sampling.

		unsigned int NSamplingAttempts{ 30 }; //>! The number of attempts for Poisson sampling.

		unsigned int NSteps{ 20 }; //>! The number of steps for the Runge-Kutta integration
		double StepSize{ 0.1 }; //>! The step size for the Runge-Kutta integration.
		bool UseAdaptiveStepSize{ true }; //>! Whether to use adaptive step size for the Runge-Kutta integration.
		double Tolerance{ 1e-6 }; //>! The tolerance for the Runge-Kutta integration.
		unsigned int MaxIterations{ 10'000 }; //>! The maximum number of iterations for the Runge-Kutta method.

		unsigned int Seed{ std::random_device{}() }; // Default seed is random
	};

	/**
	* \brief Computes a stream line from a given vector field.
	* \param vectorGrid      input grid.
	* \param settings        settings for streamlines
	*
	* \return a list of stream lines.
	*/
	[[nodiscard]] std::vector<std::vector<pmp::Point2>> CalculateStreamLines(const VectorGrid2D& vectorGrid, const StreamLineSettings& settings);

	/**
	* \brief Computes a list of stream lines from a given vector field (using RK4 Method).
	* \param vectorGrid     input grid.
	* \param seedPoints     starting points of the stream lines.
	* \param numSteps       number of steps to take.
	* \param stepSize       size of each step.
	* \param useAdaptiveStepSize    whether to use adaptive step size.
	* \param tolerance      tolerance for adaptive step size.
	* \param maxIterations  maximum number of iterations for adaptive step size.
	* 
	* \return a list of stream lines.
	*/
	[[nodiscard]] std::vector<std::vector<pmp::Point2>> CalculateStreamLines(
		const VectorGrid2D& vectorGrid,
		const std::vector<pmp::Point2>& seedPoints,
		unsigned int numSteps,
		double stepSize,
		bool useAdaptiveStepSize = false,
		double tolerance = 1e-5,
		unsigned int maxIterations = 1000);

	/**
	 * \brief A parameter container for the generation of Euler stream line paths on a vector grid
	 * \struct EulerStreamLineSettings
	 */
	struct EulerStreamLineSettings
	{
		unsigned int NSamplePts{ 10 };         // Number of seed points for Poisson disk sampling
		unsigned int NSamplingAttempts{ 30 };  // Number of attempts for Poisson sampling
		unsigned int NSteps{ 100 };            // Number of steps for path integration
		double StepSize{ 0.1 };                // Step size for path integration
		unsigned int MaxIterations{ 10'000 };    // Maximum number of iterations per path
		unsigned int Seed{ std::random_device{}() }; // Seed for reproducibility
	};

	/**
	* \brief Computes a list of stream lines from a given vector field (using Euler method).
	* \param vectorGrid      input grid.
	* \param settings        settings for stream lines
	*
	* \return a list of stream lines.
	*/
	[[nodiscard]] std::vector<std::vector<pmp::Point2>> CalculateStreamLinesEuler(const VectorGrid2D& vectorGrid, const EulerStreamLineSettings& settings);

	/**
	* \brief Computes a list of stream lines from a given vector field (using Euler method).
	* \param vectorGrid     input grid.
	* \param seedPoints     starting points of the stream lines.
	* \param numSteps       number of steps to take.
	* \param stepSize       size of each step.
	* \param maxIterations  maximum number of iterations for adaptive step size.
	*
	* \return a list of stream lines.
	*/
	[[nodiscard]] std::vector<std::vector<pmp::Point2>> CalculateStreamLinesEuler(
		const VectorGrid2D& vectorGrid,
		const std::vector<pmp::Point2>& seedPoints,
		unsigned int numSteps,
		double stepSize,
		unsigned int maxIterations = 1000);

	// ..........................................................................................

	/// \brief A wrapper for the 8 direction vectors and their corresponding grid points
	struct FieldDirectionsOctet
	{
		std::array<pmp::vec2, 8> GridPts{};
		std::array<pmp::dvec2, 8> Directions{};
	};

	/**
	* \brief Finds the closest grid point near samplePt and its 8 neighboring points where the normalized field directions are extracted.
	* \param samplePt        point whose 8 neighboring cells are going to be evaluated
	* \param vectorGrid      input grid.
	*
	* \return the evaluated FieldDirectionsOctet
	*/
	[[nodiscard]] FieldDirectionsOctet GetFieldDirectionsAroundPoint(const pmp::Point2& samplePt, const VectorGrid2D& vectorGrid);

	// ==========================================================================================
	//             Interpolation
	// ------------------------------------------------------------------------------------------

	/**
	 * \brief Retrieves the value of the nearest neighbor from the surrounding cell values of a sampled point.
	 * \param samplePt    point where the grid is sampled.
	 * \param grid        interpolated scalar grid.
	 * \return interpolated value.
	 *
	 * DISCLAIMER: For samplePt outside of grid.Box(), the values are clamped to boundary values.
	 */
	[[nodiscard]] double GetNearestNeighborScalarValue(const pmp::vec3& samplePt, const ScalarGrid& grid);

	/**
	 * \brief Retrieves the value of the nearest neighbor from the surrounding cell vector values of a sampled point.
	 * \param samplePt    point where the grid is sampled.
	 * \param grid        interpolated vector grid.
	 * \return interpolated vector value.
	 *
	 * DISCLAIMER: For samplePt outside of grid.Box(), the values are clamped to boundary values.
	 */
	[[nodiscard]] pmp::dvec3 GetNearestNeighborVectorValue(const pmp::vec3& samplePt, const VectorGrid& grid);

	/**
	 * \brief Retrieves the value of the nearest neighbor from the surrounding cell values of a sampled point.
	 * \param samplePt    point where the grid is sampled.
	 * \param grid        interpolated scalar grid.
	 * \return interpolated value.
	 *
	 * DISCLAIMER: For samplePt outside of grid.Box(), the values are clamped to boundary values.
	 */
	[[nodiscard]] double GetNearestNeighborScalarValue2D(const pmp::vec2& samplePt, const ScalarGrid2D& grid);

	/**
	 * \brief Retrieves the value of the nearest neighbor from the surrounding cell vector values of a sampled point.
	 * \param samplePt    point where the grid is sampled.
	 * \param grid        interpolated vector grid.
	 * \return interpolated vector value.
	 *
	 * DISCLAIMER: For samplePt outside of grid.Box(), the values are clamped to boundary values.
	 */
	[[nodiscard]] pmp::dvec2 GetNearestNeighborVectorValue2D(const pmp::vec2& samplePt, const VectorGrid2D& grid);

	/**
	 * \brief Trilinearly interpolates from the surrounding cell values of a sampled point.
	 * \param samplePt    point where the grid is sampled.
	 * \param grid        interpolated scalar grid.
	 * \return interpolated value.
	 *
	 * DISCLAIMER: For samplePt outside of grid.Box(), the values are clamped to boundary values.
	 */
	[[nodiscard]] double TrilinearInterpolateScalarValue(const pmp::vec3& samplePt, const ScalarGrid& grid);

	/**
	 * \brief Trilinearly interpolates from the surrounding cell vector values of a sampled point.
	 * \param samplePt    point where the grid is sampled.
	 * \param grid        interpolated vector grid.
	 * \return interpolated vector value.
	 *
	 * DISCLAIMER: For samplePt outside of grid.Box(), the values are clamped to boundary values.
	 */
	[[nodiscard]] pmp::dvec3 TrilinearInterpolateVectorValue(const pmp::vec3& samplePt, const VectorGrid& grid);

	/**
	 * \brief Bilinearly interpolates from the surrounding cell values of a sampled point.
	 * \param samplePt    point where the grid is sampled.
	 * \param grid        interpolated scalar grid.
	 * \return interpolated value.
	 *
	 * DISCLAIMER: For samplePt outside of grid.Box(), the values are clamped to boundary values.
	 */
	[[nodiscard]] double BilinearInterpolateScalarValue(const pmp::vec2& samplePt, const ScalarGrid2D& grid);

	/**
	 * \brief Bilinearly interpolates from the surrounding cell vector values of a sampled point.
	 * \param samplePt    point where the grid is sampled.
	 * \param grid        interpolated vector grid.
	 * \return interpolated vector value.
	 *
	 * DISCLAIMER: For samplePt outside of grid.Box(), the values are clamped to boundary values.
	 */
	[[nodiscard]] pmp::dvec2 BilinearInterpolateVectorValue(const pmp::vec2& samplePt, const VectorGrid2D& grid);

	// ------------------------------------------------------------------------------------------------

	/**
	 * \brief Computes sign values (-1 or 1) for each grid point using mesh normals.
	 * \param grid       ScalarGrid where values are stored.
	 * \param mesh       input mesh (needs to have normals).
	 */
	void ComputeInteriorExteriorSignFromMeshNormals(ScalarGrid& grid, const pmp::SurfaceMesh& mesh);

	/**
	 * \brief Computes distance for each grid point using mesh normals [Baerentzen et al. 2005].
	 * \param grid       ScalarGrid where values are stored.
	 * \param mesh       input mesh (needs to have normals).
	 * NOTE: Not using DistanceFieldGenerator because its pipeline is unnecessary for this approach.
	 */
	void ComputeMeshSignedDistanceFromNormals(ScalarGrid& grid, const pmp::SurfaceMesh& mesh);

	/// \brief a functor for evaluating a new grid value according to a boolean operation.
	using ScalarGridBoolOpFunction = std::function<double(const double&, const double&)>;

	//
	// ======= Bool Op Function Gallery =============
	// see Section 5.1.1 in [Sanchez, M., Continuous Signed Distance Field Representation of Polygonal Meshes, 2011]
	// https://nccastaff.bournemouth.ac.uk/jmacey/MastersProject/MSc11/Mathieu/msanchez-sdf-thesis.pdf

	/**
	 * \brief A union operator between functional representations: f1 ^ f2 = max(f1, f2).
	 * \param f1Val     value of function 1 to be used for evaluation.
	 * \param f2Val     value of function 2 to be used for evaluation.
	 * \return new value.
	 */
	[[nodiscard]] double SimpleUnion(const double& f1Val, const double& f2Val);

	/**
	 * \brief An intersection operator between functional representations: f1 v f2 = min(f1, f2).
	 * \param f1Val     value of function 1 to be used for evaluation.
	 * \param f2Val     value of function 2 to be used for evaluation.
	 * \return new value.
	 */
	[[nodiscard]] double SimpleIntersection(const double& f1Val, const double& f2Val);

	/**
	 * \brief A difference operator between functional representations: f1 \ f2 = max(f1 - f2, 0.0).
	 * \param f1Val     value of function 1 to be used for evaluation.
	 * \param f2Val     value of function 2 to be used for evaluation.
	 * \return new value.
	 */
	[[nodiscard]] double SimpleDifference(const double& f1Val, const double& f2Val);

	/**
	 * \brief A blended union operator between functional representations: f1 ^b f2 = f1 + f2 + sqrt(f1^2 + f2^2).
	 * \param f1Val     value of function 1 to be used for evaluation.
	 * \param f2Val     value of function 2 to be used for evaluation.
	 * \return new value.
	 */
	[[nodiscard]] double BlendedUnion(const double& f1Val, const double& f2Val);

	/**
	 * \brief A blended intersection operator between functional representations: f1 vb f2 = f1 + f2 - sqrt(f1^2 + f2^2).
	 * \param f1Val     value of function 1 to be used for evaluation.
	 * \param f2Val     value of function 2 to be used for evaluation.
	 * \return new value.
	 */
	[[nodiscard]] double BlendedIntersection(const double& f1Val, const double& f2Val);

	/**
	 * \brief A blended difference operator between functional representations: f1 \b f2 = f1 - f2 - sqrt(f1^2 + f2^2).
	 * \param f1Val     value of function 1 to be used for evaluation.
	 * \param f2Val     value of function 2 to be used for evaluation.
	 * \return new value.
	 */
	[[nodiscard]] double BlendedDifference(const double& f1Val, const double& f2Val);

	/**
	 * \brief A union operator between distance functional representations: f1 v f2 = min(f1, f2).
	 * \param f1Val     value of function 1 to be used for evaluation.
	 * \param f2Val     value of function 2 to be used for evaluation.
	 * \return new value.
	 */
	[[nodiscard]] double DistanceUnion(const double& f1Val, const double& f2Val);

	/**
	 * \brief A parameter container for the metaball generator for ScalarGrid.
	 * \struct MetaBallParams
	 */
	struct MetaBallParams
	{
		pmp::vec3 Center{}; //! center of the metaball object 
		pmp::Scalar Radius{ 1.0 }; //! effective radius of the metaball object's zero isosurface
		ScalarGridBoolOpFunction BoolOpFunction{ SimpleUnion }; //! a boolean functor for blending values of a metaball object.
	};

	/**
	 * \brief Applies a metaball object with given MetaBallParams to a ScalarGrid.
	 * \param grid      grid, where metaball values are to be applied.
	 * \param params    parameters of the applied metaball object.
	 */
	void ApplyMetaBallToGrid(ScalarGrid& grid, const MetaBallParams& params);

	/**
	 * \brief A parameter container for the capsule generator for ScalarGrid.
	 * \struct CapsuleParams
	 */
	struct CapsuleParams
	{
		pmp::vec3 Position{}; //! base position of the capsule object
		pmp::Scalar Height{ 2.0 }; //! height of the capsule object
		pmp::Scalar Radius{ 1.0 }; //! radius of the capsule object
		ScalarGridBoolOpFunction BoolOpFunction{ SimpleUnion }; //! a boolean functor for blending values of a capsule object.
	};

	/**
	 * \brief Applies a capsule object with given CapsuleParams to a ScalarGrid.
	 * \param grid      grid, where capsule values are to be applied.
	 * \param params    parameters of the applied capsule object.
	 */
	void ApplyCapsuleDistanceFieldToGrid(ScalarGrid& grid, const CapsuleParams& params);

	/**
	 * \brief A parameter container for the torus generator for ScalarGrid.
	 * \struct TorusParams
	 */
	struct TorusParams
	{
		pmp::vec3 Center{}; //! base position of the torus object
		pmp::Scalar RingRadius{ 1.0 }; //! radius of the torus loop
		pmp::Scalar TubeRadius{ 0.3 }; //! radius of the torus tube
		ScalarGridBoolOpFunction BoolOpFunction{ SimpleUnion }; //! a boolean functor for blending values of a torus object.
	};

	/**
	 * \brief Applies a torus object with given TorusParams to a ScalarGrid.
	 * \param grid      grid, where torus values are to be applied.
	 * \param params    parameters of the applied torus object.
	 */
	void ApplyTorusDistanceFieldToGrid(ScalarGrid& grid, const TorusParams& params);

	/**
	 * \brief A parameter container for the single-sheet quadric generator for ScalarGrid.
	 * \struct QuadricParams
	 */
	struct QuadricParams
	{
		pmp::vec3 Center{}; //! base position of the single-sheet hyperboloid object
		pmp::Scalar a{ 1.0 };
		pmp::Scalar b{ 1.0 };
		pmp::Scalar c{ 1.0 };
		pmp::vec3 HalfSize{ 0.5, 0.5, 0.5 }; //! the half-size of the infinite hyperboloid's ROI.
		ScalarGridBoolOpFunction BoolOpFunction{ SimpleUnion }; //! a boolean functor for blending values of a torus object.
	};

	/**
	 * \brief Applies a single-sheet hyperboloid object with given QuadricParams to a ScalarGrid.
	 * \param grid      grid, where hyperboloid values are to be applied.
	 * \param params    parameters of the applied torus object.
	 */
	void ApplyHyperboloidDistanceFieldToGrid(ScalarGrid& grid, const QuadricParams& params);

	/**
	 * \brief Applies a single-sheet ellipsoid object with given QuadricParams to a ScalarGrid.
	 * \param grid      grid, where ellipsoid values are to be applied.
	 * \param params    parameters of the applied torus object.
	 */
	void ApplyEllipsoidDistanceFieldToGrid(ScalarGrid& grid, const QuadricParams& params);

	// ==================================================================================================================

	/**
	 * \brief Re-samples grid data to a new grid with the same bounds and different cell size.
	 * \param newCellSize     cell size of the re-sampled grid.
	 * \param origGrid            original grid
	 * \return re-sampled grid.
	 * NOTE: Using trilinear interpolation because the new cells don't need to align with the original grid's cells.
	 */
	[[nodiscard]] ScalarGrid ExtractReSampledGrid(const pmp::Scalar& newCellSize, const ScalarGrid& origGrid);

	/**
	 * \brief Searches the cell (ix, iy) and its neighbors at given radius for a local maximum of a 2D quadratic polynomial.
	 * \return true if the maximum is found between the 9 cells.
	 */
	[[nodiscard]] bool ContainsLocalMaximumNearScalarGridCell(const ScalarGrid2D& grid, unsigned int ix, unsigned int iy, unsigned int radius = 1);

	/**
	 * \brief Searches the cell (ix, iy) and its neighbors at given radius for a local extreme of a 2D quadratic polynomial.
	 * \return true if the extreme is found between the 9 cells.
	 */
	[[nodiscard]] bool ContainsLocalExtremesNearScalarGridCell(const ScalarGrid2D& grid, unsigned int ix, unsigned int iy, unsigned int radius = 1);

	/**
	 * \brief Searches the cell (ix, iy) and its neighbors for a local maximum of a 2D quadratic polynomial.
	 * \return optional point of local maximum. std::nullopt if the maximum isn't found between the neighboring cells.
	 */
	[[nodiscard]] std::optional<pmp::Point2> FindLocalMaximumNearScalarGridCell(const ScalarGrid2D& grid, unsigned int ix, unsigned int iy, unsigned int radius = 1);

	/// \brief Verifies whether the vector field has a non-zero divergence at position (ix, iy) within a given radius.
	[[nodiscard]] bool IsConvergentOrDivergentNearCell(const VectorGrid2D& vecGrid, unsigned int ix, unsigned int iy, unsigned int radius = 1);

	/**
	 * \brief Searches the cell (ix, iy, iz) and its neighbors for a local maximum of a 3D quadratic polynomial.
	 * \return optional point of local maximum. std::nullopt if the maximum isn't found between the neighboring cells.
	 */
	[[nodiscard]] std::optional<pmp::Point> FindLocalMaximumNearScalarGridCell_EXPERIMENTAL(const ScalarGrid& grid, unsigned int ix, unsigned int iy, unsigned int iz, unsigned int radius = 1);

	/**
	 * \brief Searches the cell (ix, iy, iz) and its neighbors for a local maximum of a 2D quadratic polynomial for each z-slice, and then averaging them.
	 * \return optional point of local maximum. std::nullopt if the maximum isn't found between the neighboring cells.
	 */
	[[nodiscard]] std::optional<pmp::Point> FindLocalMaximumNearScalarGridCell(const ScalarGrid& grid, unsigned int ix, unsigned int iy, unsigned int iz, unsigned int radius = 1);

	/// \brief Verifies whether the vector field has a non-zero divergence at position (ix, iy, iz) within a given radius.
	[[nodiscard]] bool IsConvergentOrDivergentNearCell(const VectorGrid& vecGrid, unsigned int ix, unsigned int iy, unsigned int iz, unsigned int radius = 1);

	/// \brief Extracts a sub-grid from given indices
	[[nodiscard]] ScalarGrid2D ExtractSubGrid2D(const ScalarGrid2D& grid,
		unsigned int ix0, unsigned int iy0,
		unsigned int ix1, unsigned int iy1
	);

	/// \brief Extracts a sub-grid from given indices
	[[nodiscard]] ScalarGrid ExtractSubGrid(const ScalarGrid& grid,
		unsigned int ix0, unsigned int iy0, unsigned int iz0,
		unsigned int ix1, unsigned int iy1, unsigned int iz1
	);

} // namespace Geometry