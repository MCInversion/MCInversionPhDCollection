# pragma once

#include "Grid.h"

#include <functional>
#include <optional>

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
	void RepairScalarGrid(ScalarGrid& grid);

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
		float Radius{ 1.0f }; //! effective radius of the metaball object's zero isosurface
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
		float Height{ 2.0f }; //! height of the capsule object
		float Radius{ 1.0f }; //! radius of the capsule object
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
		float RingRadius{ 1.0f }; //! radius of the torus loop
		float TubeRadius{ 0.3f }; //! radius of the torus tube
		ScalarGridBoolOpFunction BoolOpFunction{ SimpleUnion }; //! a boolean functor for blending values of a torus object.
	};

	/**
	 * \brief Applies a torus object with given TorusParams to a ScalarGrid.
	 * \param grid      grid, where torus values are to be applied.
	 * \param params    parameters of the applied torus object.
	 */
	void ApplyTorusDistanceFieldToGrid(ScalarGrid& grid, const TorusParams& params);

	// ==================================================================================================================

	/**
	 * \brief Re-samples grid data to a new grid with the same bounds and different cell size.
	 * \param newCellSize     cell size of the re-sampled grid.
	 * \param origGrid            original grid
	 * \return re-sampled grid.
	 * NOTE: Using trilinear interpolation because the new cells don't need to align with the original grid's cells.
	 */
	[[nodiscard]] ScalarGrid ExtractReSampledGrid(const float& newCellSize, const ScalarGrid& origGrid);

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