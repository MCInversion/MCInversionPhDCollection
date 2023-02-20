# pragma once

#include "Grid.h"

namespace pmp
{
	class SurfaceMesh;
}

namespace Geometry
{
	/**
	 * \brief Negate the values of a scalar grid.
	 * \param grid    input scalar grid.
	 */
	inline void NegateGrid(ScalarGrid& grid)
	{
		if (!grid.IsValid())
			return;
		grid *= -1.0;
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
	 * \brief Looks for nans and infinities in the grid, if the cell neighbors have valid values, averaged value is written for an invalid cell. Otherwise the cell value is set to default init value.
	 * \param grid    input grid.
	 */
	void RepairScalarGrid(ScalarGrid& grid);

	/**
	 * \brief Computes a gradient from a given scalar grid.
	 * \param scalarGrid     input grid.
	 * \return gradient field.
	 *
	 * DISCLAIMER: This function uses central difference for approximating partial derivatives of scalarGrid values. Boundary voxels are skipped and contain default vector values.
	 */
	[[nodiscard]] VectorGrid ComputeGradient(const ScalarGrid& scalarGrid);

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
	 * \brief Computes sign values (-1 or 1) for each grid point using mesh normals.
	 * \param grid       ScalarGrid where values are stored.
	 * \param mesh       input mesh (needs to have normals).
	 */
	void ComputeInteriorExteriorSignFromMeshNormals(ScalarGrid& grid, const pmp::SurfaceMesh& mesh);

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

} // namespace Geometry