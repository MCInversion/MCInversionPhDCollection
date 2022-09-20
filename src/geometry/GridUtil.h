# pragma once

#include "Grid.h"

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

} // namespace Geometry