#pragma once

#include "geometry/Grid.h"

namespace SDF
{
	/// \brief A set of settings for the FastSweep function.
	struct SweepSolverSettings
	{
		unsigned int NSweeps{ 8 }; //>! the number of sweeps for the solver.
		double EikonalRHS{ 1.0 }; //>! the right-hand side f of the Eikonal equation ||grad(u)|| = f.
	};

	/**
	 * \brief A utility for solving the Eikonal equation using the Fast-Sweeping algorithm [Zhao, 2005].
	 * \param grid        a modifiable input grid with proper initial condition setup.
	 * \param settings    settings for this solver function.
	 */
	void FastSweep(Geometry::ScalarGrid& grid, const SweepSolverSettings& settings);

} // namespace SDF