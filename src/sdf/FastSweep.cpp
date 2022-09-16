#include "FastSweep.h"

#include "pmp/BoundingBox.h"

namespace SDF
{
	void FastSweep(Geometry::ScalarGrid& grid, const SweepSolverSettings& settings)
	{
		const double f = settings.EikonalRHS;
		assert(grid.CellSize() > 0.0);
		const auto h = static_cast<double>(grid.CellSize());
		const unsigned int NSweeps = settings.NSweeps;

		const auto& dims = grid.Dimensions();
		auto& gridValues = grid.Values();
		const auto& gridFrozen = grid.FrozenValues();

		const int Nx = static_cast<int>(dims.Nx);
		const int Ny = static_cast<int>(dims.Ny);
		const int Nz = static_cast<int>(dims.Nz);

		// sweep directions { start, end, step }
		const int dirX[8][3] = {
			{ 0, Nx - 1, 1 }, { Nx - 1, 0, -1 }, { Nx - 1, 0, -1 }, { Nx - 1, 0, -1 },
			{ Nx - 1, 0, -1 }, { 0, Nx - 1, 1 }, { 0, Nx - 1, 1 }, { 0, Nx - 1, 1 } };
		const int dirY[8][3] = {
			{ 0, Ny - 1, 1 }, { 0, Ny - 1, 1 }, { Ny - 1, 0, -1 }, { Ny - 1, 0, -1 },
			{ 0, Ny - 1, 1 }, { 0, Ny - 1, 1 }, { Ny - 1, 0, -1 }, { Ny - 1, 0, -1 } };
		const int dirZ[8][3] = {
			{ 0, Nz - 1, 1 }, { 0, Nz - 1, 1 }, { 0, Nz - 1, 1 }, { Nz - 1, 0, -1 },
			{ Nz - 1, 0, -1 }, { Nz - 1, 0, -1 }, { Nz - 1, 0, -1 }, { 0, Nz - 1, 1 } };

		int s, ix, iy, iz, gridPos;
		double aa[3], tmp, eps = 1e-6;
		double d_curr, d_new, a, b, c, D;

		for (s = 0; s < NSweeps; s++) {
			// std::cout << "sweep " << s << " ... " << std::endl;
			for (iz = dirZ[s][0]; dirZ[s][2] * iz <= dirZ[s][1]; iz += dirZ[s][2]) {
				for (iy = dirY[s][0]; dirY[s][2] * iy <= dirY[s][1]; iy += dirY[s][2]) {
					for (ix = dirX[s][0]; dirX[s][2] * ix <= dirX[s][1]; ix += dirX[s][2]) {

						gridPos = ((iz * Ny + iy) * Nx + ix);
						if (gridFrozen[gridPos])
							continue;

						// === neighboring cells (Upwind Godunov) ===
						if (iz == 0 || iz == (Nz - 1)) 
						{
							if (iz == 0)
							{
								aa[2] = gridValues[gridPos] < gridValues[((iz + 1) * Ny + iy) * Nx + ix] ? gridValues[gridPos] : gridValues[((iz + 1) * Ny + iy) * Nx + ix];
							}
							if (iz == (Nz - 1)) 
							{
								aa[2] = gridValues[((iz - 1) * Ny + iy) * Nx + ix] < gridValues[gridPos] ? gridValues[((iz - 1) * Ny + iy) * Nx + ix] : gridValues[gridPos];
							}
						}
						else 
						{
							aa[2] = gridValues[((iz - 1) * Ny + iy) * Nx + ix] < gridValues[((iz + 1) * Ny + iy) * Nx + ix] ? gridValues[((iz - 1) * Ny + iy) * Nx + ix] : gridValues[((iz + 1) * Ny + iy) * Nx + ix];
						}

						if (iy == 0 || iy == (Ny - 1)) 
						{
							if (iy == 0)
							{
								aa[1] = gridValues[gridPos] < gridValues[(iz * Ny + (iy + 1)) * Nx + ix] ? gridValues[gridPos] : gridValues[(iz * Ny + (iy + 1)) * Nx + ix];
							}
							if (iy == (Ny - 1))
							{
								aa[1] = gridValues[(iz * Ny + (iy - 1)) * Nx + ix] < gridValues[gridPos] ? gridValues[(iz * Ny + (iy - 1)) * Nx + ix] : gridValues[gridPos];
							}
						}
						else 
						{
							aa[1] = gridValues[(iz * Ny + (iy - 1)) * Nx + ix] < gridValues[(iz * Ny + (iy + 1)) * Nx + ix] ? gridValues[(iz * Ny + (iy - 1)) * Nx + ix] : gridValues[(iz * Ny + (iy + 1)) * Nx + ix];
						}

						if (ix == 0 || ix == (Nx - 1)) 
						{
							if (ix == 0) 
							{
								aa[0] = gridValues[gridPos] < gridValues[(iz * Ny + iy) * Nx + (ix + 1)] ? gridValues[gridPos] : gridValues[(iz * Ny + iy) * Nx + (ix + 1)];
							}
							if (ix == (Nx - 1))
							{
								aa[0] = gridValues[(iz * Ny + iy) * Nx + (ix - 1)] < gridValues[gridPos] ? gridValues[(iz * Ny + iy) * Nx + (ix - 1)] : gridValues[gridPos];
							}
						}
						else 
						{
							aa[0] = gridValues[(iz * Ny + iy) * Nx + (ix - 1)] < gridValues[(iz * Ny + iy) * Nx + (ix + 1)] ? gridValues[(iz * Ny + iy) * Nx + (ix - 1)] : gridValues[(iz * Ny + iy) * Nx + (ix + 1)];
						}

						// simple bubble sort
						if (aa[0] > aa[1]) { tmp = aa[0]; aa[0] = aa[1]; aa[1] = tmp; }
						if (aa[1] > aa[2]) { tmp = aa[1]; aa[1] = aa[2]; aa[2] = tmp; }
						if (aa[0] > aa[1]) { tmp = aa[0]; aa[0] = aa[1]; aa[1] = tmp; }

						d_curr = aa[0] + h * f; // first estimate

						if (d_curr <= (aa[1] + eps)) // second estimate
						{
							d_new = d_curr;
						}
						else 
						{
							// third estimate
							a = 2.0; b = -2.0 * (aa[0] + aa[1]);
							c = aa[0] * aa[0] + aa[1] * aa[1] - h * h * f * f;
							D = b * b - 4.0 * a * c;

							if (D > 0)
							{
								d_curr = (-b + sqrt(D));
							}
							else
							{
								d_curr = (-b + sqrt(-D));
							}
							d_curr /= (2.0 * a);

							if (d_curr <= (aa[2] + eps))
								d_new = d_curr;
							else 
							{
								a = 3.0;
								b = -2.0 * (aa[0] + aa[1] + aa[2]);
								c = aa[0] * aa[0] + aa[1] * aa[1] + aa[2] * aa[2] - h * h * f * f;
								D = b * b - 4.0 * a * c;

								if (D > 0)
								{
									d_new = (-b + sqrt(D));
								}
								else
								{
									d_new = (-b + sqrt(-D));
								}
								d_new /= (2.0 * a);
							}
						}

						if (gridValues[gridPos] >= d_new)
							gridValues[gridPos] = d_new;						
					}
				}
			}
		}

	}
}