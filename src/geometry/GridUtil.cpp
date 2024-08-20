#include "GridUtil.h"

#include "pmp/SurfaceMesh.h"
#include "pmp/algorithms/BarycentricCoordinates.h"
#include "pmp/algorithms/TriangleKdTree.h"

namespace Geometry
{
	/// \brief constant for kernel radius of a "narrow" kernel.
	constexpr unsigned int NARROW_KERNEL_RADIUS = 1;
	/// \brief constant for kernel radius of a "wide" kernel.
	constexpr unsigned int WIDE_KERNEL_RADIUS = 2;

	/// \brief a blur kernel value wrapper.
	struct BlurKernel
	{
		std::vector<double> KernelValues{};
		unsigned int Radius;
	};

	/**
	 * \brief Computes an averaging kernel for a voxel grid with a given radius.
	 * \param radius     number of voxels (in each axis direction) by which the kernel reaches beyond the central voxel.
	 * \return a kernel with radius and values.
	 */
	[[nodiscard]] BlurKernel GetAveragingKernel(const unsigned int& radius)
	{
		assert(radius > 0);
		const unsigned int nCellsAxis = (2 * radius + 1);
		const unsigned int size = nCellsAxis * nCellsAxis * nCellsAxis;
		const double voxWeight = 1.0 / size;
		return { std::vector(size, voxWeight), radius };
	}

	/**
	 * \brief Computes an averaging kernel for a 2D pixel grid with a given radius.
	 * \param radius     number of pixel (in each axis direction) by which the kernel reaches beyond the central pixel.
	 * \return a kernel with radius and values.
	 */
	[[nodiscard]] BlurKernel GetAveragingKernel2D(const unsigned int& radius)
	{
		assert(radius > 0);
		const unsigned int nCellsAxis = (2 * radius + 1);
		const unsigned int size = nCellsAxis * nCellsAxis;
		const double cellWeight = 1.0 / size;
		return { std::vector(size, cellWeight), radius };
	}



	/**
	 * \brief Universal internal procedure for applying a blur kernel onto a ScalarGrid.
	 * \param grid      input grid.
	 * \param kernel    blur kernel to be applied.
	 */
	void ApplyBlurKernelInternal(ScalarGrid& grid, const BlurKernel& kernel)
	{
		const auto& kernelVals = kernel.KernelValues;
		const int rad = static_cast<int>(kernel.Radius);
		const auto& values = grid.Values();
		auto resultFieldValues = values;
		const auto& dims = grid.Dimensions();
		const int Nx = static_cast<int>(dims.Nx), Ny = static_cast<int>(dims.Ny), Nz = static_cast<int>(dims.Nz);

		const int nKernelCellsAxis = 2 * rad + 1;
		const int nKernelCellsAxisSq = nKernelCellsAxis * nKernelCellsAxis;
		const auto applyKernelToVoxel = [&](const int& ix, const int& iy, const int& iz)
		{
			double weightedSum = 0.0;
			for (int i = -rad; i <= rad; i++)
			{
				for (int j = -rad; j <= rad; j++)
				{
					for (int k = -rad; k <= rad; k++)
					{
						// =========== x-indices =====================
						int iKernX = nKernelCellsAxis + i, iGridX = ix + i;
						// ...... x-index bounds adjustment .....
						if (iGridX < 0)
						{
							iGridX = 0; iKernX = 0;
						}
						else if (iGridX > Nx - 1)
						{
							iGridX = Nx - 1; iKernX = nKernelCellsAxis - 1;
						}

						// =========== y-indices =====================
						int iKernY = nKernelCellsAxis + j, iGridY = iy + j;
						// ...... y-index bounds adjustment .....
						if (iGridY < 0)
						{
							iGridY = 0; iKernY = 0;
						}
						else if (iGridY > Ny - 1)
						{
							iGridY = Ny - 1; iKernY = nKernelCellsAxis - 1;
						}

						// =========== z-indices =====================
						int iKernZ = nKernelCellsAxis + k, iGridZ = iz + k;
						// ...... y-index bounds adjustment .....
						if (iGridZ < 0)
						{
							iGridZ = 0; iKernZ = 0;
						}
						else if (iGridZ > Nz - 1)
						{
							iGridZ = Nz - 1; iKernZ = nKernelCellsAxis - 1;
						}

						const unsigned int kernelGridPos = nKernelCellsAxisSq * iKernZ + nKernelCellsAxis * iKernY + iKernX;
						const unsigned int blurredGridPos = Nx * Ny * iGridZ + Nx * iGridY + iGridX;

						weightedSum += values[blurredGridPos] * kernelVals[kernelGridPos];
					}
				}
			}
			resultFieldValues[Nx * Ny * iz + Nx * iy + ix] = weightedSum;
		};

		// apply kernel for each voxel;
		for (int iz = 0; iz < Nz; iz++)
		{
			for (int iy = 0; iy < Ny; iy++)
			{
				for (int ix = 0; ix < Nx; ix++)
				{
					applyKernelToVoxel(ix, iy, iz);
				}
			}
		}
	}

	/**
	 * \brief Universal internal procedure for applying a blur kernel onto a ScalarGrid2D.
	 * \param grid      input grid.
	 * \param kernel    blur kernel to be applied.
	 */
	void ApplyBlurKernelInternal(ScalarGrid2D& grid, const BlurKernel& kernel)
	{
		const auto& kernelVals = kernel.KernelValues;
		const int rad = static_cast<int>(kernel.Radius);
		const auto& values = grid.Values();
		auto resultFieldValues = values;
		const auto& dims = grid.Dimensions();
		const int Nx = static_cast<int>(dims.Nx), Ny = static_cast<int>(dims.Ny);

		const int nKernelCellsAxis = 2 * rad + 1;
		const auto applyKernelToGridPoint = [&](const int& ix, const int& iy)
		{
			double weightedSum = 0.0;
			for (int i = -rad; i <= rad; i++)
			{
				for (int j = -rad; j <= rad; j++)
				{
					// =========== x-indices =====================
					int iKernX = rad + i, iGridX = ix + i;
					// ...... x-index bounds adjustment .....
					if (iGridX < 0)
					{
						iGridX = 0; iKernX = 0;
					}
					else if (iGridX > Nx - 1)
					{
						iGridX = Nx - 1; iKernX = nKernelCellsAxis - 1;
					}

					// =========== y-indices =====================
					int iKernY = rad + j, iGridY = iy + j;
					// ...... y-index bounds adjustment .....
					if (iGridY < 0)
					{
						iGridY = 0; iKernY = 0;
					}
					else if (iGridY > Ny - 1)
					{
						iGridY = Ny - 1; iKernY = nKernelCellsAxis - 1;
					}

					const unsigned int kernelGridPos = nKernelCellsAxis * iKernY + iKernX;
					const unsigned int blurredGridPos = Nx * iGridY + iGridX;

					weightedSum += values[blurredGridPos] * kernelVals[kernelGridPos];
				}
			}
			resultFieldValues[Nx * iy + ix] = weightedSum;
		};

		// Apply kernel for each grid point
		for (int iy = 0; iy < Ny; iy++)
		{
			for (int ix = 0; ix < Nx; ix++)
			{
				applyKernelToGridPoint(ix, iy);
			}
		}

		grid.Values() = resultFieldValues;
	}

	// ---------------------------------------------------

	void ApplyNarrowAveragingBlur(ScalarGrid& grid)
    {
		constexpr int rad = static_cast<int>(NARROW_KERNEL_RADIUS);
		const auto kernel = GetAveragingKernel(rad);
		ApplyBlurKernelInternal(grid, kernel);
    }

    void ApplyWideAveragingBlur(ScalarGrid& grid)
    {
		constexpr int rad = static_cast<int>(WIDE_KERNEL_RADIUS);
		const auto kernel = GetAveragingKernel(rad);
		ApplyBlurKernelInternal(grid, kernel);
    }

	// ---------------------------------------------------

	void ApplyNarrowAveragingBlur2D(ScalarGrid2D& grid)
	{
		constexpr int rad = static_cast<int>(NARROW_KERNEL_RADIUS);
		const auto kernel = GetAveragingKernel2D(rad);
		ApplyBlurKernelInternal(grid, kernel);
	}

	void ApplyWideAveragingBlur2D(ScalarGrid2D& grid)
	{
		constexpr int rad = static_cast<int>(WIDE_KERNEL_RADIUS);
		const auto kernel = GetAveragingKernel2D(rad);
		ApplyBlurKernelInternal(grid, kernel);
	}

	/**
	 * \brief Computes a Gaussian kernel for a voxel grid with a given radius.
	 * \param radius     number of voxels (in each axis direction) by which the kernel reaches beyond the central voxel.
	 * \return a kernel with radius and values.
	 */
	[[nodiscard]] BlurKernel GetGaussianKernel(const unsigned int& radius)
	{
		assert(radius > 0);
		const unsigned int nCellsAxis = (2 * radius + 1);
		const unsigned int size = nCellsAxis * nCellsAxis * nCellsAxis;
		const unsigned int nCellsAxisSq = nCellsAxis * nCellsAxis;
		auto resultVals = std::vector(size, 0.0);

		const double sigma = radius / 3.0;
		double wSum = 0.0;

		// direct integration of the kernel using error function (erf)
		for (unsigned int i = 0; i <= 2 * radius; i++)
		{
			for (unsigned int j = 0; j <= 2 * radius; j++)
			{
				for (unsigned int k = 0; k <= 2 * radius; k++)
				{
					const double x0 = static_cast<double>(2 * static_cast<int>(i) - 1) * 0.5 - radius;
					const double x1 = static_cast<double>(2 * static_cast<int>(i) + 1) * 0.5 - radius;

					const double y0 = static_cast<double>(2 * static_cast<int>(j) - 1) * 0.5 - radius;
					const double y1 = static_cast<double>(2 * static_cast<int>(j) + 1) * 0.5 - radius;

					const double z0 = static_cast<double>(2 * static_cast<int>(k) - 1) * 0.5 - radius;
					const double z1 = static_cast<double>(2 * static_cast<int>(k) + 1) * 0.5 - radius;

					// 3D gaussians are independent G(x,y,z) = G(x) G(y) G(z) and so are their integrals
					const double voxIntegral = 0.25 * 
						(erf(x1 / (M_SQRT2 * sigma)) - erf(x0 / (M_SQRT2 * sigma))) *
						(erf(y1 / (M_SQRT2 * sigma)) - erf(y0 / (M_SQRT2 * sigma))) * 
						(erf(z1 / (M_SQRT2 * sigma)) - erf(z0 / (M_SQRT2 * sigma)));

					const unsigned int id = i * nCellsAxisSq + j * nCellsAxis + k;
					resultVals[id] = voxIntegral;
					wSum += voxIntegral;					
				}
			}
		}

		// normalize
		for (auto& val : resultVals)
			val /= wSum;

		return {resultVals, radius };
	}

	/**
	 * \brief Computes a Gaussian kernel for a pixel grid with a given radius.
	 * \param radius     number of pixel (in each axis direction) by which the kernel reaches beyond the central pixel.
	 * \return a kernel with radius and values.
	 */
	[[nodiscard]] BlurKernel GetGaussianKernel2D(const unsigned int& radius)
	{
		assert(radius > 0);
		const unsigned int nCellsAxis = (2 * radius + 1);
		const unsigned int size = nCellsAxis * nCellsAxis;
		auto resultVals = std::vector(size, 0.0);

		const double sigma = radius / 3.0;
		double wSum = 0.0;

		// Direct integration of the kernel using error function (erf)
		for (unsigned int i = 0; i <= 2 * radius; i++)
		{
			for (unsigned int j = 0; j <= 2 * radius; j++)
			{
				const double x0 = static_cast<double>(2 * static_cast<int>(i) - 1) * 0.5 - radius;
				const double x1 = static_cast<double>(2 * static_cast<int>(i) + 1) * 0.5 - radius;

				const double y0 = static_cast<double>(2 * static_cast<int>(j) - 1) * 0.5 - radius;
				const double y1 = static_cast<double>(2 * static_cast<int>(j) + 1) * 0.5 - radius;

				// 2D gaussians are independent G(x,y) = G(x) G(y) and so are their integrals
				const double cellIntegral = 0.25 *
					(std::erf(x1 / (M_SQRT2 * sigma)) - std::erf(x0 / (M_SQRT2 * sigma))) *
					(std::erf(y1 / (M_SQRT2 * sigma)) - std::erf(y0 / (M_SQRT2 * sigma)));

				const unsigned int id = i * nCellsAxis + j;
				resultVals[id] = cellIntegral;
				wSum += cellIntegral;
			}
		}

		// Normalize
		for (auto& val : resultVals)
			val /= wSum;

		return { resultVals, radius };
	}

	// -----------------------------------------------------

    void ApplyNarrowGaussianBlur(ScalarGrid& grid)
    {
		constexpr int rad = static_cast<int>(NARROW_KERNEL_RADIUS);
		const auto kernel = GetGaussianKernel(rad);
		ApplyBlurKernelInternal(grid, kernel);
    }

    void ApplyWideGaussianBlur(ScalarGrid& grid)
    {
		constexpr int rad = static_cast<int>(WIDE_KERNEL_RADIUS);
		const auto kernel = GetGaussianKernel(rad);
		ApplyBlurKernelInternal(grid, kernel);
    }

	// --------------------------------------------------------

	void ApplyNarrowGaussianBlur2D(ScalarGrid2D& grid)
	{
		constexpr int rad = static_cast<int>(NARROW_KERNEL_RADIUS);
		const auto kernel = GetGaussianKernel2D(rad);
		ApplyBlurKernelInternal(grid, kernel);
	}

	void ApplyWideGaussianBlur2D(ScalarGrid2D& grid)
	{
		constexpr int rad = static_cast<int>(WIDE_KERNEL_RADIUS);
		const auto kernel = GetGaussianKernel2D(rad);
		ApplyBlurKernelInternal(grid, kernel);
	}

	//
	// ============================================================================
	//

	void NegateGridSubVolume(ScalarGrid& grid, const pmp::BoundingBox& subBox)
	{
		const auto& gridBox = grid.Box();
		const auto dMin = subBox.min() - gridBox.min();
		const auto dMax = subBox.max() - gridBox.min();
		assert(dMin[0] >= 0.0f && dMin[1] >= 0.0f && dMin[2] >= 0.0f);
		assert(dMax[0] >= 0.0f && dMax[1] >= 0.0f && dMax[2] >= 0.0f);

		// compute sub-grid bounds
		const auto& cellSize = grid.CellSize();
		const unsigned int iXStart = std::floor(dMin[0] / cellSize);
		const unsigned int iYStart = std::floor(dMin[1] / cellSize);
		const unsigned int iZStart = std::floor(dMin[2] / cellSize);

		const unsigned int iXEnd = std::ceil(dMax[0] / cellSize);
		const unsigned int iYEnd = std::ceil(dMax[1] / cellSize);
		const unsigned int iZEnd = std::ceil(dMax[2] / cellSize);
		const auto& dims = grid.Dimensions();
		auto& values = grid.Values();
		
		for (unsigned int iz = iZStart; iz < iZEnd; iz++)
		{
			for (unsigned int iy = iYStart; iy < iYEnd; iy++)
			{
				for (unsigned int ix = iXStart; ix < iXEnd; ix++)
				{
					const unsigned int gridPos = dims.Nx * dims.Ny * iz + dims.Nx * iy + ix;
					values[gridPos] *= -1.0;
				}
			}
		}
	}

	/// \brief a validation helper for scalar grid values. Also increments nanCount and infCount if encountering nan or inf.
	[[nodiscard]] bool IsGridValueValid(const double& val, size_t& nanCount, size_t& infCount)
	{
		if (std::isnan(val) || !std::isnormal(val))
		{
			nanCount++;
			return false;
		}

		if (std::isinf(val) || std::abs<double>(val) >= DBL_MAX)
		{
			infCount++;
			return false;
		}

		return true;
	}

	/// \brief a basic validation helper for scalar grid values.
	[[nodiscard]] bool IsGridValueValidBasic(const double& val)
	{
		if (std::isnan(val) || !std::isnormal(val))
		{
			return false;
		}

		if (std::isinf(val) || std::abs<double>(val) >= DBL_MAX)
		{
			return false;
		}

		return true;
	}

	/// \brief a boundary validation helper for scalar grid values.
	[[nodiscard]] bool IsABoundaryCell(
		const int& gridPosPrevX, const int& gridPosNextX,
		const int& gridPosPrevY, const int& gridPosNextY,
		const int& gridPosPrevZ, const int& gridPosNextZ, const size_t& gridSize)
	{
		if (gridPosPrevX < 0)
			return true;

		if (gridPosPrevY < 0)
			return true;

		if (gridPosPrevZ < 0)
			return true;

		if (gridPosNextX >= static_cast<int>(gridSize))
			return true;

		if (gridPosNextY >= static_cast<int>(gridSize))
			return true;

		if (gridPosNextZ >= static_cast<int>(gridSize))
			return true;

		return false;
	}

	void RepairScalarGrid(ScalarGrid& grid)
	{
		size_t nanCount = 0;
		size_t infCount = 0;

		auto& values = grid.Values();
		const auto nValues = values.size();
		const auto& dim = grid.Dimensions();

		const auto Nx = static_cast<unsigned int>(dim.Nx);
		const auto Ny = static_cast<unsigned int>(dim.Ny);
		const auto Nz = static_cast<unsigned int>(dim.Nz);

		std::cout << "----------------------------------------------------------\n";
		std::cout << "RepairScalarGrid: repairing scalar grid...\n";

		for (unsigned int iz = 0; iz < Nz; iz++) 
		{
			for (unsigned int iy = 0; iy < Ny; iy++) 
			{
				for (unsigned int ix = 0; ix < Nx; ix++)
				{
					const unsigned int gridPos = Nx * Ny * iz + Nx * iy + ix;

					if (IsGridValueValid(values[gridPos], nanCount, infCount))
						continue;

					// check neighbors
					const int gridPosPrevX = Nx * Ny * iz + Nx * iy + (ix - 1);
					const int gridPosNextX = Nx * Ny * iz + Nx * iy + (ix + 1);

					const int gridPosPrevY = Nx * Ny * iz + Nx * (iy - 1) + ix;
					const int gridPosNextY = Nx * Ny * iz + Nx * (iy + 1) + ix;

					const int gridPosPrevZ = Nx * Ny * (iz - 1) + Nx * iy + ix;
					const int gridPosNextZ = Nx * Ny * (iz + 1) + Nx * iy + ix;

					if (IsABoundaryCell(gridPosPrevX, gridPosNextX, gridPosPrevY, gridPosNextY, gridPosPrevZ, gridPosNextZ, nValues))
					{
						values[gridPos] = DEFAULT_SCALAR_GRID_INIT_VAL;
						continue;
					}

					if (!IsGridValueValidBasic(values[gridPosPrevX]) || !IsGridValueValidBasic(values[gridPosNextX]) ||
						!IsGridValueValidBasic(values[gridPosPrevY]) || !IsGridValueValidBasic(values[gridPosNextY]) ||
						!IsGridValueValidBasic(values[gridPosPrevZ]) || !IsGridValueValidBasic(values[gridPosNextZ]))
					{
						values[gridPos] = DEFAULT_SCALAR_GRID_INIT_VAL;
						continue;
					}

					// all neighbors are valid. Computing average value:
					values[gridPos] = 
						(values[gridPosPrevX] + values[gridPosNextX] +
						 values[gridPosPrevY] + values[gridPosNextY] +
						 values[gridPosPrevZ] + values[gridPosNextZ]) / 6.0;
				}
			}
		}

		if (nanCount > 0 || infCount > 0)
		{
			std::cout << "RepairScalarGrid: [WARNING]: Encountered & repaired invalid cell values! nanCount: " << nanCount << ", infCount: " << infCount << ".\n";
			std::cout << "----------------------------------------------------------\n";
			return;
		}

		std::cout << "RepairScalarGrid: All grid values are valid.\n";
		std::cout << "----------------------------------------------------------\n";
	}

	void NormalizeScalarGridValues(ScalarGrid& grid)
	{
		// find max absolute value
		double maxAbsVal = 0;
		for (const auto& val : grid.Values())
		{
			if (std::abs(val) > maxAbsVal) maxAbsVal = std::abs(val);
		}
		if (maxAbsVal < DBL_EPSILON)
		{
			std::cerr << "NormalizeScalarGridValues: maxAbsVal < DBL_EPSILON!\n";
			return;
		}

		// normalize
		for (auto& val : grid.Values())
		{
			val = val / maxAbsVal; // Normalize to [-1, 1]
		}
	}

	/// \brief if true, additional checks will be performed before writing gradient values.
#define EXPECT_INVALID_VALUES false

	VectorGrid ComputeGradient(const ScalarGrid& scalarGrid)
	{
		if (!scalarGrid.IsValid())
			throw std::invalid_argument("ComputeGradient: scalarGrid to be processed is invalid!\n");

		VectorGrid result(scalarGrid);
		const auto cellSize = static_cast<double>(result.CellSize());
		const auto& dim = result.Dimensions();

		const auto Nx = static_cast<unsigned int>(dim.Nx);
		const auto Ny = static_cast<unsigned int>(dim.Ny);
		const auto Nz = static_cast<unsigned int>(dim.Nz);

		const auto& gridValues = scalarGrid.Values();

		auto& gradValsX = result.ValuesX();
		auto& gradValsY = result.ValuesY();
		auto& gradValsZ = result.ValuesZ();

		for (unsigned int iz = 1; iz < Nz - 1; iz++) {
			for (unsigned int iy = 1; iy < Ny - 1; iy++) {
				for (unsigned int ix = 1; ix < Nx - 1; ix++) {

					const unsigned int gridPosPrevX = Nx * Ny * iz + Nx * iy + (ix - 1);
					const unsigned int gridPosNextX = Nx * Ny * iz + Nx * iy + (ix + 1);

					const unsigned int gridPosPrevY = Nx * Ny * iz + Nx * (iy - 1) + ix;
					const unsigned int gridPosNextY = Nx * Ny * iz + Nx * (iy + 1) + ix;

					const unsigned int gridPosPrevZ = Nx * Ny * (iz - 1) + Nx * iy + ix;
					const unsigned int gridPosNextZ = Nx * Ny * (iz + 1) + Nx * iy + ix;

					const unsigned int gradPos = Nx * Ny * iz + Nx * iy + ix;

					// central difference for non-boundary voxels
					const double grad_x = (gridValues[gridPosNextX] - gridValues[gridPosPrevX]) / (2.0 * cellSize);
					const double grad_y = (gridValues[gridPosNextY] - gridValues[gridPosPrevY]) / (2.0 * cellSize);
					const double grad_z = (gridValues[gridPosNextZ] - gridValues[gridPosPrevZ]) / (2.0 * cellSize);

#if EXPECT_INVALID_VALUES
					if (std::isnan(grad_x) || std::isinf(grad_x) ||
						std::isnan(grad_y) || std::isinf(grad_y) ||
						std::isnan(grad_y) || std::isinf(grad_z))
					{
						const std::string msg = "ComputeGradient: nans or infs encountered for cell " + std::to_string(gradPos) + "! Setting value to zero.\n";
						assert(false);
						std::cerr << msg;
						continue;
					}
#endif

					gradValsX[gradPos] = grad_x;
					gradValsY[gradPos] = grad_y;
					gradValsZ[gradPos] = grad_z;
				}
			}
		}

		return result;
	}

	VectorGrid2D ComputeGradient(const ScalarGrid2D& scalarGrid)
	{
		if (!scalarGrid.IsValid())
			throw std::invalid_argument("ComputeGradient: scalarGrid to be processed is invalid!\n");

		VectorGrid2D result(scalarGrid);
		const auto cellSize = static_cast<double>(result.CellSize());
		const auto& dim = result.Dimensions();

		const auto Nx = static_cast<unsigned int>(dim.Nx);
		const auto Ny = static_cast<unsigned int>(dim.Ny);

		const auto& gridValues = scalarGrid.Values();

		auto& gradValsX = result.ValuesX();
		auto& gradValsY = result.ValuesY();

			for (unsigned int iy = 1; iy < Ny - 1; iy++) {
				for (unsigned int ix = 1; ix < Nx - 1; ix++) {

					const unsigned int gridPosPrevX = Nx * iy + (ix - 1);
					const unsigned int gridPosNextX = Nx * iy + (ix + 1);

					const unsigned int gridPosPrevY = Nx * (iy - 1) + ix;
					const unsigned int gridPosNextY = Nx * (iy + 1) + ix;

					const unsigned int gradPos = Nx * iy + ix;

					// central difference for non-boundary pixels
					const double grad_x = (gridValues[gridPosNextX] - gridValues[gridPosPrevX]) / (2.0 * cellSize);
					const double grad_y = (gridValues[gridPosNextY] - gridValues[gridPosPrevY]) / (2.0 * cellSize);

#if EXPECT_INVALID_VALUES
					if (std::isnan(grad_x) || std::isinf(grad_x) ||
						std::isnan(grad_y) || std::isinf(grad_y))
					{
						const std::string msg = "ComputeGradient: nans or infs encountered for cell " + std::to_string(gradPos) + "! Setting value to zero.\n";
						assert(false);
						std::cerr << msg;
						continue;
					}
#endif

					gradValsX[gradPos] = grad_x;
					gradValsY[gradPos] = grad_y;
				}
			}

		return result;
	}

	//! tolerance for gradient vector norms
	constexpr double NORM_EPSILON = 1e-6;

	VectorGrid ComputeNormalizedGradient(const ScalarGrid& scalarGrid)
	{
		if (!scalarGrid.IsValid())
			throw std::invalid_argument("ComputeGradient: scalarGrid to be processed is invalid!\n");

		VectorGrid result(scalarGrid);
		const auto cellSize = static_cast<double>(result.CellSize());
		const auto& dim = result.Dimensions();

		const auto Nx = static_cast<unsigned int>(dim.Nx);
		const auto Ny = static_cast<unsigned int>(dim.Ny);
		const auto Nz = static_cast<unsigned int>(dim.Nz);

		const auto& gridValues = scalarGrid.Values();

		auto& gradValsX = result.ValuesX();
		auto& gradValsY = result.ValuesY();
		auto& gradValsZ = result.ValuesZ();

		for (unsigned int iz = 1; iz < Nz - 1; iz++) {
			for (unsigned int iy = 1; iy < Ny - 1; iy++) {
				for (unsigned int ix = 1; ix < Nx - 1; ix++) {

					const unsigned int gridPosPrevX = Nx * Ny * iz + Nx * iy + (ix - 1);
					const unsigned int gridPosNextX = Nx * Ny * iz + Nx * iy + (ix + 1);

					const unsigned int gridPosPrevY = Nx * Ny * iz + Nx * (iy - 1) + ix;
					const unsigned int gridPosNextY = Nx * Ny * iz + Nx * (iy + 1) + ix;

					const unsigned int gridPosPrevZ = Nx * Ny * (iz - 1) + Nx * iy + ix;
					const unsigned int gridPosNextZ = Nx * Ny * (iz + 1) + Nx * iy + ix;

					const unsigned int gradPos = Nx * Ny * iz + Nx * iy + ix;

					// central difference for non-boundary voxels
					const double grad_x = (gridValues[gridPosNextX] - gridValues[gridPosPrevX]) / (2.0 * cellSize);
					const double grad_y = (gridValues[gridPosNextY] - gridValues[gridPosPrevY]) / (2.0 * cellSize);
					const double grad_z = (gridValues[gridPosNextZ] - gridValues[gridPosPrevZ]) / (2.0 * cellSize);

#if EXPECT_INVALID_VALUES
					if (std::isnan(grad_x) || std::isinf(grad_x) ||
						std::isnan(grad_y) || std::isinf(grad_y) ||
						std::isnan(grad_y) || std::isinf(grad_z))
					{
						const std::string msg = "ComputeNormalizedGradient: nans or infs encountered for cell " + std::to_string(gradPos) + "! Setting value to zero.\n";
						assert(false);
						std::cerr << msg;
						continue;
					}
#endif

					const double norm = sqrt(grad_x * grad_x + grad_y * grad_y + grad_z * grad_z);

					assert(!std::isnan(norm) && !std::isinf(norm));

					if (norm < NORM_EPSILON)
						continue; // keep zero init val

					gradValsX[gradPos] = grad_x / norm;
					gradValsY[gradPos] = grad_y / norm;
					gradValsZ[gradPos] = grad_z / norm;
				}
			}
		}

		return result;
	}

	VectorGrid ComputeNormalizedNegativeGradient(const ScalarGrid& scalarGrid)
	{
		if (!scalarGrid.IsValid())
			throw std::invalid_argument("ComputeNormalizedNegativeGradient: scalarGrid to be processed is invalid!\n");

		VectorGrid result(scalarGrid);
		const auto cellSize = static_cast<double>(result.CellSize());
		const auto& dim = result.Dimensions();

		const auto Nx = static_cast<unsigned int>(dim.Nx);
		const auto Ny = static_cast<unsigned int>(dim.Ny);
		const auto Nz = static_cast<unsigned int>(dim.Nz);

		const auto& gridValues = scalarGrid.Values();

		auto& gradValsX = result.ValuesX();
		auto& gradValsY = result.ValuesY();
		auto& gradValsZ = result.ValuesZ();

		for (unsigned int iz = 1; iz < Nz - 1; iz++) {
			for (unsigned int iy = 1; iy < Ny - 1; iy++) {
				for (unsigned int ix = 1; ix < Nx - 1; ix++) {

					const unsigned int gridPosPrevX = Nx * Ny * iz + Nx * iy + (ix - 1);
					const unsigned int gridPosNextX = Nx * Ny * iz + Nx * iy + (ix + 1);

					const unsigned int gridPosPrevY = Nx * Ny * iz + Nx * (iy - 1) + ix;
					const unsigned int gridPosNextY = Nx * Ny * iz + Nx * (iy + 1) + ix;

					const unsigned int gridPosPrevZ = Nx * Ny * (iz - 1) + Nx * iy + ix;
					const unsigned int gridPosNextZ = Nx * Ny * (iz + 1) + Nx * iy + ix;

					const unsigned int gradPos = Nx * Ny * iz + Nx * iy + ix;

					// central difference for non-boundary voxels
					const double grad_x = (gridValues[gridPosNextX] - gridValues[gridPosPrevX]) / (2.0 * cellSize);
					const double grad_y = (gridValues[gridPosNextY] - gridValues[gridPosPrevY]) / (2.0 * cellSize);
					const double grad_z = (gridValues[gridPosNextZ] - gridValues[gridPosPrevZ]) / (2.0 * cellSize);

#if EXPECT_INVALID_VALUES
					if (std::isnan(grad_x) || std::isinf(grad_x) ||
						std::isnan(grad_y) || std::isinf(grad_y) ||
						std::isnan(grad_y) || std::isinf(grad_z))
					{
						const std::string msg = "ComputeNormalizedNegativeGradient: nans or infs encountered for cell " + std::to_string(gradPos) + "! Setting value to zero.\n";
						assert(false);
						std::cerr << msg;
						continue;
					}
#endif

					const double norm = -1.0 * sqrt(grad_x * grad_x + grad_y * grad_y + grad_z * grad_z);

					assert(!std::isnan(norm) && !std::isinf(norm));

					if (norm > -NORM_EPSILON)
						continue; // keep zero init val

					gradValsX[gradPos] = grad_x / norm;
					gradValsY[gradPos] = grad_y / norm;
					gradValsZ[gradPos] = grad_z / norm;
				}
			}
		}

		return result;
	}

	VectorGrid2D ComputeNormalizedNegativeGradient(const ScalarGrid2D& scalarGrid)
	{
		if (!scalarGrid.IsValid())
			throw std::invalid_argument("ComputeNormalizedNegativeGradient: scalarGrid to be processed is invalid!\n");

		VectorGrid2D result(scalarGrid);
		const auto cellSize = static_cast<double>(result.CellSize());
		const auto& dim = result.Dimensions();

		const auto Nx = static_cast<unsigned int>(dim.Nx);
		const auto Ny = static_cast<unsigned int>(dim.Ny);

		const auto& gridValues = scalarGrid.Values();

		auto& gradValsX = result.ValuesX();
		auto& gradValsY = result.ValuesY();

		for (unsigned int iy = 1; iy < Ny - 1; iy++) {
			for (unsigned int ix = 1; ix < Nx - 1; ix++) {

				const unsigned int gridPosPrevX = Nx * iy + (ix - 1);
				const unsigned int gridPosNextX = Nx * iy + (ix + 1);

				const unsigned int gridPosPrevY = Nx * (iy - 1) + ix;
				const unsigned int gridPosNextY = Nx * (iy + 1) + ix;

				const unsigned int gradPos =  Nx * iy + ix;

				// central difference for non-boundary voxels
				const double grad_x = (gridValues[gridPosNextX] - gridValues[gridPosPrevX]) / (2.0 * cellSize);
				const double grad_y = (gridValues[gridPosNextY] - gridValues[gridPosPrevY]) / (2.0 * cellSize);

#if EXPECT_INVALID_VALUES
				if (std::isnan(grad_x) || std::isinf(grad_x) ||
					std::isnan(grad_y) || std::isinf(grad_y))
				{
					const std::string msg = "ComputeNormalizedNegativeGradient: nans or infs encountered for cell " + std::to_string(gradPos) + "! Setting value to zero.\n";
					assert(false);
					std::cerr << msg;
					continue;
				}
#endif

				const double norm = -1.0 * sqrt(grad_x * grad_x + grad_y * grad_y);

				assert(!std::isnan(norm) && !std::isinf(norm));

				if (norm > -NORM_EPSILON)
					continue; // keep zero init val

				gradValsX[gradPos] = grad_x / norm;
				gradValsY[gradPos] = grad_y / norm;
			}
		}

		return result;
	}

	// ==================================================================================================
	//                                     Interpolation utils 
	// --------------------------------------------------------------------------------------------------

	/// \brief zero index with respect to axis-aligned grid dimension.
	constexpr size_t ZERO_CELL_ID = 0;

	double GetNearestNeighborScalarValue(const pmp::vec3& samplePt, const ScalarGrid& grid)
	{
		const auto& [Nx, Ny, Nz] = grid.Dimensions();
		const auto& boxMin = grid.Box().min();
		const float cellSize = grid.CellSize();

		const auto ix = std::max(std::min(static_cast<size_t>(std::round((samplePt[0] - boxMin[0]) / cellSize)), Nx - 1), ZERO_CELL_ID);
		const auto iy = std::max(std::min(static_cast<size_t>(std::round((samplePt[1] - boxMin[1]) / cellSize)), Ny - 1), ZERO_CELL_ID);
		const auto iz = std::max(std::min(static_cast<size_t>(std::round((samplePt[2] - boxMin[2]) / cellSize)), Nz - 1), ZERO_CELL_ID);

		const auto i = Nx * Ny * iz + Nx * iy + ix;
		return grid.Values()[i];
	}

	pmp::dvec3 GetNearestNeighborVectorValue(const pmp::vec3& samplePt, const VectorGrid& grid)
	{
		const auto& [Nx, Ny, Nz] = grid.Dimensions();
		const auto& boxMin = grid.Box().min();
		const float cellSize = grid.CellSize();

		const auto ix = std::max(std::min(static_cast<size_t>(std::round((samplePt[0] - boxMin[0]) / cellSize)), Nx - 1), ZERO_CELL_ID);
		const auto iy = std::max(std::min(static_cast<size_t>(std::round((samplePt[1] - boxMin[1]) / cellSize)), Ny - 1), ZERO_CELL_ID);
		const auto iz = std::max(std::min(static_cast<size_t>(std::round((samplePt[2] - boxMin[2]) / cellSize)), Nz - 1), ZERO_CELL_ID);

		const auto i = Nx * Ny * iz + Nx * iy + ix;
		return pmp::dvec3(grid.ValuesX()[i], grid.ValuesY()[i], grid.ValuesZ()[i]);
	}

	// ------------------------------------------------

	double GetNearestNeighborScalarValue2D(const pmp::vec2& samplePt, const ScalarGrid2D& grid)
	{
		const auto& [Nx, Ny] = grid.Dimensions();
		const auto& boxMin = grid.Box().min();
		const float cellSize = grid.CellSize();

		const auto ix = std::max(std::min(static_cast<size_t>(std::round((samplePt[0] - boxMin[0]) / cellSize)), Nx - 1), ZERO_CELL_ID);
		const auto iy = std::max(std::min(static_cast<size_t>(std::round((samplePt[1] - boxMin[1]) / cellSize)), Ny - 1), ZERO_CELL_ID);

		const auto i = Nx * iy + ix;
		return grid.Values()[i];
	}

	pmp::dvec2 GetNearestNeighborVectorValue2D(const pmp::vec2& samplePt, const VectorGrid2D& grid)
	{
		const auto& [Nx, Ny] = grid.Dimensions();
		const auto& boxMin = grid.Box().min();
		const float cellSize = grid.CellSize();

		const auto ix = std::max(std::min(static_cast<size_t>(std::round((samplePt[0] - boxMin[0]) / cellSize)), Nx - 1), ZERO_CELL_ID);
		const auto iy = std::max(std::min(static_cast<size_t>(std::round((samplePt[1] - boxMin[1]) / cellSize)), Ny - 1), ZERO_CELL_ID);

		const auto i = Nx * iy + ix;
		return pmp::dvec2(grid.ValuesX()[i], grid.ValuesY()[i]);
	}

	// --------------------------------------------------------------------------------------
	//

	/// \brief A wrapper for grid interpolation scalar values
	struct SurroundingScalarCells
	{
		pmp::vec3 MinPt{};
		pmp::vec3 MaxPt{};
		// cell point values:
		double Val000{ 0.0 };
		double Val100{ 0.0 };
		double Val010{ 0.0 };
		double Val110{ 0.0 };
		double Val001{ 0.0 };
		double Val101{ 0.0 };
		double Val011{ 0.0 };
		double Val111{ 0.0 };
	};

	/**
	 * \brief Gets surrounding cell values of a sampled point from a scalar grid.
	 * \param samplePt    point where the grid is sampled.
	 * \param grid        interpolated grid.
	 * \return surrounding scalar cells.
	 *
	 * DISCLAIMER: For samplePt outside of grid.Box(), the values are clamped to boundary values.
	 */
	[[nodiscard]] SurroundingScalarCells GetSurroundingCells(const pmp::vec3& samplePt, const ScalarGrid& grid)
	{
		const auto& [Nx, Ny, Nz] = grid.Dimensions();
		const auto& boxMin = grid.Box().min();
		const float cellSize = grid.CellSize();

		const auto ix = std::max(std::min(static_cast<size_t>(std::floor((samplePt[0] - boxMin[0]) / cellSize)), Nx - 1), ZERO_CELL_ID);
		const auto iy = std::max(std::min(static_cast<size_t>(std::floor((samplePt[1] - boxMin[1]) / cellSize)), Ny - 1), ZERO_CELL_ID);
		const auto iz = std::max(std::min(static_cast<size_t>(std::floor((samplePt[2] - boxMin[2]) / cellSize)), Nz - 1), ZERO_CELL_ID);

		const auto ix1 = std::min(ix + 1, Nx - 1);
		const auto iy1 = std::min(iy + 1, Ny - 1);
		const auto iz1 = std::min(iz + 1, Nz - 1);

		const auto i000 = Nx * Ny * iz + Nx * iy + ix;
		const auto i100 = Nx * Ny * iz + Nx * iy + ix1;
		const auto i010 = Nx * Ny * iz + Nx * iy1 + ix;
		const auto i110 = Nx * Ny * iz + Nx * iy1 + ix1;
		const auto i001 = Nx * Ny * iz1 + Nx * iy + ix;
		const auto i101 = Nx * Ny * iz1 + Nx * iy + ix1;
		const auto i011 = Nx * Ny * iz1 + Nx * iy1 + ix;
		const auto i111 = Nx * Ny * iz1 + Nx * iy1 + ix1;

		const auto& values = grid.Values();
		return {
			pmp::vec3{
				boxMin[0] + static_cast<float>(ix) * cellSize,
				boxMin[1] + static_cast<float>(iy) * cellSize,
				boxMin[2] + static_cast<float>(iz) * cellSize
			},
			pmp::vec3{
				boxMin[0] + static_cast<float>(ix + 1) * cellSize,
				boxMin[1] + static_cast<float>(iy + 1) * cellSize,
				boxMin[2] + static_cast<float>(iz + 1) * cellSize
			},
			values[i000], values[i100], values[i010], values[i110],
			values[i001], values[i101], values[i011], values[i111]
		};
	}

	double TrilinearInterpolateScalarValue(const pmp::vec3& samplePt, const ScalarGrid& grid)
	{
		assert(grid.IsValid() && grid.CellSize() > 0.0f);
		const auto surrCellVals = GetSurroundingCells(samplePt, grid);
		const auto x = static_cast<double>(samplePt[0]), y = static_cast<double>(samplePt[1]), z = static_cast<double>(samplePt[2]);
		// cell min
		const auto x0 = static_cast<double>(surrCellVals.MinPt[0]), y0 = static_cast<double>(surrCellVals.MinPt[1]), z0 = static_cast<double>(surrCellVals.MinPt[2]);
		// cell max
		const auto x1 = static_cast<double>(surrCellVals.MaxPt[0]), y1 = static_cast<double>(surrCellVals.MaxPt[1]), z1 = static_cast<double>(surrCellVals.MaxPt[2]);

		// cell values
		const double c000 = surrCellVals.Val000, c100 = surrCellVals.Val100, c010 = surrCellVals.Val010, c110 = surrCellVals.Val110;
		const double c001 = surrCellVals.Val001, c101 = surrCellVals.Val101, c011 = surrCellVals.Val011, c111 = surrCellVals.Val111;

		const auto det = (x0 - x1) * (y0 - y1) * (z0 - z1);
		const double a0 =
			(c111 * x0 * y0 * z0 - c011 * x1 * y0 * z0 - c101 * x0 * y1 * z0 + c001 * x1 * y1 * z0 -
				c110 * x0 * y0 * z1 + c010 * x1 * y0 * z1 + c100 * x0 * y1 * z1 - c000 * x1 * y1 * z1) / det;

		const double a1 =
			(c011 * y0 * z0 - c111 * y0 * z0 - c001 * y1 * z0 + c101 * y1 * z0 - c010 * y0 * z1 +
				c110 * y0 * z1 + c000 * y1 * z1 - c100 * y1 * z1) / det;

		const double a2 =
			(c101 * x0 * z0 - c111 * x0 * z0 - c001 * x1 * z0 + c011 * x1 * z0 - c100 * x0 * z1 +
				c110 * x0 * z1 + c000 * x1 * z1 - c010 * x1 * z1) / det;

		const double a3 =
			(c110 * x0 * y0 - c111 * x0 * y0 - c010 * x1 * y0 + c011 * x1 * y0 - c100 * x0 * y1 +
				c101 * x0 * y1 + c000 * x1 * y1 - c001 * x1 * y1) / det;

		const double a4 =
			(c001 * z0 - c011 * z0 - c101 * z0 + c111 * z0 - c000 * z1 + c010 * z1 + c100 * z1 -
				c110 * z1) / det;

		const double a5 =
			(c010 * y0 - c011 * y0 - c110 * y0 + c111 * y0 - c000 * y1 + c001 * y1 + c100 * y1 -
				c101 * y1) / det;

		const double a6 =
			(c100 * x0 - c101 * x0 - c110 * x0 + c111 * x0 - c000 * x1 + c001 * x1 + c010 * x1 -
				c011 * x1) / det;

		const double a7 =
			(c000 - c001 - c010 + c011 - c100 + c101 + c110 - c111) / det;

		return a0 + a1 * x + a2 * y + a3 * z + a4 * x * y + a5 * x * z + a6 * y * z + a7 * x * y * z;
	}

	/// \brief A wrapper for grid interpolation vector values
	struct SurroundingVectorCells
	{
		pmp::vec3 MinPt{};
		pmp::vec3 MaxPt{};
		// cell point values:
		pmp::dvec3 Val000{ 0.0 };
		pmp::dvec3 Val100{ 0.0 };
		pmp::dvec3 Val010{ 0.0 };
		pmp::dvec3 Val110{ 0.0 };
		pmp::dvec3 Val001{ 0.0 };
		pmp::dvec3 Val101{ 0.0 };
		pmp::dvec3 Val011{ 0.0 };
		pmp::dvec3 Val111{ 0.0 };
	};

	/**
	 * \brief Gets surrounding vector cell values of a sampled point from a vector grid.
	 * \param samplePt    point where the grid is sampled.
	 * \param grid        interpolated grid.
	 * \return surrounding vector cells.
	 *
	 * DISCLAIMER: For samplePt outside of grid.Box(), the values are clamped to boundary values.
	 */
	[[nodiscard]] SurroundingVectorCells GetSurroundingCells(const pmp::vec3& samplePt, const VectorGrid& grid)
	{
		const auto& [Nx, Ny, Nz] = grid.Dimensions();
		const auto& boxMin = grid.Box().min();
		const float cellSize = grid.CellSize();

		const auto ix = std::max(std::min(static_cast<size_t>(std::floor((samplePt[0] - boxMin[0]) / cellSize)), Nx - 1), ZERO_CELL_ID);
		const auto iy = std::max(std::min(static_cast<size_t>(std::floor((samplePt[1] - boxMin[1]) / cellSize)), Ny - 1), ZERO_CELL_ID);
		const auto iz = std::max(std::min(static_cast<size_t>(std::floor((samplePt[2] - boxMin[2]) / cellSize)), Nz - 1), ZERO_CELL_ID);

		const auto ix1 = std::min(ix + 1, Nx - 1);
		const auto iy1 = std::min(iy + 1, Ny - 1);
		const auto iz1 = std::min(iz + 1, Nz - 1);

		const auto i000 = Nx * Ny * iz + Nx * iy + ix;
		const auto i100 = Nx * Ny * iz + Nx * iy + ix1;
		const auto i010 = Nx * Ny * iz + Nx * iy1 + ix;
		const auto i110 = Nx * Ny * iz + Nx * iy1 + ix1;
		const auto i001 = Nx * Ny * iz1 + Nx * iy + ix;
		const auto i101 = Nx * Ny * iz1 + Nx * iy + ix1;
		const auto i011 = Nx * Ny * iz1 + Nx * iy1 + ix;
		const auto i111 = Nx * Ny * iz1 + Nx * iy1 + ix1;

		const auto& valuesX = grid.ValuesX();
		const auto& valuesY = grid.ValuesY();
		const auto& valuesZ = grid.ValuesZ();
		return {
			pmp::vec3{
				boxMin[0] + static_cast<float>(ix) * cellSize,
				boxMin[1] + static_cast<float>(iy) * cellSize,
				boxMin[2] + static_cast<float>(iz) * cellSize
			},
			pmp::vec3{
				boxMin[0] + static_cast<float>(ix + 1) * cellSize,
				boxMin[1] + static_cast<float>(iy + 1) * cellSize,
				boxMin[2] + static_cast<float>(iz + 1) * cellSize
			},
			pmp::dvec3{valuesX[i000], valuesY[i000], valuesZ[i000]},
			pmp::dvec3{valuesX[i100], valuesY[i100], valuesZ[i100]},
			pmp::dvec3{valuesX[i010], valuesY[i010], valuesZ[i010]},
			pmp::dvec3{valuesX[i110], valuesY[i110], valuesZ[i110]},

			pmp::dvec3{valuesX[i001], valuesY[i001], valuesZ[i001]},
			pmp::dvec3{valuesX[i101], valuesY[i101], valuesZ[i101]},
			pmp::dvec3{valuesX[i011], valuesY[i011], valuesZ[i011]},
			pmp::dvec3{valuesX[i111], valuesY[i111], valuesZ[i111]}
		};
	}

	pmp::dvec3 TrilinearInterpolateVectorValue(const pmp::vec3& samplePt, const VectorGrid& grid)
	{
		assert(grid.IsValid() && grid.CellSize() > 0.0f);
		const auto surrCellVals = GetSurroundingCells(samplePt, grid);
		const auto x = static_cast<double>(samplePt[0]), y = static_cast<double>(samplePt[1]), z = static_cast<double>(samplePt[2]);
		// cell min
		const auto x0 = static_cast<double>(surrCellVals.MinPt[0]), y0 = static_cast<double>(surrCellVals.MinPt[1]), z0 = static_cast<double>(surrCellVals.MinPt[2]);
		// cell max
		const auto x1 = static_cast<double>(surrCellVals.MaxPt[0]), y1 = static_cast<double>(surrCellVals.MaxPt[1]), z1 = static_cast<double>(surrCellVals.MaxPt[2]);

		const auto interpolateValueFromCoordId = [&](const unsigned int& i)
		{
			assert(i <= 2);
			// cell values
			const double c000 = surrCellVals.Val000[i], c100 = surrCellVals.Val100[i], c010 = surrCellVals.Val010[i], c110 = surrCellVals.Val110[i];
			const double c001 = surrCellVals.Val001[i], c101 = surrCellVals.Val101[i], c011 = surrCellVals.Val011[i], c111 = surrCellVals.Val111[i];

			const auto det = (x0 - x1) * (y0 - y1) * (z0 - z1);
			const double a0 =
				(c111 * x0 * y0 * z0 - c011 * x1 * y0 * z0 - c101 * x0 * y1 * z0 + c001 * x1 * y1 * z0 -
					c110 * x0 * y0 * z1 + c010 * x1 * y0 * z1 + c100 * x0 * y1 * z1 - c000 * x1 * y1 * z1) / det;

			const double a1 =
				(c011 * y0 * z0 - c111 * y0 * z0 - c001 * y1 * z0 + c101 * y1 * z0 - c010 * y0 * z1 +
					c110 * y0 * z1 + c000 * y1 * z1 - c100 * y1 * z1) / det;

			const double a2 =
				(c101 * x0 * z0 - c111 * x0 * z0 - c001 * x1 * z0 + c011 * x1 * z0 - c100 * x0 * z1 +
					c110 * x0 * z1 + c000 * x1 * z1 - c010 * x1 * z1) / det;

			const double a3 =
				(c110 * x0 * y0 - c111 * x0 * y0 - c010 * x1 * y0 + c011 * x1 * y0 - c100 * x0 * y1 +
					c101 * x0 * y1 + c000 * x1 * y1 - c001 * x1 * y1) / det;

			const double a4 =
				(c001 * z0 - c011 * z0 - c101 * z0 + c111 * z0 - c000 * z1 + c010 * z1 + c100 * z1 -
					c110 * z1) / det;

			const double a5 =
				(c010 * y0 - c011 * y0 - c110 * y0 + c111 * y0 - c000 * y1 + c001 * y1 + c100 * y1 -
					c101 * y1) / det;

			const double a6 =
				(c100 * x0 - c101 * x0 - c110 * x0 + c111 * x0 - c000 * x1 + c001 * x1 + c010 * x1 -
					c011 * x1) / det;

			const double a7 =
				(c000 - c001 - c010 + c011 - c100 + c101 + c110 - c111) / det;

			return a0 + a1 * x + a2 * y + a3 * z + a4 * x * y + a5 * x * z + a6 * y * z + a7 * x * y * z;
		};

		return pmp::dvec3(
			interpolateValueFromCoordId(0),
			interpolateValueFromCoordId(1),
			interpolateValueFromCoordId(2)
		);
	}

	// ------------------------------------------------------------------------------------------

	/// \brief A wrapper for grid interpolation scalar values
	struct SurroundingScalarCells2D
	{
		pmp::vec2 MinPt{};
		pmp::vec2 MaxPt{};
		double Val00{ 0.0 };
		double Val10{ 0.0 };
		double Val01{ 0.0 };
		double Val11{ 0.0 };
	};

	SurroundingScalarCells2D GetSurroundingCells2D(const pmp::vec2& samplePt, const ScalarGrid2D& grid)
	{
		const auto& [Nx, Ny] = grid.Dimensions();
		const auto& boxMin = grid.Box().min();
		const float cellSize = grid.CellSize();

		const auto ix = std::max(std::min(static_cast<size_t>(std::floor((samplePt[0] - boxMin[0]) / cellSize)), Nx - 1), ZERO_CELL_ID);
		const auto iy = std::max(std::min(static_cast<size_t>(std::floor((samplePt[1] - boxMin[1]) / cellSize)), Ny - 1), ZERO_CELL_ID);

		const auto ix1 = std::min(ix + 1, Nx - 1);
		const auto iy1 = std::min(iy + 1, Ny - 1);

		const auto i00 = Nx * iy + ix;
		const auto i10 = Nx * iy + ix1;
		const auto i01 = Nx * iy1 + ix;
		const auto i11 = Nx * iy1 + ix1;

		const auto& values = grid.Values();
		return {
			pmp::vec2{boxMin[0] + static_cast<float>(ix) * cellSize, boxMin[1] + static_cast<float>(iy) * cellSize},
			pmp::vec2{boxMin[0] + static_cast<float>(ix + 1) * cellSize, boxMin[1] + static_cast<float>(iy + 1) * cellSize},
			values[i00], values[i10], values[i01], values[i11]
		};
	}

	double BilinearInterpolateScalarValue(const pmp::vec2& samplePt, const ScalarGrid2D& grid)
	{
		assert(grid.IsValid() && grid.CellSize() > 0.0f);
		const auto surrCellVals = GetSurroundingCells2D(samplePt, grid);
		const auto x = static_cast<double>(samplePt[0]), y = static_cast<double>(samplePt[1]);

		const auto x0 = static_cast<double>(surrCellVals.MinPt[0]), y0 = static_cast<double>(surrCellVals.MinPt[1]);
		const auto x1 = static_cast<double>(surrCellVals.MaxPt[0]), y1 = static_cast<double>(surrCellVals.MaxPt[1]);

		const auto dx = (x - x0) / (x1 - x0);
		const auto dy = (y - y0) / (y1 - y0);

		return (1 - dx) * (1 - dy) * surrCellVals.Val00 +
			dx * (1 - dy) * surrCellVals.Val10 +
			(1 - dx) * dy * surrCellVals.Val01 +
			dx * dy * surrCellVals.Val11;
	}

	// --------------------------------------------------------------------------------------------

	/// \brief A wrapper for grid interpolation vector values
	struct SurroundingVectorCells2D
	{
		pmp::vec2 MinPt{};
		pmp::vec2 MaxPt{};
		pmp::dvec2 Val00{ 0.0 };
		pmp::dvec2 Val10{ 0.0 };
		pmp::dvec2 Val01{ 0.0 };
		pmp::dvec2 Val11{ 0.0 };
	};

	SurroundingVectorCells2D GetSurroundingCells2D(const pmp::vec2& samplePt, const VectorGrid2D& grid)
	{
		const auto& [Nx, Ny] = grid.Dimensions();
		const auto& boxMin = grid.Box().min();
		const float cellSize = grid.CellSize();

		const auto ix = std::max(std::min(static_cast<size_t>(std::floor((samplePt[0] - boxMin[0]) / cellSize)), Nx - 1), ZERO_CELL_ID);
		const auto iy = std::max(std::min(static_cast<size_t>(std::floor((samplePt[1] - boxMin[1]) / cellSize)), Ny - 1), ZERO_CELL_ID);

		const auto ix1 = std::min(ix + 1, Nx - 1);
		const auto iy1 = std::min(iy + 1, Ny - 1);

		const auto i00 = Nx * iy + ix;
		const auto i10 = Nx * iy + ix1;
		const auto i01 = Nx * iy1 + ix;
		const auto i11 = Nx * iy1 + ix1;

		const auto& valuesX = grid.ValuesX();
		const auto& valuesY = grid.ValuesY();
		return {
			pmp::vec2{boxMin[0] + static_cast<float>(ix) * cellSize, boxMin[1] + static_cast<float>(iy) * cellSize},
			pmp::vec2{boxMin[0] + static_cast<float>(ix + 1) * cellSize, boxMin[1] + static_cast<float>(iy + 1) * cellSize},
			pmp::dvec2{valuesX[i00], valuesY[i00]},
			pmp::dvec2{valuesX[i10], valuesY[i10]},
			pmp::dvec2{valuesX[i01], valuesY[i01]},
			pmp::dvec2{valuesX[i11], valuesY[i11]}
		};
	}

	pmp::dvec2 BilinearInterpolateVectorValue(const pmp::vec2& samplePt, const VectorGrid2D& grid)
	{
		assert(grid.IsValid() && grid.CellSize() > 0.0f);
		const auto surrCellVals = GetSurroundingCells2D(samplePt, grid);
		const auto x = static_cast<double>(samplePt[0]), y = static_cast<double>(samplePt[1]);

		const auto x0 = static_cast<double>(surrCellVals.MinPt[0]), y0 = static_cast<double>(surrCellVals.MinPt[1]);
		const auto x1 = static_cast<double>(surrCellVals.MaxPt[0]), y1 = static_cast<double>(surrCellVals.MaxPt[1]);

		const auto dx = (x - x0) / (x1 - x0);
		const auto dy = (y - y0) / (y1 - y0);

		return
			(1 - dx) * (1 - dy) * surrCellVals.Val00 +
			dx * (1 - dy) * surrCellVals.Val10 +
			(1 - dx) * dy * surrCellVals.Val01 +
			dx * dy * surrCellVals.Val11;
	}

	// --------------------------------------------------------------------------------------------

	void ComputeInteriorExteriorSignFromMeshNormals(ScalarGrid& grid, const pmp::SurfaceMesh& mesh)
	{
		if (!mesh.has_vertex_property("v:normal"))
		{
			std::cerr << "Geometry::ComputeInteriorExteriorSignFromMeshNormals: Input mesh has no normals!\n";
			return; // nothing to compute from
		}

		if (!mesh.is_triangle_mesh())
		{
			std::cerr << "Geometry::ComputeInteriorExteriorSignFromMeshNormals: Must be a triangle mesh! Triangulate before processing!\n";
			return; // must be a triangle mesh
		}

		const auto vNormalProp = mesh.get_vertex_property<pmp::Normal>("v:normal");
		const auto ptrMeshKDTree = std::make_unique<pmp::TriangleKdTree>(mesh, 0);

		auto& values = grid.Values();

		const auto& dim = grid.Dimensions();
		const auto& orig = grid.Box().min();
		const float cellSize = grid.CellSize();

		const auto Nx = static_cast<unsigned int>(dim.Nx);
		const auto Ny = static_cast<unsigned int>(dim.Ny);
		const auto Nz = static_cast<unsigned int>(dim.Nz);

		const unsigned int progressStep = Nz / 10;

		for (unsigned int iz = 0; iz < Nz; iz++)
		{
			// ----------------------------------
			if (iz % progressStep == 0)
			{
				const float progress = static_cast<float>(iz) / static_cast<float>(Nz);
				std::cout << "Geometry::ComputeInteriorExteriorSignFromMeshNormals: " << progress << " %\n";
			}
			// ----------------------------------

			for (unsigned int iy = 0; iy < Ny; iy++)
			{
				for (unsigned int ix = 0; ix < Nx; ix++)
				{
					const auto gridPt = pmp::Point{
						orig[0] + static_cast<float>(ix) * cellSize,
						orig[1] + static_cast<float>(iy) * cellSize,
						orig[2] + static_cast<float>(iz) * cellSize
					};

					const auto nearestNeighbor = ptrMeshKDTree->nearest(gridPt);

					const auto nearestPt = nearestNeighbor.nearest;
					const auto nearestFace = nearestNeighbor.face;

					std::vector<pmp::Vertex> vertices{};
					vertices.reserve(3);
					for (const auto& v : mesh.vertices(nearestFace))
					{
						vertices.push_back(v);
					}
					assert(vertices.size() == 3);

					const auto v0 = mesh.position(vertices[0]);
					const auto v1 = mesh.position(vertices[1]);
					const auto v2 = mesh.position(vertices[2]);

					const auto v0Normal = vNormalProp[vertices[0]];
					const auto v1Normal = vNormalProp[vertices[1]];
					const auto v2Normal = vNormalProp[vertices[2]];

					// get barycentric coordinates
					pmp::Point b = barycentric_coordinates(nearestPt, v0, v1, v2);

					// interpolate normal
					pmp::Point n;
					n = (v0Normal * b[0]);
					n += (v1Normal * b[1]);
					n += (v2Normal * b[2]);
					n.normalize();
					assert(!std::isnan(n[0]));
					const auto dotProd = static_cast<double>(dot(n, gridPt - nearestPt));

					const unsigned int gridPos = Nx * Ny * iz + Nx * iy + ix;
					values[gridPos] = (dotProd > 0.0 ? 1.0 : -1.0);
				}
			}
		}
	}

	/// \brief sign func.
	template <typename T> int sgn(T val) {
		return (T(0) < val) - (val < T(0));
	}

	void ComputeMeshSignedDistanceFromNormals(ScalarGrid& grid, const pmp::SurfaceMesh& mesh)
	{
		if (!mesh.has_vertex_property("v:normal"))
		{
			std::cerr << "Geometry::ComputeMeshSignedDistanceFromNormals: Input mesh has no normals!\n";
			return; // nothing to compute from
		}

		if (!mesh.is_triangle_mesh())
		{
			std::cerr << "Geometry::ComputeMeshSignedDistanceFromNormals: Must be a triangle mesh! Triangulate before processing!\n";
			return; // must be a triangle mesh
		}

		const auto vNormalProp = mesh.get_vertex_property<pmp::Normal>("v:normal");
		const auto ptrMeshKDTree = std::make_unique<pmp::TriangleKdTree>(mesh, 0);

		auto& values = grid.Values();

		const auto& dim = grid.Dimensions();
		const auto& orig = grid.Box().min();
		const float cellSize = grid.CellSize();

		const auto Nx = static_cast<unsigned int>(dim.Nx);
		const auto Ny = static_cast<unsigned int>(dim.Ny);
		const auto Nz = static_cast<unsigned int>(dim.Nz);

		pmp::Point gridPt{};		
		pmp::Point n; // interpolated normal
		pmp::Point b; // closest pt on triangle

		const unsigned int progressStep = Nz / 10;

		for (unsigned int iz = 0; iz < Nz; iz++)
		{
			// ----------------------------------
			if (iz % progressStep == 0)
			{
				const float progress = static_cast<float>(iz) / static_cast<float>(Nz);
				std::cout << "Geometry::ComputeMeshSignedDistanceFromNormals: " << progress << " %\n";
			}
			// ----------------------------------

			for (unsigned int iy = 0; iy < Ny; iy++)
			{
				for (unsigned int ix = 0; ix < Nx; ix++)
				{
					gridPt[0] = orig[0] + static_cast<float>(ix) * cellSize;
					gridPt[1] = orig[1] + static_cast<float>(iy) * cellSize;
					gridPt[2] = orig[2] + static_cast<float>(iz) * cellSize;

					const auto nearestNeighbor = ptrMeshKDTree->nearest(gridPt);

					const auto& nearestPt = nearestNeighbor.nearest;
					const auto& nearestFace = nearestNeighbor.face;

					std::vector<pmp::Vertex> vertices{};
					vertices.reserve(3);
					for (const auto& v : mesh.vertices(nearestFace))
					{
						vertices.push_back(v);
					}
					assert(vertices.size() == 3);

					const auto& v0 = mesh.position(vertices[0]);
					const auto& v1 = mesh.position(vertices[1]);
					const auto& v2 = mesh.position(vertices[2]);

					const auto& v0Normal = vNormalProp[vertices[0]];
					const auto& v1Normal = vNormalProp[vertices[1]];
					const auto& v2Normal = vNormalProp[vertices[2]];

					b = barycentric_coordinates(nearestPt, v0, v1, v2);

					n = (v0Normal * b[0]);
					n += (v1Normal * b[1]);
					n += (v2Normal * b[2]);
					n.normalize();
					assert(!std::isnan(n[0]));
					const auto vecToMesh = gridPt - nearestPt;
					const auto dotProd = static_cast<double>(dot(n, vecToMesh));

					const unsigned int gridPos = Nx * Ny * iz + Nx * iy + ix;
					values[gridPos] = sgn(dotProd) * norm(vecToMesh);
				}
			}
		}
	}

	double SimpleUnion(const double& f1Val, const double& f2Val)
	{
		return std::max(f1Val, f2Val);
	}

	double SimpleIntersection(const double& f1Val, const double& f2Val)
	{
		return std::min(f1Val, f2Val);
	}

	double SimpleDifference(const double& f1Val, const double& f2Val)
	{
		return std::max(f1Val - f2Val, 0.0);
	}

	double BlendedUnion(const double& f1Val, const double& f2Val)
	{
		return f1Val + f2Val + sqrt(f1Val * f1Val + f2Val * f2Val);
	}

	double BlendedIntersection(const double& f1Val, const double& f2Val)
	{
		return f1Val + f2Val - sqrt(f1Val * f1Val + f2Val * f2Val);
	}

	double BlendedDifference(const double& f1Val, const double& f2Val)
	{
		return f1Val - f2Val - sqrt(f1Val * f1Val + f2Val * f2Val);
	}

	double DistanceUnion(const double& f1Val, const double& f2Val)
	{
		return std::min(f1Val, f2Val);
	}

	// Zero of function g(r) = r^4 - r^2 + 0.25 in interval [0, 1], [Section III. from Ryan Geiss http://www.geisswerks.com/ryan/BLOBS/blobs.html]
	constexpr double DECAY_POLYNOMIAL_ZERO_LVL_SQUARED = 0.49;

	void ApplyMetaBallToGrid(ScalarGrid& grid, const MetaBallParams& params)
	{
		// parameter check
		if (params.Radius < FLT_EPSILON)
		{
			std::cerr << "ApplyMetaBallToGrid: Invalid parameter. params.Radius <= 0.0f!\n";
			return;
		}

		// grid
		auto& values = grid.Values();

		const auto& dim = grid.Dimensions();
		const auto& orig = grid.Box().min();
		const float cellSize = grid.CellSize();

		// metaball ROI (region of influence (box))
		const float radius = params.Radius;
		const auto& center = params.Center;
		const auto radiusVec = pmp::vec3(radius, radius, radius);
		const float radiusSq = radius * radius;
		const auto roi = pmp::BoundingBox(center - radiusVec, center + radiusVec);

		if (!grid.Box().Intersects(roi))
		{
			return; // nothing happens
		}

		const auto trueROI = grid.Box().Intersect(roi);

		const auto ixMin = static_cast<unsigned int>(std::floor((trueROI.min()[0] - orig[0]) / cellSize));
		const auto iyMin = static_cast<unsigned int>(std::floor((trueROI.min()[1] - orig[1]) / cellSize));
		const auto izMin = static_cast<unsigned int>(std::floor((trueROI.min()[2] - orig[2]) / cellSize));

		const auto ixMax = static_cast<unsigned int>(std::floor((trueROI.max()[0] - orig[0]) / cellSize));
		const auto iyMax = static_cast<unsigned int>(std::floor((trueROI.max()[1] - orig[1]) / cellSize));
		const auto izMax = static_cast<unsigned int>(std::floor((trueROI.max()[2] - orig[2]) / cellSize));
		assert(ixMax <= dim.Nx); assert(iyMax <= dim.Ny); assert(izMax <= dim.Nz);

		const auto Nx = static_cast<unsigned int>(dim.Nx);
		const auto Ny = static_cast<unsigned int>(dim.Ny);

		for (unsigned int iz = izMin; iz < izMax; iz++)
		{
			for (unsigned int iy = iyMin; iy < iyMax; iy++)
			{
				for (unsigned int ix = ixMin; ix < ixMax; ix++)
				{
					const auto gridPt = pmp::Point{
						orig[0] + static_cast<float>(ix) * cellSize,
						orig[1] + static_cast<float>(iy) * cellSize,
						orig[2] + static_cast<float>(iz) * cellSize
					};
					const auto posVect = gridPt - center;
					const auto distSq = static_cast<double>(dot(posVect, posVect) / radiusSq);
					// Source: [Section III. from Ryan Geiss http://www.geisswerks.com/ryan/BLOBS/blobs.html]
					const double val = (distSq < DECAY_POLYNOMIAL_ZERO_LVL_SQUARED ? (distSq * distSq - distSq + 0.25): 0.0);

					const unsigned int gridPos = Nx * Ny * iz + Nx * iy + ix;
					values[gridPos] = params.BoolOpFunction(values[gridPos], val);
				}
			}
		}
	}

	void ApplyCapsuleDistanceFieldToGrid(ScalarGrid& grid, const CapsuleParams& params)
	{
		// parameter check
		if (params.Height < FLT_EPSILON)
		{
			std::cerr << "ApplyCapsuleDistanceFieldToGrid: Invalid parameter. params.Height <= 0.0f!\n";
			return;
		}
		if (params.Radius < FLT_EPSILON)
		{
			std::cerr << "ApplyCapsuleDistanceFieldToGrid: Invalid parameter. params.Radius <= 0.0f!\n";
			return;
		}

		// grid
		auto& values = grid.Values();

		const auto& dim = grid.Dimensions();
		const auto& orig = grid.Box().min();
		const float cellSize = grid.CellSize();

		// capsule ROI (region of influence (box))
		const float radius = params.Radius;
		const float height = params.Height;
		const auto& position = params.Position;
		const auto radiusVec = pmp::vec3(radius, radius, radius);
		const auto heightVec = pmp::vec3(0.0f, 0.0f, height);
		const auto roi = pmp::BoundingBox(position - radiusVec, position + heightVec + radiusVec);

		if (!grid.Box().Intersects(roi))
		{
			return; // nothing happens
		}

		const auto Nx = static_cast<unsigned int>(dim.Nx);
		const auto Ny = static_cast<unsigned int>(dim.Ny);
		const auto Nz = static_cast<unsigned int>(dim.Nz);

		for (unsigned int iz = 0; iz < Nz; iz++)
		{
			for (unsigned int iy = 0; iy < Ny; iy++)
			{
				for (unsigned int ix = 0; ix < Nx; ix++)
				{
					const auto gridPt = pmp::Point{
						orig[0] + static_cast<float>(ix) * cellSize,
						orig[1] + static_cast<float>(iy) * cellSize,
						orig[2] + static_cast<float>(iz) * cellSize
					};
					const auto posVect = gridPt - position;
					const auto posVectClamped = pmp::vec3{ posVect[0], posVect[1], posVect[2] - std::clamp<float>(posVect[2], 0.0f, height + radius) };
					const auto dist = static_cast<double>(norm(posVectClamped)) - static_cast<double>(radius);
					const unsigned int gridPos = Nx * Ny * iz + Nx * iy + ix;
					values[gridPos] = params.BoolOpFunction(values[gridPos], dist);
				}
			}
		}
	}

	//
	// ========================================================================================
	//

	void ApplyTorusDistanceFieldToGrid(ScalarGrid& grid, const TorusParams& params)
	{
		// parameter check
		if (params.RingRadius < FLT_EPSILON) {
			std::cerr << "ApplyTorusDistanceFieldToGrid: Invalid parameter. params.RingRadius <= 0.0f!\n";
			return;
		}
		if (params.TubeRadius < FLT_EPSILON) {
			std::cerr << "ApplyTorusDistanceFieldToGrid: Invalid parameter. params.TubeRadius <= 0.0f!\n";
			return;
		}

		// grid
		auto& values = grid.Values();
		const auto& dim = grid.Dimensions();
		const auto& orig = grid.Box().min();
		const float cellSize = grid.CellSize();

		// torus ROI (region of influence (box))
		const float ringRadius = params.RingRadius;
		const float tubeRadius = params.TubeRadius;
		const auto& center = params.Center;

		const auto roi = pmp::BoundingBox{
			pmp::Point{center[0] - ringRadius - tubeRadius, center[1] - ringRadius - tubeRadius, center[2] - tubeRadius},
			pmp::Point{center[0] + ringRadius + tubeRadius, center[1] + ringRadius + tubeRadius, center[2] + tubeRadius}
		};

		if (!grid.Box().Intersects(roi)) {
			return; // nothing happens
		}

		const auto Nx = static_cast<unsigned int>(dim.Nx);
		const auto Ny = static_cast<unsigned int>(dim.Ny);
		const auto Nz = static_cast<unsigned int>(dim.Nz);

		for (unsigned int iz = 0; iz < Nz; iz++) {
			for (unsigned int iy = 0; iy < Ny; iy++) {
				for (unsigned int ix = 0; ix < Nx; ix++) {
					const auto gridPt = pmp::Point{
						orig[0] + static_cast<float>(ix) * cellSize,
						orig[1] + static_cast<float>(iy) * cellSize,
						orig[2] + static_cast<float>(iz) * cellSize
					};
					const auto posVect = gridPt - center;

					// Calculate distance from grid point to the torus
					const float len = sqrt(posVect[0] * posVect[0] + posVect[1] * posVect[1]) - ringRadius;
					const auto dist = static_cast<double>(sqrt(len * len + posVect[2] * posVect[2]) - tubeRadius);

					const unsigned int gridPos = Nx * Ny * iz + Nx * iy + ix;
					values[gridPos] = params.BoolOpFunction(values[gridPos], dist);
				}
			}
		}
	}

	ScalarGrid ExtractReSampledGrid(const float& newCellSize, const ScalarGrid& origGrid)
	{
		if (std::abs(newCellSize - origGrid.CellSize()) < FLT_EPSILON)
			return origGrid;

		ScalarGrid result(newCellSize, origGrid.Box());
		auto& values = result.Values();

		const auto& dim = result.Dimensions();
		const auto& newOrigin = result.Box().min();

		const auto Nx = static_cast<unsigned int>(dim.Nx);
		const auto Ny = static_cast<unsigned int>(dim.Ny);
		const auto Nz = static_cast<unsigned int>(dim.Nz);

		for (unsigned int iz = 0; iz < Nz; iz++)
		{
			for (unsigned int iy = 0; iy < Ny; iy++)
			{
				for (unsigned int ix = 0; ix < Nx; ix++)
				{
					const auto newGridPt = pmp::Point{
						newOrigin[0] + static_cast<float>(ix) * newCellSize,
						newOrigin[1] + static_cast<float>(iy) * newCellSize,
						newOrigin[2] + static_cast<float>(iz) * newCellSize
					};

					const unsigned int newGridPos = Nx * Ny * iz + Nx * iy + ix;
					values[newGridPos] = TrilinearInterpolateScalarValue(newGridPt, origGrid);
				}
			}
		}

		return result;
	}


	bool ContainsLocalMaximumNearScalarGridCell(const ScalarGrid2D& grid, unsigned int ix, unsigned int iy, unsigned int radius)
	{
		using namespace Eigen;

		const auto& values = grid.Values();
		const auto Nx = static_cast<unsigned int>(grid.Dimensions().Nx);
		const auto Ny = static_cast<unsigned int>(grid.Dimensions().Ny);
		assert(ix < Nx);
		assert(iy < Ny);

		// Collect the values in the 3x3 neighborhood
		std::array<pmp::Scalar, 9> neighborhood_values;
		int idx = 0;
		const auto r = static_cast<int>(radius);
		for (int di = -r; di <= r; di += r)
		{
			for (int dj = -r; dj <= r; dj += r)
			{
				const int adjusted_ix = std::clamp(static_cast<int>(ix) + di, 0, static_cast<int>(Nx) - 1);
				const int adjusted_iy = std::clamp(static_cast<int>(iy) + dj, 0, static_cast<int>(Ny) - 1);

				neighborhood_values[idx++] = values[Nx * adjusted_iy + adjusted_ix];
			}
		}

		// Fit a quadratic function f(x, y) = a00 + a10*x + a01*y + a20*x^2 + a11*x*y + a02*y^2
		Matrix<float, 9, 6> A;
		Vector<float, 9> b;

		// normalized cell coordinates
		A << 1, -1, -1, 1, 1, 1,
			1, 0, -1, 0, 0, 1,
			1, 1, -1, 1, -1, 1,
			1, -1, 0, 1, 0, 0,
			1, 0, 0, 0, 0, 0,
			1, 1, 0, 1, 0, 0,
			1, -1, 1, 1, -1, 1,
			1, 0, 1, 0, 0, 1,
			1, 1, 1, 1, 1, 1;

		for (int i = 0; i < 9; ++i)
		{
			b[i] = neighborhood_values[i];
		}

		const Vector<float, 6> coeffs = A.colPivHouseholderQr().solve(b);
		const float a10 = coeffs[1];
		const float a01 = coeffs[2];
		const float a20 = coeffs[3];
		const float a11 = coeffs[4];
		const float a02 = coeffs[5];

		// Solve for the critical point
		Matrix2f H;
		H << 2 * a20, a11, a11, 2 * a02;
		Vector2f grad;
		grad << a10, a01;

		SelfAdjointEigenSolver<Matrix2f> solver(H);
		if (solver.eigenvalues()(0) > FLT_EPSILON || solver.eigenvalues()(1) > FLT_EPSILON)
		{
			return false; // Hessian is not negative definite
		}

		Vector2f critical_point = -H.inverse() * grad;

		if (critical_point[0] < -1 || critical_point[0] > 1 || critical_point[1] < -1 || critical_point[1] > 1)
		{
			return false; // Critical point outside the 9 grid cells
		}

		return true;
	}

	bool ContainsLocalExtremesNearScalarGridCell(const ScalarGrid2D& grid, unsigned int ix, unsigned int iy, unsigned int radius)
	{
		using namespace Eigen;

		const auto& values = grid.Values();
		const auto Nx = static_cast<unsigned int>(grid.Dimensions().Nx);
		const auto Ny = static_cast<unsigned int>(grid.Dimensions().Ny);
		assert(ix < Nx);
		assert(iy < Ny);

		// Collect the values in the 3x3 neighborhood
		std::array<pmp::Scalar, 9> neighborhood_values;
		int idx = 0;
		const auto r = static_cast<int>(radius);
		for (int di = -r; di <= r; di += r)
		{
			for (int dj = -r; dj <= r; dj += r)
			{
				// Ensure indices are within grid bounds
				const int clamped_ix = std::clamp(static_cast<int>(ix) + di, 0, static_cast<int>(Nx) - 1);
				const int clamped_iy = std::clamp(static_cast<int>(iy) + dj, 0, static_cast<int>(Ny) - 1);
				neighborhood_values[idx++] = values[Nx * clamped_iy + clamped_ix];
			}
		}

		// Fit a quadratic function f(x, y) = a00 + a10*x + a01*y + a20*x^2 + a11*x*y + a02*y^2
		Matrix<float, 9, 6> A;
		Vector<float, 9> b;

		// normalized cell coordinates
		A << 1, -1, -1, 1, 1, 1,
			1, 0, -1, 0, 0, 1,
			1, 1, -1, 1, -1, 1,
			1, -1, 0, 1, 0, 0,
			1, 0, 0, 0, 0, 0,
			1, 1, 0, 1, 0, 0,
			1, -1, 1, 1, -1, 1,
			1, 0, 1, 0, 0, 1,
			1, 1, 1, 1, 1, 1;

		for (int i = 0; i < 9; ++i)
		{
			b[i] = neighborhood_values[i];
		}

		const Vector<float, 6> coeffs = A.colPivHouseholderQr().solve(b);
		const float a10 = coeffs[1];
		const float a01 = coeffs[2];
		const float a20 = coeffs[3];
		const float a11 = coeffs[4];
		const float a02 = coeffs[5];

		// Solve for the critical point
		Matrix2f H;
		H << 2 * a20, a11, a11, 2 * a02;
		Vector2f grad;
		grad << a10, a01;

		Vector2f critical_point = -H.inverse() * grad;

		// Check if the critical point lies within the [-1, 1] x [-1, 1] range of the neighborhood
		if (critical_point[0] >= -1 && critical_point[0] <= 1 && critical_point[1] >= -1 && critical_point[1] <= 1)
		{
			return true;  // A local extreme (max or min) is present
		}

		return false;
	}

	std::optional<pmp::Point2> FindLocalMaximumNearScalarGridCell(const ScalarGrid2D& grid, unsigned int ix, unsigned int iy, unsigned int radius)
	{
		using namespace Eigen;

		const auto& values = grid.Values();
		const auto Nx = static_cast<unsigned int>(grid.Dimensions().Nx);
		const auto Ny = static_cast<unsigned int>(grid.Dimensions().Ny);
		assert(ix < Nx);
		assert(iy < Ny);
		const auto cellSize = grid.CellSize();
		const auto& orig = grid.Box().min();

		// Collect the values in the 3x3 neighborhood
		std::array<pmp::Scalar, 9> neighborhood_values;
		int idx = 0;
		const auto r = static_cast<int>(radius);
		for (int di = -r; di <= r; di += r)
		{
			for (int dj = -r; dj <= r; dj += r)
			{
				const int adjusted_ix = std::clamp(static_cast<int>(ix) + di, 0, static_cast<int>(Nx) - 1);
				const int adjusted_iy = std::clamp(static_cast<int>(iy) + dj, 0, static_cast<int>(Ny) - 1);

				neighborhood_values[idx++] = values[Nx * adjusted_iy + adjusted_ix];
			}
		}

		// Fit a quadratic function f(x, y) = a00 + a10*x + a01*y + a20*x^2 + a11*x*y + a02*y^2
		Matrix<float, 9, 6> A;
		Vector<float, 9> b;

		// normalized cell coordinates
		A << 1, -1, -1, 1, 1, 1,
			1, 0, -1, 0, 0, 1,
			1, 1, -1, 1, -1, 1,
			1, -1, 0, 1, 0, 0,
			1, 0, 0, 0, 0, 0,
			1, 1, 0, 1, 0, 0,
			1, -1, 1, 1, -1, 1,
			1, 0, 1, 0, 0, 1,
			1, 1, 1, 1, 1, 1;

		for (int i = 0; i < 9; ++i)
		{
			b[i] = neighborhood_values[i];
		}

		const Vector<float, 6> coeffs = A.colPivHouseholderQr().solve(b);
		const float a10 = coeffs[1];
		const float a01 = coeffs[2];
		const float a20 = coeffs[3];
		const float a11 = coeffs[4];
		const float a02 = coeffs[5];

		// Solve for the critical point
		Matrix2f H;
		H << 2 * a20, a11, a11, 2 * a02;
		Vector2f grad;
		grad << a10, a01;

		SelfAdjointEigenSolver<Matrix2f> solver(H);
		if (solver.eigenvalues()(0) > FLT_EPSILON || solver.eigenvalues()(1) > FLT_EPSILON)
		{
			return {}; // Hessian is not negative definite
		}

		Vector2f critical_point = -H.inverse() * grad;

		if (critical_point[0] < -1 || critical_point[0] > 1 || critical_point[1] < -1 || critical_point[1] > 1)
		{
			return {}; // Critical point outside the 9 grid cells
		}

		return orig + pmp::Point2(
			(ix + radius * critical_point[0]) * cellSize,
			(iy + radius * critical_point[1]) * cellSize);
	}

	bool IsConvergentOrDivergentNearCell(const VectorGrid2D& vecGrid, unsigned int ix, unsigned int iy, unsigned int radius)
	{
		const auto Nx = static_cast<unsigned int>(vecGrid.Dimensions().Nx);
		const auto Ny = static_cast<unsigned int>(vecGrid.Dimensions().Ny);

		const auto& valuesX = vecGrid.ValuesX();
		const auto& valuesY = vecGrid.ValuesY();

		float divergenceSum = 0.0f;
		const int r = static_cast<int>(radius);

		// Iterate over the 4 key points in the neighborhood of the cell (ix, iy)
		const std::array<std::pair<int, int>, 4> keyPoints = { {
			{ -r, 0 }, // left
			{ r, 0 },  // right
			{ 0, -r }, // bottom
			{ 0, r }   // top
		} };

		for (const auto& [di, dj] : keyPoints)
		{
			// Ensure we don't go out of bounds
			const int ni = std::clamp(static_cast<int>(ix) + di, 1, static_cast<int>(Nx) - 2);
			const int nj = std::clamp(static_cast<int>(iy) + dj, 1, static_cast<int>(Ny) - 2);

			// Approximate divergence: calculating partial derivatives of vector components with respect to x and y
			float divX = 0.0f;
			float divY = 0.0f;

			// Compute divergence for X component (partial derivative with respect to x)
			if (ni > 0 && ni < Nx - 1)
			{
				divX = (valuesX[Nx * nj + ni + 1] - valuesX[Nx * nj + ni - 1]) / 2.0f;
			}

			// Compute divergence for Y component (partial derivative with respect to y)
			if (nj > 0 && nj < Ny - 1)
			{
				divY = (valuesY[Nx * (nj + 1) + ni] - valuesY[Nx * (nj - 1) + ni]) / 2.0f;
			}

			// Sum up the divergence for the selected points
			divergenceSum += divX + divY;
		}

		// Average divergence over the neighborhood
		const float avgDivergence = divergenceSum / 9.0f;

		// We consider this point convergent or divergent if the average divergence is significantly non-zero
		constexpr float divergenceThreshold = 1e-3f;
		return std::abs(avgDivergence) > divergenceThreshold;
	}

} // namespace Geometry
