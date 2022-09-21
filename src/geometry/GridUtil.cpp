#include "GridUtil.h"

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
					const double x0 = static_cast<double>(2 * i - 1) * 0.5 - radius;
					const double x1 = static_cast<double>(2 * i + 1) * 0.5 - radius;

					const double y0 = static_cast<double>(2 * j - 1) * 0.5 - radius;
					const double y1 = static_cast<double>(2 * j + 1) * 0.5 - radius;

					const double z0 = static_cast<double>(2 * k - 1) * 0.5 - radius;
					const double z1 = static_cast<double>(2 * k + 1) * 0.5 - radius;

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

} // namespace Geometry
