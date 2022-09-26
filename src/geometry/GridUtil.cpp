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
					gradValsX[gradPos] = (gridValues[gridPosNextX] - gridValues[gridPosPrevX]) / (2.0 * cellSize);
					gradValsY[gradPos] = (gridValues[gridPosNextY] - gridValues[gridPosPrevY]) / (2.0 * cellSize);
					gradValsZ[gradPos] = (gridValues[gridPosNextZ] - gridValues[gridPosPrevZ]) / (2.0 * cellSize);
				}
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

					const double norm = sqrt(grad_x * grad_x + grad_y * grad_y + grad_z * grad_z);

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

					const double norm = -1.0 * sqrt(grad_x * grad_x + grad_y * grad_y + grad_z * grad_z);

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

	/// \brief zero index with respect to axis-aligned grid dimension.
	constexpr size_t ZERO_CELL_ID = 0;

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

} // namespace Geometry
