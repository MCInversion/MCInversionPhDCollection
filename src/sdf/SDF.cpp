#include "SDF.h"

#include "CollisionKdTree.h"
#include "OctreeVoxelizer.h"
#include "geometry/GeometryUtil.h"

namespace SDF
{
	// this one's painfully slow!
	/*void PreprocessGrid(Geometry::ScalarGrid& grid, const pmp::SurfaceMesh& inputMesh, const SplitFunction& spltFunc)
	{
		auto& gridVals = grid.Values();
		auto& gridFrozenVals = grid.FrozenValues();

		const auto kdTree = CollisionKdTree(inputMesh, spltFunc);
		const auto& vertexPositions = kdTree.VertexPositions();
		const auto& triangles = kdTree.TriVertexIds();

		const auto& dims = grid.Dimensions();
		const auto& orig = grid.Box().min();
		const float cellSize = grid.CellSize();

		std::vector triangle{ pmp::vec3(), pmp::vec3(), pmp::vec3() };

		for (unsigned int iz = 0; iz < dims.Nz; iz++)
		{
			for (unsigned int iy = 0; iy < dims.Ny; iy++)
			{
				for (unsigned int ix = 0; ix < dims.Nx; ix++)
				{
					const pmp::vec3 voxelCenter{
						orig[0] + static_cast<float>(ix) * cellSize,
						orig[1] + static_cast<float>(iy) * cellSize,
						orig[2] + static_cast<float>(iz) * cellSize
					};

					const pmp::vec3 voxelMin{
						voxelCenter[0] - 0.5f * cellSize,
						voxelCenter[1] - 0.5f * cellSize,
						voxelCenter[2] - 0.5f * cellSize
					};

					const pmp::vec3 voxelMax{
						voxelCenter[0] + 0.5f * cellSize,
						voxelCenter[1] + 0.5f * cellSize,
						voxelCenter[2] + 0.5f * cellSize
					};

					const pmp::BoundingBox voxelBox{ voxelMin , voxelMax };
					std::vector<unsigned int> voxelTriangleIds{};
					kdTree.GetTrianglesInABox(voxelBox, voxelTriangleIds);

					if (voxelTriangleIds.empty())
						continue; // no triangles found

					double distToTriSq = DBL_MAX;

					for (const auto& triId : voxelTriangleIds)
					{
						triangle[0] = vertexPositions[triangles[triId].v0Id];
						triangle[1] = vertexPositions[triangles[triId].v1Id];
						triangle[2] = vertexPositions[triangles[triId].v2Id];

						const double currentTriDistSq = Geometry::GetDistanceToTriangleSq(triangle, voxelCenter);

						if (currentTriDistSq < distToTriSq)
							distToTriSq = currentTriDistSq;
					}

					assert(distToTriSq < DBL_MAX);

					const unsigned int gridPos = dims.Nx * dims.Ny * iz + dims.Nx * iy + ix;
					gridVals[gridPos] = sqrt(distToTriSq);
					gridFrozenVals[gridPos] = true;
				}
			}
		}
	}*/

	void PreprocessGrid(Geometry::ScalarGrid& grid, const pmp::SurfaceMesh& inputMesh, const SplitFunction& spltFunc)
	{
		auto& gridVals = grid.Values();
		auto& gridFrozenVals = grid.FrozenValues();

		const auto kdTree = CollisionKdTree(inputMesh, spltFunc);
		const float cellSize = grid.CellSize();
		const auto octreeVox = OctreeVoxelizer(kdTree, grid.Box(), cellSize);

		// extract leaf boxes and distance values from their centroids
		std::vector<pmp::BoundingBox*> boxBuffer{};
		std::vector<double> valueBuffer{};
		octreeVox.GetLeafBoxesAndValues(boxBuffer, valueBuffer);
		const size_t nOutlineVoxels = boxBuffer.size();
		const auto& oBox = octreeVox.GetRootBox();
		const float oBoxMinX = oBox.min()[0]; const float oBoxMinY = oBox.min()[1];	const float oBoxMinZ = oBox.min()[2];

		const auto& dims = grid.Dimensions();
		const size_t Nx = dims.Nx; const size_t Ny = dims.Ny; const size_t Nz = dims.Nz;

		unsigned int ix, iy, iz, gridPos;

		for (size_t i = 0; i < nOutlineVoxels; i++)
		{
			// transform from real space to grid index space
			ix = static_cast<unsigned int>(std::floor(0.5f * (boxBuffer[i]->min()[0] + boxBuffer[i]->max()[0]) - oBoxMinX) / cellSize);
			iy = static_cast<unsigned int>(std::floor(0.5f * (boxBuffer[i]->min()[1] + boxBuffer[i]->max()[1]) - oBoxMinY) / cellSize);
			iz = static_cast<unsigned int>(std::floor(0.5f * (boxBuffer[i]->min()[2] + boxBuffer[i]->max()[2]) - oBoxMinZ) / cellSize);

			gridPos = Nx * Ny * iz + Nx * iy + ix;
			gridVals[gridPos] = valueBuffer[i];
			gridFrozenVals[gridPos] = true; // freeze initial condition for FastSweep
		}
	}

	/**
	 * \brief Provides a split functor according to the given setting.
	 * \param splitType      split type identifier.
	 * \return the split function identified by splitType.
	 */
	[[nodiscard]] SplitFunction GetSplitFunction(const KDTreeSplitType& splitType)
	{
		if (splitType == KDTreeSplitType::Center)
			return CenterSplitFunction;

		return AdaptiveSplitFunction;
	}

	//
	// ===============================================================================================
	//


	Geometry::ScalarGrid ComputeDistanceField(
		const pmp::SurfaceMesh& inputMesh, const float& cellSize, const KDTreeSplitType& splitType,
		const float& volumeExpansion, const bool& computeSign)
	{
		assert(cellSize > 0.0f);
		assert(volumeExpansion >= 0.0f);

		auto sdfBBox = inputMesh.bounds();
		if (volumeExpansion > 0.0f)
		{
			const auto size = sdfBBox.max() - sdfBBox.min();
			const float minSize = std::min({ size[0], size[1], size[2] });
			const float expansion = volumeExpansion * minSize;
			sdfBBox.expand(expansion, expansion, expansion);
		}

		Geometry::ScalarGrid resultGrid(cellSize, sdfBBox);
		PreprocessGrid(resultGrid, inputMesh, GetSplitFunction(splitType));

		return resultGrid;
	}
} // namespace SDF