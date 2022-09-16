#include "SDF.h"

#include "CollisionKdTree.h"
#include "FastSweep.h"
#include "OctreeVoxelizer.h"
#include "geometry/GeometryUtil.h"

namespace SDF
{
	/**
	 * \brief A preprocessing approach for distance grid using CollisionKdTree to create "voxel outline" of inputMesh.
	 * \param grid         modifiable input grid.
	 * \param inputMesh    input mesh.
	 * \param spltFunc     split function used in CollisionKdTree.
	 */
	void PreprocessGridNoOctree(Geometry::ScalarGrid& grid, const pmp::SurfaceMesh& inputMesh, const SplitFunction& spltFunc)
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
		pmp::vec3 voxelCenter, voxelMin, voxelMax;

		for (unsigned int iz = 0; iz < dims.Nz; iz++)
		{
			for (unsigned int iy = 0; iy < dims.Ny; iy++)
			{
				for (unsigned int ix = 0; ix < dims.Nx; ix++)
				{
					voxelCenter[0] = orig[0] + (static_cast<float>(ix) + 0.5f) * cellSize;
					voxelCenter[1] = orig[1] + (static_cast<float>(iy) + 0.5f) * cellSize;
					voxelCenter[2] = orig[2] + (static_cast<float>(iz) + 0.5f) * cellSize;

					voxelMin[0] = voxelCenter[0] - 0.5f * cellSize;
					voxelMin[1] = voxelCenter[1] - 0.5f * cellSize;
					voxelMin[2] = voxelCenter[2] - 0.5f * cellSize;

					voxelMax[0] = voxelCenter[0] + 0.5f * cellSize;
					voxelMax[1] = voxelCenter[1] + 0.5f * cellSize;
					voxelMax[2] = voxelCenter[2] + 0.5f * cellSize;

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
	}

	/**
	 * \brief A preprocessing approach for distance grid using CollisionKdTree and OctreeVoxelizer to create "voxel outline" of inputMesh.
	 * \param grid         modifiable input grid.
	 * \param inputMesh    input mesh.
	 * \param spltFunc     split function used in CollisionKdTree.
	 */
	void PreprocessGridWithOctree(Geometry::ScalarGrid& grid, const pmp::SurfaceMesh& inputMesh, const SplitFunction& spltFunc)
	{
		auto& gridVals = grid.Values();
		auto& gridFrozenVals = grid.FrozenValues();

		const auto kdTree = CollisionKdTree(inputMesh, spltFunc);
		const float cellSize = grid.CellSize();
		const auto& gridBox = grid.Box();
		const auto octreeVox = OctreeVoxelizer(kdTree, gridBox, cellSize);

		// extract leaf boxes and distance values from their centroids
		std::vector<pmp::BoundingBox*> boxBuffer{};
		std::vector<double> valueBuffer{};
		octreeVox.GetLeafBoxesAndValues(boxBuffer, valueBuffer);
		const size_t nOutlineVoxels = boxBuffer.size();
		const float gBoxMinX = gridBox.min()[0]; const float gBoxMinY = gridBox.min()[1];	const float gBoxMinZ = gridBox.min()[2];

		const auto& dims = grid.Dimensions();
		const unsigned int Nx = dims.Nx; const unsigned int Ny = dims.Ny;

		unsigned int ix, iy, iz, gridPos;

		for (size_t i = 0; i < nOutlineVoxels; i++)
		{
			// transform from real space to grid index space
			ix = static_cast<unsigned int>(std::floor((0.5f * (boxBuffer[i]->min()[0] + boxBuffer[i]->max()[0]) - gBoxMinX) / cellSize));
			iy = static_cast<unsigned int>(std::floor((0.5f * (boxBuffer[i]->min()[1] + boxBuffer[i]->max()[1]) - gBoxMinY) / cellSize));
			iz = static_cast<unsigned int>(std::floor((0.5f * (boxBuffer[i]->min()[2] + boxBuffer[i]->max()[2]) - gBoxMinZ) / cellSize));

			gridPos = Nx * Ny * iz + Nx * iy + ix;
			gridVals[gridPos] = valueBuffer[i];
			gridFrozenVals[gridPos] = true; // freeze initial condition for FastSweep
		}
	}

	using PreprocessingFunction = std::function<void(Geometry::ScalarGrid&, const pmp::SurfaceMesh&, const SplitFunction&)>;

	/**
	 * \brief Provides a preprocessing functor according to the given setting.
	 * \param preprocType      preprocessing type identifier.
	 * \return the preprocessing function identified by preprocType.
	 */
	[[nodiscard]] PreprocessingFunction GetPreprocessingFunction(const PreprocessingType& preprocType)
	{
		if (preprocType == PreprocessingType::Octree)
			return PreprocessGridWithOctree;

		return PreprocessGridNoOctree;
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
		const pmp::SurfaceMesh& inputMesh, const DistanceFieldSettings& settings)
	{
		assert(settings.CellSize > 0.0f);
		assert(settings.VolumeExpansionFactor >= 0.0f);

		auto sdfBBox = inputMesh.bounds();
		const auto size = sdfBBox.max() - sdfBBox.min();
		const float minSize = std::min({ size[0], size[1], size[2] });

		if (settings.VolumeExpansionFactor > 0.0f)
		{
			const float expansion = settings.VolumeExpansionFactor * minSize;
			sdfBBox.expand(expansion, expansion, expansion);
		}

		// percentage of the minimum half-size of the mesh's bounding box.
		const double truncationValue = (settings.TruncationFactor < DBL_MAX ? settings.TruncationFactor * (static_cast<double>(minSize) / 2.0) : DBL_MAX);
		Geometry::ScalarGrid resultGrid(settings.CellSize, sdfBBox, truncationValue);

		const auto preprocessGrid = GetPreprocessingFunction(settings.PreprocType);
		preprocessGrid(resultGrid, inputMesh, GetSplitFunction(settings.KDTreeSplit));

		if (truncationValue > 0.0)
		{
			constexpr SweepSolverSettings fsSettings{};
			FastSweep(resultGrid, fsSettings);			
		}

		return resultGrid;
	}
} // namespace SDF