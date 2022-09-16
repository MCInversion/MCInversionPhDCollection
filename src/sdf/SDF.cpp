#include "SDF.h"


#include "geometry/GeometryUtil.h"
#include "geometry/GridUtil.h"

#include "CollisionKdTree.h"
#include "FastSweep.h"
#include "OctreeVoxelizer.h"
#include "pmp/algorithms/HoleFilling.h"

#include <stack>
#include <nmmintrin.h>

namespace SDF
{
	void DistanceFieldGenerator::PreprocessGridNoOctree(Geometry::ScalarGrid& grid, const pmp::SurfaceMesh& inputMesh, const SplitFunction& spltFunc)
	{
		assert(m_KdTree);
		auto& gridVals = grid.Values();
		auto& gridFrozenVals = grid.FrozenValues();
		const auto& vertexPositions = m_KdTree->VertexPositions();
		const auto& triangles = m_KdTree->TriVertexIds();

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
					m_KdTree->GetTrianglesInABox(voxelBox, voxelTriangleIds);

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

	void DistanceFieldGenerator::PreprocessGridWithOctree(Geometry::ScalarGrid& grid, const pmp::SurfaceMesh& inputMesh, const SplitFunction& spltFunc)
	{
		assert(m_KdTree);
		auto& gridVals = grid.Values();
		auto& gridFrozenVals = grid.FrozenValues();
		const float cellSize = grid.CellSize();
		const auto& gridBox = grid.Box();
		const auto octreeVox = OctreeVoxelizer(*m_KdTree, gridBox, cellSize);

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

	PreprocessingFunction DistanceFieldGenerator::GetPreprocessingFunction(const PreprocessingType& preprocType)
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

	/// \brief a blur function for postprocessing the distance field.
	using BlurFunction = std::function<void(Geometry::ScalarGrid&)>;

	/**
	 * \brief Provides a blur functor according to the given setting.
	 * \param blurType      blur type identifier.
	 * \return the blur function identified by splitType.
	 */
	[[nodiscard]] BlurFunction GetBlurFunction(const BlurPostprocessingType& blurType)
	{
		if (blurType == BlurPostprocessingType::ThreeCubedVoxelAveraging)
			return Geometry::ApplyNarrowAveragingBlur;

		if (blurType == BlurPostprocessingType::FiveCubedVoxelAveraging)
			return Geometry::ApplyWideAveragingBlur;

		if (blurType == BlurPostprocessingType::ThreeCubedVoxelGaussian)
			return Geometry::ApplyNarrowGaussianBlur;

		if (blurType == BlurPostprocessingType::FiveCubedVoxelGaussian)
			return Geometry::ApplyWideGaussianBlur;

		return {}; // empty blur function
	}

	SignFunction DistanceFieldGenerator::GetSignFunction(const SignComputation& signCompType)
	{
		if (signCompType == SignComputation::VoxelFloodFill)
			return ComputeSignUsingFloodFill; // [&](auto& g) { ComputeSignUsingFloodFill(g); };

		if (signCompType == SignComputation::RayFromAHoleFilledMesh)
			return ComputeSignUsingRays; //[&](auto& g) { ComputeSignUsingRays(g); };

		return {}; // empty sign function
	}

	/**
	 * \brief Fill all mesh holes.
	 * \param mesh        mesh to have its holes filled.
	 */
	void FillMeshHoles(pmp::SurfaceMesh& mesh)
	{
		pmp::HoleFilling hf(mesh);
		for (const auto& h : mesh.halfedges())
		{
			if (!mesh.is_boundary(h) || !mesh.is_manifold(mesh.to_vertex(h)))
				continue;
			hf.fill_hole(h);
		}
	}

	//
	// ===============================================================================================
	//

	Geometry::ScalarGrid DistanceFieldGenerator::Generate(const pmp::SurfaceMesh& inputMesh, const DistanceFieldSettings& settings)
	{
		assert(settings.CellSize > 0.0f);
		assert(settings.VolumeExpansionFactor >= 0.0f);

		m_Mesh = inputMesh;
		if (settings.SignMethod != SignComputation::None)
		{
			FillMeshHoles(m_Mesh); // make mesh watertight
		}

		auto sdfBBox = m_Mesh.bounds();
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

		m_KdTree = std::make_unique<CollisionKdTree>(m_Mesh, GetSplitFunction(settings.KDTreeSplit));
		const auto preprocessGrid = GetPreprocessingFunction(settings.PreprocType);
		preprocessGrid(resultGrid, m_Mesh, GetSplitFunction(settings.KDTreeSplit));

		if (truncationValue > 0.0)
		{
			constexpr SweepSolverSettings fsSettings{};
			FastSweep(resultGrid, fsSettings);
		}

		if (settings.SignMethod != SignComputation::None)
		{
			const auto signFunction = GetSignFunction(settings.SignMethod);
			signFunction(resultGrid);
		}

		// blur postprocessing
		if (settings.BlurType != BlurPostprocessingType::None)
		{
			const auto blurFunction = GetBlurFunction(settings.BlurType);
			blurFunction(resultGrid);
		}
		return resultGrid;
	}

	void DistanceFieldGenerator::ComputeSignUsingFloodFill(Geometry::ScalarGrid& grid)
	{
		const auto origFrozenFlags = grid.FrozenValues();
		Geometry::NegateGrid(grid);
		auto& gridVals = grid.Values();
		auto& gridFrozenVals = grid.FrozenValues();
		const auto& dim = grid.Dimensions();
		unsigned int nx = dim.Nx - 1;
		unsigned int ny = dim.Ny - 1;
		unsigned int nz = dim.Nz - 1;

		double val; unsigned int gridPos;
		unsigned int iz = 0, iy = 0, ix = 0;

		union { __m128i idsTriple; unsigned int ids[3]; };
		union { __m128i idsMask; unsigned int imask[3]; };
		std::stack<__m128i> stack = {};

		idsTriple = _mm_setr_epi32(ix, iy, iz, INT_MAX);

		// find the first unfrozen cell
		gridPos = 0;
		while (gridFrozenVals[gridPos]) {
			idsMask = _mm_cmplt_epi32(idsTriple, _mm_setr_epi32(nx, ny, nz, INT_MAX));
			idsTriple = _mm_add_epi32(
				idsTriple, _mm_setr_epi32(
					(imask[0] > 0) * 1, (imask[1] > 0) * 1, (imask[2] > 0) * 1, 0)
			);
			ix = ids[0]; iy = ids[1]; iz = ids[2];
			gridPos = dim.Nx * dim.Ny * iz + dim.Nx * iy + ix;
		}

		ids[0] = ix; ids[1] = iy; ids[2] = iz;
		stack.push(idsTriple);

		// a simple voxel flood
		while (!stack.empty()) {
			idsTriple = stack.top();
			stack.pop();

			ix = ids[0]; iy = ids[1]; iz = ids[2];
			gridPos = dim.Nx * dim.Ny * iz + dim.Nx * iy + ix;

			if (!gridFrozenVals[gridPos]) {
				val = -1.0 * gridVals[gridPos];
				gridVals[gridPos] = val;
				gridFrozenVals[gridPos] = true; // freeze cell when done

				idsMask = _mm_cmpgt_epi32(idsTriple, _mm_set1_epi32(0)); // lower bounds

				if (imask[0] > 0) {
					stack.push(_mm_setr_epi32(ix - 1, iy, iz, INT_MAX));
				}
				if (imask[1] > 0) {
					stack.push(_mm_setr_epi32(ix, iy - 1, iz, INT_MAX));
				}
				if (imask[2] > 0) {
					stack.push(_mm_setr_epi32(ix, iy, iz - 1, INT_MAX));
				}

				idsMask = _mm_cmplt_epi32(idsTriple, _mm_setr_epi32(nx, dim.Ny, nz, INT_MAX)); // upper bounds

				if (imask[0] > 0) {
					stack.push(_mm_setr_epi32(ix + 1, iy, iz, INT_MAX));
				}
				if (imask[1] > 0) {
					stack.push(_mm_setr_epi32(ix, iy + 1, iz, INT_MAX));
				}
				if (imask[2] > 0) {
					stack.push(_mm_setr_epi32(ix, iy, iz + 1, INT_MAX));
				}
			}
		}
		grid.FrozenValues() = origFrozenFlags;
	}

	void DistanceFieldGenerator::ComputeSignUsingRays(Geometry::ScalarGrid& grid)
	{
		auto& gridVals = grid.Values();
		const auto& dims = grid.Dimensions();
		const auto& orig = grid.Box().min();
		const float cellSize = grid.CellSize();

		Ray ray{};
		std::vector triangle{ pmp::vec3(), pmp::vec3(), pmp::vec3() };
		pmp::vec3 gridPt;

		for (unsigned int iz = 0; iz < dims.Nz; iz++)
		{
			for (unsigned int iy = 0; iy < dims.Ny; iy++)
			{
				for (unsigned int ix = 0; ix < dims.Nx; ix++) 
				{
					gridPt[0] = orig[0] + ix * cellSize;
					gridPt[1] = orig[1] + iy * cellSize;
					gridPt[2] = orig[2] + iz * cellSize;

					ray.StartPt = gridPt;
					if (!m_KdTree->RayIntersectsATriangle(ray))
						continue;

					const unsigned int gridPos = dims.Nx * dims.Ny * iz + dims.Nx * iy + ix;
					gridVals[gridPos] *= -1.0;
				}
			}
		}
	}

} // namespace SDF