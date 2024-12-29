#include "SDF.h"

#include "geometry/GeometryUtil.h"
#include "geometry/GridUtil.h"

#include "FastSweep.h"
#include "Voxelizer.h"
#include "pmp/algorithms/HoleFilling.h"

#include <stack>
#include <nmmintrin.h>

// TODO: remove
#include <fstream>

namespace SDF
{
	void DistanceFieldGenerator::PreprocessGridNoOctree(Geometry::ScalarGrid& grid)
	{
		assert(m_KdTree);
		auto& gridVals = grid.Values();
		auto& gridFrozenVals = grid.FrozenValues();
		const auto& vertexPositions = m_KdTree->VertexPositions();
		const auto& triangles = m_KdTree->TriVertexIds();

		const auto& dims = grid.Dimensions();
		const auto& orig = grid.Box().min();
		const pmp::Scalar cellSize = grid.CellSize();

		std::vector triangle{ pmp::vec3(), pmp::vec3(), pmp::vec3() };
		pmp::vec3 voxelCenter, voxelMin, voxelMax;

		for (unsigned int iz = 0; iz < dims.Nz; iz++)
		{
			for (unsigned int iy = 0; iy < dims.Ny; iy++)
			{
				for (unsigned int ix = 0; ix < dims.Nx; ix++)
				{
					voxelCenter[0] = orig[0] + (static_cast<pmp::Scalar>(ix) + 0.5) * cellSize;
					voxelCenter[1] = orig[1] + (static_cast<pmp::Scalar>(iy) + 0.5) * cellSize;
					voxelCenter[2] = orig[2] + (static_cast<pmp::Scalar>(iz) + 0.5) * cellSize;

					voxelMin[0] = voxelCenter[0] - 0.5 * cellSize;
					voxelMin[1] = voxelCenter[1] - 0.5 * cellSize;
					voxelMin[2] = voxelCenter[2] - 0.5 * cellSize;

					voxelMax[0] = voxelCenter[0] + 0.5 * cellSize;
					voxelMax[1] = voxelCenter[1] + 0.5 * cellSize;
					voxelMax[2] = voxelCenter[2] + 0.5 * cellSize;

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

	void DistanceFieldGenerator::PreprocessGridWithOctree(Geometry::ScalarGrid& grid)
	{
		assert(m_KdTree);
		auto& gridVals = grid.Values();
		auto& gridFrozenVals = grid.FrozenValues();
		const pmp::Scalar cellSize = grid.CellSize();
		const auto& gridBox = grid.Box();
		const auto octreeVox = OctreeVoxelizer(*m_KdTree, gridBox, cellSize);

		// extract leaf boxes and distance values from their centroids
		std::vector<pmp::BoundingBox*> boxBuffer{};
		std::vector<double> valueBuffer{};
		octreeVox.GetLeafBoxesAndValues(boxBuffer, valueBuffer);
		const size_t nOutlineVoxels = boxBuffer.size();
		const pmp::Scalar gBoxMinX = gridBox.min()[0];
		const pmp::Scalar gBoxMinY = gridBox.min()[1];
		const pmp::Scalar gBoxMinZ = gridBox.min()[2];

		const auto& dims = grid.Dimensions();
		const auto Nx = static_cast<unsigned int>(dims.Nx);
		const auto Ny = static_cast<unsigned int>(dims.Ny);
		unsigned int ix, iy, iz, gridPos;

		for (size_t i = 0; i < nOutlineVoxels; i++)
		{
			// transform from real space to grid index space
			ix = static_cast<unsigned int>(std::floor((0.5 * (boxBuffer[i]->min()[0] + boxBuffer[i]->max()[0]) - gBoxMinX) / cellSize));
			iy = static_cast<unsigned int>(std::floor((0.5 * (boxBuffer[i]->min()[1] + boxBuffer[i]->max()[1]) - gBoxMinY) / cellSize));
			iz = static_cast<unsigned int>(std::floor((0.5 * (boxBuffer[i]->min()[2] + boxBuffer[i]->max()[2]) - gBoxMinZ) / cellSize));

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
	[[nodiscard]] Geometry::SplitFunction GetSplitFunction(const KDTreeSplitType& splitType)
	{
		if (splitType == KDTreeSplitType::Center)
			return Geometry::CenterSplitFunction;

		return Geometry::AdaptiveSplitFunction;
	}

	///**
	// * \brief Provides a split functor according to the given setting.
	// * \param splitType      split type identifier.
	// * \return the split function identified by splitType.
	// */
	//[[nodiscard]] Geometry::SplitFunction2D GetSplitFunction2D(const KDTreeSplitType& splitType)
	//{
	//	if (splitType == KDTreeSplitType::Center)
	//		return Geometry::CenterSplitFunction2D;

	//	return Geometry::AdaptiveSplitFunction2D;
	//}

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
			return ComputeSignUsingFloodFill;

		if (signCompType == SignComputation::RayFromAHoleFilledMesh)
			return ComputeSignUsingRays;

		return {}; // empty sign function
	}

	/**
	 * \brief Fill all mesh holes.
	 * \param meshAdapter       adapter for the mesh to have its holes filled.
	 */
	void FillMeshHoles(Geometry::MeshAdapter& meshAdapter)
	{
		if (const auto pmpAdapter = dynamic_cast<Geometry::PMPSurfaceMeshAdapter*>(&meshAdapter))
		{
			pmp::HoleFilling hf(pmpAdapter->GetMesh());
			for (const auto& h : pmpAdapter->GetMesh().halfedges())
			{
				if (!pmpAdapter->GetMesh().is_boundary(h) || !pmpAdapter->GetMesh().is_manifold(pmpAdapter->GetMesh().to_vertex(h)))
					continue;
				hf.fill_hole(h);
			}
		}
		//throw std::runtime_error("SDF::FillMeshHoles: meshAdapter not supported!\n");
	}

#define REPORT_SDF_STEPS false // Note: may affect performance

	[[nodiscard]] std::string PrintKDTreeSplitType(const KDTreeSplitType& type)
	{
		if (type == KDTreeSplitType::Adaptive)
			return "KDTreeSplitType::Adaptive";

		return "KDTreeSplitType::Center";
	}

	[[nodiscard]] std::string PrintSignMethod(const SignComputation& type)
	{
		if (type == SignComputation::VoxelFloodFill)
			return "SignComputation::VoxelFloodFill";
		if (type == SignComputation::RayFromAHoleFilledMesh)
			return "SignComputation::RayFromAHoleFilledMesh";

		return "SignComputation::None";
	}

	[[nodiscard]] std::string PrintBlurType(const BlurPostprocessingType& type)
	{
		if (type == BlurPostprocessingType::ThreeCubedVoxelAveraging)
			return "BlurPostprocessingType::ThreeCubedVoxelAveraging";
		if (type == BlurPostprocessingType::FiveCubedVoxelAveraging)
			return "BlurPostprocessingType::FiveCubedVoxelAveraging";
		if (type == BlurPostprocessingType::ThreeCubedVoxelGaussian)
			return "BlurPostprocessingType::ThreeCubedVoxelGaussian";
		if (type == BlurPostprocessingType::FiveCubedVoxelGaussian)
			return "BlurPostprocessingType::FiveCubedVoxelGaussian";

		return "BlurPostprocessingType::None";
	}

	[[nodiscard]] std::string PrintPreprocessingType(const PreprocessingType& type)
	{
		if (type == PreprocessingType::Octree)
			return "PreprocessingType::Octree";

		return "PreprocessingType::NoOctree";
	}

	void ReportInput(const pmp::SurfaceMesh& inputMesh, const DistanceFieldSettings& settings, std::ostream& os)
	{
		const auto bbox = inputMesh.bounds();
		os << "======================================================================\n";
		os << "> > > > > Initiating DistanceFieldGenerator::Generate for: < < < < < <\n";
		os << "inputMesh: " << inputMesh.name() << "\n";
		os << "           NVertices: " << inputMesh.n_vertices() << ",\n";
		os << "           NFaces: " << inputMesh.n_faces() << ",\n";
		os << "           Min: {" << bbox.min()[0] << ", " << bbox.min()[1] << ", " << bbox.min()[2] << "},\n";
		os << "           Max: {" << bbox.max()[0] << ", " << bbox.max()[1] << ", " << bbox.max()[2] << "},\n";
		os << "           Size: {"
			<< (bbox.max()[0] - bbox.min()[0]) << ", "
			<< (bbox.max()[1] - bbox.min()[1]) << ", "
			<< (bbox.max()[2] - bbox.min()[2]) << "},\n";
		os << "           Center: {"
			<< 0.5 * (bbox.max()[0] + bbox.min()[0]) << ", "
		    << 0.5 * (bbox.max()[1] + bbox.min()[1]) << ", "
		    << 0.5 * (bbox.max()[2] + bbox.min()[2]) << "},\n";
		os << "----------------------------------------------------------------------\n";
		os << "CellSize: " << settings.CellSize << "\n";
		os << "VolumeExpansionFactor: " << settings.VolumeExpansionFactor << "\n";
		os << "TruncationFactor: " << (settings.TruncationFactor < DBL_MAX ? std::to_string(settings.TruncationFactor) : "DBL_MAX") << "\n";
		os << "......................................................................\n";
		os << "PreprocessingType: " << PrintPreprocessingType(settings.PreprocType) << "\n";
		os << "KDTreeSplit: " << PrintKDTreeSplitType(settings.KDTreeSplit) << "\n";
		os << "SignMethod: " << PrintSignMethod(settings.SignMethod) << "\n";
		os << "BlurType: " << PrintBlurType(settings.BlurType) << "\n";
		os << "----------------------------------------------------------------------\n";
	}

	void ReportOutput(const Geometry::ScalarGrid& grid, std::ostream& os)
	{
		os << "----------------------------------------------------------------------\n";
		os << "> > > > > > > > > > > > Exporting ScalarGrid < < < < < < < < < < < < <\n";
		const auto& dim = grid.Dimensions();
		os << "Dimensions: " << dim.Nx << " x " << dim.Ny << " x " << dim.Nz << "\n";
		os << "Cell Size: " << grid.CellSize() << "\n";
		double maxVal = -DBL_MAX, minVal = DBL_MAX;
		for (const auto& val : grid.Values())
		{
			if (val > maxVal) maxVal = val;
			if (val < minVal) minVal = val;
		}
		os << "Min value: " << minVal << ", Max value: " << maxVal << "\n";
		const auto& box = grid.Box();
		os << "Bounds:    Min: {" << box.min()[0] << ", " << box.min()[1] << ", " << box.min()[2] << "},\n";
		os << "           Max: {" << box.max()[0] << ", " << box.max()[1] << ", " << box.max()[2] << "},\n";
		os << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
	}

	//
	// ===============================================================================================
	//

	Geometry::ScalarGrid DistanceFieldGenerator::Generate(
		const Geometry::MeshAdapter& inputMesh, 
		const DistanceFieldSettings& settings,
		const std::optional<pmp::BoundingBox>& customFieldBox)
	{
		assert(settings.CellSize > 0.0);
		assert(settings.VolumeExpansionFactor >= 0.0);
		assert(settings.TruncationFactor > 0);

		m_Mesh = inputMesh.Clone();
		if (m_Mesh->IsEmpty())
		{
			throw std::invalid_argument("DistanceFieldGenerator::Generate: m_Mesh->IsEmpty()! inputMesh contains no spatial data!\n");
		}

		if (customFieldBox.has_value() && customFieldBox->is_empty())
		{
			throw std::invalid_argument("DistanceFieldGenerator::Generate: customFieldBox is not valid!\n");
		}

		if (customFieldBox.has_value())
		{
			for (const auto& p : m_Mesh->GetVertices())
			{
				if (!customFieldBox->Contains(p))
				{
					throw std::invalid_argument("DistanceFieldGenerator::Generate: customFieldBox does not contain all mesh vertices!\n");
				}
			}
		}

		if (settings.SignMethod != SignComputation::None)
		{
#if REPORT_SDF_STEPS
			std::cout << "FillMeshHoles ... ";
#endif
			FillMeshHoles(*m_Mesh); // make mesh watertight
#if REPORT_SDF_STEPS
			std::cout << "done\n";
#endif
		}

		auto sdfBBox = customFieldBox.value_or(m_Mesh->GetBounds());
		const auto size = sdfBBox.max() - sdfBBox.min();
		const pmp::Scalar minSize = std::min({ size[0], size[1], size[2] });

		if (settings.VolumeExpansionFactor > 0.0 && !customFieldBox.has_value())
		{
			const pmp::Scalar expansion = settings.VolumeExpansionFactor * minSize;
			sdfBBox.expand(expansion, expansion, expansion);
		}

		// percentage of the minimum half-size of the mesh's bounding box.
		const double truncationValue = (settings.TruncationFactor < Geometry::DEFAULT_SCALAR_GRID_INIT_VAL ? settings.TruncationFactor * static_cast<double>(minSize) : Geometry::DEFAULT_SCALAR_GRID_INIT_VAL);
		Geometry::ScalarGrid resultGrid(settings.CellSize, sdfBBox, truncationValue);
#if REPORT_SDF_STEPS
		std::cout << "truncationValue: " << truncationValue << "\n";
		std::cout << "CollisionKdTree ... ";
#endif
		m_KdTree = std::make_unique<Geometry::CollisionKdTree>(*m_Mesh, GetSplitFunction(settings.KDTreeSplit));
#if REPORT_SDF_STEPS
		std::cout << "done\n";
#endif
#if REPORT_SDF_STEPS
		std::cout << "preprocessGrid ... ";
#endif
		const auto preprocessGrid = GetPreprocessingFunction(settings.PreprocType);
		preprocessGrid(resultGrid);
#if REPORT_SDF_STEPS
		std::cout << "done\n";
#endif

		if (truncationValue > 0.0)
		{
#if REPORT_SDF_STEPS
			std::cout << "FastSweep ... ";
#endif
			constexpr SweepSolverSettings fsSettings{};
			FastSweep(resultGrid, fsSettings);
#if REPORT_SDF_STEPS
			std::cout << "done\n";
#endif
		}

		if (settings.SignMethod != SignComputation::None)
		{
#if REPORT_SDF_STEPS
			std::cout << "signFunction ... ";
#endif
			const auto signFunction = GetSignFunction(settings.SignMethod);
			signFunction(resultGrid);
#if REPORT_SDF_STEPS
			std::cout << "done\n";
#endif
		}

		// blur postprocessing
		if (settings.BlurType != BlurPostprocessingType::None)
		{
#if REPORT_SDF_STEPS
			std::cout << "blurFunction ... ";
#endif
			const auto blurFunction = GetBlurFunction(settings.BlurType);
			blurFunction(resultGrid);
#if REPORT_SDF_STEPS
			std::cout << "done\n";
#endif
		}
		return resultGrid;
	}

	//! \brief if true __m128i avx buffers will be used for iterating through the field.
#define USE_INTRINSICS false

	void DistanceFieldGenerator::ComputeSignUsingFloodFill(Geometry::ScalarGrid& grid)
	{
		const auto origFrozenFlags = grid.FrozenValues();
		Geometry::NegateGrid(grid);
		auto& gridVals = grid.Values();
		auto& gridFrozenVals = grid.FrozenValues();
		const auto& dim = grid.Dimensions();

		const auto Nx = static_cast<unsigned int>(dim.Nx);
		const auto Ny = static_cast<unsigned int>(dim.Ny);
		const int nx = static_cast<int>(dim.Nx) - 1;
		const int ny = static_cast<int>(dim.Ny) - 1;
		const int nz = static_cast<int>(dim.Nz) - 1;

		unsigned int gridPos;
		int iz = 0, iy = 0, ix = 0;

#if USE_INTRINSICS
		// TODO: some distance fields remain negated (when using intrinsics) even though they're closed.
		union { __m128i idsTriple; int ids[3]; };
		union { __m128i idsMask; int imask[3]; };
		std::stack<__m128i> stack = {};

		idsTriple = _mm_setr_epi32(ix, iy, iz, INT_MAX);

		// find the first unfrozen cell
		gridPos = 0;
		while (gridFrozenVals[gridPos]) 
		{
			idsMask = _mm_cmplt_epi32(idsTriple, _mm_setr_epi32(nx, ny, nz, INT_MAX));
			idsTriple = _mm_add_epi32(
				idsTriple, _mm_setr_epi32(
					(imask[0] > 0) * 1, (imask[1] > 0) * 1, (imask[2] > 0) * 1, 0)
			);
			ix = ids[0]; iy = ids[1]; iz = ids[2];
			gridPos = Nx * Ny * iz + Nx * iy + ix;
		}

		ids[0] = ix;
		ids[1] = iy;
		ids[2] = iz;
		stack.push(idsTriple);

		// a simple voxel flood
		while (!stack.empty()) 
		{
			idsTriple = stack.top();
			stack.pop();

			ix = ids[0]; iy = ids[1]; iz = ids[2];
			gridPos = Nx * Ny * iz + Nx * iy + ix;

			if (!gridFrozenVals[gridPos])
			{
				gridVals[gridPos] = -1.0 * gridVals[gridPos];
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

				idsMask = _mm_cmplt_epi32(idsTriple, _mm_setr_epi32(nx, ny, nz, INT_MAX)); // upper bounds

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
#else
		std::stack<std::tuple<int, int, int>> stack = {};
		std::tuple<int, int, int> idsTriple{};
		// find the first unfrozen cell
		gridPos = 0;
		while (gridFrozenVals[gridPos])
		{
			ix += (ix < nx ? 1 : 0);
			iy += (iy < ny ? 1 : 0);
			iz += (iz < nz ? 1 : 0);
			gridPos = Nx * Ny * iz + Nx * iy + ix;
		}

		stack.push({ ix, iy, iz });
		// a simple voxel flood
		while (!stack.empty())
		{
			idsTriple = stack.top();
			stack.pop();
			ix = std::get<0>(idsTriple);
			iy = std::get<1>(idsTriple);
			iz = std::get<2>(idsTriple);
			gridPos = Nx * Ny * iz + Nx * iy + ix;
			if (!gridFrozenVals[gridPos])
			{
				gridVals[gridPos] = -1.0 * gridVals[gridPos];
				gridFrozenVals[gridPos] = true; // freeze cell when done

				if (ix > 0) {
					stack.push({ ix - 1, iy, iz });
				}
				if (ix < nx) {
					stack.push({ ix + 1, iy, iz });
				}
				if (iy > 0) {
					stack.push({ ix, iy - 1, iz });
				}
				if (iy < ny) {
					stack.push({ ix, iy + 1, iz });
				}
				if (iz > 0) {
					stack.push({ ix, iy, iz - 1 });
				}
				if (iz < nz) {
					stack.push({ ix, iy, iz + 1 });
				}
			}
		}
#endif

		grid.FrozenValues() = origFrozenFlags;
	}

	void DistanceFieldGenerator::ComputeSignUsingRays(Geometry::ScalarGrid& grid)
	{
		const auto origFrozenFlags = grid.FrozenValues();
		//Geometry::NegateGrid(grid);
		const auto origMeshBox = m_Mesh->GetBounds();
		Geometry::NegateGridSubVolume(grid, origMeshBox);
		const auto& gridBox = grid.Box();

		// grid sub-bounds
		const auto dMin = origMeshBox.min() - gridBox.min();
		const auto dMax = origMeshBox.max() - gridBox.min();
		assert(dMin[0] >= 0.0 && dMin[1] >= 0.0 && dMin[2] >= 0.0);
		assert(dMax[0] >= 0.0 && dMax[1] >= 0.0 && dMax[2] >= 0.0);

		// compute sub-grid bounds
		const auto& cellSize = grid.CellSize();
		const unsigned int iXStart = std::floor(dMin[0] / cellSize);
		const unsigned int iYStart = std::floor(dMin[1] / cellSize);
		const unsigned int iZStart = std::floor(dMin[2] / cellSize);

		const unsigned int iXEnd = std::ceil(dMax[0] / cellSize);
		const unsigned int iYEnd = std::ceil(dMax[1] / cellSize);
		const unsigned int iZEnd = std::ceil(dMax[2] / cellSize);

		auto& gridVals = grid.Values();
		const auto& dims = grid.Dimensions();
		const auto& orig = grid.Box().min();

		pmp::vec3 gridPt;
		const pmp::vec3 rayXDir{ 1.0, 0.0, 0.0 };
		Geometry::Ray ray{ gridPt, rayXDir };

		for (unsigned int iz = iZStart; iz < iZEnd; iz++)
		{
			for (unsigned int iy = iYStart; iy < iYEnd; iy++)
			{
				for (unsigned int ix = iXStart; ix < iXEnd; ix++)
				{
					gridPt[0] = orig[0] + ix * cellSize;
					gridPt[1] = orig[1] + iy * cellSize;
					gridPt[2] = orig[2] + iz * cellSize;

					ray.StartPt = gridPt;
					ray.HitParam = FLT_MAX;
					if (m_KdTree->IsRayStartPointInsideTriangleMesh(ray))
						continue; // skip negated values of interior grid points

					// negate negated values for exterior grid points
					const unsigned int gridPos = dims.Nx * dims.Ny * iz + dims.Nx * iy + ix;
					gridVals[gridPos] *= -1.0;
				}
			}
		}
		grid.FrozenValues() = origFrozenFlags;
	}

	Geometry::ScalarGrid PointCloudDistanceFieldGenerator::Generate(const std::vector<pmp::vec3>& inputPoints, const PointCloudDistanceFieldSettings& settings)
	{
		assert(settings.CellSize > 0.0);
		assert(settings.VolumeExpansionFactor >= 0.0);
		assert(settings.TruncationFactor > 0);
		if (inputPoints.empty())
		{
			throw std::invalid_argument("PointCloudDistanceFieldGenerator::Generate: inputPoints.empty()! inputPoints contains no spatial data!\n");
		}

		pmp::BoundingBox dfBBox(inputPoints);
		const auto size = dfBBox.max() - dfBBox.min();
		const pmp::Scalar minSize = std::min({ size[0], size[1], size[2] });

		if (settings.VolumeExpansionFactor > 0.0)
		{
			const pmp::Scalar expansion = settings.VolumeExpansionFactor * minSize;
			dfBBox.expand(expansion, expansion, expansion);
		}

		// percentage of the minimum half-size of the mesh's bounding box.
		const double truncationValue = (settings.TruncationFactor < Geometry::DEFAULT_SCALAR_GRID_INIT_VAL ? settings.TruncationFactor * static_cast<double>(minSize) : Geometry::DEFAULT_SCALAR_GRID_INIT_VAL);
		Geometry::ScalarGrid resultGrid(settings.CellSize, dfBBox, truncationValue);

		m_Points = inputPoints;
		PreprocessGridFromPoints(resultGrid);

		if (truncationValue > 0.0)
		{
#if REPORT_SDF_STEPS
			std::cout << "FastSweep ... ";
#endif
			constexpr SweepSolverSettings fsSettings{};
			FastSweep(resultGrid, fsSettings);
#if REPORT_SDF_STEPS
			std::cout << "done\n";
#endif
		}

		// blur postprocessing
		if (settings.BlurType != BlurPostprocessingType::None)
		{
#if REPORT_SDF_STEPS
			std::cout << "blurFunction ... ";
#endif
			const auto blurFunction = GetBlurFunction(settings.BlurType);
			blurFunction(resultGrid);
#if REPORT_SDF_STEPS
			std::cout << "done\n";
#endif
		}
		return resultGrid;
	}

	void PointCloudDistanceFieldGenerator::PreprocessGridFromPoints(Geometry::ScalarGrid& grid)
	{
		if (m_Points.empty())
		{
			std::cerr << "PointCloudDistanceFieldGenerator::PreprocessGridFromPoints: m_Points.empty()!\n";
			return;
		}
		auto& gridVals = grid.Values();
		auto& gridFrozenVals = grid.FrozenValues();
		const pmp::Scalar cellSize = grid.CellSize();

		const auto& gridBox = grid.Box();
		const pmp::Scalar gBoxMinX = gridBox.min()[0];
		const pmp::Scalar gBoxMinY = gridBox.min()[1];
		const pmp::Scalar gBoxMinZ = gridBox.min()[2];

		const auto& dims = grid.Dimensions();
		const auto Nx = static_cast<unsigned int>(dims.Nx);
		const auto Ny = static_cast<unsigned int>(dims.Ny);
		unsigned int ix, iy, iz, gridPos;
		pmp::vec3 gridPt;

		for (const auto& p : m_Points)
		{
			// transform from real space to grid index space
			ix = static_cast<unsigned int>(std::floor((p[0] - gBoxMinX) / cellSize));
			iy = static_cast<unsigned int>(std::floor((p[1] - gBoxMinY) / cellSize));
			iz = static_cast<unsigned int>(std::floor((p[2] - gBoxMinZ) / cellSize));

			gridPt[0] = gBoxMinX + ix * cellSize;
			gridPt[1] = gBoxMinY + iy * cellSize;
			gridPt[2] = gBoxMinZ + iz * cellSize;

			gridPos = Nx * Ny * iz + Nx * iy + ix;
			assert(gridPos < gridVals.size());
			gridVals[gridPos] = norm(gridPt - p);
			gridFrozenVals[gridPos] = true; // freeze initial condition for FastSweep
		}
	}

	Geometry::ScalarGrid2D PlanarDistanceFieldGenerator::Generate(
		const Geometry::CurveAdapter& inputCurve, 
		const DistanceField2DSettings& settings, 
		const std::optional<pmp::BoundingBox2>& customFieldBox)
	{
		assert(settings.CellSize > 0.0);
		assert(settings.AreaExpansionFactor >= 0.0);
		assert(settings.TruncationFactor > 0);

		m_Curve = inputCurve.Clone();
		if (m_Curve->IsEmpty())
		{
			throw std::invalid_argument("PlanarDistanceFieldGenerator::Generate: m_Curve->IsEmpty()! inputCurve contains no spatial data!\n");
		}

		if (customFieldBox.has_value() && customFieldBox->is_empty())
		{
			throw std::invalid_argument("PlanarDistanceFieldGenerator::Generate: customFieldBox is not valid!\n");
		}

		if (customFieldBox.has_value())
		{
			for (const auto& p : m_Curve->GetVertices())
			{
				if (!customFieldBox->Contains(p))
				{
					throw std::invalid_argument("PlanarDistanceFieldGenerator::Generate: customFieldBox does not contain all curve vertices!\n");
				}
			}
		}

//		if (settings.SignMethod != SignComputation2D::None)
//		{
//#if REPORT_SDF_STEPS
//			std::cout << "FillCurveHoles ... ";
//#endif
//			FillCurveHoles(*m_Curve); // make the curve watertight
//#if REPORT_SDF_STEPS
//			std::cout << "done\n";
//#endif
//		}

		auto sdfBBox = customFieldBox.value_or(m_Curve->GetBounds());
		const auto size = sdfBBox.max() - sdfBBox.min();
		const pmp::Scalar minSize = std::min( size[0], size[1] );

		if (settings.AreaExpansionFactor > 0.0 && !customFieldBox.has_value())
		{
			const pmp::Scalar expansion = settings.AreaExpansionFactor * minSize;
			sdfBBox.expand(expansion, expansion);
		}

		// percentage of the minimum half-size of the mesh's bounding box.
		const double truncationValue = (settings.TruncationFactor < Geometry::DEFAULT_SCALAR_GRID_INIT_VAL ? settings.TruncationFactor * static_cast<double>(minSize) : Geometry::DEFAULT_SCALAR_GRID_INIT_VAL);
		Geometry::ScalarGrid2D resultGrid(settings.CellSize, sdfBBox, truncationValue);
#if REPORT_SDF_STEPS
		std::cout << "truncationValue: " << truncationValue << "\n";
		std::cout << "CollisionKdTree ... ";
#endif
		m_KdTree = std::make_unique<Geometry::Collision2DTree>(*m_Curve, Geometry::CenterSplitFunction2D);
#if REPORT_SDF_STEPS
		std::cout << "done\n";
#endif
#if REPORT_SDF_STEPS
		std::cout << "preprocessGrid ... ";
#endif
		const auto preprocessGrid = GetPreprocessingFunction(settings.PreprocType);
		preprocessGrid(resultGrid);
#if REPORT_SDF_STEPS
		std::cout << "done\n";
#endif

		if (truncationValue > 0.0)
		{
#if REPORT_SDF_STEPS
			std::cout << "FastSweep ... ";
#endif
			constexpr SweepSolver2DSettings fsSettings{};
			FastSweep2D(resultGrid, fsSettings);
#if REPORT_SDF_STEPS
			std::cout << "done\n";
#endif
		}

		if (settings.SignMethod != SignComputation2D::None)
		{
#if REPORT_SDF_STEPS
			std::cout << "signFunction ... ";
#endif
			const auto signFunction = GetSignFunction(settings.SignMethod);
			signFunction(resultGrid);
#if REPORT_SDF_STEPS
			std::cout << "done\n";
#endif
		}

		// blur postprocessing
//		if (settings.BlurType != BlurPostprocessingType::None)
//		{
//#if REPORT_SDF_STEPS
//			std::cout << "blurFunction ... ";
//#endif
//			const auto blurFunction = GetBlurFunction(settings.BlurType);
//			blurFunction(resultGrid);
//#if REPORT_SDF_STEPS
//			std::cout << "done\n";
//#endif
//		}
		return resultGrid;
	}

	SignFunction2D PlanarDistanceFieldGenerator::GetSignFunction(const SignComputation2D& signCompType)
	{
		if (signCompType == SignComputation2D::PixelFloodFill)
			return ComputeSignUsingFloodFill;

		return {}; // empty sign function
	}

	PreprocessingFunction2D PlanarDistanceFieldGenerator::GetPreprocessingFunction(const PreprocessingType2D& preprocType)
	{
		if (preprocType == PreprocessingType2D::Quadtree)
			return PreprocessGridWithQuadtree;

		return PreprocessGridNoQuadtree;
	}

	void PlanarDistanceFieldGenerator::PreprocessGridNoQuadtree(Geometry::ScalarGrid2D& grid)
	{
		assert(m_KdTree);
		auto& gridVals = grid.Values();
		auto& gridFrozenVals = grid.FrozenValues();
		const auto& vertexPositions = m_KdTree->VertexPositions();
		const auto& edges = m_KdTree->EdgeVertexIds();

		const auto& dims = grid.Dimensions();
		const auto& orig = grid.Box().min();
		const pmp::Scalar cellSize = grid.CellSize();

		std::vector edge{ pmp::vec2(), pmp::vec2() };
		pmp::vec2 pixelCenter, pixelMin, pixelMax;

		for (unsigned int iy = 0; iy < dims.Ny; iy++)
		{
			for (unsigned int ix = 0; ix < dims.Nx; ix++)
			{
				pixelCenter[0] = orig[0] + (static_cast<pmp::Scalar>(ix) + 0.5) * cellSize;
				pixelCenter[1] = orig[1] + (static_cast<pmp::Scalar>(iy) + 0.5) * cellSize;

				pixelMin[0] = pixelCenter[0] - 0.5 * cellSize;
				pixelMin[1] = pixelCenter[1] - 0.5 * cellSize;

				pixelMax[0] = pixelCenter[0] + 0.5 * cellSize;
				pixelMax[1] = pixelCenter[1] + 0.5 * cellSize;

				const pmp::BoundingBox2 pixelBox{ pixelMin , pixelMax };
				std::vector<unsigned int> pixelEdgeIds{};
				m_KdTree->GetEdgesInABox(pixelBox, pixelEdgeIds);

				if (pixelEdgeIds.empty())
					continue; // no edges found

				double distToEdgeSq = DBL_MAX;

				for (const auto& eId : pixelEdgeIds)
				{
					edge[0] = vertexPositions[edges[eId].v0Id];
					edge[1] = vertexPositions[edges[eId].v1Id];

					const double currentEdgeDistSq = Geometry::GetDistanceToLine2DSq(edge, pixelCenter);

					if (currentEdgeDistSq < distToEdgeSq)
						distToEdgeSq = currentEdgeDistSq;
				}

				assert(distToEdgeSq < DBL_MAX);

				const unsigned int gridPos = dims.Nx * iy + ix;
				gridVals[gridPos] = sqrt(distToEdgeSq);
				gridFrozenVals[gridPos] = true;
			}
		}
	}

	void PlanarDistanceFieldGenerator::PreprocessGridWithQuadtree(Geometry::ScalarGrid2D& grid)
	{
		assert(m_KdTree);
		auto& gridVals = grid.Values();
		auto& gridFrozenVals = grid.FrozenValues();
		const pmp::Scalar cellSize = grid.CellSize();
		auto gridBox = grid.Box();
		const auto gridBoxSize = gridBox.max() - gridBox.min();
		gridBox.expand(Geometry::BOX_INFLATION * gridBoxSize[0], Geometry::BOX_INFLATION * gridBoxSize[1]);
		const auto quadtreeVox = QuadtreeVoxelizer(*m_KdTree, gridBox, cellSize);

		// extract leaf boxes and distance values from their centroids
		std::vector<pmp::BoundingBox2*> boxBuffer{};
		std::vector<double> valueBuffer{};
		quadtreeVox.GetLeafBoxesAndValues(boxBuffer, valueBuffer);
		const size_t nOutlinePixels = boxBuffer.size();
		const pmp::Scalar gBoxMinX = gridBox.min()[0];
		const pmp::Scalar gBoxMinY = gridBox.min()[1];

		const auto& dims = grid.Dimensions();
		const auto Nx = static_cast<unsigned int>(dims.Nx);
		unsigned int ix, iy, gridPos;

		//std::ofstream file("C:\\Users\\Martin\\source\\repos\\MCInversionPhDCollection\\output\\sdf_tests\\boxCoords.txt");
		//file << "outlinePts = [\n";

		for (size_t i = 0; i < nOutlinePixels; i++)
		{
			// transform from real space to grid index space
			//ix = static_cast<unsigned int>(std::floor((0.5 * (boxBuffer[i]->min()[0] + boxBuffer[i]->max()[0]) - gBoxMinX) / cellSize));
			//iy = static_cast<unsigned int>(std::floor((0.5 * (boxBuffer[i]->min()[1] + boxBuffer[i]->max()[1]) - gBoxMinY) / cellSize));

			//const auto outlinePt = boxBuffer[i]->min();
			//const auto outlinePt = boxBuffer[i]->center();
			//file << "\t[" << outlinePt[0] << ", " << outlinePt[1] << "],\n";

			const double param = 0.5;
			ix = static_cast<unsigned int>(std::floor((((1.0 - param) * boxBuffer[i]->min()[0] + param * boxBuffer[i]->max()[0]) - gBoxMinX) / cellSize));
			iy = static_cast<unsigned int>(std::floor((((1.0 - param) * boxBuffer[i]->min()[1] + param * boxBuffer[i]->max()[1]) - gBoxMinY) / cellSize));

			gridPos = Nx * iy + ix;
			gridVals[gridPos] = valueBuffer[i];
			gridFrozenVals[gridPos] = true; // freeze initial condition for FastSweep
		}

		//file << "]";
		//file.close();
	}

	void PlanarDistanceFieldGenerator::ComputeSignUsingFloodFill(Geometry::ScalarGrid2D& grid)
	{
		const auto origFrozenFlags = grid.FrozenValues();
		Geometry::NegateGrid(grid);
		auto& gridVals = grid.Values();
		auto& gridFrozenVals = grid.FrozenValues();
		const auto& dim = grid.Dimensions();

		const auto Nx = static_cast<unsigned int>(dim.Nx);
		const int nx = static_cast<int>(dim.Nx) - 1;
		const int ny = static_cast<int>(dim.Ny) - 1;

		unsigned int gridPos;
		int iy = 0, ix = 0;

		std::stack<std::tuple<int, int>> stack = {};
		std::tuple<int, int> idsPair;
		// find the first unfrozen cell
		gridPos = 0;
		while (gridFrozenVals[gridPos])
		{
			ix += (ix < nx ? 1 : 0);
			iy += (iy < ny ? 1 : 0);
			gridPos = Nx * iy + ix;
		}
		stack.push({ ix, iy });
		// a simple pixel flood
		while (!stack.empty())
		{
			idsPair = stack.top();
			stack.pop();
			ix = std::get<0>(idsPair);
			iy = std::get<1>(idsPair);
			gridPos = Nx * iy + ix;
			if (!gridFrozenVals[gridPos]) 
			{
				gridVals[gridPos] = -1.0 * gridVals[gridPos];
				gridFrozenVals[gridPos] = true; // freeze cell when done
				if (ix > 0) 
				{
					stack.push({ ix - 1, iy });
				}
				if (ix < nx)
				{
					stack.push({ ix + 1, iy });
				}
				if (iy > 0) 
				{
					stack.push({ ix, iy - 1 });
				}
				if (iy < ny)
				{
					stack.push({ ix, iy + 1 });
				}
			}
		}

		grid.FrozenValues() = origFrozenFlags;
	}

	Geometry::ScalarGrid2D PlanarPointCloudDistanceFieldGenerator::Generate(const std::vector<pmp::Point2>& inputPoints, const PointCloudDistanceField2DSettings& settings)
	{
		assert(settings.CellSize > 0.0);
		assert(settings.AreaExpansionFactor >= 0.0);
		assert(settings.TruncationFactor > 0);
		if (inputPoints.empty())
		{
			throw std::invalid_argument("PlanarPointCloudDistanceFieldGenerator::Generate: inputPoints.empty()! inputPoints contains no spatial data!\n");
		}

		pmp::BoundingBox2 dfBBox(inputPoints);
		const auto size = dfBBox.max() - dfBBox.min();
		const pmp::Scalar minSize = std::min(size[0], size[1]);

		if (settings.AreaExpansionFactor > 0.0)
		{
			const pmp::Scalar expansion = settings.AreaExpansionFactor * minSize;
			dfBBox.expand(expansion, expansion);
		}

		// shift by half-pixel
		dfBBox.min() += pmp::vec2{ settings.CellSize, settings.CellSize } * (pmp::Scalar)0.5;
		dfBBox.max() += pmp::vec2{ settings.CellSize, settings.CellSize } * (pmp::Scalar)0.5;

		// percentage of the minimum half-size of the mesh's bounding box.
		const double truncationValue = (settings.TruncationFactor < Geometry::DEFAULT_SCALAR_GRID_INIT_VAL ? settings.TruncationFactor * static_cast<double>(minSize) : Geometry::DEFAULT_SCALAR_GRID_INIT_VAL);
		Geometry::ScalarGrid2D resultGrid(settings.CellSize, dfBBox, truncationValue);

		m_Points = inputPoints;
		PreprocessGridFromPoints(resultGrid);

		if (truncationValue > 0.0)
		{
#if REPORT_SDF_STEPS
			std::cout << "FastSweep ... ";
#endif
			constexpr SweepSolver2DSettings fsSettings{};
			FastSweep2D(resultGrid, fsSettings);
#if REPORT_SDF_STEPS
			std::cout << "done\n";
#endif
		}

//		// blur postprocessing
//		if (settings.BlurType != BlurPostprocessingType::None)
//		{
//#if REPORT_SDF_STEPS
//			std::cout << "blurFunction ... ";
//#endif
//			const auto blurFunction = GetBlurFunction(settings.BlurType);
//			blurFunction(resultGrid);
//#if REPORT_SDF_STEPS
//			std::cout << "done\n";
//#endif
//		}
		return resultGrid;
	}

	void PlanarPointCloudDistanceFieldGenerator::PreprocessGridFromPoints(Geometry::ScalarGrid2D& grid)
	{
		if (m_Points.empty())
		{
			std::cerr << "PlanarPointCloudDistanceFieldGenerator::PreprocessGridFromPoints: m_Points.empty()!\n";
			return;
		}
		auto& gridVals = grid.Values();
		auto& gridFrozenVals = grid.FrozenValues();
		const pmp::Scalar cellSize = grid.CellSize();

		const auto& gridBox = grid.Box();
		const pmp::Scalar gBoxMinX = gridBox.min()[0];
		const pmp::Scalar gBoxMinY = gridBox.min()[1];

		const auto& dims = grid.Dimensions();
		const auto Nx = static_cast<unsigned int>(dims.Nx);
		unsigned int ix, iy, gridPos;
		pmp::vec2 gridPt;

		for (const auto& p : m_Points)
		{
			// transform from real space to grid index space
			ix = static_cast<unsigned int>(std::floor((p[0] - gBoxMinX) / cellSize));
			iy = static_cast<unsigned int>(std::floor((p[1] - gBoxMinY) / cellSize));

			gridPt[0] = gBoxMinX + ix * cellSize;
			gridPt[1] = gBoxMinY + iy * cellSize;

			gridPos = Nx * iy + ix;
			assert(gridPos < gridVals.size());
			gridVals[gridPos] = norm(gridPt - p);
			gridFrozenVals[gridPos] = true; // freeze initial condition for FastSweep
		}
	}

} // namespace SDF