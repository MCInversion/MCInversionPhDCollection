#pragma once

#include "geometry/Grid.h"
#include "pmp/SurfaceMesh.h"

namespace SDF
{
	/// \brief enumerator for grid preprocessing function type for distance field.
	enum class [[nodiscard]] PreprocessingType
	{
		Octree = 0,  //>! preprocess ScalarGrid using an OctreeVoxelizer.
		NoOctree = 1 //>! preprocess ScalarGrid without using an OctreeVoxelizer.
	};

	/// \brief enumerator for split function type for Kd tree.
	enum class [[nodiscard]] KDTreeSplitType
	{
		Center = 0, // simplest split function is chosen, evaluates the split position to box center
		Adaptive = 1 // a more robust split function, adaptively re-samples the kd-node's box according to triangle distribution within.
	};

	/// \brief enumerator for the approach to the computation of sign of the distance field to a mesh (negative inside, positive outside).
	enum class [[nodiscard]] SignComputation
	{
		None = 0, //>! no sign for the distance field is to be computed.
		VoxelFloodFill = 1, //>! negate and apply recursive flood-fill algorithm for non-frozen voxels.
		RayFromAHoleFilledMesh = 2 //>! compute distance field to a mesh after applying a pmp::HoleFilling, then all voxels whose outgoing ray intersects the mesh are interior.
	};

	/// \brief enumerator for which type of blur filter to apply after the distance field is completed.
	enum class [[nodiscard]] BlurPostprocessingType
	{
		None = 0, //>! no blur is to be applied.
		ThreeCubedVoxelAveraging = 1, //>! a 3x3x3 averaging kernel is applied.
		FiveCubedVoxelAveraging = 2,  //>! a 5x5x5 averaging kernel is applied.
		ThreeCubedVoxelGaussian = 3,  //>! a 3x3x3 Gaussian kernel is applied.
		FiveCubedVoxelGaussian = 4 //>! a 5x5x5 Gaussian kernel is applied.
	};

	/// \brief A wrapper for input settings for computing distance field.
	struct DistanceFieldSettings
	{
		float CellSize{ 1.0f }; //>! size of a single distance voxel.
		KDTreeSplitType KDTreeSplit{ KDTreeSplitType::Center }; //>! the choice of a split function for KD-tree.
		float VolumeExpansionFactor{ 1.0f }; //>! expansion factor (how many times the minimum dimension of mesh's bounding box) for the resulting scalar grid volume.
		double TruncationFactor{ 0.1 }; //>! factor by which the minimum half-dimension of mesh's bounding box gives rise to a truncation (cutoff) value for the distance field.
		SignComputation SignMethod{ SignComputation::None }; //>! method by which the sign of the distance field should be computed.
		BlurPostprocessingType BlurType{ BlurPostprocessingType::None }; //>! type of blur filter to be used for post-processing.
		PreprocessingType PreprocType{ PreprocessingType::Octree }; //>! function type for the preprocessing of distance field scalar grid.
	};

	/**
	 * \brief Computes the signed distance field of given input mesh.
	 * \param inputMesh               evaluated mesh.
	 * \param settings                settings for the distance field.
	 * \return the computed distance field's ScalarGrid.
	 */
	[[nodiscard]] Geometry::ScalarGrid ComputeDistanceField(
		const pmp::SurfaceMesh& inputMesh, const DistanceFieldSettings& settings);

} // namespace SDF
