#pragma once

#include "geometry/Grid.h"
#include "pmp/SurfaceMesh.h"

namespace SDF
{
	/// \brief enumerator for split function type for Kd tree.
	enum class [[nodiscard]] KDTreeSplitType
	{
		Center = 0, // simplest split function is chosen, evaluates the split position to box center
		Adaptive = 1 // a more robust split function, adaptively re-samples the kd-node's box according to triangle distribution within.
	};

	/**
	 * \brief Computes the signed distance field of given input mesh.
	 * \param inputMesh               evaluated mesh.
	 * \param cellSize                size of a single distance voxel.
	 * \param splitType               the choice of a split function for KD-tree.
	 * \param volumeExpansion         expansion factor (how many times the minimum dimension of mesh's bounding box) for the resulting scalar grid volume.
	 * \param computeSign             if true sign of the distance field should be computed.
	 * \return the computed distance field's ScalarGrid.
	 */
	[[nodiscard]] Geometry::ScalarGrid ComputeDistanceField(
		const pmp::SurfaceMesh& inputMesh, const float& cellSize, const KDTreeSplitType& splitType = KDTreeSplitType::Center, const float& volumeExpansion = 1.0f, const bool& computeSign = true);
}
