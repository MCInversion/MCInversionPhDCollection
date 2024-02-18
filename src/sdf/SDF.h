#pragma once

#include "geometry/CollisionKdTree.h"
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
		float VolumeExpansionFactor{ 1.0f }; //>! expansion factor (how many times the minimum dimension of mesh's bounding box) for the resulting scalar grid volume.
		double TruncationFactor{ 0.1 }; //>! factor by which the minimum half-dimension of mesh's bounding box gives rise to a truncation (cutoff) value for the distance field.
		KDTreeSplitType KDTreeSplit{ KDTreeSplitType::Center }; //>! the choice of a split function for KD-tree.
		SignComputation SignMethod{ SignComputation::None }; //>! method by which the sign of the distance field should be computed.
		BlurPostprocessingType BlurType{ BlurPostprocessingType::None }; //>! type of blur filter to be used for post-processing.
		PreprocessingType PreprocType{ PreprocessingType::Octree }; //>! function type for the preprocessing of distance field scalar grid.
	};

	/// \brief a functor for computing the sign of the distance field.
	using SignFunction = std::function<void(Geometry::ScalarGrid&)> const;

	/// \brief a functor for preprocessing the scalar grid for distance field.
	using PreprocessingFunction = std::function<void(Geometry::ScalarGrid&)>;

	/// \brief A singleton object for computing distance fields to triangle meshes.
	class DistanceFieldGenerator
	{
	public:
		/// \brief a deleted default constructor because m_KdTree is not default-constructable.
		DistanceFieldGenerator() = delete;

		/**
		 * \brief Compute the signed distance field of given input mesh.
		 * \param inputMesh               evaluated mesh.
		 * \param settings                settings for the distance field.
		 * \return the computed distance field's ScalarGrid.
		 */
		static [[nodiscard]] Geometry::ScalarGrid Generate(const pmp::SurfaceMesh& inputMesh, const DistanceFieldSettings& settings);
		
	private:
		inline static pmp::SurfaceMesh m_Mesh; //>! mesh to be (pre)processed.
		inline static std::unique_ptr<Geometry::CollisionKdTree> m_KdTree{ nullptr }; //>! mesh kd tree.

		/**
		 * \brief provides the SignFunction, a function from this generator's private interface that computes the sign of the distance field.
		 * \param signCompType      sign computation function type identifier.
		 * \return the sign function identified by signCompType.
		 */
		static [[nodiscard]] SignFunction GetSignFunction(const SignComputation& signCompType);

		/**
		 * \brief Provides a preprocessing functor according to the given setting.
		 * \param preprocType      preprocessing type identifier.
		 * \return the preprocessing function identified by preprocType.
		 */
		static [[nodiscard]] PreprocessingFunction GetPreprocessingFunction(const PreprocessingType& preprocType);

		/**
		 * \brief A preprocessing approach for distance grid using CollisionKdTree to create "voxel outline" of inputMesh.
		 * \param grid         modifiable input grid.
		 */
		static void PreprocessGridNoOctree(Geometry::ScalarGrid& grid);

		/**
		 * \brief A preprocessing approach for distance grid using CollisionKdTree and OctreeVoxelizer to create "voxel outline" of inputMesh.
		 * \param grid         modifiable input grid.
		 */
		static void PreprocessGridWithOctree(Geometry::ScalarGrid& grid);

		/// \brief computes sign of the distance field by negating and applying a recursive flood-fill algorithm for non-frozen voxels.
		static void ComputeSignUsingFloodFill(Geometry::ScalarGrid& grid);

		/// \brief after applying a pmp::HoleFilling, then all voxels whose outgoing ray intersects the mesh are interior.
		static void ComputeSignUsingRays(Geometry::ScalarGrid& grid);
	};

	/// \brief A wrapper for input settings for computing distance field.
	struct PointCloudDistanceFieldSettings
	{
		float CellSize{ 1.0f }; //>! size of a single distance voxel.
		float VolumeExpansionFactor{ 1.0f }; //>! expansion factor (how many times the minimum dimension of the point cloud's bounding box) for the resulting scalar grid volume.
		double TruncationFactor{ 0.1 }; //>! factor by which the minimum half-dimension of point cloud's bounding box gives rise to a truncation (cutoff) value for the distance field.
		BlurPostprocessingType BlurType{ BlurPostprocessingType::None }; //>! type of blur filter to be used for post-processing.
	};

	/// \brief A singleton object for computing distance fields to point clouds.
	class PointCloudDistanceFieldGenerator
	{
	public:
		/// \brief a deleted default constructor because m_KdTree is not default-constructable.
		PointCloudDistanceFieldGenerator() = delete;

		/**
		 * \brief Compute the signed distance field of given input point cloud.
		 * \param inputPoints             evaluated point cloud.
		 * \param settings                settings for the distance field.
		 * \return the computed distance field's ScalarGrid.
		 */
		static [[nodiscard]] Geometry::ScalarGrid Generate(const std::vector<pmp::vec3>& inputPoints, const PointCloudDistanceFieldSettings& settings);

	private:

		/**
		 * \brief A preprocessing approach for distance grid using nearest neighbor approximation
		 * \param grid         modifiable input grid.
		 */
		static void PreprocessGridFromPoints(Geometry::ScalarGrid& grid);

		//
		// ===================================================
		//

		inline static std::vector<pmp::vec3> m_Points{}; //>! the input point cloud
	};

	/**
	 * \brief Reports DistanceFieldGenerator's input to a given stream.
	 * \param inputMesh   input mesh for DistanceFieldGenerator.
	 * \param settings    input settings for DistanceFieldGenerator.
	 * \param os          output stream.
	 */
	void ReportInput(const pmp::SurfaceMesh& inputMesh, const DistanceFieldSettings& settings, std::ostream& os);

	/**
	 * \brief Reports DistanceFieldGenerator's output to a given stream.
	 * \param grid    resulting grid.
	 * \param os      output stream.
	 */
	void ReportOutput(const Geometry::ScalarGrid& grid, std::ostream& os);

} // namespace SDF
