#pragma once

#include "geometry/CollisionKdTree.h"
#include "geometry/Grid.h"

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
		pmp::Scalar CellSize{ 1.0 }; //>! size of a single distance voxel.
		pmp::Scalar VolumeExpansionFactor{ 1.0 }; //>! expansion factor (how many times the minimum dimension of mesh's bounding box) for the resulting scalar grid volume.
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
		 * \param inputMesh               an adapter for the evaluated mesh.
		 * \param settings                settings for the distance field.
		 * \param customFieldBox          optional custom bounding box for the distance field.
		 * \return the computed distance field's ScalarGrid.
		 * \throw std::invalid_argument if inputMesh contains no spatial data!
		 */
		static [[nodiscard]] Geometry::ScalarGrid Generate(
			const Geometry::MeshAdapter& inputMesh, 
			const DistanceFieldSettings& settings,
			const std::optional<pmp::BoundingBox>& customFieldBox = std::nullopt);
		
	private:
		inline static std::unique_ptr<Geometry::MeshAdapter> m_Mesh{ nullptr }; //>! mesh to be (pre)processed.
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
		pmp::Scalar CellSize{ 1.0 }; //>! size of a single distance voxel.
		pmp::Scalar VolumeExpansionFactor{ 1.0 }; //>! expansion factor (how many times the minimum dimension of the point cloud's bounding box) for the resulting scalar grid volume.
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
		 * \throw std::invalid_argument if inputMesh contains no spatial data!
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

	//
	// ================================================================================
	//         2D Distance field generators
	// --------------------------------------------------------------------------------
	//

	/// \brief enumerator for grid preprocessing function type for distance field.
	enum class [[nodiscard]] PreprocessingType2D
	{
		Quadtree = 0,  //>! preprocess ScalarGrid2D using an QuadtreeVoxelizer.
		NoQuadtree = 1 //>! preprocess ScalarGrid2D without using an QuadtreeVoxelizer.
	};

	/// \brief enumerator for the approach to the computation of sign of the distance field to a curve (negative inside, positive outside).
	enum class [[nodiscard]] SignComputation2D
	{
		None = 0, //>! no sign for the distance field is to be computed.
		PixelFloodFill = 1, //>! negate and apply recursive flood-fill algorithm for non-frozen voxels.
	};

	// TODO: 2D Blur enum types (if useful)

	/// \brief A wrapper for input settings for computing distance field.
	struct DistanceField2DSettings
	{
		pmp::Scalar CellSize{ 1.0 }; //>! size of a single distance voxel.
		pmp::Scalar AreaExpansionFactor{ 1.0 }; //>! expansion factor (how many times the minimum dimension of mesh's bounding box) for the resulting scalar grid volume.
		double TruncationFactor{ 0.1 }; //>! factor by which the minimum half-dimension of mesh's bounding box gives rise to a truncation (cutoff) value for the distance field.
		KDTreeSplitType KDTreeSplit{ KDTreeSplitType::Center }; //>! the choice of a split function for KD-tree.
		SignComputation2D SignMethod{ SignComputation2D::None }; //>! method by which the sign of the distance field should be computed.
		PreprocessingType2D PreprocType{ PreprocessingType2D::Quadtree }; //>! function type for the preprocessing of distance field scalar grid.
	};

	/// \brief a functor for computing the sign of the distance field.
	using SignFunction2D = std::function<void(Geometry::ScalarGrid2D&)> const;

	/// \brief a functor for preprocessing the scalar grid for distance field.
	using PreprocessingFunction2D = std::function<void(Geometry::ScalarGrid2D&)>;

	/// \brief A singleton object for computing distance fields to 2D curves.
	class PlanarDistanceFieldGenerator
	{
	public:
		/// \brief a deleted default constructor because m_KdTree is not default-constructable.
		PlanarDistanceFieldGenerator() = delete;

		/**
		 * \brief Compute the signed distance field of given input curve.
		 * \param inputCurve              an adapter for the evaluated curve.
		 * \param settings                settings for the distance field.
		 * \param customFieldBox          optional custom bounding box for the distance field.
		 * \return the computed distance field's ScalarGrid2D.
		 * \throw std::invalid_argument if inputMesh contains no spatial data!
		 */
		static [[nodiscard]] Geometry::ScalarGrid2D Generate(
			const Geometry::CurveAdapter& inputCurve, 
			const DistanceField2DSettings& settings, 
			const std::optional<pmp::BoundingBox2>& customFieldBox = std::nullopt);

	private:
		inline static std::unique_ptr<Geometry::CurveAdapter> m_Curve{ nullptr }; //>! curve to be (pre)processed.
		inline static std::unique_ptr<Geometry::Collision2DTree> m_KdTree{ nullptr }; //>! curve kd tree.

		/**
		 * \brief provides the SignFunction, a function from this generator's private interface that computes the sign of the distance field.
		 * \param signCompType      sign computation function type identifier.
		 * \return the sign function identified by signCompType.
		 */
		static [[nodiscard]] SignFunction2D GetSignFunction(const SignComputation2D& signCompType);

		/**
		 * \brief Provides a preprocessing functor according to the given setting.
		 * \param preprocType      preprocessing type identifier.
		 * \return the preprocessing function identified by preprocType.
		 */
		static [[nodiscard]] PreprocessingFunction2D GetPreprocessingFunction(const PreprocessingType2D& preprocType);

		/**
		 * \brief A preprocessing approach for distance grid using Collision2DTree to create "pixel outline" of inputCurve.
		 * \param grid         modifiable input grid.
		 */
		static void PreprocessGridNoQuadtree(Geometry::ScalarGrid2D& grid);

		/**
		 * \brief A preprocessing approach for distance grid using Collision2DTree and QuadtreeVoxelizer to create "pixel outline" of inputMesh.
		 * \param grid         modifiable input grid.
		 */
		static void PreprocessGridWithQuadtree(Geometry::ScalarGrid2D& grid);

		/// \brief computes sign of the distance field by negating and applying a recursive flood-fill algorithm for non-frozen pixels.
		static void ComputeSignUsingFloodFill(Geometry::ScalarGrid2D& grid);
	};

	/// \brief A wrapper for input settings for computing distance field.
	struct PointCloudDistanceField2DSettings
	{
		pmp::Scalar CellSize{ 1.0 }; //>! size of a single distance voxel.
		pmp::Scalar AreaExpansionFactor{ 1.0 }; //>! expansion factor (how many times the minimum dimension of the point cloud's bounding box) for the resulting scalar grid volume.
		double TruncationFactor{ 0.1 }; //>! factor by which the minimum half-dimension of point cloud's bounding box gives rise to a truncation (cutoff) value for the distance field.
	};

	/// \brief A singleton object for computing distance fields to 2D point clouds.
	class PlanarPointCloudDistanceFieldGenerator
	{
	public:
		/// \brief a deleted default constructor because m_KdTree is not default-constructable.
		PlanarPointCloudDistanceFieldGenerator() = delete;

		/**
		 * \brief Compute the signed distance field of given input point cloud.
		 * \param inputPoints             evaluated point cloud.
		 * \param settings                settings for the distance field.
		 * \return the computed distance field's ScalarGrid.
		 * \throw std::invalid_argument if inputMesh contains no spatial data!
		 */
		static [[nodiscard]] Geometry::ScalarGrid2D Generate(const std::vector<pmp::Point2>& inputPoints, const PointCloudDistanceField2DSettings& settings);

	private:

		/**
		 * \brief A preprocessing approach for distance grid using nearest neighbor approximation
		 * \param grid         modifiable input grid.
		 */
		static void PreprocessGridFromPoints(Geometry::ScalarGrid2D& grid);

		//
		// ===================================================
		//

		inline static std::vector<pmp::Point2> m_Points{}; //>! the input point cloud
	};

	/// \brief A wrapper for input settings for computing distance field.
	struct ImageDistanceField2DSettings
	{
		pmp::Scalar CellSize{ 1.0 }; //>! size of a single distance voxel.
	};

	/// \brief A singleton object for computing distance fields to a chosen zero level of the input image.
	class ImageDistanceFieldGenerator
	{
	public:
		/// \brief a deleted default constructor because m_KdTree is not default-constructable.
		ImageDistanceFieldGenerator() = delete;

		/**
		 * \brief Compute the signed distance field of given input point cloud.
		 * \param absFileName             image filename.
		 * \param settings                settings for the distance field.
		 * \return the computed distance field's ScalarGrid.
		 * \throw std::invalid_argument if absFileName couldn't be processed
		 */
		static [[nodiscard]] Geometry::ScalarGrid2D Generate(const std::string& absFileName, const ImageDistanceField2DSettings& settings);

	private:
	};

} // namespace SDF
