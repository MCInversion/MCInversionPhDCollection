#pragma once

#include "CollisionKdTree.h"

#include <pmp/SurfaceMesh.h>

#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace Geometry
{
	/** ===============================================================================================
	 *
	 * \brief A smart container which gathers the important self-intersection spatial data (points, orientation, bounding volume)
	 *        for working with the intersections of neighboring faces.
	 * \class MeshSelfIntersectionBucket
	 *
	 =============================================================================================== */
	class MeshSelfIntersectionBucket
	{
	public:

		/// \brief Maps the data onto the given face index and updates the bucket's bounding volume.
		void AddFacePoints(const pmp::SurfaceMesh& mesh, const pmp::Face& face, const std::vector<pmp::vec3>& data);

		/// \brief Updates the bucket's orientation to account for orientation deviation of neighboring faces.
		void UpdateBucketOrientation(const pmp::Normal& normal)
		{
			m_Orientation += normal;
			m_Orientation = normalize(m_Orientation);
		}

		/// \brief Intersection points getter for a given face Id.
		[[nodiscard]] const std::vector<pmp::vec3>& GetFacePoints(const pmp::Face& face) const
		{
			return m_Data.at(face.idx());
		}

		/// \brief A getter for this bucket's orientation.
		[[nodiscard]] const pmp::Normal& GetOrientation() const
		{
			return m_Orientation;
		}

		/// \brief A getter for this bucket's bounding box.
		[[nodiscard]] const pmp::BoundingBox& GetBoundingBox() const
		{
			return m_BBox;
		}

	private:

		pmp::BoundingBox m_BBox{}; //> Bounds all the point and face vertex data of this bucket.
		pmp::Normal m_Orientation{}; //> Represents the face orientation of this bucket.
		std::unordered_map<unsigned int, std::vector<pmp::vec3>> m_Data{}; //> Maps the face id to the points intersecting the given face.
	};

	/// \brief Maps the index of a face to the indices of all faces that intersect it.
	using FaceIntersectionMap = std::unordered_multimap<unsigned int, unsigned int>;

	/** ===============================================================================================
	 *
	 * \brief An object for collecting buckets of intersection spatial data (points) from triangle-triangle intersections of the given surface mesh.
	 * \class MeshSelfIntersectionBucketCollector
	 *
	 =============================================================================================== */
	class MeshSelfIntersectionBucketCollector
	{
	public:
		/// \brief Constructor.
		explicit MeshSelfIntersectionBucketCollector(const pmp::SurfaceMesh& mesh)
			: m_Mesh(mesh) {}

		/// \brief The main functionality of this collector object. Extracts the necessary data, and fills the buckets with intersection points.
		[[nodiscard]] std::vector<MeshSelfIntersectionBucket> Retrieve(const bool& checkMesh = true);

	private:
		/**
		 * \brief Uses the CollisionKdTree to accelerate spatial search for triangle-triangle intersections which will be indexed in m_FaceIntersections.
		 */
		void ExtractFaceIntersectionMap();

		/**
		 * \brief Proceeds to the next face multimap entry. Because the mesh is supposedly watertight and manifold, if we iterate by picking the first available neighbor, we run in a loop.
		 * \param f   Current face id.
		 * \return pair { next iterator if valid, whether to terminate current bucket }
		 */
		[[nodiscard]] std::pair<FaceIntersectionMap::iterator, bool> Proceed(const pmp::Face& f);

		pmp::SurfaceMesh m_Mesh; //> the processed mesh instance.
		std::unique_ptr<CollisionKdTree> m_ptrKdTree{ nullptr }; //> A kD-tree instance for accelerated computing of box-triangle queries.

		std::unique_ptr<MeshSelfIntersectionBucket> m_ptrCurrentBucket{nullptr}; //> an instance of the bucket that is being filled by neighborhood querying.
		FaceIntersectionMap m_FaceIntersections; //> a multimap mapping fId -> { ids of all faces intersecting f }

		std::unordered_set<unsigned int> m_NeighborIds; //> a set for quick access to neighboring face ids.
	};

	
} // namespace Geometry