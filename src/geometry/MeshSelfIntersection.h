#pragma once

#include <pmp/SurfaceMesh.h>

#include <vector>
#include <unordered_map>

namespace Geometry
{
	/// \brief A utility for welding close-enough points on a triangle.
	void RemoveDuplicates(const pmp::SurfaceMesh& mesh, const pmp::Face& face, std::vector<pmp::vec3>& data);

	class MeshSelfIntersectionBucket
	{
	public:

		void AddFacePoints(const pmp::SurfaceMesh& mesh, const pmp::Face& face, const std::vector<pmp::vec3>& data)
		{
			m_Data[face.idx()] = data;
			RemoveDuplicates(mesh, face, m_Data[face.idx()]);
		}

		[[nodiscard]] const std::vector<pmp::vec3>& GetFacePoints(const unsigned int& faceId)
		{
			return m_Data[faceId];
		}

	private:

		std::unordered_map<unsigned int, std::vector<pmp::vec3>> m_Data{};
	};

	using FaceIntersectionMap = std::unordered_multimap<unsigned int, unsigned int>;

	class MeshSelfIntersectionBucketCollector
	{
	public:

		explicit MeshSelfIntersectionBucketCollector(const pmp::SurfaceMesh& mesh)
			: m_Mesh(mesh) {}

		[[nodiscard]] std::vector<MeshSelfIntersectionBucket> Retrieve();

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

		pmp::SurfaceMesh m_Mesh;

		FaceIntersectionMap m_FaceIntersections; //> a multimap mapping fId -> { ids of all faces intersecting f }

	};

	
} // namespace Geometry