#include "MeshSelfIntersection.h"

#include "CollisionKdTree.h"
#include "GeometryUtil.h"

#include "pmp/algorithms/BarycentricCoordinates.h"
#include "pmp/algorithms/Normals.h"

constexpr float POLYLINE_END_DISTANCE_TOLERANCE = 1e-6f;

namespace
{
	void FillFaceVertices(const pmp::SurfaceMesh& mesh, const pmp::Face& face, std::vector<pmp::vec3>& result)
	{
		result.clear();
		result.reserve(3);
		for (const auto v : mesh.vertices(face))
		{
			result.push_back(mesh.position(v));
		}
	}

	void FillNeighboringFaceIds(const pmp::SurfaceMesh& mesh, const pmp::Face& face, const Geometry::FaceIntersectionMap& faceIntersections, std::unordered_set<unsigned int>& neighboringFaceIds)
	{
		neighboringFaceIds.reserve(3);
		for (const auto he : mesh.halfedges(face))
		{
			const auto oppHe = mesh.opposite_halfedge(he);
			if (!oppHe.is_valid())
				continue;

			const auto newNeighborId = mesh.face(oppHe).idx();
			if (neighboringFaceIds.contains(newNeighborId))
				continue; // Skip duplicates

			if (!faceIntersections.contains(newNeighborId))
				continue; // Skip non-intersecting faces

			neighboringFaceIds.insert(mesh.face(oppHe).idx());
		}
	}

	/// \brief A utility for obtaining vertex and neighboring face data for self-intersection mapping evaluation
	void FillFaceCombinatorialAndSpatialData(const pmp::SurfaceMesh& mesh, const pmp::Face& face,
		std::vector<pmp::vec3>& vertices, 
		std::unordered_set<unsigned int>& neighboringFaceIds,
		pmp::BoundingBox& fBBox)
	{
		vertices.clear();
		neighboringFaceIds.clear();
		vertices.reserve(3);
		// an average vertex has valence 6 which yields 12 neighboring faces per 3 vertices sharing a triangle
		neighboringFaceIds.reserve(12);
		for (const auto v : mesh.vertices(face))
		{
			// we need to consider also vertex-adjacent faces
			for (const auto nf : mesh.faces(v))
			{
				neighboringFaceIds.insert(nf.idx());
			}
			const auto& vPos = mesh.position(v);
			vertices.push_back(vPos);
			fBBox += vPos;
		}
	}

	/// \brief A utility for welding close-enough points on a triangle.
	void RemoveDuplicates(const pmp::SurfaceMesh& mesh, const pmp::Face& face, std::vector<pmp::vec3>& data)
	{
		std::vector<pmp::vec3> vertices;
		FillFaceVertices(mesh, face, vertices);
		std::vector<std::pair<pmp::Point, pmp::vec3>> sortedPts;
		sortedPts.reserve(data.size());
		for (const auto& pt : data)
		{
			sortedPts.emplace_back(
				barycentric_coordinates(pt, vertices[0], vertices[1], vertices[2]),
				pt
			);
		}
		std::ranges::sort(sortedPts, [](const auto& a, const auto& b) {
			if (std::abs(a.first[0] - b.first[0]) > POLYLINE_END_DISTANCE_TOLERANCE) return a.first[0] < b.first[0];
			return a.first[1] < b.first[1];
			});
		data.clear();
		data.reserve(sortedPts.size());
		data.push_back(sortedPts[0].second);
		for (int i = 1; i < sortedPts.size(); ++i)
		{
			if (norm(sortedPts[i].second - sortedPts[i - 1].second) < POLYLINE_END_DISTANCE_TOLERANCE)
				continue;

			data.push_back(sortedPts[i].second);
		}
	}
} // namespace

/// \brief A simple verification utility for the mesh. An empty initializer list (implicitly converted to an empty vector) is returned.
#define VERIFY_MESH(msg, check)          \
if (check && !m_Mesh.is_triangle_mesh()) \
{                                        \
	std::cerr << msg;                    \
	return {};                           \
}

namespace Geometry
{
	std::vector<MeshSelfIntersectionBucket> MeshSelfIntersectionBucketCollector::Retrieve(const bool& checkMesh)
	{
		VERIFY_MESH("MeshSelfIntersectionBucketCollector::Retrieve: non-triangle SurfaceMesh not supported for this function!\n", checkMesh)

		std::vector<MeshSelfIntersectionBucket> buckets;
		ExtractFaceIntersectionMap();

		m_ptrCurrentBucket = std::make_unique<MeshSelfIntersectionBucket>();
		auto fIt = m_FaceIntersections.begin();
		while (fIt != m_FaceIntersections.end())
		{
			auto intersectingFaceRange = m_FaceIntersections.equal_range(fIt->first);
			const auto baseFace = pmp::Face(fIt->first);
			std::vector<pmp::vec3> vertices0;
			FillFaceVertices(m_Mesh, baseFace, vertices0);
			FillNeighboringFaceIds(m_Mesh, baseFace, m_FaceIntersections, m_NeighborIds);

			// process all faces intersecting baseFace
			std::vector<pmp::vec3> bucketPointData;
			for (auto rangeIt = intersectingFaceRange.first; rangeIt != intersectingFaceRange.second; ++rangeIt)
			{
				if (!m_FaceIntersections.contains(rangeIt->second))
					continue;

				const auto intersectingFace = pmp::Face(rangeIt->second);
				std::vector<pmp::vec3> vertices1;
				FillFaceVertices(m_Mesh, intersectingFace, vertices1);

				const auto intersectionLineOpt = ComputeTriangleTriangleIntersectionLine(vertices0, vertices1);
				if (!intersectionLineOpt.has_value())
					continue;

				const auto& intersectionLine = intersectionLineOpt.value();
				if (norm(intersectionLine.second - intersectionLine.first) < POLYLINE_END_DISTANCE_TOLERANCE)
					continue; // degenerate intersection line

				bucketPointData.push_back(intersectionLine.first);
				bucketPointData.push_back(intersectionLine.second);
			}
			if (!bucketPointData.empty())
			{
				RemoveDuplicates(m_Mesh, baseFace, bucketPointData);
				m_ptrCurrentBucket->AddFacePoints(m_Mesh, baseFace, bucketPointData);
				const auto faceNormal = pmp::Normals::compute_face_normal(m_Mesh, baseFace);
				m_ptrCurrentBucket->UpdateBucketOrientation(faceNormal);
			}
			m_FaceIntersections.erase(fIt->first); // erasing all entries for this face

			bool shouldTerminateBucket = false;
			std::tie(fIt, shouldTerminateBucket) = Proceed(baseFace);
			if (shouldTerminateBucket && !m_ptrCurrentBucket->Empty())
			{
				buckets.push_back(std::move(*m_ptrCurrentBucket));
				m_ptrCurrentBucket = std::make_unique<MeshSelfIntersectionBucket>();
			}
		}

		return buckets;
	}

	void MeshSelfIntersectionBucketCollector::ExtractFaceIntersectionMap()
	{
		m_FaceIntersections.clear();
		m_ptrKdTree = std::make_unique<CollisionKdTree>(m_Mesh, CenterSplitFunction);

		for (const auto f : m_Mesh.faces())
		{
			pmp::BoundingBox fBBox;
			std::vector<pmp::Point> vertices0;
			std::unordered_set<unsigned int> neighboringFaceIds;
			FillFaceCombinatorialAndSpatialData(m_Mesh, f, vertices0, neighboringFaceIds, fBBox);

			std::vector<unsigned int> candidateIds;
			m_ptrKdTree->GetTrianglesInABox(fBBox, candidateIds);

			for (const auto& ci : candidateIds)
			{
				if (ci == f.idx() || neighboringFaceIds.contains(ci))
				{
					continue; // Skip self and neighboring faces
				}

				const auto cf = pmp::Face(ci);
				std::vector<pmp::Point> vertices1;
				FillFaceVertices(m_Mesh, cf, vertices1);

				if (TriangleIntersectsTriangle(vertices0, vertices1))
				{
					m_FaceIntersections.emplace(f.idx(), ci);
				}
			}
		}
	}

	std::pair<FaceIntersectionMap::iterator, bool> MeshSelfIntersectionBucketCollector::Proceed(const pmp::Face& f)
	{
		if (!m_NeighborIds.empty())
		{
			//  Because the mesh is supposedly watertight and manifold,
			// if we iterate across neighboring faces, we fill the bucket with all the correct faces.
			const auto nId = *m_NeighborIds.begin();
			m_NeighborIds.erase(m_NeighborIds.begin());

			const auto fIt = m_FaceIntersections.find(nId);
			if (fIt != m_FaceIntersections.cend())
			{
				return { fIt, false };
			}
		}

		if (!m_FaceIntersections.empty() && m_ptrCurrentBucket)
		{
			// If there are still faces in the map, but no more neighbors to process, we should terminate the bucket
			return { m_FaceIntersections.begin(), true };
		}

		return { m_FaceIntersections.end(), true };
	}

	void MeshSelfIntersectionBucket::AddFacePoints(const pmp::SurfaceMesh& mesh, const pmp::Face& face, const std::vector<pmp::vec3>& data)
	{
		m_Data[face.idx()] = data;

		// The bounding box must be expanded by face vertices because it's later used for neighborhood querying
		// the convex hull of the face already contains the inserted point data
		std::vector<pmp::vec3> vertices;
		FillFaceVertices(mesh, face, vertices);
		m_BBox += vertices;
	}

	std::vector<std::pair<unsigned int, pmp::vec3>> ExtractFaceData(const MeshSelfIntersectionBucket& bucket, pmp::Point& barycenter)
	{
		barycenter = pmp::Point(0, 0, 0);
		std::vector<std::pair<unsigned int, pmp::vec3>> faceData;
		size_t nPts = 0;
		for (const auto& [faceId, points] : bucket)
		{
			for (const auto& pt : points) 
			{
				barycenter += pt;
				faceData.emplace_back(faceId, pt);
				++nPts;
			}
		}
		barycenter /= static_cast<pmp::Scalar>(nPts);
		return faceData;
	}

} // namespace Geometry