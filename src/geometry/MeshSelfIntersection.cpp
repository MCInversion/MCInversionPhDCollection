#include "MeshSelfIntersection.h"

#include "CollisionKdTree.h"
#include "GeometryUtil.h"

#include "pmp/algorithms/BarycentricCoordinates.h"

#include <unordered_set>

namespace
{
	void FillFaceVertices(const pmp::SurfaceMesh& mesh, const pmp::Face& face, std::vector<pmp::vec3>& result)
	{
		std::vector<pmp::vec3> vertices;
		vertices.reserve(3);
		for (const auto v : mesh.vertices(face))
		{
			vertices.push_back(mesh.position(v));
		}
	}

	void FillNeighboringFaceIds(const pmp::SurfaceMesh& mesh, const pmp::Face& face, std::unordered_set<unsigned int>& neighboringFaceIds)
	{
		for (const auto v : mesh.vertices(face))
		{
			for (const auto nf : mesh.faces(v))
			{
				neighboringFaceIds.insert(nf.idx());
			}
		}
	}

	void FillFaceCombinatorialAndSpatialData(const pmp::SurfaceMesh& mesh, const pmp::Face& face,
		std::vector<pmp::vec3>& vertices, 
		std::unordered_set<unsigned int>& neighboringFaceIds,
		pmp::BoundingBox& fBBox)
	{
		vertices.reserve(3);
		for (const auto v : mesh.vertices(face))
		{
			for (const auto nf : mesh.faces(v))
			{
				neighboringFaceIds.insert(nf.idx());
			}
			const auto& vPos = mesh.position(v);
			vertices.push_back(vPos);
			fBBox += vPos;
		}
	}


}

namespace Geometry
{
	constexpr float POLYLINE_END_DISTANCE_TOLERANCE = 1e-6f;

	std::vector<MeshSelfIntersectionBucket> MeshSelfIntersectionBucketCollector::Retrieve()
	{
		if (!m_Mesh.is_triangle_mesh())
		{
			std::cerr << "MeshSelfIntersectionBucketCollector::Retrieve: non-triangle SurfaceMesh not supported for this function!\n";
			return {};
		}

		std::vector<MeshSelfIntersectionBucket> buckets;
		ExtractFaceIntersectionMap();

		MeshSelfIntersectionBucket currentBucket;
		auto fIt = m_FaceIntersections.begin();
		while (!m_FaceIntersections.empty())
		{
			auto intersectingFaceRange = m_FaceIntersections.equal_range(fIt->first);
			const auto baseFace = pmp::Face(fIt->first);
			std::vector<pmp::vec3> vertices0;
			FillFaceVertices(m_Mesh, baseFace, vertices0);

			// process all faces intersecting baseFace
			std::vector<pmp::vec3> bucketPointData;
			for (auto rangeIt = intersectingFaceRange.first; rangeIt != intersectingFaceRange.second; ++rangeIt)
			{
				const auto intersectingFace = pmp::Face(rangeIt->second);
				std::vector<pmp::vec3> vertices1;
				FillFaceVertices(m_Mesh, intersectingFace, vertices1);

				auto intersectionLineOpt = ComputeTriangleTriangleIntersectionLine(vertices0, vertices1);
				if (!intersectionLineOpt.has_value())
					continue;

				if (norm(intersectionLineOpt->second - intersectionLineOpt->first) < POLYLINE_END_DISTANCE_TOLERANCE)
					continue; // degenerate intersection line

				bucketPointData.push_back(intersectionLineOpt->first);
				bucketPointData.push_back(intersectionLineOpt->second);
			}
			currentBucket.AddFacePoints(m_Mesh, baseFace, bucketPointData);

			m_FaceIntersections.erase(fIt->first); // erasing all entries for this face

			bool shouldTerminateBucket = false;
			std::tie(fIt, shouldTerminateBucket) = Proceed(baseFace);
			if (shouldTerminateBucket)
			{
				buckets.push_back(currentBucket);
				currentBucket = MeshSelfIntersectionBucket();
			}
		}

		return buckets;
	}

	void MeshSelfIntersectionBucketCollector::ExtractFaceIntersectionMap()
	{
		std::unordered_multimap<unsigned int, unsigned int> faceToIntersectingFaces;
		const auto ptrMeshCollisionKdTree = std::make_unique<CollisionKdTree>(m_Mesh, CenterSplitFunction);

		for (const auto f : m_Mesh.faces())
		{
			pmp::BoundingBox fBBox;
			std::vector<pmp::Point> vertices0;
			std::unordered_set<unsigned int> neighboringFaceIds;
			FillFaceCombinatorialAndSpatialData(m_Mesh, f, vertices0, neighboringFaceIds, fBBox);

			std::vector<unsigned int> candidateIds;
			ptrMeshCollisionKdTree->GetTrianglesInABox(fBBox, candidateIds);

			for (const auto& ci : candidateIds)
			{
				if (ci == f.idx() || neighboringFaceIds.contains(ci))
				{
					continue; // Skip self and neighboring faces
				}

				std::vector<pmp::Point> vertices1;
				vertices1.reserve(3);
				const auto cf = pmp::Face(ci);
				for (const auto cv : m_Mesh.vertices(cf))
				{
					vertices1.push_back(m_Mesh.position(cv));
				}

				if (TriangleIntersectsTriangle(vertices0, vertices1))
				{
					faceToIntersectingFaces.emplace(f.idx(), ci);
				}
			}
		}
	}

	std::pair<FaceIntersectionMap::iterator, bool> MeshSelfIntersectionBucketCollector::Proceed(const pmp::Face& f)
	{
		std::unordered_set<unsigned int> neighborIds;
		FillNeighboringFaceIds(m_Mesh, f, neighborIds);
		for (auto neighborId : neighborIds)
		{
			const auto fIt = m_FaceIntersections.find(neighborId);
			if (fIt != m_FaceIntersections.cend())
			{
				//  Because the mesh is supposedly watertight and manifold, if we iterate by picking the first available neighbor, we run in a loop.
				return { fIt, false };
			}
		}

		if (!m_FaceIntersections.empty())
		{
			return { m_FaceIntersections.begin(), true };
		}

		return { m_FaceIntersections.end(), true };
	}

	void RemoveDuplicates(const pmp::Face& face, const pmp::SurfaceMesh& mesh, std::vector<pmp::vec3>& data)
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

} // namespace Geometry