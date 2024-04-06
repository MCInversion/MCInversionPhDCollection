#include "ConvexHull.h"


#include "pmp/algorithms/Normals.h"

#include <numeric>
#include <queue>
#include <ranges>
#include <unordered_set>

namespace
{
	using VisitedPoints = std::unordered_set<unsigned int>;

	[[nodiscard]] std::tuple<pmp::SurfaceMesh, VisitedPoints, std::vector<unsigned int>> FindInitialTetrahedron(const std::vector<pmp::Point>& points)
	{
		struct InitialTetrahedronIds
		{
			unsigned int i0{ UINT_MAX }, i1{ UINT_MAX }, i2{ UINT_MAX }, i3{ UINT_MAX };
		};

		// find the 6 points with min/max coordinate values
		std::vector<unsigned int> minMaxIds(6, 0); // Initialize with the first point
		for (unsigned int i = 0; i < points.size(); ++i)
		{
			if (points[i][0] < points[minMaxIds[0]][0]) minMaxIds[0] = i;
			if (points[i][1] < points[minMaxIds[1]][1]) minMaxIds[1] = i;
			if (points[i][2] < points[minMaxIds[2]][2]) minMaxIds[2] = i;
			if (points[i][0] > points[minMaxIds[3]][0]) minMaxIds[3] = i;
			if (points[i][1] > points[minMaxIds[4]][1]) minMaxIds[4] = i;
			if (points[i][2] > points[minMaxIds[5]][2]) minMaxIds[5] = i;
		}

		VisitedPoints visitedTetraPtIds; // tracking visited points

		// farthest two points of the tetrahedron base
		pmp::Scalar maxDistSq = 0.0f;
		InitialTetrahedronIds tetra;
		for (int i = 0; i < 6; ++i)
		{
			for (int j = i + 1; j < 6; ++j)
			{
				const pmp::Scalar distSq = sqrnorm(points[minMaxIds[i]] - points[minMaxIds[j]]);
				if (distSq > maxDistSq)
				{
					maxDistSq = distSq;
					tetra.i0 = minMaxIds[i];
					tetra.i1 = minMaxIds[j];
				}
			}
		}
		visitedTetraPtIds.insert({ tetra.i0, tetra.i1 });

		// third point to create the tetrahedron base
		maxDistSq = 0.0f;
		for (unsigned int i = 0; i < points.size(); ++i)
		{
			if (visitedTetraPtIds.contains(i))
				continue;

			const pmp::Scalar area = norm(cross(points[tetra.i1] - points[tetra.i0], points[i] - points[tetra.i0]));
			if (area > maxDistSq)
			{
				maxDistSq = area;
				tetra.i2 = i;
			}
		}
		visitedTetraPtIds.insert(tetra.i2);

		// find the fourth point of the tetrahedron
		maxDistSq = -FLT_MAX;
		const pmp::Normal planeNormal = normalize(cross(points[tetra.i1] - points[tetra.i0], points[tetra.i2] - points[tetra.i0]));
		const pmp::Scalar planeDist = dot(planeNormal, points[tetra.i0]);
		for (unsigned int i = 0; i < points.size(); ++i)
		{
			if (visitedTetraPtIds.contains(i))
				continue;

			const pmp::Scalar dist = std::abs(dot(planeNormal, points[i]) - planeDist);
			if (dist > maxDistSq)
			{
				maxDistSq = dist;
				tetra.i3 = i;
			}
		}
		visitedTetraPtIds.insert(tetra.i3);

		std::vector<unsigned int> remainingPoints;
		remainingPoints.reserve(points.size());
		for (unsigned int i = 0; i < points.size(); ++i)
		{
			if (visitedTetraPtIds.contains(i)) continue;
			remainingPoints.push_back(i);
		}
		remainingPoints.shrink_to_fit();

		// initialize the initial convex hull as a pmp::SurfaceMesh
		pmp::SurfaceMesh convexHullMesh;

		// Add points to the mesh
		const auto v0 = convexHullMesh.add_vertex(points[tetra.i0]);
		const auto v1 = convexHullMesh.add_vertex(points[tetra.i1]);
		const auto v2 = convexHullMesh.add_vertex(points[tetra.i2]);
		const auto v3 = convexHullMesh.add_vertex(points[tetra.i3]);

		// Add faces to form a tetrahedron
		convexHullMesh.add_triangle(v0, v1, v2);
		convexHullMesh.add_triangle(v0, v2, v3);
		convexHullMesh.add_triangle(v0, v3, v1);
		convexHullMesh.add_triangle(v1, v3, v2);

		return { convexHullMesh, visitedTetraPtIds, remainingPoints };
	}


	using FaceVisibilityLists = std::vector<std::vector<unsigned int>>;
	using FarthestVertexPerFaceList = std::vector<std::pair<unsigned int, pmp::Scalar>>;

	[[nodiscard]] std::pair<FaceVisibilityLists, FarthestVertexPerFaceList> CollectVisibilityListsAndFarthestPointsPerFace(
		const pmp::SurfaceMesh& convexHull, 
		const std::vector<unsigned int>& remainingPoints, 
		const std::vector<pmp::Point>& points,
		const pmp::FaceProperty<pmp::Normal>& fNormalsProp,
		const pmp::Scalar& distTolerance)
	{
		std::vector<std::vector<unsigned int>> pointsPerConvexHullFace(convexHull.n_faces());
		std::vector farthestVertexPerFace(convexHull.n_faces(), std::make_pair(UINT_MAX, 0.0f));
		for (const auto& remPtId : remainingPoints)
		{
			for (const auto f : convexHull.faces())
			{
				auto fVit = convexHull.vertices(f).begin();
				const pmp::Scalar dist = dot(points[remPtId] - convexHull.position(*fVit), fNormalsProp[f]);
				if (dist < distTolerance)
					continue;

				pointsPerConvexHullFace[f.idx()].push_back(remPtId);
				if (dist > farthestVertexPerFace[f.idx()].second)
				{
					farthestVertexPerFace[f.idx()].second = dist;
					farthestVertexPerFace[f.idx()].first = remPtId;
				}
			}
		}

		return { pointsPerConvexHullFace , farthestVertexPerFace };
	}

	using PointConstIterator = std::vector<pmp::Point>::const_iterator;

	using VertexToFaceDistance = std::pair<pmp::Vertex, pmp::Scalar>;

	//[[nodiscard]] std::pair<PointConstIterator, pmp::Scalar> FindFarthestPoint(const std::vector<pmp::Point>& points, const pmp::SurfaceMesh& mesh)
	//{
	//	PointConstIterator farthestPointIt = points.cend();
	//	pmp::Scalar maxDistance = -FLT_MAX;

	//	for (const auto f : mesh.faces())
	//	{
	//		std::vector<pmp::Point> facePoints;
	//		for (const auto v : mesh.vertices(f))
	//		{
	//			facePoints.push_back(mesh.position(v));
	//		}

	//		if (facePoints.size() >= 3)
	//		{
	//			const pmp::Point& v0 = facePoints[0];
	//			const pmp::Point normal = normalize(cross(facePoints[1] - v0, facePoints[2] - v0));
	//			const pmp::Scalar d = dot(normal, v0);

	//			for (auto it = points.cbegin(); it != points.cend(); ++it)
	//			{
	//				const pmp::Scalar distance = std::abs(dot(normal, *it) - d);
	//				if (distance > maxDistance)
	//				{
	//					farthestPointIt = it;
	//					maxDistance = distance;
	//				}
	//			}
	//		}
	//	}

	//	return { farthestPointIt, maxDistance };
	//}

	//void UpdateVisibilityData(const std::vector<pmp::Point>& inputPoints,
	//	pmp::SurfaceMesh& convexHull,
	//	const pmp::Scalar& distTolerance,
	//	std::vector<std::vector<unsigned int>>& listVertexPerFace,
	//	std::vector<std::pair<unsigned int, pmp::Scalar>>& farthestVertexPerFace,
	//	const std::unordered_set<unsigned int>& visited)
	//{
	//	// Reinitialize the vectors based on the current number of faces in convexHull
	//	listVertexPerFace.clear();
	//	listVertexPerFace.resize(convexHull.n_faces());

	//	farthestVertexPerFace.clear();
	//	farthestVertexPerFace.resize(convexHull.n_faces(), std::make_pair(0, -std::numeric_limits<pmp::Scalar>::max()));

	//	unsigned int faceIdx = 0;
	//	for (const auto f : convexHull.faces())
	//	{
	//		auto normal = pmp::Normals::compute_face_normal(convexHull, f); // Assume this function exists or is implemented
	//		const auto he = convexHull.halfedge(f);
	//		auto refPoint = convexHull.position(convexHull.to_vertex(he)); // Reference point on the face

	//		for (unsigned int i = 0; i < inputPoints.size(); ++i)
	//		{
	//			if (visited.contains(i))
	//				continue;// Skip if visited
	//			pmp::Scalar dist = dot(inputPoints[i] - refPoint, normal);
	//			if (dist < distTolerance)
	//				continue;

	//			listVertexPerFace[faceIdx].push_back(i);
	//			if (dist > farthestVertexPerFace[faceIdx].second)
	//			{
	//				farthestVertexPerFace[faceIdx] = std::make_pair(i, dist);
	//			}
	//		}
	//		++faceIdx;
	//	}
	//}

	//[[nodiscard]] std::pair<std::set<pmp::Face>, std::vector<pmp::Halfedge>> IdentifyVisibleFacesAndTheirBoundary(pmp::SurfaceMesh& mesh, const pmp::Point& newPoint)
	//{
	//	std::set<pmp::Face> visibleFacesSet;
	//	std::unordered_map<unsigned int, unsigned int> halfedgeVisibilityCount; // Counts visibility of edges

	//	// Identify visible faces and count visibility of each halfedge
	//	for (const auto f : mesh.faces())
	//	{
	//		auto fVit = mesh.vertices(f).begin(); // Get a vertex iterator of the face
	//		auto normal = pmp::Normals::compute_face_normal(mesh, f);
	//		if (dot(normal, newPoint - mesh.position(*fVit)) < 0)
	//			continue;

	//		// If dot product is positive, the face is visible
	//		visibleFacesSet.insert(f);
	//		for (const auto h : mesh.halfedges(f))
	//		{
	//			++halfedgeVisibilityCount[h.idx()];
	//		}
	//	}

	//	// Collect boundary halfedges of the visible region
	//	std::vector<pmp::Halfedge> boundaryHalfedges;
	//	for (const auto& [h, nVisited] : halfedgeVisibilityCount)
	//	{
	//		if (nVisited != 1)
	//			continue;

	//		// Halfedge is part of the boundary if only one of its faces is visible
	//		boundaryHalfedges.emplace_back(pmp::Halfedge(h));
	//	}

	//	return { visibleFacesSet, boundaryHalfedges };
	//}

	//void DeleteVisibleFacesAndInternalVertices(pmp::SurfaceMesh& mesh, const std::set<pmp::Face>& facesToDelete)
	//{
	//	// Step 1: Delete faces without altering mesh connectivity or invoking cleanup
	//	for (const auto f : facesToDelete)
	//	{
	//		mesh.delete_face_raw(f); // Marks the face as deleted
	//	}

	//	// Step 2: Identify vertices that should be deleted
	//	std::unordered_set<pmp::IndexType> vertexIndicesToDelete;
	//	for (const auto f : facesToDelete)
	//	{
	//		for (const auto v : mesh.vertices(f))
	//		{
	//			vertexIndicesToDelete.insert(v.idx()); // Initially mark all vertices of deleted faces
	//		}
	//	}

	//	// Step 3: Exclude vertices that are part of non-deleted faces (boundary vertices)
	//	for (const auto f : mesh.faces())
	//	{
	//		if (facesToDelete.contains(f))
	//			continue;

	//		// Face is not marked for deletion
	//		for (const auto v : mesh.vertices(f))
	//		{
	//			vertexIndicesToDelete.erase(v.idx()); // Preserve vertices part of non-deleted faces
	//		}
	//	}

	//	// Step 4: Delete marked vertices, ensuring connectivity is maintained
	//	for (const auto idx : vertexIndicesToDelete)
	//	{
	//		const pmp::Vertex v(idx);
	//		if (!mesh.is_valid(v))
	//			continue;

	//		mesh.mark_vertex_and_connections_for_deletion(v);
	//	}

	//	// Note: Garbage collection is deferred to a later stage, post re-triangulation
	//}


	//void UpdateConvexHull(pmp::SurfaceMesh& mesh, const pmp::Point& newPoint)
	//{
	//	const auto vNew = mesh.add_vertex(newPoint);
	//	const auto [visibleFaces, visibleFacesBoundaryHalfEdges] = IdentifyVisibleFacesAndTheirBoundary(mesh, newPoint);

	//	DeleteVisibleFacesAndInternalVertices(mesh, visibleFaces);

	//	// Triangulate the Hole with the new vertex
	//	for (const auto bdHe : visibleFacesBoundaryHalfEdges)
	//	{
	//		const auto start = mesh.from_vertex(bdHe);
	//		const auto end = mesh.to_vertex(bdHe);
	//		mesh.add_triangle(vNew, start, end); // add a new triangle without topological checks
	//	}

	//	mesh.garbage_collection();
	//}

} // anonymous namespace

namespace Geometry
{
	//void FillConvexHullWithPoints(pmp::SurfaceMesh& convexHull, const std::vector<pmp::Point>& points, const float& distTolerance)
	//{
	//	std::unordered_set<unsigned int> visited; // To track visited points

	//	// Step 1: Find the initial simplex (tetrahedron) from the point cloud
	//	auto [initialSimplex, remainingPoints] = FindInitialTetrahedron(points);
	//	InitializeFromTetrahedron(convexHull, initialSimplex); // Add initial simplex to the convexHull

	//	// Initialize visibility data for the initial simplex
	//	std::vector<std::vector<unsigned int>> listVertexPerFace;
	//	std::vector<std::pair<unsigned int, pmp::Scalar>> farthestVertexPerFace;
	//	UpdateVisibilityData(remainingPoints, convexHull, distTolerance, listVertexPerFace, farthestVertexPerFace, visited);

	//	// Step 2: Iteratively add points to the convex hull
	//	while (!remainingPoints.empty())
	//	{
	//		auto [farthestPointIt, distance] = FindFarthestPoint(remainingPoints, convexHull);
	//		if (distance <= distTolerance || farthestPointIt == remainingPoints.cend())
	//		{
	//			// All points are within the hull or closer than the tolerance
	//			break;
	//		}

	//		// Mark the point as visited
	//		visited.insert(static_cast<unsigned int>(std::distance(remainingPoints.cbegin(), farthestPointIt)));

	//		UpdateConvexHull(convexHull, *farthestPointIt);
	//		remainingPoints.erase(farthestPointIt);
	//		UpdateVisibilityData(remainingPoints, convexHull, distTolerance, listVertexPerFace, farthestVertexPerFace, visited);
	//	}
	//}

	pmp::SurfaceMesh ComputeConvexHullFromPoints(const std::vector<pmp::Point>& points, const float& distTolerance)
	{
		auto [convexHull, visitedPts, remainingPoints] = FindInitialTetrahedron(points);

		pmp::Normals::compute_face_normals(convexHull);
		auto fNormalsProp = convexHull.face_property<pmp::Normal>("f:normal");

		// Build a list of visible vertices for each face of the initial convex hull and find the farthest point
		auto [pointsPerConvexHullFace, farthestVertexPerFace] = 
			CollectVisibilityListsAndFarthestPointsPerFace(convexHull, remainingPoints, points, fNormalsProp, distTolerance);

		for (unsigned int i = 0; i < pointsPerConvexHullFace.size(); i++)
		{
			if (pointsPerConvexHullFace[i].empty())
				continue;

			const auto& newVertId = farthestVertexPerFace[i].first; // ID of the farthest point for this face
			std::queue<unsigned int> faceQueue; // Queue to process faces in breadth-first search manner
			std::unordered_set<unsigned int> processedFaces; // Track processed faces to avoid repetition
			std::vector<unsigned int> boundaryFaceIds; // Track boundary faces for hole filling

			faceQueue.push(i); // Start with the current face
			while (!faceQueue.empty())
			{
				unsigned int fId = faceQueue.front();
				faceQueue.pop();
				
				if (!processedFaces.insert(fId).second) 
					continue; // Skip if this face has already been processed

				for (const auto h : convexHull.halfedges(pmp::Face(fId)))
				{
					const auto oppH = convexHull.opposite_halfedge(h);
					if (convexHull.is_boundary(oppH)) 
						continue; // Skip if opposite halfedge is a boundary (no neighboring face)

					const auto neighborF = convexHull.face(oppH);
					if (!neighborF.is_valid() || processedFaces.contains(neighborF.idx()))
						continue;

					const auto baseV = convexHull.from_vertex(oppH);
					pmp::Scalar dist = dot(points[newVertId] - convexHull.position(baseV), fNormalsProp[neighborF]);

					if (dist > distTolerance) 
					{
						// This neighboring face is visible to the new vertex and should be marked for deletion
						processedFaces.insert(neighborF.idx());
						faceQueue.push(neighborF.idx());
					}
					else
					{
						// This face is on the boundary relative to the new vertex
						boundaryFaceIds.push_back(neighborF.idx());
					}
				}
			}

			// Atomically delete all processed faces (includes connectivity cleanups)
			for (const auto& fDelId : processedFaces)
				convexHull.delete_face(pmp::Face(fDelId));

			// Insert a new vertex if there are faces to connect it to.
			pmp::Vertex newConvexHullVertHandle{};
			if (!boundaryFaceIds.empty())
			{
				newConvexHullVertHandle = convexHull.add_vertex(points[newVertId]);
				visitedPts.insert(newVertId);
				remainingPoints.erase(remainingPoints.begin() + newVertId);
			}

			// Re-triangulate the hole with the new vertex
			for (const auto& boundaryFaceIdx : boundaryFaceIds)
			{
				for (const auto h : convexHull.halfedges(pmp::Face(boundaryFaceIdx)))
				{
					// iterate along the boundary face half-edges until finding a half-edge whose opposite belongs to a face-to-be-deleted.
					auto oppH = convexHull.opposite_halfedge(h);
					auto neighboringFace = convexHull.face(oppH);
					if (!processedFaces.contains(neighboringFace.idx()))
						continue;

					// Create a new face using the new vertex and the edge vertices of the current boundary face
					auto fromV = convexHull.from_vertex(h);
					auto toV = convexHull.to_vertex(h);
					auto newFace = convexHull.add_triangle(newConvexHullVertHandle, fromV, toV);

					fNormalsProp[newFace] = pmp::Normals::compute_face_normal(convexHull, newFace);

					// collect the visibility data for the new face
					std::vector<unsigned int> tempPtIds;
					std::vector<unsigned int> pointsToTest(pointsPerConvexHullFace[boundaryFaceIdx].size() + pointsPerConvexHullFace[neighboringFace.idx()].size());
					std::ranges::set_union(pointsPerConvexHullFace[boundaryFaceIdx],pointsPerConvexHullFace[neighboringFace.idx()], std::back_inserter(pointsToTest));

					auto newInfo = std::make_pair(UINT_MAX, 0.0f);
					for (const auto& ptToTestId : pointsToTest)
					{
						if (visitedPts.contains(ptToTestId))
							continue;

						auto fVit = convexHull.vertices(newFace).begin();
						const pmp::Scalar dist = dot(points[ptToTestId] - convexHull.position(*fVit), fNormalsProp[newFace]);
						if (dist < distTolerance)
							continue;

						tempPtIds.push_back(ptToTestId);
						pointsPerConvexHullFace[newFace.idx()].push_back(ptToTestId);
						if (dist > newInfo.second)
						{
							newInfo.second = dist;
							newInfo.first = ptToTestId;
						}
					}
					pointsPerConvexHullFace.push_back(tempPtIds);
					farthestVertexPerFace.push_back(newInfo);
				}
			}
		}
		convexHull.garbage_collection();

		return convexHull;
	}

} // namespace Geometry