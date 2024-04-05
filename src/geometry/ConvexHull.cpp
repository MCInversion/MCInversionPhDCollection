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
		visitedTetraPtIds.insert(tetra.i0, tetra.i1);

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
		auto [pointsPerConvexHullFace, farthestVertexPerFace] = CollectVisibilityListsAndFarthestPointsPerFace(convexHull, remainingPoints, points, fNormalsProp, distTolerance);

		std::vector<std::unordered_set<unsigned int>> facesToDelete(pointsPerConvexHullFace.size());
		std::vector<std::unordered_set<unsigned int>> boundaryFacesPerPoint(pointsPerConvexHullFace.size());
		for (unsigned int i = 0; i < pointsPerConvexHullFace.size(); i++)
		{
			if (pointsPerConvexHullFace[i].empty())
				continue;

			const auto& newVertId = farthestVertexPerFace[i].first; // ID of the farthest point for this face
			std::queue<unsigned int> faceQueue; // Queue to process faces in breadth-first search manner
			std::unordered_set<unsigned int> processedFaces; // Track processed faces to avoid repetition

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
						facesToDelete[i].insert(neighborF.idx());
						faceQueue.push(neighborF.idx());
					}
					else
					{
						// This face is on the boundary relative to the new vertex
						boundaryFacesPerPoint[i].insert(neighborF.idx());
					}
				}
			}

			pmp::Vertex newConvexHullVertHandle{};
			if (!boundaryFacesPerPoint[i].empty())
			{
				newConvexHullVertHandle = convexHull.add_vertex(points[newVertId]);
				visitedPts.insert(newVertId);
			}

			// Add a new face for each border
			std::unordered_map<unsigned int, std::vector<unsigned int>> newFacesMap;

			// TODO: I'm not sure this will work
			// Assuming `boundaryFacesPerPoint[i]` contains indices of boundary faces for the current point
			for (const auto& boundaryFaceIdx : boundaryFacesPerPoint[i]) 
			{
				for (const auto h : convexHull.halfedges(pmp::Face(boundaryFaceIdx)))
				{
					auto oppH = convexHull.opposite_halfedge(h);
					auto fromV = convexHull.from_vertex(h);
					auto toV = convexHull.to_vertex(h);

					// Create a new face using the new vertex and the edge vertices of the current boundary face
					auto newFace = convexHull.add_triangle(newConvexHullVertHandle, fromV, toV);

					// Store the new face information for potential further topology updates
					newFacesMap[boundaryFaceIdx].push_back(newFace.idx());
				}
			}

			// TODO: more work

			std::unordered_map< CHVertexPointer, std::pair<int, char> > fanMap;
			for (const auto& bdFaceId : bdFaceIds)
			{
				CHFacePointer f = &convexHull.face[bdFaceId];
				for (int j = 0; j < 3; j++)
				{
					if (f->IsB(j))
					{
						f->ClearB(j);
						//Add new face
						CHFaceIterator fi = vcg::tri::Allocator<CHMesh>::AddFace(convexHull, &convexHull.vert.back(), f->V1(j), f->V0(j));
						(*fi).N() = vcg::NormalizedTriangleNormal(*fi);
						f = &convexHull.face[indexFace];
						int newFace = vcg::tri::Index(convexHull, *fi);
						//Update convex hull FF topology
						CHVertexPointer vp[] = { f->V1(j), f->V0(j) };
						for (int ii = 0; ii < 2; ii++)
						{
							int indexE = ii * 2;
							typename std::unordered_map< CHVertexPointer, std::pair<int, char> >::iterator vIter = fanMap.find(vp[ii]);
							if (vIter != fanMap.end())
							{
								CHFacePointer f2 = &convexHull.face[(*vIter).second.first];
								char edgeIndex = (*vIter).second.second;
								f2->FFp(edgeIndex) = &convexHull.face.back();
								f2->FFi(edgeIndex) = indexE;
								fi->FFp(indexE) = f2;
								fi->FFi(indexE) = edgeIndex;
							}
							else
							{
								fanMap[vp[ii]] = std::make_pair(newFace, indexE);
							}
						}
						//Build the visibility list for the new face
						std::vector<InputVertexPointer> tempVect;
						int indices[2] = { indexFace, int(vcg::tri::Index(convexHull, f->FFp(j))) };
						std::vector<InputVertexPointer> vertexToTest(pointsPerConvexHullFace[indices[0]].size() + pointsPerConvexHullFace[indices[1]].size());
						typename std::vector<InputVertexPointer>::iterator tempIt = std::set_union(pointsPerConvexHullFace[indices[0]].begin(), pointsPerConvexHullFace[indices[0]].end(), pointsPerConvexHullFace[indices[1]].begin(), pointsPerConvexHullFace[indices[1]].end(), vertexToTest.begin());
						vertexToTest.resize(tempIt - vertexToTest.begin());

						Pair newInfo = std::make_pair((InputVertexPointer)NULL, 0.0f);
						for (size_t ii = 0; ii < vertexToTest.size(); ii++)
						{
							if (!(*vertexToTest[ii]).IsV())
							{
								float dist = ((*vertexToTest[ii]).P() - (*fi).P(0)).dot((*fi).N());
								if (dist > distTolerance)
								{
									tempVect.push_back(vertexToTest[ii]);
									if (dist > newInfo.second)
									{
										newInfo.second = dist;
										newInfo.first = vertexToTest[ii];
									}
								}
							}
						}
						pointsPerConvexHullFace.push_back(tempVect);
						farthestVertexPerFace.push_back(newInfo);
						//Update topology of the new face
						CHFacePointer ffp = f->FFp(j);
						int ffi = f->FFi(j);
						ffp->FFp(ffi) = ffp;
						ffp->FFi(ffi) = ffi;
						f->FFp(j) = &convexHull.face.back();
						f->FFi(j) = 1;
						fi->FFp(1) = f;
						fi->FFi(1) = j;
					}
				}
			}
			//Delete the faces inside the updated convex hull
			for (size_t j = 0; j < visFace.size(); j++)
			{
				if (!convexHull.face[visFace[j]].IsD())
				{
					std::vector<InputVertexPointer> emptyVec;
					vcg::tri::Allocator<CHMesh>::DeleteFace(convexHull, convexHull.face[visFace[j]]);
					pointsPerConvexHullFace[visFace[j]].swap(emptyVec);
				}
			}
		}

		return convexHull;
	}

} // namespace Geometry