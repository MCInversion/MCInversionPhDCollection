#include "TriangulationUtils.h"

#include <unordered_set>

#include "Fade_2D.h"
#include "poly2tri/Poly2Tri_CDT.h"

namespace
{
    /// \brief Function to compute the normal of a triangle given its vertices
    pmp::Point ComputeNormal(const pmp::Point& p1, const pmp::Point& p2, const pmp::Point& p3)
    {
        pmp::Point u = p2 - p1;
        pmp::Point v = p3 - p1;
        pmp::Point normal = pmp::cross(u, v);
        normal /= pmp::norm(normal);
        return normal;
    }

    //[[nodiscard]] std::vector<pmp::Point> GetSubdividedOffsetPolyline(const std::vector<pmp::Point>& polyline, pmp::Scalar maxSegmentLength, pmp::Scalar offsetDistance)
    //{
    //    std::vector<pmp::Point> subdivided;

    //    auto offsetPoint = [offsetDistance](const pmp::Point& p1, const pmp::Point& p2) -> pmp::Point {
    //        pmp::Point direction = p2 - p1;
    //        pmp::Point perpendicular(-direction[1], direction[0], 0.0); // Perpendicular to the segment in the xy plane
    //        perpendicular /= std::sqrt(perpendicular[0] * perpendicular[0] + perpendicular[1] * perpendicular[1]); // Normalize
    //        return p1 + offsetDistance * perpendicular;
    //    };

    //    for (size_t i = 0; i < polyline.size(); ++i)
    //    {
    //        pmp::Point p1 = polyline[i];
    //        pmp::Point p2 = polyline[(i + 1) % polyline.size()]; // Wrap around to the start for closed loop
    //        subdivided.push_back(offsetPoint(p1, p2)); // Offset the first point of the segment
    //        const pmp::Scalar segmentLength = norm(p2 - p1);
    //        if (segmentLength > maxSegmentLength)
    //        {
    //            const int numSubdivisions = static_cast<int>(std::ceil(segmentLength / maxSegmentLength));
    //            for (int j = 1; j < numSubdivisions; ++j)
    //            {
    //                pmp::Scalar t = static_cast<pmp::Scalar>(j) / numSubdivisions;
    //                pmp::Point interpolated = (1.0 - t) * p1 + t * p2;
    //                subdivided.push_back(offsetPoint(interpolated, p2)); // Offset the interpolated points
    //            }
    //        }
    //        subdivided.push_back(offsetPoint(p2, p1)); // Offset the last point of the segment
    //    }
    //    return subdivided;
    //}

	/// \brief Function to get subdivided offset polyline indices considering only 2D points
    [[nodiscard]] std::vector<int> GetSubdividedOffsetPolylineIndices2D(const std::vector<pmp::Point>& polyline, const std::vector<pmp::Point>& pointSet, pmp::Scalar maxSegmentLength)
    {
        std::vector<pmp::Point> subdivided;
        std::unordered_set<int> uniqueIndices;

        auto subdivideSegment = [&](const pmp::Point& p1, const pmp::Point& p2) {
            subdivided.push_back(p1);
            const pmp::Scalar segmentLength = norm(p2 - p1);
            if (segmentLength > maxSegmentLength)
            {
                const int numSubdivisions = static_cast<int>(std::ceil(segmentLength / maxSegmentLength));
                for (int j = 1; j < numSubdivisions; ++j)
                {
                    pmp::Scalar t = static_cast<pmp::Scalar>(j) / numSubdivisions;
                    pmp::Point interpolated = (1.0 - t) * p1 + t * p2;
                    subdivided.push_back(interpolated);
                }
            }
            subdivided.push_back(p2);
        };

        for (size_t i = 0; i < polyline.size(); ++i)
        {
            std::cout << "Subdividing segment " << i + 1 << " of " << polyline.size() << "\n";
            pmp::Point p1 = polyline[i];
            pmp::Point p2 = polyline[(i + 1) % polyline.size()]; // Wrap around to the start for closed loop
            subdivideSegment(p1, p2);
        }

        std::cout << "Subdivided polyline has " << subdivided.size() << " points\n";

        // Find unique indices of the closest points
        std::vector<int> indices;
        indices.reserve(subdivided.size());
        for (const auto& point : subdivided)
        {
            int closestIndex = Geometry::GetClosestPointIndex2D(pointSet, point);
            std::cout << "Closest index: " << closestIndex << "\n";
            if (closestIndex != -1 && !uniqueIndices.contains(closestIndex))
            {
                indices.push_back(closestIndex);
                uniqueIndices.insert(closestIndex);
            }
        }

        return indices;
    }
}

namespace Geometry
{
    void TriangulateWithFade2D(BaseMeshGeometryData& data, const std::vector<pmp::Point>& constraintPolyline)
    {
        if (data.Vertices.empty())
        {
            std::cerr << "Geometry::TriangulateWithFade2D: nothing to triangulate!\n";
            return;
        }

        std::cout << "Geometry::TriangulateWithFade2D: begin\n";

        // Extract the contour
        const auto meanDistance = ComputeNearestNeighborMeanInterVertexDistance(data.Vertices, 6);
        std::cout << "Geometry::TriangulateWithFade2D: meanDistance = " << meanDistance << "\n";
        const auto contourIds = GetSubdividedOffsetPolylineIndices2D(constraintPolyline, data.Vertices, 0.7 * meanDistance);
        std::cout << "Geometry::TriangulateWithFade2D: " << contourIds.size() << " contourIds extracted\n";

        // Create a Fade_2D triangulation object
        GEOM_FADE2D::Fade_2D dt;

        // Add points to the triangulation and map them to their indices
        std::vector<GEOM_FADE2D::Point2> fadePts;
        for (const auto& vertex : data.Vertices)
        {
            fadePts.push_back(GEOM_FADE2D::Point2(vertex[0], vertex[1]));
        }
        // Triangulation takes place
        std::vector<GEOM_FADE2D::Point2*> fadePtPtrs;
        dt.insert(fadePts, fadePtPtrs);

        std::cout << "Geometry::TriangulateWithFade2D: Fade2D triangulation complete, applying constraint zone.\n";

        // Add constraint polyline as segments
        std::vector<GEOM_FADE2D::Segment2> constraintSegments;
        for (size_t i = 0; i < contourIds.size(); ++i)
        {
            const size_t next_i = (i + 1) % contourIds.size(); // Ensure the last point connects back to the first point if closed
            constraintSegments.emplace_back(
                GEOM_FADE2D::Point2(data.Vertices[contourIds[i]][0], data.Vertices[contourIds[i]][1]),
                GEOM_FADE2D::Point2(data.Vertices[contourIds[next_i]][0], data.Vertices[contourIds[next_i]][1])
            );
        }
        GEOM_FADE2D::ConstraintGraph2* pTriangleBordersGraph = dt.createConstraint(constraintSegments, GEOM_FADE2D::CIS_CONSTRAINED_DELAUNAY, true);

        // Extract all dt triangles
        std::vector<GEOM_FADE2D::Triangle2*> triangles;
        dt.getTrianglePointers(triangles);

        // Create an inside seed point and grow a zone from there
        GEOM_FADE2D::Point2& p0(*triangles[0]->getCorner(0));
        GEOM_FADE2D::Point2& p1(*triangles[0]->getCorner(1));
        GEOM_FADE2D::Point2& p2(*triangles[0]->getCorner(2));
        GEOM_FADE2D::Point2 seedPoint((p0.x() + p1.x() + p2.x()) / 3.0, (p0.y() + p1.y() + p2.y()) / 3.0);
        GEOM_FADE2D::Zone2* pTriangleZone = dt.createZone(pTriangleBordersGraph, GEOM_FADE2D::ZL_GROW, seedPoint, true);

        // Retrieve the triangles and add them to your result
        std::vector<GEOM_FADE2D::Triangle2*> zoneTriangles;
        pTriangleZone->getTriangles(zoneTriangles);

        std::cout << "Geometry::TriangulateWithFade2D: Constraint zone extracted.\n";

        // Clear any existing PolyIndices
        data.PolyIndices.clear();

        // Extract the indices of the vertices of the triangles
        std::unordered_map<const GEOM_FADE2D::Point2*, unsigned int> vertexIndices;
        for (unsigned int count = 0; const auto* p : fadePtPtrs)
        {
            vertexIndices[p] = count++;
        }
        // Convert Fade2D triangles to PolyIndices
        for (const auto* tri : zoneTriangles)
        {
            data.PolyIndices.push_back({
				vertexIndices[tri->getCorner(0)],
				vertexIndices[tri->getCorner(1)],
				vertexIndices[tri->getCorner(2)]
			});
		}

        std::cout << "Geometry::TriangulateWithFade2D: Vertex indices extracted.\n";
    }

	constexpr pmp::Scalar MAGIC_RADIUS_MULTIPLIER = 2.5;

	void TriangulateWithVCGBPA(BaseMeshGeometryData& data)
	{
		if (data.Vertices.empty())
		{
			std::cerr << "Geometry::TriangulateWithVCGBPA: nothing to triangulate!\n";
			return;
		}

		// Compute an appropriate radius based on point distribution
		const auto meanDistance = Geometry::ComputeNearestNeighborMeanInterVertexDistance(data.Vertices, 6);
		const auto ballRadius = meanDistance * MAGIC_RADIUS_MULTIPLIER;
		constexpr auto clusteringPercentage = (1.0 / MAGIC_RADIUS_MULTIPLIER) * 100.0;

		const auto meshDataOpt = ComputeBallPivotingMeshFromPoints(data.Vertices, ballRadius, clusteringPercentage);
		if (!meshDataOpt.has_value())
		{
			std::cerr << "Geometry::TriangulateWithVCGBPA: Internal algorithm error!\n";
			return;
		}
		const auto& meshData = meshDataOpt.value();
        data.PolyIndices = meshData.PolyIndices;

        // Reorient all triangles upwards
        for (auto& triangle : data.PolyIndices)
        {
            const pmp::Point& p1 = data.Vertices[triangle[0]];
            const pmp::Point& p2 = data.Vertices[triangle[1]];
            const pmp::Point& p3 = data.Vertices[triangle[2]];

            pmp::Point normal = ComputeNormal(p1, p2, p3);
            if (normal[2] < 0) // Normal is pointing downwards
            {
                std::swap(triangle[1], triangle[2]); // Flip the triangle
            }
        }
	}

    constexpr unsigned int TRIANG_MAX_TRIES = 42;

    void TriangulateWithPoly2Tri(BaseMeshGeometryData& data, const std::vector<pmp::Point>& constraintPolyline)
    {
        if (data.Vertices.empty())
        {
            std::cerr << "Geometry::TriangulateWithPoly2Tri: nothing to triangulate!\n";
            return;
        }

        // Extract the contour
        const auto meanDistance = ComputeNearestNeighborMeanInterVertexDistance(data.Vertices, 6);
        const auto contourIds = GetSubdividedOffsetPolylineIndices2D(constraintPolyline, data.Vertices, 0.7 * meanDistance);

        // Project the 3D points onto a 2D plane (assuming projection along Z-axis)
        std::vector<Poly2Tri::Point*> projections;
        projections.reserve(data.Vertices.size());
        for (const auto& vertex : data.Vertices)
        {
            projections.push_back(new Poly2Tri::Point(vertex[0], vertex[1]));
        }
        Poly2Tri::CDT* constrainedDelaunay = nullptr;

        // Seed the random number generator
        std::srand(static_cast<unsigned int>(std::time(nullptr)));

        // Attempt triangulation with jitter to avoid collinear points and duplicates
        for (unsigned int triesUtilized = 0; triesUtilized < TRIANG_MAX_TRIES; triesUtilized++)
        {
            std::vector<Poly2Tri::Point*> points;
            points.reserve(projections.size());
            for (unsigned int i = 0; i < projections.size(); i++)
            {
                auto pt = new Poly2Tri::Point(
                    projections[i]->x + (triesUtilized * (static_cast<double>(std::rand()) / RAND_MAX * 0.001 - 0.0005)),
                    projections[i]->y + (triesUtilized * (static_cast<double>(std::rand()) / RAND_MAX * 0.001 - 0.0005))
                );
                pt->_vIdx = i;
                points.emplace_back(pt);
            }

            try
            {
                std::vector<Poly2Tri::Point*> contour;
                contour.reserve(contourIds.size());
                for (const auto& id : contourIds)
                {
                    contour.push_back(new Poly2Tri::Point(data.Vertices[id][0], data.Vertices[id][1]));
                }
                constrainedDelaunay = new Poly2Tri::CDT(contour);

                // Add points as Steiner points (interior points)
                for (const auto& pt : points)
                {
                    constrainedDelaunay->AddPoint(pt);
                }

                if (constrainedDelaunay == nullptr)
                {
                    throw(std::runtime_error("TriangulateWithPoly2Tri: Failed triangulation! constrainedDelaunay == nullptr\n"));
                }
                constrainedDelaunay->Triangulate();
                for (const auto* pt : contour)
                {
                    delete pt;
                }
                break; // Assuming it all goes well
            }
            catch (const std::runtime_error& err)
            {
                std::cerr << err.what() << "\nTriangulateWithPoly2Tri: Triangulation was not able to finish with "
                    << TRIANG_MAX_TRIES << " tries!\n";
            }
        }

        if (constrainedDelaunay == nullptr)
        {
            throw(std::runtime_error("TriangulateWithPoly2Tri: Failed triangulation! constrainedDelaunay == nullptr\n"));
        }

        const auto triangles = constrainedDelaunay->GetTriangles();

        // Clear any existing PolyIndices
        data.PolyIndices.clear();
        data.PolyIndices.reserve(triangles.size());

        // Convert Poly2Tri triangles to PolyIndices
        for (const auto& tri : triangles)
        {
            data.PolyIndices.push_back({
                tri->GetPoint(0)->_vIdx,
                tri->GetPoint(1)->_vIdx,
                tri->GetPoint(2)->_vIdx
                });
        }

        // Cleanup
        for (const auto* pt : projections)
        {
            delete pt;
        }
        delete constrainedDelaunay;
    }
} // namespace Geometry