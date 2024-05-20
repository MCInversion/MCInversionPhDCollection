#include "TriangulationUtils.h"

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

    [[nodiscard]] std::vector<pmp::Point> GetSubdividedOffsetPolyline(const std::vector<pmp::Point>& polyline, float maxSegmentLength, float offsetDistance)
    {
        std::vector<pmp::Point> subdivided;

        auto offsetPoint = [offsetDistance](const pmp::Point& p1, const pmp::Point& p2) -> pmp::Point {
            pmp::Point direction = p2 - p1;
            pmp::Point perpendicular(-direction[1], direction[0], 0.0f); // Perpendicular to the segment in the xy plane
            perpendicular /= std::sqrt(perpendicular[0] * perpendicular[0] + perpendicular[1] * perpendicular[1]); // Normalize
            return p1 + offsetDistance * perpendicular;
        };

        for (size_t i = 0; i < polyline.size(); ++i)
        {
            pmp::Point p1 = polyline[i];
            pmp::Point p2 = polyline[(i + 1) % polyline.size()]; // Wrap around to the start for closed loop
            subdivided.push_back(offsetPoint(p1, p2)); // Offset the first point of the segment
            float segmentLength = norm(p2 - p1);
            if (segmentLength > maxSegmentLength)
            {
                int numSubdivisions = static_cast<int>(std::ceil(segmentLength / maxSegmentLength));
                for (int j = 1; j < numSubdivisions; ++j)
                {
                    float t = static_cast<float>(j) / numSubdivisions;
                    pmp::Point interpolated = (1.0f - t) * p1 + t * p2;
                    subdivided.push_back(offsetPoint(interpolated, p2)); // Offset the interpolated points
                }
            }
            subdivided.push_back(offsetPoint(p2, p1)); // Offset the last point of the segment
        }
        return subdivided;
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

        // Create a Fade_2D triangulation object
        GEOM_FADE2D::Fade_2D dt;

        // Add points to the triangulation and map them to their indices
        std::vector<GEOM_FADE2D::Point2*> pointPtrs;
        for (const auto& vertex : data.Vertices)
        {
            pointPtrs.push_back(dt.insert(GEOM_FADE2D::Point2(vertex[0], vertex[1])));
        }

        // Extract the dense offset contour
        const auto meanDistance = Geometry::ComputeNearestNeighborMeanInterVertexDistance(data.Vertices, 6);
        const auto denseOffsetContour = GetSubdividedOffsetPolyline(constraintPolyline, meanDistance, meanDistance);

        // Add constraint polyline as segments
        std::vector<GEOM_FADE2D::Segment2> constraintSegments;
        for (size_t i = 0; i < denseOffsetContour.size(); ++i)
        {
            const size_t next_i = (i + 1) % denseOffsetContour.size(); // Ensure the last point connects back to the first point if closed
            constraintSegments.emplace_back(
                GEOM_FADE2D::Point2(denseOffsetContour[i][0], denseOffsetContour[i][1]),
                GEOM_FADE2D::Point2(denseOffsetContour[next_i][0], denseOffsetContour[next_i][1])
            );
        }

        dt.createConstraint(constraintSegments, GEOM_FADE2D::CIS_CONSTRAINED_DELAUNAY);

        // Perform the triangulation
        dt.applyConstraintsAndZones();

        // Retrieve the triangles and add them to your result
        std::vector<GEOM_FADE2D::Triangle2*> triangles;
        dt.getTrianglePointers(triangles);

        // Clear any existing PolyIndices
        data.PolyIndices.clear();

        // Map from Point2* to index in data.Vertices
        std::map<GEOM_FADE2D::Point2*, unsigned int> pointToIndexMap;
        for (unsigned int i = 0; i < pointPtrs.size(); ++i)
        {
            pointToIndexMap[pointPtrs[i]] = i;
        }

        // Convert Fade2D triangles to PolyIndices
        for (const auto& triangle : triangles)
        {
            unsigned int idx0 = pointToIndexMap[triangle->getCorner(0)];
            unsigned int idx1 = pointToIndexMap[triangle->getCorner(1)];
            unsigned int idx2 = pointToIndexMap[triangle->getCorner(2)];

            data.PolyIndices.push_back({ idx0, idx1, idx2 });
        }
    }

	constexpr float MAGIC_RADIUS_MULTIPLIER = 2.5f;

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
		constexpr auto clusteringPercentage = (1.0f / MAGIC_RADIUS_MULTIPLIER) * 100.0f;

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

        // Extract the dense offset contour
        const auto meanDistance = Geometry::ComputeNearestNeighborMeanInterVertexDistance(data.Vertices, 6);
        const auto denseOffsetContour = GetSubdividedOffsetPolyline(constraintPolyline, meanDistance, meanDistance);

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
                contour.reserve(denseOffsetContour.size());
                for (const auto& point : denseOffsetContour)
                {
                    //points.push_back(new Poly2Tri::Point(point[0], point[1]));
                    contour.push_back(new Poly2Tri::Point(point[0], point[1], points.size() + contour.size() - 1));
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