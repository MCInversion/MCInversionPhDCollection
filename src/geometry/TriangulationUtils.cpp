#include "TriangulationUtils.h"

#include "Fade_2D.h"

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

        // Add constraint polyline as segments
        std::vector<GEOM_FADE2D::Segment2> constraintSegments;
        for (size_t i = 0; i < constraintPolyline.size(); ++i)
        {
            const size_t next_i = (i + 1) % constraintPolyline.size(); // Ensure the last point connects back to the first point if closed
            constraintSegments.emplace_back(
                GEOM_FADE2D::Point2(constraintPolyline[i][0], constraintPolyline[i][1]),
                GEOM_FADE2D::Point2(constraintPolyline[next_i][0], constraintPolyline[next_i][1])
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
} // namespace Geometry