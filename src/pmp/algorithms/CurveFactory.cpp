#include "CurveFactory.h"
// Copyright 2011-2021 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "CurveFactory.h"


namespace pmp
{
	namespace
	{
	    void DeformCircleWithSineWave(ManifoldCurve2D& curve, float amplitude, float freq)
	    {
	        for (const auto v : curve.vertices())
	        {
	            Point2 p = curve.position(v) - Point2(0, 0);
	            const float angle = atan2(p[1], p[0]);
	            vec2 direction = normalize(p);
	            curve.position(v) += direction * amplitude * (sin(freq * angle) + 1.0f);
	        }
	    }
	} // anonymous namespace

	ManifoldCurve2D CurveFactory::circle(
		const Point2& center, 
		Scalar radius, 
		size_t nSegments, 
		Scalar startAngle, 
		Scalar endAngle)
	{
        ManifoldCurve2D curve;
        const bool isFullCircle = std::fabs(endAngle - startAngle - 2.0f * M_PI) < 1e-6f;
        const size_t nVerts = (isFullCircle ? nSegments : nSegments + 1);
        curve.reserve(nVerts, nSegments);
        std::vector<Vertex> vertices;
        vertices.reserve(nVerts);
        const Scalar angleIncrement = (endAngle - startAngle) / static_cast<Scalar>(nSegments);
        for (size_t i = 0; i < nSegments; ++i)
        {
            const Scalar angle = startAngle + i * angleIncrement;
            Point2 position(
                center[0] + radius * std::cos(angle),
                center[1] + radius * std::sin(angle));
            vertices.push_back(curve.add_vertex(position));
        }

        // If it's a full circle, ensure the last vertex matches the first vertex exactly
        if (isFullCircle)
        {
            vertices.push_back(vertices[0]);
        }
        else
        {
            const Scalar angle = endAngle;
            const Point2 position(
                center[0] + radius * std::cos(angle),
                center[1] + radius * std::sin(angle));
            vertices.push_back(curve.add_vertex(position));
        }

        // Create edges between consecutive vertices
        for (size_t i = 0; i < vertices.size() - 1; ++i)
        {
            curve.new_edge(vertices[i], vertices[i + 1]);
        }

        return curve;
	}

    ManifoldCurve2D CurveFactory::sine_deformed_circle(const Point2& center, Scalar radius, size_t nSegments, float amplitude, float freq, Scalar startAngle, Scalar endAngle)
    {
        auto result = circle(center, radius, nSegments, startAngle, endAngle);
        DeformCircleWithSineWave(result, amplitude, freq);
        return result;
    }

    ManifoldCurve2D CurveFactory::rectangle(const Point2& center, Scalar sideX, Scalar sideY, size_t nSegments, bool chamferCorners)
    {
        ManifoldCurve2D curve;

        // Half lengths of the sides
        const Scalar halfSideX = sideX / 2.0;
        const Scalar halfSideY = sideY / 2.0;

        // Rectangle corner points (before chamfering)
        std::vector<Point2> corners{
            Point2{center[0] - halfSideX, center[1] - halfSideY},  // Bottom-left
            Point2{center[0] + halfSideX, center[1] - halfSideY},  // Bottom-right
            Point2{center[0] + halfSideX, center[1] + halfSideY},  // Top-right
            Point2{center[0] - halfSideX, center[1] + halfSideY}   // Top-left
        };

        // To store the vertices
        std::vector<Vertex> vertices;
        vertices.reserve(nSegments);

        // Iterate over each side of the rectangle
        for (size_t i = 0; i < 4; ++i)
        {
            Point2 currentCorner = corners[i];
            Point2 nextCorner = corners[(i + 1) % 4];

            // Add straight segments between two corners
            for (size_t j = 0; j <= nSegments; ++j)
            {
                Scalar t = static_cast<Scalar>(j) / static_cast<Scalar>(nSegments);
                Point2 position = (1 - t) * currentCorner + t * nextCorner;

                // If chamfering is enabled, skip adding the exact corner point
                if (chamferCorners && (j == nSegments || j == 0))
                {
                    continue;  // Skip the middle point (corner point)
                }

                // Add the position as a vertex to the curve
                vertices.push_back(curve.add_vertex(position));
            }
        }

        // Create edges between consecutive vertices
        for (size_t i = 0; i < vertices.size() - 1; ++i)
        {
            curve.new_edge(vertices[i], vertices[i + 1]);
        }

        // Close the loop by connecting the last vertex with the first
        curve.new_edge(vertices.back(), vertices[0]);

        return curve;
    }

    ManifoldCurve2D CurveFactory::sampled_polygon(const std::vector<Point2>& polyVertices, size_t nSegments, bool chamferCorners)
    {
        ManifoldCurve2D curve;
        if (polyVertices.size() < 3)
        {
            std::cerr << "Invalid polygon: at least 3 vertices required.\n";
            return curve;
        }

        // TODO: fix this nonsense

        // Calculate total length of the polyline
        Scalar totalLength = 0.0;
        for (size_t i = 0; i < polyVertices.size(); ++i)
        {
            totalLength += norm(polyVertices[(i + 1) % polyVertices.size()] - polyVertices[i]);
        }
        const Scalar meanSegmentLength = totalLength / static_cast<Scalar>(nSegments);

        for (size_t i = 0; i < polyVertices.size(); ++i)
        {
            Point2 startVertex = polyVertices[i];
            Point2 endVertex = polyVertices[(i + 1) % polyVertices.size()];
            const auto edgeLength = norm(endVertex - startVertex);

            for (size_t j = 0; j <= nSegments; ++j)
            {
                const Scalar t = static_cast<Scalar>(j) / static_cast<Scalar>(nSegments);
                Point2 sampledPoint = (1 - t) * startVertex + t * endVertex;
            }
        }

        // If chamfering, calculate chamfer length (proportional to the average segment length)
        const Scalar chamferLength = chamferCorners ? segmentLength : 0.0f;

        // Prepare to store vertices
        std::vector<Vertex> vertices;
        vertices.reserve(nSegments + (chamferCorners ? polyVertices.size() : 0));  // Reserve space for vertices

        // Loop over each segment of the polygon
        for (size_t i = 0; i < polyVertices.size(); ++i)
        {
            Point2 currentVertex = polyVertices[i];
            Point2 nextVertex = polyVertices[(i + 1) % polyVertices.size()];  // Loop around

            // If chamfering, skip exact corners and add chamfer points
            if (chamferCorners)
            {
                Point2 prevVertex = polyVertices[(i == 0) ? polyVertices.size() - 1 : i - 1];
                Point2 chamferStart = currentVertex + chamferLength * normalize(prevVertex - currentVertex);
                Point2 chamferEnd = currentVertex + chamferLength * normalize(nextVertex - currentVertex);

                vertices.push_back(curve.add_vertex(chamferStart));  // Start of chamfer
                vertices.push_back(curve.add_vertex(chamferEnd));    // End of chamfer

                // Sample remaining points between chamfer end and next vertex
                size_t remainingSegments = nSegments / polyVertices.size();  // Segments per edge
                for (size_t j = 1; j <= remainingSegments; ++j)
                {
                    Scalar t = static_cast<Scalar>(j) / static_cast<Scalar>(remainingSegments);
                    Point2 sampledPoint = (1 - t) * chamferEnd + t * nextVertex;
                    vertices.push_back(curve.add_vertex(sampledPoint));
                }
            }
            else
            {
                // Sample along the edge between the current and next vertex
                size_t remainingSegments = nSegments / polyVertices.size();  // Segments per edge
                for (size_t j = 0; j <= remainingSegments; ++j)
                {
                    Scalar t = static_cast<Scalar>(j) / static_cast<Scalar>(remainingSegments);
                    Point2 sampledPoint = (1 - t) * currentVertex + t * nextVertex;
                    vertices.push_back(curve.add_vertex(sampledPoint));
                }
            }
        }

        // Add edges between consecutive vertices
        for (size_t i = 0; i < vertices.size() - 1; ++i)
        {
            curve.new_edge(vertices[i], vertices[i + 1]);
        }

        // Close the loop by connecting the last vertex with the first
        curve.new_edge(vertices.back(), vertices[0]);

        return curve;
    }


} // namespace pmp