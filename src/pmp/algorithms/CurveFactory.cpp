#include "CurveFactory.h"
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

        [[nodiscard]] Scalar CalculatePolyPointsSignedArea(const std::vector<Point2>& polyVertices)
        {
            if (polyVertices.size() < 3)
            {
                throw std::invalid_argument("CalculatePolyPointsSignedArea: At least 3 vertices are required.");
            }

            Scalar signedArea = 0.0;
            for (size_t i = 0; i < polyVertices.size(); ++i)
            {
                const auto& current = polyVertices[i];
                const auto& next = polyVertices[(i + 1) % polyVertices.size()]; // Wrap around to first vertex
                signedArea += (current[0] * next[1] - next[0] * current[1]);
            }

            return 0.5 * signedArea;
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

    ManifoldCurve2D CurveFactory::sampled_polygon(const std::vector<Point2>& polyVertices, size_t nSegments, bool chamferCorners, bool closeLoop, bool checkOrientation)
    {
        ManifoldCurve2D curve;
        if (polyVertices.size() < 3)
        {
            std::cerr << "Invalid polygon: at least 3 vertices required.\n";
            return curve;
        }

        std::vector<Point2> polyVerticesCopy{ polyVertices };
        if (checkOrientation && CalculatePolyPointsSignedArea(polyVerticesCopy) < 0.0f)
        {
            std::ranges::reverse(polyVerticesCopy);
        }

        // Step 1: Calculate total length of the polygon edges and the mean segment length
        Scalar totalLength = 0.0;
        std::vector<Scalar> polyEdgeLengths(polyVerticesCopy.size());
        for (size_t i = 0; i < polyVerticesCopy.size(); ++i)
        {
            polyEdgeLengths[i] = norm(polyVerticesCopy[(i + 1) % polyVerticesCopy.size()] - polyVerticesCopy[i]);
            totalLength += polyEdgeLengths[i];
        }
        const Scalar meanSegmentLength = totalLength / static_cast<Scalar>(nSegments);

        // Step 2: Calculate the number of segments per edge
        std::vector<size_t> nSegmentsPerEdge;
        nSegmentsPerEdge.reserve(polyVerticesCopy.size());
        size_t totalAssignedSegments = 0;
        for (size_t i = 0; i < polyEdgeLengths.size(); ++i)
        {
            // Compute the number of segments proportional to the edge length
            size_t segments = static_cast<size_t>(std::round(polyEdgeLengths[i] / meanSegmentLength));
            nSegmentsPerEdge.push_back(segments);
            totalAssignedSegments += segments;
        }
        // If there's a difference in total assigned segments (due to rounding), adjust
        if (totalAssignedSegments > nSegments)
        {
            // Reduce the number of segments by removing them from the longest edges
            while (totalAssignedSegments > nSegments)
            {
                auto maxIt = std::ranges::max_element(nSegmentsPerEdge);
                if (*maxIt > 1)
                {
                    (*maxIt)--;
                    totalAssignedSegments--;
                }
            }
        }
        else if (totalAssignedSegments < nSegments)
        {
            // Increase the number of segments by adding them to the longest edges
            while (totalAssignedSegments < nSegments)
            {
                auto maxIt = std::ranges::max_element(nSegmentsPerEdge);
                (*maxIt)++;
                totalAssignedSegments++;
            }
        }

        // Step 3: Tessellate the polygon
        std::vector<Vertex> vertices;
        vertices.reserve(nSegments + (chamferCorners ? polyVerticesCopy.size() : 0));
        const size_t mainPolySegments = polyVerticesCopy.size() - (closeLoop ? 0 : 1);
        for (size_t i = 0; i < mainPolySegments; ++i)
        {
            Point2 startVertex = polyVerticesCopy[i];
            Point2 endVertex = polyVerticesCopy[(i + 1) % polyVerticesCopy.size()];
            size_t segmentsOnEdge = nSegmentsPerEdge[i];

            for (size_t j = 0; j <= segmentsOnEdge; ++j)
            {
                const Scalar t = static_cast<Scalar>(j) / static_cast<Scalar>(segmentsOnEdge + 1);
                Point2 sampledPoint = (1 - t) * startVertex + t * endVertex;

                // Chamfer logic: Skip corner points
                if (chamferCorners && (j == segmentsOnEdge || j == 0))
                {
                    continue;  // Skip the exact corner point
                }

                // Add the sampled point as a vertex
                vertices.push_back(curve.add_vertex(sampledPoint));
            }
        }

        // Step 4: Create edges between consecutive vertices
        for (size_t i = 0; i < vertices.size() - 1; ++i)
        {
            curve.new_edge(vertices[i], vertices[i + 1]);
        }
        if (closeLoop)
            curve.new_edge(vertices.back(), vertices[0]);

        return curve;
    }

    ManifoldCurve2D CurveFactory::hyper_ellipse(
        const Point2& center,
        Scalar radiusX,
        Scalar radiusY,
        size_t degree,
        size_t nSegments,
        Scalar startAngle,
        Scalar endAngle)
    {
        ManifoldCurve2D curve;
        const bool isFullEllipse = std::fabs(endAngle - startAngle - 2.0f * M_PI) < 1e-6f;
        const size_t nVerts = (isFullEllipse ? nSegments : nSegments + 1);
        curve.reserve(nVerts, nSegments);
        std::vector<Vertex> vertices;
        vertices.reserve(nVerts);

        const Scalar angleIncrement = (endAngle - startAngle) / static_cast<Scalar>(nSegments);

        // Loop over the number of segments to generate points along the hyperellipse
        for (size_t i = 0; i < nSegments; ++i)
        {
            const Scalar angle = startAngle + i * angleIncrement;
            const Scalar cosAngle = std::cos(angle);
            const Scalar sinAngle = std::sin(angle);

            // Use the degree to modify the cos and sin values, creating a superellipse shape
            Point2 position(
                center[0] + radiusX * std::pow(std::abs(cosAngle), Scalar(2) / degree) * (cosAngle >= 0 ? 1 : -1),
                center[1] + radiusY * std::pow(std::abs(sinAngle), Scalar(2) / degree) * (sinAngle >= 0 ? 1 : -1)
            );
            vertices.push_back(curve.add_vertex(position));
        }

        // Handle the last vertex
        if (isFullEllipse)
        {
            vertices.push_back(vertices[0]);
        }
        else
        {
            const Scalar angle = endAngle;
            const Scalar cosAngle = std::cos(angle);
            const Scalar sinAngle = std::sin(angle);
            const Point2 position(
                center[0] + radiusX * std::pow(std::abs(cosAngle), Scalar(2) / degree) * (cosAngle >= 0 ? 1 : -1),
                center[1] + radiusY * std::pow(std::abs(sinAngle), Scalar(2) / degree) * (sinAngle >= 0 ? 1 : -1)
            );
            vertices.push_back(curve.add_vertex(position));
        }

        // Create edges between consecutive vertices
        for (size_t i = 0; i < vertices.size() - 1; ++i)
        {
            curve.new_edge(vertices[i], vertices[i + 1]);
        }

        return curve;
    }


} // namespace pmp