// Copyright 2011-2021 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "CurveFactory.h"

namespace pmp
{
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
} // namespace pmp