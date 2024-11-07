#include "geometry/GridUtil.h"
#include "geometry/GeometryUtil.h"

#include "CharacteristicsBuilder.h"

#include <limits>
#include <cmath>

using namespace Geometry;

//
// ================================================================================================================
//

std::vector<Ray2D> PlanarPointCloudCharacteristicsBuilder::GenerateInitialRays(const VectorGrid2D& gradDF) const
{
    std::vector<Ray2D> rays;

    for (const auto& point : m_Points)
    {
        // Get the field directions around the current point
        FieldDirectionsOctet fieldDirections = GetFieldDirectionsAroundPoint(point, gradDF);

        // Generate rays from each of the 8 directions
        for (size_t i = 0; i < 8; ++i)
        {
            const auto& direction = fieldDirections.Directions[i];
            if (norm(direction) > std::numeric_limits<double>::epsilon()) // Check for valid direction
            {
                // Create a ray with the current point as the start and the normalized direction
                rays.emplace_back(point, pmp::vec2(static_cast<float>(direction[0]), static_cast<float>(direction[1])));
            }
        }
    }

    return rays;
}

void PlanarPointCloudCharacteristicsBuilder::CullRays(std::vector<Ray2D>& rays) const
{
    // Iterate over each ray and evaluate intersections with other rays to update the HitParam and ParamMax
    for (size_t i = 0; i < rays.size(); ++i)
    {
        for (size_t j = 0; j < rays.size(); ++j)
        {
            if (i == j) continue; // Skip self-comparison

            // Use the operator+= to evaluate and update the parametric distances
            rays[i] += rays[j];
        }
    }
}

std::vector<std::vector<pmp::Point2>> PlanarPointCloudCharacteristicsBuilder::ConvertRaysToPolylines(const std::vector<Ray2D>& rays) const
{
    std::vector<std::vector<pmp::Point2>> polylines;

    for (const auto& ray : rays)
    {
        // Calculate the end point of the ray using the HitParam
        pmp::Point2 endPt = ray.StartPt + ray.HitParam * ray.Direction;

        // Create a polyline containing the start and end points
        std::vector<pmp::Point2> polyline = { ray.StartPt, endPt };
        polylines.push_back(polyline);
    }

    return polylines;
}

//
// ================================================================================================================
//

std::vector<std::vector<pmp::Point2>> PlanarPointCloudCharacteristicsBuilder::Build()
{
    if (m_Points.empty())
    {
        std::cerr << "PlanarPointCloudCharacteristicsBuilder::Build: Cannot calculate characteristics from an empty point cloud!\n";
        return {};
    }

    // Compute fields
    const auto df = SDF::PlanarPointCloudDistanceFieldGenerator::Generate(m_Points, m_Settings);
    const auto gradDF = ComputeNormalizedGradient(df);

    // Compute rays from field directions around m_Points
    auto rays = GenerateInitialRays(gradDF);

    // Process intersections
    CullRays(rays);

    // Convert rays to polylines and return
    return ConvertRaysToPolylines(rays);
}
