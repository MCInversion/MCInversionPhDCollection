
#include "pmp/algorithms/DifferentialGeometry.h"
#include "pmp/algorithms/Normals.h"

#include "geometry/GridUtil.h"
#include "geometry/GeometryUtil.h"

#include "CharacteristicsBuilder.h"

#include <limits>

using namespace Geometry;

namespace
{
    constexpr float DEFAULT_FAN_STEP_ANGLE{ M_PI / 12.0f };
    constexpr float ANGLE_STEP_FACTOR{ 4.0 };

    /// \brief Calculates an angle step in radians considering the input curve is a uniformly sampled circle centered in the middle of the input box
    [[nodiscard]] float CalculateAngularFanSegmentSize(const pmp::BoundingBox2& box, const pmp::ManifoldCurve2D& curve)
    {
        if (curve.n_edges() == 0)
            return DEFAULT_FAN_STEP_ANGLE;

        float meanEdgeLength = 0.0f;
        for (const auto e : curve.edges())
        {
            meanEdgeLength += curve.edge_length(e);
        }
        meanEdgeLength /= curve.n_edges();
        if (meanEdgeLength < FLT_EPSILON)
            return DEFAULT_FAN_STEP_ANGLE;

        const auto size = box.max() - box.min();
        const float boxCircumference = 2.0f * (size[0] + size[0]);
        const size_t nBoxBoundaryIntersections = std::round(boxCircumference / meanEdgeLength);

        return ANGLE_STEP_FACTOR * 2.0f * M_PI / nBoxBoundaryIntersections;
    }
	
} // anonymous namespace

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

void PlanarCharacteristicsBuilder::CullRays(std::vector<Ray2D>& rays, const pmp::BoundingBox2& clipBox)
{
    // Iterate over each ray and evaluate intersections with other rays to update the HitParam and ParamMax
    //for (size_t i = 0; i < rays.size(); ++i)
    //{
    //    for (size_t j = 0; j < rays.size(); ++j)
    //    {
    //        if (i == j) continue; // Skip self-comparison

    //        // Use the operator+= to evaluate and update the parametric distances
    //        rays[i] += rays[j];
    //    }
    //}

    // Clip rays by clipBox
    for (auto& ray : rays)
    {
        // Calculate the intersection of the ray with the bounding box
        float tMin, tMax;
        if (RayBoxIntersection2D(ray.StartPt, ray.Direction, clipBox, tMin, tMax))
        {
            // Update ParamMax if the box intersection parameter is closer than the current value
            ray.ParamMax = std::min(ray.ParamMax, tMax);
            ray.HitParam = std::min(ray.HitParam, ray.ParamMax);
        }
        else
        {
            // If the ray does not intersect the bounding box, mark it as invalid (e.g., set HitParam to 0)
            ray.HitParam = 0.0f;
        }
    }
}

std::vector<std::vector<pmp::Point2>> PlanarCharacteristicsBuilder::ConvertRaysToPolylines(const std::vector<Ray2D>& rays) const
{
    std::vector<std::vector<pmp::Point2>> polylines;

    for (const auto& ray : rays)
    {
        // Calculate the end point of the ray using the HitParam
        pmp::Point2 endPt = ray.StartPt + ray.HitParam * ray.Direction;

        // Create a polyline containing the start and end points
        std::vector polyline = { ray.StartPt, endPt };
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
    CullRays(rays, gradDF.Box());

    // Convert rays to polylines and return
    return ConvertRaysToPolylines(rays);
}

//
// =========================================================================================================
//

std::vector<Ray2D> PlanarManifoldCurveCharacteristicsBuilder::GenerateInitialRays(const float& fanAngleStep) const
{
    std::vector<Ray2D> rays;

    for (const auto v : m_Curve.vertices())
    {
        const auto vCurvature = vertex_curvature(m_Curve, v);
    	const auto pos = m_Curve.position(v);
        const auto [eTo, eFrom] = m_Curve.edges(v);

        const auto eToNormal = pmp::Normals2::compute_edge_normal(m_Curve, eTo);
        const auto eFromNormal = pmp::Normals2::compute_edge_normal(m_Curve, eFrom);
        const auto angleBetweenNormals = angle(eToNormal, eFromNormal);

        if (angleBetweenNormals > FLT_EPSILON && vCurvature > FLT_EPSILON)
        {
	        // normals of adjacent edges diverge, we need to create a fan of characteristics
            const size_t nAngleSegments = std::round(angleBetweenNormals / fanAngleStep);
            if (nAngleSegments < 2)
                continue;

            // Calculate rotation axis (assuming vectors are in 2D, rotate counterclockwise)
            const float startAngle = std::atan2(eToNormal[1], eToNormal[0]);
            float endAngle = std::atan2(eFromNormal[1], eFromNormal[0]);

            // Adjust the angles to handle wrapping around correctly
            if (endAngle < startAngle)
                endAngle += 2.0f * static_cast<float>(M_PI);

            for (size_t i = 1; i < nAngleSegments; ++i)
            {
                const float currentAngle = startAngle + i * fanAngleStep;

                // Ensure we do not overshoot the end angle
                if (currentAngle > endAngle)
                    break;

                // Convert the angle back to a unit vector (cosine and sine)
                pmp::vec2 currentNormal = pmp::vec2(std::cos(currentAngle), std::sin(currentAngle));

                // Create a ray with the current normal direction
                rays.emplace_back(pos, currentNormal);
            }
        }
        else
        {
            // outward pointing normals of adjacent edges converge, one characteristic suffices
            const auto normal = (eToNormal + eFromNormal) * 0.5f;
	        rays.emplace_back(pos, normal);
            if (m_ConstructInwardCharacteristics)
            {
				rays.emplace_back(pos, -normal);	
            }
        }

        if (m_ConstructInwardCharacteristics && angleBetweenNormals > FLT_EPSILON && vCurvature < -FLT_EPSILON)
        {
            // Inward characteristics for concave vertices
            const size_t nAngleSegments = std::round(fanAngleStep / std::abs(angleBetweenNormals));
            if (nAngleSegments < 2)
                continue;

            const float startAngle = std::atan2(-eToNormal[1], -eToNormal[0]);
            float endAngle = std::atan2(-eFromNormal[1], -eFromNormal[0]);

            if (endAngle < startAngle)
                endAngle += 2.0f * static_cast<float>(M_PI);

            for (size_t i = 1; i < nAngleSegments; ++i)
            {
                const float currentAngle = startAngle + i * fanAngleStep;
                if (currentAngle > endAngle)
                    break;

                pmp::vec2 currentNormal = pmp::vec2(std::cos(currentAngle), std::sin(currentAngle));
                rays.emplace_back(pos, currentNormal);
            }
        }
    }

    return rays;
}

//
// ---------------------------------------------------------------------------------------------------------
//

std::vector<std::vector<pmp::Point2>> PlanarManifoldCurveCharacteristicsBuilder::Build()
{
    if (m_Curve.is_empty())
    {
        std::cerr << "PlanarManifoldCurveCharacteristicsBuilder::Build: Cannot calculate characteristics from an empty ManifoldCurve2D!\n";
        return {};
    }

    // construct bounding volume
    auto bbox = m_Curve.bounds();
    const auto size = bbox.max() - bbox.min();
    const float minSize = std::min(size[0], size[1]);
    if (m_ExpansionFactor > 0.0f)
    {
        const float expansion = m_ExpansionFactor * minSize;
        bbox.expand(expansion, expansion);
    }
    const auto fanAngleStep = CalculateAngularFanSegmentSize(bbox, m_Curve);

    // Generate rays
    auto rays = GenerateInitialRays(fanAngleStep);

    // Process intersections
    CullRays(rays, bbox);

    // Convert rays to polylines and return
    return ConvertRaysToPolylines(rays);
}

//
// =========================================================================================================
//
