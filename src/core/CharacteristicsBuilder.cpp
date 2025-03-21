
#include "pmp/algorithms/DifferentialGeometry.h"
#include "pmp/algorithms/Normals.h"

#include "geometry/GridUtil.h"
#include "geometry/GeometryUtil.h"

#include "CharacteristicsBuilder.h"

#include <limits>

using namespace Geometry;

namespace
{
    constexpr pmp::Scalar DEFAULT_FAN_STEP_ANGLE{ M_PI / 12.0 };
    constexpr pmp::Scalar ANGLE_STEP_FACTOR{ 7.0 };

    /// \brief Calculates an angle step in radians considering the input curve is a uniformly sampled circle centered in the middle of the input box
    [[nodiscard]] pmp::Scalar CalculateAngularFanSegmentSize(const pmp::BoundingBox2& box, const pmp::ManifoldCurve2D& curve)
    {
        if (curve.n_edges() == 0)
            return DEFAULT_FAN_STEP_ANGLE;

        pmp::Scalar meanEdgeLength = 0.0;
        for (const auto e : curve.edges())
        {
            meanEdgeLength += curve.edge_length(e);
        }
        meanEdgeLength /= curve.n_edges();
        if (meanEdgeLength < FLT_EPSILON)
            return DEFAULT_FAN_STEP_ANGLE;

        const auto size = box.max() - box.min();
        const pmp::Scalar boxCircumference = 2.0 * (size[0] + size[0]);
        const size_t nBoxBoundaryIntersections = std::round(boxCircumference / meanEdgeLength);

        return ANGLE_STEP_FACTOR * 2.0 * M_PI / nBoxBoundaryIntersections;
    }

    void ClipRayByDivergenceField(Ray2D& ray, const pmp::Scalar& divFieldThreshold, const ScalarGrid2D& divField)
    {
        const auto [Nx, Ny] = divField.Dimensions();
        const auto cellSize = divField.CellSize();
        const auto& gridBox = divField.Box();

        // Convert the ray's start point to grid coordinates
        int ix = static_cast<int>((ray.StartPt[0] - gridBox.min()[0]) / cellSize);
        int iy = static_cast<int>((ray.StartPt[1] - gridBox.min()[1]) / cellSize);

        // Check if the starting point is within grid bounds
        if (ix < 0 || ix >= Nx || iy < 0 || iy >= Ny) return;

        // Get ray direction components
        pmp::Scalar dx = ray.Direction[0];
        pmp::Scalar dy = ray.Direction[1];

        // Determine step direction for x and y
        int stepX = (dx > 0) ? 1 : (dx < 0) ? -1 : 0;
        int stepY = (dy > 0) ? 1 : (dy < 0) ? -1 : 0;

        // Calculate initial tMax and tDelta values
        pmp::Scalar tMaxX, tMaxY;
        if (dx != 0) {
            pmp::Scalar nextBoundaryX = gridBox.min()[0] + (ix + (stepX > 0 ? 1 : 0)) * cellSize;
            tMaxX = (nextBoundaryX - ray.StartPt[0]) / dx;
        }
        else {
            tMaxX = std::numeric_limits<pmp::Scalar>::infinity();
        }

        if (dy != 0) {
            pmp::Scalar nextBoundaryY = gridBox.min()[1] + (iy + (stepY > 0 ? 1 : 0)) * cellSize;
            tMaxY = (nextBoundaryY - ray.StartPt[1]) / dy;
        }
        else {
            tMaxY = std::numeric_limits<pmp::Scalar>::infinity();
        }

        pmp::Scalar tDeltaX = (dx != 0) ? cellSize / std::abs(dx) : std::numeric_limits<pmp::Scalar>::infinity();
        pmp::Scalar tDeltaY = (dy != 0) ? cellSize / std::abs(dy) : std::numeric_limits<pmp::Scalar>::infinity();

        // Traverse the grid cells along the ray path
        while (ix >= 0 && ix < Nx && iy >= 0 && iy < Ny)
        {
            // Check divergence value at the current cell
            const size_t cellIndex = ix + iy * Nx;
            if (divField.Values()[cellIndex] < divFieldThreshold)
            {
                // Calculate the world coordinates of the critical grid point
                pmp::Point2 criticalPoint{
                    gridBox.min()[0] + ix * cellSize,
                    gridBox.min()[1] + iy * cellSize
                };

                // Project the critical grid point onto the ray to compute HitParam
                ray.HitParam = dot(criticalPoint - ray.StartPt, ray.Direction);
                return;
            }

            // Move to the next cell in the grid
            if (tMaxX < tMaxY)
            {
                tMaxX += tDeltaX;
                ix += stepX;
            }
            else
            {
                tMaxY += tDeltaY;
                iy += stepY;
            }
        }
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
                rays.emplace_back(point, pmp::vec2(static_cast<pmp::Scalar>(direction[0]), static_cast<pmp::Scalar>(direction[1])));
            }
        }
    }

    return rays;
}

void PlanarCharacteristicsBuilder::CullRays(std::vector<Ray2D>& rays, const pmp::BoundingBox2& clipBox)
{
    // Clip rays by clipBox
    for (auto& ray : rays)
    {
        // Calculate the intersection of the ray with the bounding box
        pmp::Scalar tMin, tMax;
        if (RayBoxIntersection2D(ray.StartPt, ray.Direction, clipBox, tMin, tMax))
        {
            // Update ParamMax if the box intersection parameter is closer than the current value
            ray.ParamMax = std::min(ray.ParamMax, tMax);
            ray.HitParam = std::min(ray.HitParam, ray.ParamMax);
        }
        else
        {
            // If the ray does not intersect the bounding box, mark it as invalid (e.g., set HitParam to 0)
            ray.HitParam = 0.0;
        }
    }

    if (!m_DivergenceField)
    {
        std::cerr << "PlanarCharacteristicsBuilder::CullRays: m_DivergenceField == nullptr!\n";
        return;
    }

    // Clip rays by divergence field
    for (auto& ray : rays)
    {
        ClipRayByDivergenceField(ray, m_Settings.DivFieldThresholdFactor, *m_DivergenceField);
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
    const auto df = SDF::PlanarPointCloudDistanceFieldGenerator::Generate(m_Points, m_Settings.DFSettings);
    auto blurredDf = df;
    ApplyWideGaussianBlur2D(blurredDf);
    const auto gradDF = ComputeNormalizedGradient(blurredDf);
    m_DivergenceField = std::make_shared<ScalarGrid2D>(ComputeDivergenceField(gradDF));

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

std::vector<Ray2D> PlanarManifoldCurveCharacteristicsBuilder::GenerateInitialRays(const pmp::Scalar& fanAngleStep) const
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
            if (!m_Settings.ConstructOutwardCharacteristics)
                continue;

	        // normals of adjacent edges diverge, we need to create a fan of characteristics
            const size_t nAngleSegments = std::round(angleBetweenNormals / fanAngleStep);
            if (nAngleSegments < 2)
                continue;

            // Calculate rotation axis (assuming vectors are in 2D, rotate counterclockwise)
            const pmp::Scalar startAngle = std::atan2(eToNormal[1], eToNormal[0]);
            pmp::Scalar endAngle = std::atan2(eFromNormal[1], eFromNormal[0]);

            // Adjust the angles to handle wrapping around correctly
            if (endAngle < startAngle)
                endAngle += 2.0 * static_cast<pmp::Scalar>(M_PI);

            for (size_t i = 0; i <= nAngleSegments; ++i)
            {
                const pmp::Scalar currentAngle = startAngle + i * fanAngleStep;

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
            const auto normal = (eToNormal + eFromNormal) * 0.5;
            if (m_Settings.ConstructOutwardCharacteristics)
            {
				rays.emplace_back(pos, normal);	            
            }
            if (m_Settings.ConstructInwardCharacteristics)
            {
				rays.emplace_back(pos, -normal);	
            }
        }

        if (m_Settings.ConstructInwardCharacteristics && angleBetweenNormals > FLT_EPSILON && vCurvature < -FLT_EPSILON)
        {
            // Inward characteristics for concave vertices
            const size_t nAngleSegments = std::round(angleBetweenNormals / fanAngleStep);
            if (nAngleSegments < 2)
                continue;

            const pmp::Scalar startAngle = std::atan2(-eFromNormal[1], -eFromNormal[0]);
            pmp::Scalar endAngle = std::atan2(-eToNormal[1], -eToNormal[0]);

            if (endAngle < startAngle)
                endAngle += 2.0 * static_cast<pmp::Scalar>(M_PI);

            for (size_t i = 0; i <= nAngleSegments; ++i)
            {
                const pmp::Scalar currentAngle = startAngle + i * fanAngleStep;
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

    // Compute divergence field
    const auto df = SDF::PlanarPointCloudDistanceFieldGenerator::Generate(m_Curve.positions(), m_Settings.DFSettings);
    auto blurredDf = df;
    ApplyWideGaussianBlur2D(blurredDf);
    const auto gradDF = ComputeNormalizedGradient(blurredDf);
    m_DivergenceField = std::make_shared<ScalarGrid2D>(ComputeDivergenceField(gradDF));

    // calculate fan angle step
    const auto fanAngleStep = CalculateAngularFanSegmentSize(gradDF.Box(), m_Curve);

    // Generate rays
    auto rays = GenerateInitialRays(fanAngleStep);

    // Process intersections
    CullRays(rays, gradDF.Box());

    // Convert rays to polylines and return
    return ConvertRaysToPolylines(rays);
}

//
// =========================================================================================================
//
