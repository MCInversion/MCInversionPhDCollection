#pragma once

#include "pmp/Types.h"
#include "sdf/SDF.h"

#include <vector>

// Forward declarations
namespace Geometry
{
    class VectorGrid2D;
    struct Ray2D;
}

/**
 * \brief A builder for generating distance field characteristics from a planar point cloud using the ray-casting method.
 * \class PlanarPointCloudCharacteristicsBuilder
 */
class PlanarPointCloudCharacteristicsBuilder
{
public:
    /// \brief Constructor
    explicit PlanarPointCloudCharacteristicsBuilder(const std::vector<pmp::Point2>& points, const SDF::PointCloudDistanceField2DSettings& settings)
        : m_Points(points), m_Settings(settings)
    {}

    /// \brief Main method of the builder
    [[nodiscard]] std::vector<std::vector<pmp::Point2>> Build();

private:
    /// \brief Generates initial rays from each point in the point cloud using field directions.
    [[nodiscard]] std::vector<Geometry::Ray2D> GenerateInitialRays(const Geometry::VectorGrid2D& gradDF) const;

    /// \brief Culls rays by limiting their parametric length based on intersections with other rays.
    void CullRays(std::vector<Geometry::Ray2D>& rays) const;

    /// \brief Converts rays to polylines.
    [[nodiscard]] std::vector<std::vector<pmp::Point2>> ConvertRaysToPolylines(const std::vector<Geometry::Ray2D>& rays) const;

    // -----------------------------------------------------------------------------------

    std::vector<pmp::Point2> m_Points; //>! The input point cloud emanating the characteristics
    SDF::PointCloudDistanceField2DSettings m_Settings; //>! Settings for distance field generation
};
