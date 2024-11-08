#pragma once

#include "pmp/Types.h"
#include "pmp/ManifoldCurve2D.h"

#include "sdf/SDF.h"

#include <vector>

// Forward declarations
namespace Geometry
{
    class VectorGrid2D;
    struct Ray2D;
}

/**
 * \brief A base builder for generating distance field characteristics from a given geometry using the ray-casting method.
 * \class PlanarPointCloudCharacteristicsBuilder
 */
class PlanarCharacteristicsBuilder
{
public:
    /// \brief Main method of the builder
    virtual [[nodiscard]] std::vector<std::vector<pmp::Point2>> Build() = 0;

protected:

    /// \brief Culls rays by limiting their parametric length based on intersections with other rays.
    static void CullRays(std::vector<Geometry::Ray2D>& rays, const pmp::BoundingBox2& clipBox);

    /// \brief Converts rays to polylines.
    [[nodiscard]] std::vector<std::vector<pmp::Point2>> ConvertRaysToPolylines(const std::vector<Geometry::Ray2D>& rays) const;
};

/**
 * \brief A builder for generating distance field characteristics from a planar point cloud using the ray-casting method.
 * \class PlanarPointCloudCharacteristicsBuilder
 */
class PlanarPointCloudCharacteristicsBuilder : public PlanarCharacteristicsBuilder
{
public:
    /// \brief Constructor
    explicit PlanarPointCloudCharacteristicsBuilder(const std::vector<pmp::Point2>& points, const SDF::PointCloudDistanceField2DSettings& settings)
        : m_Points(points), m_Settings(settings)
    {}

    /// \brief Main method of the builder
    [[nodiscard]] std::vector<std::vector<pmp::Point2>> Build() override;

private:
    /// \brief Generates initial rays from each point in the point cloud using field directions.
    [[nodiscard]] std::vector<Geometry::Ray2D> GenerateInitialRays(const Geometry::VectorGrid2D& gradDF) const;

    // -----------------------------------------------------------------------------------

    std::vector<pmp::Point2> m_Points; //>! The input point cloud emanating the characteristics
    SDF::PointCloudDistanceField2DSettings m_Settings; //>! Settings for distance field generation
};

/**
 * \brief A builder for generating distance field characteristics from a ManifoldCurve2D using the ray-casting method.
 * \class PlanarManifoldCurveCharacteristicsBuilder
 */
class PlanarManifoldCurveCharacteristicsBuilder : public PlanarCharacteristicsBuilder
{
public:

    explicit PlanarManifoldCurveCharacteristicsBuilder(
        pmp::ManifoldCurve2D& curve, 
        const bool& constructInwardCharacteristics = false,
        const bool& constructOutwardCharacteristics = true,
        const float& boxExpansionFactor = 1.0f)
	    : m_Curve(curve),
	      m_ConstructInwardCharacteristics(constructInwardCharacteristics),
		  m_ConstructOutwardCharacteristics(constructOutwardCharacteristics),
	      m_ExpansionFactor(boxExpansionFactor)
    {
    }

    /// \brief Main method of the builder
    [[nodiscard]] std::vector<std::vector<pmp::Point2>> Build() override;

private:

    /// \brief Generates initial rays from points on the curve with estimated normals.
    [[nodiscard]] std::vector<Geometry::Ray2D> GenerateInitialRays(const float& fanAngleStep) const;


    // -------------------------------------------------------------------------------------

    pmp::ManifoldCurve2D& m_Curve;
    bool m_ConstructInwardCharacteristics;
    bool m_ConstructOutwardCharacteristics;
    float m_ExpansionFactor;
};
