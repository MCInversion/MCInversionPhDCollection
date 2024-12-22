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

struct CharacteristicsBuilderSettings
{
    SDF::PointCloudDistanceField2DSettings DFSettings;
    bool ConstructInwardCharacteristics{ false };
    bool ConstructOutwardCharacteristics{ true };
    pmp::Scalar DivFieldThresholdFactor{ 0.01 }; //>! during clipping of rays by divergence field this threshold is used for valuating whether a ray should be clipped.
};

/**
 * \brief A base builder for generating distance field characteristics from a given geometry using the ray-casting method.
 * \class PlanarPointCloudCharacteristicsBuilder
 */
class PlanarCharacteristicsBuilder
{
public:

    explicit PlanarCharacteristicsBuilder(const CharacteristicsBuilderSettings& settings)
        : m_Settings(settings) {}

    /// \brief Main method of the builder
    virtual [[nodiscard]] std::vector<std::vector<pmp::Point2>> Build() = 0;

protected:

    /// \brief Culls rays by limiting their parametric length based on intersections with other rays.
    void CullRays(std::vector<Geometry::Ray2D>& rays, const pmp::BoundingBox2& clipBox);

    /// \brief Converts rays to polylines.
    [[nodiscard]] std::vector<std::vector<pmp::Point2>> ConvertRaysToPolylines(const std::vector<Geometry::Ray2D>& rays) const;

    // -----------------------------------------------------------------------------------

    CharacteristicsBuilderSettings m_Settings;
    std::shared_ptr<Geometry::ScalarGrid2D> m_DivergenceField{ nullptr };
};

/**
 * \brief A builder for generating distance field characteristics from a planar point cloud using the ray-casting method.
 * \class PlanarPointCloudCharacteristicsBuilder
 */
class PlanarPointCloudCharacteristicsBuilder : public PlanarCharacteristicsBuilder
{
public:
    /// \brief Constructor
    explicit PlanarPointCloudCharacteristicsBuilder(const std::vector<pmp::Point2>& points, const CharacteristicsBuilderSettings& settings)
        : PlanarCharacteristicsBuilder(settings), m_Points(points)
    {}

    /// \brief Main method of the builder
    [[nodiscard]] std::vector<std::vector<pmp::Point2>> Build() override;

private:
    /// \brief Generates initial rays from each point in the point cloud using field directions.
    [[nodiscard]] std::vector<Geometry::Ray2D> GenerateInitialRays(const Geometry::VectorGrid2D& gradDF) const;

    // -----------------------------------------------------------------------------------

    std::vector<pmp::Point2> m_Points; //>! The input point cloud emanating the characteristics
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
        const CharacteristicsBuilderSettings& settings)
	    : PlanarCharacteristicsBuilder(settings), m_Curve(curve)
    {
    }

    /// \brief Main method of the builder
    [[nodiscard]] std::vector<std::vector<pmp::Point2>> Build() override;

private:

    /// \brief Generates initial rays from points on the curve with estimated normals.
    [[nodiscard]] std::vector<Geometry::Ray2D> GenerateInitialRays(const pmp::Scalar& fanAngleStep) const;


    // -------------------------------------------------------------------------------------

    pmp::ManifoldCurve2D& m_Curve;
};
