#include "gtest/gtest.h"
#include "pmp/ManifoldCurve2D.h"
#include "pmp/algorithms/Remeshing.h"
#include "pmp/algorithms/CurveFactory.h"

#include <numeric>

using namespace pmp;

namespace
{
    [[nodiscard]] std::map<Vertex, float> ComputeLocalCurvature(const ManifoldCurve2D& curve)
    {
        std::map<Vertex, float> curvatures;
        for (auto v : curve.vertices())
        {
            auto [vPrev, vNext] = curve.vertices(v);
            if (!vPrev.is_valid() || !vNext.is_valid()) continue;

            Point2 p = curve.position(v);
            Point2 pPrev = curve.position(vPrev);
            Point2 pNext = curve.position(vNext);

            const float curvature = std::abs((pNext[0] - p[0]) * (p[1] - pPrev[1]) -
                (pNext[1] - p[1]) * (p[0] - pPrev[0]));
            curvatures[v] = curvature;
        }
        return curvatures;
    }

    [[nodiscard]] std::map<Vertex, float> ComputeLocalDensities(const ManifoldCurve2D& curve, float radius)
    {
        std::map<Vertex, float> densities;
        for (auto v : curve.vertices())
        {
            int count = 0;
            Point2 p = curve.position(v);
            for (auto u : curve.vertices())
            {
                if (u == v) continue;
                Point2 q = curve.position(u);
                if (distance(p, q) < radius)
                {
                    count++;
                }
            }
            densities[v] = static_cast<float>(count);
        }
        return densities;
    }

    [[nodiscard]] float ComputeMeanDensity(const std::map<Vertex, float>& densities)
    {
        float totalDensity = 0.0f;
        for (const auto& [vertex, density] : densities)
        {
            totalDensity += density;
        }
        return totalDensity / densities.size();
    }

    [[nodiscard]] float ComputeMeanCurvature(const std::map<Vertex, float>& curvatures)
    {
        float totalCurvature = 0.0f;
        for (const auto& [vertex, curvature] : curvatures)
        {
            totalCurvature += curvature;
        }
        return totalCurvature / curvatures.size();
    }

    [[nodiscard]] float ComputeDensityVariance(const std::map<Vertex, float>& densities, float mean_density)
    {
        float variance = 0.0f;
        for (const auto& [vertex, density] : densities)
        {
            variance += (density - mean_density) * (density - mean_density);
        }
        return variance / densities.size();
    }

    [[nodiscard]] float ComputeCurvatureVariance(const std::map<Vertex, float>& curvatures, float mean_curvature)
    {
        float variance = 0.0f;
        for (const auto& [vertex, curvature] : curvatures)
        {
            variance += (curvature - mean_curvature) * (curvature - mean_curvature);
        }
        return variance / curvatures.size();
    }

} // anonymous namespace

class CurveRemeshingTests_DeformedCircle : public ::testing::Test
{
protected:
    ManifoldCurve2D curve;

    void SetUp() override
    {
        curve = CurveFactory::circle(Point2(0, 0), 1.0, 32);
        DeformCircleWithSineWave(M_PI_2, 4);
    }

    void DeformCircleWithSineWave(float amplitude, float period)
    {
        for (const auto v : curve.vertices())
        {
            Point2 p = curve.position(v);
            const float angle = atan2(p[1], p[0]);
            const float radius = 1.0f + amplitude * sin(period * angle);
            curve.position(v) = Point2(radius * cos(angle), radius * sin(angle));
        }
    }
};

TEST_F(CurveRemeshingTests_DeformedCircle, UniformRemeshing_VertexDensity)
{
    // Arrange
    constexpr float edgeLength = 0.2f;
    constexpr unsigned int iterations = 10;
    CurveRemeshing remesher(curve);

    // Act
    remesher.uniform_remeshing(edgeLength, iterations);

    // Assert
    const std::map<Vertex, float> densities = ComputeLocalDensities(curve, edgeLength);
    const float meanDensity = ComputeMeanDensity(densities);
    const float densityVariance = ComputeDensityVariance(densities, meanDensity);
    EXPECT_NEAR(meanDensity, 1.0f / edgeLength, 0.1f);
    EXPECT_LT(densityVariance, 0.1f);
}

TEST_F(CurveRemeshingTests_DeformedCircle, AdaptiveRemeshing_VertexDensity)
{
    // Arrange
    constexpr float edgeLength = 0.2f;
    constexpr unsigned int iterations = 10;
    CurveRemeshing remesher(curve);
    AdaptiveRemeshingSettings settings;
    settings.MinEdgeLength = edgeLength;
    settings.MaxEdgeLength = 2.0f * edgeLength;
    settings.ApproxError = 1.5f * edgeLength;
    settings.NRemeshingIterations = iterations;
    settings.NTangentialSmoothingIters = 6;
    settings.UseProjection = true;

    // Act
    remesher.adaptive_remeshing(settings);

    // Assert
    const std::map<Vertex, float> curvatures = ComputeLocalCurvature(curve);
    const std::map<Vertex, float> densities = ComputeLocalDensities(curve, edgeLength);
    const float meanCurvature = ComputeMeanCurvature(curvatures);
    const float curvatureVariance = ComputeCurvatureVariance(curvatures, meanCurvature);

    const float meanDensity = ComputeMeanDensity(densities);
    const float densityVariance = ComputeDensityVariance(densities, meanDensity);

    EXPECT_LT(curvatureVariance, 0.1f);
    EXPECT_LT(densityVariance, 0.1f);
}