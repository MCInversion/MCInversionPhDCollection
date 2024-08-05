#include "gtest/gtest.h"
#include "pmp/ManifoldCurve2D.h"
#include "pmp/algorithms/Remeshing.h"
#include "pmp/algorithms/CurveFactory.h"

#include <numeric>

using namespace pmp;

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
        for (auto v : curve.vertices())
        {
            Point2 p = curve.position(v);
            float angle = atan2(p[1], p[0]);
            float radius = 1.0f + amplitude * sin(period * angle);
            curve.position(v) = Point2(radius * cos(angle), radius * sin(angle));
        }
    }

    [[nodiscard]] std::map<Vertex, float> ComputeLocalCurvature(const ManifoldCurve2D& curve)
    {
        std::map<Vertex, float> curvatures;
        for (auto v : curve.vertices())
        {
            auto [v_prev, v_next] = curve.vertices(v);
            if (!v_prev.is_valid() || !v_next.is_valid()) continue;

            Point2 p = curve.position(v);
            Point2 p_prev = curve.position(v_prev);
            Point2 p_next = curve.position(v_next);

            float curvature = std::abs((p_next[0] - p[0]) * (p[1] - p_prev[1]) -
                (p_next[1] - p[1]) * (p[0] - p_prev[0]));
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
};

TEST_F(CurveRemeshingTests_DeformedCircle, UniformRemeshing_VertexDensity)
{
    // Arrange
    float edgeLength = 0.2f;
    unsigned int iterations = 10;
    CurveRemeshing remesher(curve);

    // Act
    remesher.uniform_remeshing(edgeLength, iterations);

    // Assert
    std::map<Vertex, float> densities = ComputeLocalDensities(curve, edgeLength);
    float meanDensity = ComputeMeanDensity(densities);
    float densityVariance = ComputeDensityVariance(densities, meanDensity);
    EXPECT_NEAR(meanDensity, 1.0f / edgeLength, 0.1f);
    EXPECT_LT(densityVariance, 0.1f);
}

TEST_F(CurveRemeshingTests_DeformedCircle, AdaptiveRemeshing_VertexDensity)
{
    // Arrange
    float edgeLength = 0.2f;
    unsigned int iterations = 10;
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
    std::map<Vertex, float> curvatures = ComputeLocalCurvature(curve);
    std::map<Vertex, float> densities = ComputeLocalDensities(curve, edgeLength);
    float meanCurvature = ComputeMeanCurvature(curvatures);
    float curvatureVariance = ComputeCurvatureVariance(curvatures, meanCurvature);

    float meanDensity = ComputeMeanDensity(densities);
    float densityVariance = ComputeDensityVariance(densities, meanDensity);

    EXPECT_LT(curvatureVariance, 0.1f);
    EXPECT_LT(densityVariance, 0.1f);
}