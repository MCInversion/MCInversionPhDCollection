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
        deform_circle_with_sine_wave(M_PI_2, 4);
    }

    void deform_circle_with_sine_wave(float amplitude, float period)
    {
        for (auto v : curve.vertices())
        {
            Point2 p = curve.position(v);
            float angle = atan2(p[1], p[0]);
            float radius = 1.0f + amplitude * sin(period * angle);
            curve.position(v) = Point2(radius * cos(angle), radius * sin(angle));
        }
    }

    [[nodiscard]] std::map<Vertex, float> compute_local_curvatures(const ManifoldCurve2D& curve)
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

    [[nodiscard]] std::map<Vertex, float> compute_local_densities(const ManifoldCurve2D& curve, float radius)
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

    [[nodiscard]] float compute_mean_density(const std::map<Vertex, float>& densities)
    {
        float total_density = 0.0f;
        for (const auto& [vertex, density] : densities)
        {
            total_density += density;
        }
        return total_density / densities.size();
    }

    [[nodiscard]] float compute_mean_curvature(const std::map<Vertex, float>& curvatures)
    {
        float total_curvature = 0.0f;
        for (const auto& [vertex, curvature] : curvatures)
        {
            total_curvature += curvature;
        }
        return total_curvature / curvatures.size();
    }

    [[nodiscard]] float compute_density_variance(const std::map<Vertex, float>& densities, float mean_density)
    {
        float variance = 0.0f;
        for (const auto& [vertex, density] : densities)
        {
            variance += (density - mean_density) * (density - mean_density);
        }
        return variance / densities.size();
    }

    [[nodiscard]] float compute_curvature_variance(const std::map<Vertex, float>& curvatures, float mean_curvature)
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
    float edge_length = 0.2f;
    unsigned int iterations = 10;
    CurveRemeshing remesher(curve);

    // Act
    remesher.uniform_remeshing(edge_length, iterations);

    // Assert
    std::map<Vertex, float> densities = compute_local_densities(curve, edge_length);
    float mean_density = compute_mean_density(densities);
    float density_variance = compute_density_variance(densities, mean_density);
    EXPECT_NEAR(mean_density, 1.0f / edge_length, 0.1f);
    EXPECT_LT(density_variance, 0.1f);
}

TEST_F(CurveRemeshingTests_DeformedCircle, AdaptiveRemeshing_VertexDensity)
{
    // Arrange
    float edge_length = 0.2f;
    unsigned int iterations = 10;
    CurveRemeshing remesher(curve);
    AdaptiveRemeshingSettings settings;
    settings.MinEdgeLength = edge_length;
    settings.MaxEdgeLength = 2.0f * edge_length;
    settings.ApproxError = 1.5f * edge_length;
    settings.NRemeshingIterations = iterations;
    settings.NTangentialSmoothingIters = 6;
    settings.UseProjection = true;

    // Act
    remesher.adaptive_remeshing(settings);

    // Assert
    std::map<Vertex, float> curvatures = compute_local_curvatures(curve);
    std::map<Vertex, float> densities = compute_local_densities(curve, edge_length);
    float mean_curvature = compute_mean_curvature(curvatures);
    float curvature_variance = compute_curvature_variance(curvatures, mean_curvature);

    float mean_density = compute_mean_density(densities);
    float density_variance = compute_density_variance(densities, mean_density);

    EXPECT_LT(curvature_variance, 0.1f);
    EXPECT_LT(density_variance, 0.1f);
}