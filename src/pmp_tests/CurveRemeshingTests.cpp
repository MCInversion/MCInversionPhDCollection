#include "gtest/gtest.h"
#include "pmp/ManifoldCurve2D.h"
#include "pmp/algorithms/Remeshing.h"
#include "pmp/algorithms/CurveFactory.h"

#include <numeric>

using namespace pmp;

#define EXPORT_TEST_BEFORE_AFTER_CURVES true

#if EXPORT_TEST_BEFORE_AFTER_CURVES
#include <filesystem>
// set up root directory
const std::filesystem::path fsRootPath = DROOT_DIR;
const auto fsDataDirPath = fsRootPath / "data\\";
const auto fsDataOutPath = fsRootPath / "output\\";
const std::string dataDirPath = fsDataDirPath.string();
const std::string dataOutPath = fsDataOutPath.string();
#define EXPORT_BEFORE(curve, name) \
    do { \
        const auto fileName = dataOutPath + "/pmp_tests/" + (name) + "_before.ply"; \
        if (!write_to_ply((curve), fileName)) { \
            EXPECT_TRUE(false); \
        } \
    } while (0)

#define EXPORT_AFTER(curve, name) \
    do { \
        const auto fileName = dataOutPath + "/pmp_tests/" + (name) + "_after.ply"; \
        if (!write_to_ply((curve), fileName)) { \
            EXPECT_TRUE(false); \
        } \
    } while (0)

#else

#define EXPORT_BEFORE(curve, name) do {} while (0)
#define EXPORT_AFTER(curve, name) do {} while (0)

#endif

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
        DeformCircleWithSineWave(1.0, 4);
    }

    void DeformCircleWithSineWave(float amplitude, float freq)
    {
        for (const auto v : curve.vertices())
        {
            Point2 p = curve.position(v) - Point2(0, 0);
            const float angle = atan2(p[1], p[0]);
            vec2 direction = normalize(p);
            curve.position(v) += direction * amplitude * (sin(freq * angle) + 1.0f);
        }
    }
};

TEST_F(CurveRemeshingTests_DeformedCircle, UniformRemeshing_VertexDensity)
{
    // Arrange
    constexpr float edgeLength = 0.2f;
    constexpr unsigned int iterations = 10;
    CurveRemeshing remesher(curve);

    EXPORT_BEFORE(curve, "CurveRemeshingTests_DeformedCircle_UniformRemeshing");

    // Act
    remesher.uniform_remeshing(edgeLength, iterations);

    EXPORT_AFTER(curve, "CurveRemeshingTests_DeformedCircle_UniformRemeshing");

    // Assert
    const std::map<Vertex, float> densities = ComputeLocalDensities(curve, edgeLength);
    const float meanDensity = ComputeMeanDensity(densities);
    const float densityVariance = ComputeDensityVariance(densities, meanDensity);
    EXPECT_NEAR(meanDensity, 1.0f / edgeLength, 1e-3f);
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
    settings.MaxEdgeLength = 1.5f * edgeLength;
    settings.ApproxError = 0.05f * edgeLength;
    settings.NRemeshingIterations = iterations;
    settings.NTangentialSmoothingIters = 6;
    settings.UseProjection = true;

    // Act
    remesher.adaptive_remeshing(settings);

    EXPORT_AFTER(curve, "CurveRemeshingTests_DeformedCircle_AdaptiveRemeshing");

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