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
    [[nodiscard]] std::map<Vertex, Scalar> ComputeLocalCurvature(const ManifoldCurve2D& curve)
    {
        std::map<Vertex, Scalar> curvatures;
        for (auto v : curve.vertices())
        {
            auto [vPrev, vNext] = curve.vertices(v);
            if (!vPrev.is_valid() || !vNext.is_valid()) continue;

            Point2 p = curve.position(v);
            Point2 pPrev = curve.position(vPrev);
            Point2 pNext = curve.position(vNext);

            const Scalar curvature = std::abs((pNext[0] - p[0]) * (p[1] - pPrev[1]) -
                (pNext[1] - p[1]) * (p[0] - pPrev[0]));
            curvatures[v] = curvature;
        }
        return curvatures;
    }

    [[nodiscard]] std::map<Vertex, Scalar> ComputeLocalDensities(const ManifoldCurve2D& curve, Scalar edgeLength)
    {
        std::map<Vertex, Scalar> densities;

        for (auto v : curve.vertices())
        {
            Scalar localDensity = 0.0;
            int count = 0;

            // Get the outgoing and incoming edges
            auto [e_out, e_in] = curve.edges(v);

            // Sum the lengths of adjacent edges
            if (e_out.is_valid())
            {
                localDensity += curve.edge_length(e_out);
                count++;
            }

            if (e_in.is_valid())
            {
                localDensity += curve.edge_length(e_in);
                count++;
            }

            // If there are adjacent edges, compute the local density
            if (count > 0)
            {
                localDensity /= count;
                densities[v] = 1.0 / localDensity; // Local density is the inverse of average edge length
            }
            else
            {
                densities[v] = 1.0 / edgeLength; // Default to edgeLength for isolated vertices
            }
        }

        return densities;
    }

    [[nodiscard]] Scalar ComputeMeanDensity(const std::map<Vertex, Scalar>& densities)
    {
        const Scalar totalDensity = std::accumulate(densities.begin(), densities.end(), 0.0,
            [](Scalar sum, const std::pair<Vertex, Scalar>& p) { return sum + p.second; });

        return totalDensity / densities.size();
    }

    [[nodiscard]] Scalar ComputeMeanCurvature(const std::map<Vertex, Scalar>& curvatures)
    {
        Scalar totalCurvature = 0.0;
        for (const auto& [vertex, curvature] : curvatures)
        {
            totalCurvature += curvature;
        }
        return totalCurvature / curvatures.size();
    }

    [[nodiscard]] Scalar ComputeDensityVariance(const std::map<Vertex, Scalar>& densities, Scalar meanDensity)
    {
        Scalar variance = 0.0;
        for (const auto& [v, density] : densities)
        {
            variance += (density - meanDensity) * (density - meanDensity);
        }

        return variance / densities.size();
    }

    [[nodiscard]] Scalar ComputeCurvatureVariance(const std::map<Vertex, Scalar>& curvatures, Scalar meanCurvature)
    {
        Scalar variance = 0.0;
        for (const auto& [vertex, curvature] : curvatures)
        {
            variance += (curvature - meanCurvature) * (curvature - meanCurvature);
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

    void DeformCircleWithSineWave(Scalar amplitude, Scalar freq)
    {
        for (const auto v : curve.vertices())
        {
            Point2 p = curve.position(v) - Point2(0, 0);
            const Scalar angle = atan2(p[1], p[0]);
            vec2 direction = normalize(p);
            curve.position(v) += direction * amplitude * (sin(freq * angle) + 1.0);
        }
    }
};

TEST_F(CurveRemeshingTests_DeformedCircle, UniformRemeshing_VertexDensity)
{
    // Arrange
    constexpr Scalar edgeLength = 0.2;
    constexpr unsigned int iterations = 10;
    CurveRemeshing remesher(curve);

    EXPORT_BEFORE(curve, "CurveRemeshingTests_DeformedCircle_UniformRemeshing");

    // Act
    remesher.uniform_remeshing(edgeLength, iterations);

    EXPORT_AFTER(curve, "CurveRemeshingTests_DeformedCircle_UniformRemeshing");

    // Assert
    const std::map<Vertex, Scalar> densities = ComputeLocalDensities(curve, edgeLength);
    const Scalar meanDensity = ComputeMeanDensity(densities);
    const Scalar densityVariance = ComputeDensityVariance(densities, meanDensity);
    EXPECT_NEAR(1.0 / meanDensity, edgeLength, 0.05);
    EXPECT_LT(densityVariance, 0.1);
}

TEST_F(CurveRemeshingTests_DeformedCircle, AdaptiveRemeshing_VertexDensity)
{
    // Arrange
    constexpr Scalar edgeLength = 0.2;
    constexpr unsigned int iterations = 10;
    CurveRemeshing remesher(curve);
    AdaptiveRemeshingSettings settings;
    settings.MinEdgeLength = edgeLength;
    settings.MaxEdgeLength = 1.5 * edgeLength;
    settings.ApproxError = 0.05 * edgeLength;
    settings.NRemeshingIterations = iterations;
    settings.NTangentialSmoothingIters = 6;
    settings.UseProjection = true;

    // Act
    remesher.adaptive_remeshing(settings);

    EXPORT_AFTER(curve, "CurveRemeshingTests_DeformedCircle_AdaptiveRemeshing");

    // Assert
    const std::map<Vertex, Scalar> curvatures = ComputeLocalCurvature(curve);
    const std::map<Vertex, Scalar> densities = ComputeLocalDensities(curve, edgeLength);
    const Scalar meanCurvature = ComputeMeanCurvature(curvatures);
    const Scalar curvatureVariance = ComputeCurvatureVariance(curvatures, meanCurvature);

    const Scalar meanDensity = ComputeMeanDensity(densities);
    const Scalar densityVariance = ComputeDensityVariance(densities, meanDensity);

    EXPECT_LT(curvatureVariance, 0.1);
    EXPECT_LT(densityVariance, 0.1);
}