#include "gtest/gtest.h"

#include "pmp/algorithms/CurveFactory.h"
#include "pmp/ManifoldCurve2D.h"

#include "core/ManifoldEvolver.h"
#include "geometry/IcoSphereBuilder.h"

// Suite: ManifoldEvolverTests_ManifoldCurveSuite

TEST(ManifoldEvolverTests_ManifoldCurveSuite, IntersectingInnerAndOuterCircles_InvalidArgThrown)
{
    // Arrange
    pmp::ManifoldCurve2D outerCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 1.0f, 32);
    std::vector<pmp::ManifoldCurve2D> innerCurves;
    innerCurves.push_back(pmp::CurveFactory::circle(pmp::Point2(-0.7f, 0.0f), 0.5f, 32));

    ManifoldEvolutionSettings strategySettings;
    auto customStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
        strategySettings, outerCurve, innerCurves);
    GlobalManifoldEvolutionSettings globalSettings;
    ManifoldEvolver evolver(globalSettings, std::move(customStrategy));

    // Act & Assert
    EXPECT_THROW(evolver.Evolve(), std::invalid_argument);
}

TEST(ManifoldEvolverTests_ManifoldCurveSuite, ShrinkingAndExpandingCircle_NoRemeshing)
{
    // Arrange
    pmp::ManifoldCurve2D outerCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 1.0f, 32);
    std::vector<pmp::ManifoldCurve2D> innerCurves;
    innerCurves.push_back(pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 0.5f, 32));

    ManifoldEvolutionSettings strategySettings;
    strategySettings.TimeStep = 0.01;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 10;
    globalSettings.DoRemeshing = false;

    auto customStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
        strategySettings, outerCurve, innerCurves);

    ManifoldEvolver evolver(globalSettings, std::move(customStrategy));

    // Act
    evolver.Evolve();

    // Assert
    auto strategy = dynamic_cast<CustomManifoldCurveEvolutionStrategy*>(evolver.GetStrategy().get());
    auto resultOuterCurve = strategy->GetOuterCurveInOrigScale();
    auto resultInnerCurves = strategy->GetInnerCurvesInOrigScale();

    ASSERT_TRUE(resultOuterCurve != nullptr);
    ASSERT_FALSE(resultInnerCurves.empty());

    for (int i = 0; i < 32; ++i)
    {
        EXPECT_NEAR(norm(resultOuterCurve->position(pmp::Vertex(i))), std::sqrt(1.0f - 2.0f * (globalSettings.NSteps * strategySettings.TimeStep)), 1e-3f);
        EXPECT_NEAR(norm(resultInnerCurves[0]->position(pmp::Vertex(i))), std::sqrt(0.5f + 2.0f * (globalSettings.NSteps * strategySettings.TimeStep)), 1e-3f);
    }
}

TEST(ManifoldEvolverTests_ManifoldCurveSuite, ShrinkWrappingACirclePointCloud_NoInnerCurveNoRemeshing)
{
    // Arrange
    pmp::ManifoldCurve2D targetCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 0.75f, 32);
    const auto targetPts = targetCurve.positions();

    ManifoldEvolutionSettings strategySettings;
    strategySettings.UseInnerManifolds = false;
    strategySettings.Epsilon = STANDARD_EPSILON;
    strategySettings.Eta = STANDARD_ETA;
    strategySettings.TimeStep = 0.01;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 10;
    globalSettings.DoRemeshing = false;

    auto curveStrategy = std::make_shared<ManifoldCurveEvolutionStrategy>(
        strategySettings, std::make_shared<std::vector<pmp::Point2>>(targetPts));

    ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

    // Act
    evolver.Evolve();

    // Assert
    auto strategy = dynamic_cast<ManifoldCurveEvolutionStrategy*>(evolver.GetStrategy().get());
    auto resultOuterCurve = strategy->GetOuterCurveInOrigScale();
    auto resultInnerCurves = strategy->GetInnerCurvesInOrigScale();

    ASSERT_TRUE(resultOuterCurve != nullptr);
    ASSERT_TRUE(resultInnerCurves.empty());

    for (const auto vPos : resultOuterCurve->positions())
    {
        EXPECT_LT(norm(vPos), 1.0f);
        EXPECT_GT(norm(vPos), 0.75f);
    }
}

// Suite: ManifoldEvolverTests_ManifoldSurfaceSuite

namespace
{
	[[nodiscard]] pmp::SurfaceMesh ConstructIcoSphere(const pmp::Point& center, const pmp::Scalar& radius, const unsigned int& subdiv)
	{
		Geometry::IcoSphereBuilder icoBuilder({ subdiv, radius });
        icoBuilder.BuildBaseData();
        icoBuilder.BuildPMPSurfaceMesh();
        if (center == pmp::Point(0, 0, 0))
            return icoBuilder.GetPMPSurfaceMeshResult();

        auto mesh = icoBuilder.GetPMPSurfaceMeshResult();
        const auto translationMatrix = translation_matrix(center);
        mesh *= translationMatrix;
        return mesh;
	}

} // anonymous namespace

TEST(ManifoldEvolverTests_ManifoldSurfaceSuite, IntersectingInnerAndOuterSpheres_InvalidArgThrown)
{
    // Arrange
    auto outerSurface = ConstructIcoSphere(pmp::Point(0.0f, 0.0f, 0.0f), 1.0f, 2);
    std::vector<pmp::SurfaceMesh> innerSurfaces;
    innerSurfaces.push_back(ConstructIcoSphere(pmp::Point(-0.7f, 0.0f, 0.0f), 0.5f, 2));

    ManifoldEvolutionSettings strategySettings;
    auto customStrategy = std::make_shared<CustomManifoldSurfaceEvolutionStrategy>(
        strategySettings, MeshLaplacian::Barycentric, outerSurface, innerSurfaces);
    GlobalManifoldEvolutionSettings globalSettings;
    ManifoldEvolver evolver(globalSettings, std::move(customStrategy));

    // Act & Assert
    EXPECT_THROW(evolver.Evolve(), std::invalid_argument);
}