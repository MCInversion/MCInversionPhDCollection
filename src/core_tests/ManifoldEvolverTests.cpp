#include "gtest/gtest.h"

#include "pmp/algorithms/CurveFactory.h"
#include "pmp/ManifoldCurve2D.h"

#include "core/ManifoldEvolver.h"

TEST(ManifoldEvolverTests_ManifoldCurveSuite, ShrinkingAndExpandingCircle_NoRemeshing)
{
    // Arrange
    pmp::ManifoldCurve2D outerCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 1.0f, 32);
    std::vector<pmp::ManifoldCurve2D> innerCurves;
    innerCurves.push_back(pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 0.5f, 32));

    ManifoldEvolutionSettings strategySettings;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 10;
    globalSettings.TimeStep = 0.01;
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
        EXPECT_NEAR(norm(resultOuterCurve->position(pmp::Vertex(i))), std::sqrt(1.0f - 2.0f * (globalSettings.NSteps * globalSettings.TimeStep)), 1e-3f);
        EXPECT_NEAR(norm(resultInnerCurves[0]->position(pmp::Vertex(i))), std::sqrt(0.5f + 2.0f * (globalSettings.NSteps * globalSettings.TimeStep)), 1e-3f);
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
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 10;
    globalSettings.TimeStep = 0.01;
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
