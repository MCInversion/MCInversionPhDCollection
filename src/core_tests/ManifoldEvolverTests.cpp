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
    auto customStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(std::move(outerCurve), std::move(innerCurves));

    ManifoldEvolutionSettings settings;
    settings.NSteps = 10;
    settings.TimeStep = 0.01;
    settings.DoRemeshing = false;
    ManifoldEvolver evolver(settings, std::move(customStrategy));

    // Act
    evolver.Evolve();

    // Assert
    auto strategy = dynamic_cast<CustomManifoldCurveEvolutionStrategy*>(evolver.GetStrategy().get());
    auto& resultOuterCurve = strategy->GetOuterCurve();
    auto& resultInnerCurves = strategy->GetInnerCurves();
    ASSERT_TRUE(resultOuterCurve != nullptr);
    ASSERT_FALSE(resultInnerCurves.empty());
    for (int i = 0; i < 32; ++i)
    {
        EXPECT_NEAR(norm(resultOuterCurve->position(pmp::Vertex(i))), std::sqrt(1.0f - 2.0f * (settings.NSteps * settings.TimeStep)), 1e-3f);
        EXPECT_NEAR(norm(resultInnerCurves[0]->position(pmp::Vertex(i))), std::sqrt(0.5f + 2.0f * (settings.NSteps * settings.TimeStep)), 1e-3f);
    }
}

TEST(ManifoldEvolverTests_ManifoldCurveSuite, ShrinkWrappingACirclePointCloud_NoInnerCurveNoRemeshing)
{
    // Arrange
    pmp::ManifoldCurve2D targetCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 0.75f, 32);
    const auto targetPts = targetCurve.positions();
    auto curveStrategy = std::make_shared<ManifoldCurveEvolutionStrategy>(std::make_shared<std::vector<pmp::Point2>>(targetPts));

    ManifoldEvolutionSettings settings;
    settings.NSteps = 10;
    settings.TimeStep = 0.01;
    settings.DoRemeshing = false;
    ManifoldEvolver evolver(settings, std::move(curveStrategy));

    // Act
    evolver.Evolve();

    // Assert
    auto strategy = dynamic_cast<ManifoldCurveEvolutionStrategy*>(evolver.GetStrategy().get());
    auto& resultOuterCurve = strategy->GetOuterCurve();
    auto& resultInnerCurves = strategy->GetInnerCurves();
    ASSERT_TRUE(resultOuterCurve != nullptr);
    ASSERT_TRUE(resultInnerCurves.empty());
    for (int i = 0; i < 32; ++i)
    {
        EXPECT_LT(norm(resultOuterCurve->position(pmp::Vertex(i))), 1.0f);
        EXPECT_GT(norm(resultOuterCurve->position(pmp::Vertex(i))), 0.75f);
    }
}
