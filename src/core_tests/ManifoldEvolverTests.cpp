#include "gtest/gtest.h"

#include "pmp/algorithms/CurveFactory.h"

#include "core/ManifoldEvolver.h"

TEST(ManifoldEvolverTests, CustomCurveInitialization)
{
    // Arrange
    pmp::ManifoldCurve2D outerCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 1.0f, 32);
    std::vector<pmp::ManifoldCurve2D> innerCurves;
    innerCurves.push_back(pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 0.5f, 32));
    auto customStrategy = std::make_unique<CustomManifoldCurveEvolutionStrategy>(std::move(outerCurve), std::move(innerCurves));

    ManifoldEvolutionSettings settings;
    settings.NSteps = 10;
    ManifoldEvolver evolver(settings, std::move(customStrategy));

    // Act
    evolver.Evolve();

    // Assert
    auto resultOuterCurve = customStrategy->GetResultOuterCurve();
    auto resultInnerCurves = customStrategy->GetResultInnerCurves();
    //EXPECT_LT(resultOuterCurve->ComputeRadius(), 1.0f);
    //EXPECT_GT(resultInnerCurves[0]->ComputeRadius(), 0.5f);
}
