#include "gtest/gtest.h"

#include "pmp/algorithms/CurveFactory.h"

#include "core/ManifoldEvolver.h"

TEST(ManifoldEvolverTests, CustomCurveInitialization)
{
    // Create a custom outer curve (e.g., a circle)
    pmp::ManifoldCurve2D outerCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 1.0f, 32);

    // Create custom inner curves (e.g., smaller circles)
    std::vector<pmp::ManifoldCurve2D> innerCurves;
    innerCurves.push_back(pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 0.5f, 32));

    // Create the custom strategy with these curves
    auto customStrategy = std::make_unique<CustomManifoldCurveEvolutionStrategy>(std::move(outerCurve), std::move(innerCurves));

    // Create settings
    ManifoldEvolutionSettings settings;
    settings.NSteps = 10;

    // Create the evolver
    ManifoldEvolver evolver(settings, std::move(customStrategy));

    // Run the evolution
    evolver.Evolve();

    // Get the results from the strategy
    auto resultOuterCurve = customStrategy->GetResultOuterCurve();
    auto resultInnerCurves = customStrategy->GetResultInnerCurves();

    // Assert the expected results (e.g., shrinking outer curve and expanding inner curves)
    //EXPECT_LT(resultOuterCurve->ComputeRadius(), 1.0f);
    //EXPECT_GT(resultInnerCurves[0]->ComputeRadius(), 0.5f);
}
