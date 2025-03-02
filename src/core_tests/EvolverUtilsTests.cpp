#include "gtest/gtest.h"

#include "pmp/algorithms/CurveFactory.h"
#include "pmp/ManifoldCurve2D.h"

#include "core/EvolverUtilsCommon.h"

#include <filesystem>

// set up root directory
const std::filesystem::path fsRootPath = DROOT_DIR;
const auto fsDataDirPath = fsRootPath / "data\\";
const auto fsDataOutPath = fsRootPath / "output\\";
const std::string dataDirPath = fsDataDirPath.string();
const std::string dataOutPath = fsDataOutPath.string();

TEST(EvolverUtilsTests_InteractionDistanceCollectorSuite, QuadricBlendWithZeroRadiusBetweenValues_SameAsPlainMinimum)
{
	// Arrange 
	const auto minBlendStrategy = GetDistanceBlendStrategy<pmp::dvec2>(DistanceSelectionType::PlainMinimum);
	const auto quadBlendStrategy = GetDistanceBlendStrategy<pmp::dvec2>(DistanceSelectionType::QuadricBlend, 0.0);

	InteractionDistanceCollector<pmp::dvec2> interaction1{ *minBlendStrategy };
	InteractionDistanceCollector<pmp::dvec2> interaction2{ *quadBlendStrategy };

	// Act
	interaction1 << InteractionDistanceRhs<pmp::dvec2>{1.2, pmp::dvec2{ 1.0, 0.0 }};
	interaction2 << InteractionDistanceRhs<pmp::dvec2>{1.2, pmp::dvec2{ 1.0, 0.0 }};
	interaction1 << InteractionDistanceRhs<pmp::dvec2>{0.8, pmp::dvec2{ 0.0, 1.0 }};
	interaction2 << InteractionDistanceRhs<pmp::dvec2>{0.8, pmp::dvec2{ 0.0, 1.0 }};

	// Assert
	EXPECT_EQ(interaction1.Distance, interaction2.Distance);
	EXPECT_EQ(interaction1.NegGradient[0], interaction2.NegGradient[0]);
	EXPECT_EQ(interaction1.NegGradient[1], interaction2.NegGradient[1]);
}

TEST(EvolverUtilsTests_InteractionDistanceCollectorSuite, QuadricBlendWithRadiusBetweenValues_BoundedGradientAndDistanceDifference)
{
	// Arrange 
	const auto minBlendStrategy = GetDistanceBlendStrategy<pmp::dvec2>(DistanceSelectionType::PlainMinimum);
	const auto quadBlendStrategy = GetDistanceBlendStrategy<pmp::dvec2>(DistanceSelectionType::QuadricBlend, 0.7);

	InteractionDistanceCollector<pmp::dvec2> interaction1{ *minBlendStrategy };
	InteractionDistanceCollector<pmp::dvec2> interaction2{ *quadBlendStrategy };

	// Act
	interaction1 << InteractionDistanceRhs<pmp::dvec2>{1.2, pmp::dvec2{ 1.0, 0.0 }};
	interaction2 << InteractionDistanceRhs<pmp::dvec2>{1.2, pmp::dvec2{ 1.0, 0.0 }};
	const auto dist21 = interaction2.Distance;
	interaction1 << InteractionDistanceRhs<pmp::dvec2>{0.8, pmp::dvec2{ 0.0, 1.0 }};
	interaction2 << InteractionDistanceRhs<pmp::dvec2>{0.8, pmp::dvec2{ 0.0, 1.0 }};

	// Assert
	pmp::dvec2 expectedMaxGradDiff{ -1.0, 1.0 };
	constexpr double expectedMaxDistDiff{ 0.4 };
	EXPECT_LT(interaction2.Distance - dist21, expectedMaxDistDiff);
	EXPECT_LT(std::abs(interaction2.NegGradient[0] - interaction1.NegGradient[0]), std::abs(expectedMaxGradDiff[0]));
	EXPECT_LT(std::abs(interaction2.NegGradient[1] - interaction1.NegGradient[1]), std::abs(expectedMaxGradDiff[1]));
}

TEST(EvolverUtilsTests_InteractionDistanceCollectorSuite, QuadricBlendWithOutsideTransitionRegion_SameAsPlainMinimum)
{
	// Arrange 
	const auto minBlendStrategy = GetDistanceBlendStrategy<pmp::dvec2>(DistanceSelectionType::PlainMinimum);
	const auto quadBlendStrategy = GetDistanceBlendStrategy<pmp::dvec2>(DistanceSelectionType::QuadricBlend, 0.7);

	InteractionDistanceCollector<pmp::dvec2> interaction1{ *minBlendStrategy };
	InteractionDistanceCollector<pmp::dvec2> interaction2{ *quadBlendStrategy };

	// Act
	interaction1 << InteractionDistanceRhs<pmp::dvec2>{1.6, pmp::dvec2{ 1.0, 0.0 }};
	interaction2 << InteractionDistanceRhs<pmp::dvec2>{1.6, pmp::dvec2{ 1.0, 0.0 }};
	interaction1 << InteractionDistanceRhs<pmp::dvec2>{0.7, pmp::dvec2{ 0.0, 1.0 }};
	interaction2 << InteractionDistanceRhs<pmp::dvec2>{0.7, pmp::dvec2{ 0.0, 1.0 }};

	// Assert
	EXPECT_EQ(interaction1.Distance, interaction2.Distance);
	EXPECT_EQ(interaction1.NegGradient[0], interaction2.NegGradient[0]);
	EXPECT_EQ(interaction1.NegGradient[1], interaction2.NegGradient[1]);
}