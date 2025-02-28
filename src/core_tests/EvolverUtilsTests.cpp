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


}