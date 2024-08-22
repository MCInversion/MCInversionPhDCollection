
#include "gtest/gtest.h"

#include "pmp/algorithms/CurveFactory.h"
#include "pmp/ManifoldCurve2D.h"

#include "core/ManifoldEvolver.h"
#include "geometry/IcoSphereBuilder.h"

#include <filesystem>

// set up root directory
const std::filesystem::path fsRootPath = DROOT_DIR;
const auto fsDataDirPath = fsRootPath / "data\\";
const auto fsDataOutPath = fsRootPath / "output\\";
const std::string dataDirPath = fsDataDirPath.string();
const std::string dataOutPath = fsDataOutPath.string();

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

TEST(ManifoldEvolverTests_ManifoldCurveSuite, ShrinkingCircle_SemiImplicitNoRemeshing)
{
    // Arrange
    pmp::ManifoldCurve2D outerCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 1.0f, 32);
    std::vector<pmp::ManifoldCurve2D> innerCurves;

    ManifoldEvolutionSettings strategySettings;
    strategySettings.TimeStep = 0.01;
    strategySettings.UseInnerManifolds = false;
    strategySettings.TangentialVelocityWeight = 0.0;
    strategySettings.UseStabilizationViaScaling = false;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 10;
    globalSettings.DoRemeshing = false;
    globalSettings.ExportPerTimeStep = true;
    globalSettings.ProcedureName = "ShrinkingCircle_SemiImplicitNoRemeshing";
    globalSettings.OutputPath = dataOutPath + "core_tests\\";
    globalSettings.ExportResult = false;

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
    ASSERT_TRUE(resultInnerCurves.empty());

    for (int i = 0; i < 32; ++i)
    {
        EXPECT_NEAR(norm(resultOuterCurve->position(pmp::Vertex(i))), std::sqrt(1.0f - 2.0f * (globalSettings.NSteps * strategySettings.TimeStep)), 1e-2f);
    }
}

TEST(ManifoldEvolverTests_ManifoldCurveSuite, ShrinkingCircle_ExplicitNoRemeshing)
{
    // Arrange
    pmp::ManifoldCurve2D outerCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 1.0f, 32);
    std::vector<pmp::ManifoldCurve2D> innerCurves;

    ManifoldEvolutionSettings strategySettings;
    strategySettings.TimeStep = 0.005;  // smaller time step size due to instability
    strategySettings.UseInnerManifolds = false;
    strategySettings.TangentialVelocityWeight = 0.0;
    strategySettings.UseStabilizationViaScaling = false;
    strategySettings.UseSemiImplicit = false;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 6;  // fewer time steps due to instability
    globalSettings.DoRemeshing = false;
    globalSettings.ExportPerTimeStep = true;
    globalSettings.ProcedureName = "ShrinkingCircle_ExplicitNoRemeshing";
    globalSettings.OutputPath = dataOutPath + "core_tests\\";
    globalSettings.ExportResult = false;

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
    ASSERT_TRUE(resultInnerCurves.empty());

    for (int i = 0; i < 32; ++i)
    {
        EXPECT_NEAR(norm(resultOuterCurve->position(pmp::Vertex(i))), std::sqrt(1.0f - 2.0f * (globalSettings.NSteps * strategySettings.TimeStep)), 1e-2f);
    }
}

TEST(ManifoldEvolverTests_ManifoldCurveSuite, ShrinkingAndExpandingCircle_NoRemeshing)
{
    // Arrange
    pmp::ManifoldCurve2D outerCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 1.0f, 32);
    std::vector<pmp::ManifoldCurve2D> innerCurves;
    innerCurves.push_back(pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 0.5f, 16));

    ManifoldEvolutionSettings strategySettings;
    strategySettings.TimeStep = 0.01;
    strategySettings.OuterManifoldEpsilon = STANDARD_EPSILON;
    strategySettings.OuterManifoldEta = STANDARD_ETA;
    strategySettings.InnerManifoldEpsilon = strategySettings.OuterManifoldEpsilon;
    strategySettings.InnerManifoldEta = strategySettings.OuterManifoldEta;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 10;
    globalSettings.DoRemeshing = false;
    globalSettings.ExportPerTimeStep = true;
    globalSettings.ProcedureName = "ShrinkingAndExpandingCircle_NoRemeshing";
    globalSettings.OutputPath = dataOutPath + "core_tests\\";
    globalSettings.ExportResult = false;

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
        EXPECT_LT(norm(resultOuterCurve->position(pmp::Vertex(i))), 1.0f);
    }
    for (int i = 0; i < 16; ++i)
    {
	    EXPECT_GT(norm(resultInnerCurves[0]->position(pmp::Vertex(i))), 0.75f);
    }
}

TEST(ManifoldEvolverTests_ManifoldCurveSuite, ShrinkingAndExpandingCircle_SlowerInnerEtaNoRemeshing)
{
    // Arrange
    pmp::ManifoldCurve2D outerCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 1.0f, 32);
    std::vector<pmp::ManifoldCurve2D> innerCurves;
    innerCurves.push_back(pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 0.5f, 16));

    ManifoldEvolutionSettings strategySettings;
    strategySettings.TimeStep = 0.01;
    strategySettings.OuterManifoldEpsilon = STANDARD_EPSILON;
    strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
    {
        if (distance >= Geometry::DEFAULT_SCALAR_GRID_INIT_VAL)
            return 0.0;
        return 2.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
    };
    strategySettings.InnerManifoldEpsilon = strategySettings.OuterManifoldEpsilon;
    strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
    {
        if (distance >= Geometry::DEFAULT_SCALAR_GRID_INIT_VAL)
            return 0.0;
        return 0.5 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
    };
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 10;
    globalSettings.DoRemeshing = false;
    globalSettings.ExportPerTimeStep = true;
    globalSettings.ProcedureName = "ShrinkingAndExpandingCircle_SlowerInnerEtaNoRemeshing";
    globalSettings.OutputPath = dataOutPath + "core_tests\\";
    globalSettings.ExportResult = false;

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
        EXPECT_LT(norm(resultOuterCurve->position(pmp::Vertex(i))), 1.0f);
    }
    for (int i = 0; i < 16; ++i)
    {
        EXPECT_GT(norm(resultInnerCurves[0]->position(pmp::Vertex(i))), 0.75f);
    }
}

TEST(ManifoldEvolverTests_ManifoldCurveSuite, ShrinkWrappingACirclePointCloud_NoInnerCurveNoRemeshing)
{
    // Arrange
    pmp::ManifoldCurve2D targetCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 0.75f, 32);
    const auto targetPts = targetCurve.positions();

    ManifoldEvolutionSettings strategySettings;
    strategySettings.UseInnerManifolds = false;
    strategySettings.OuterManifoldEpsilon = STANDARD_EPSILON;
    strategySettings.OuterManifoldEta = STANDARD_ETA;
    strategySettings.TimeStep = 0.02;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 10;
    globalSettings.DoRemeshing = false;
    globalSettings.ExportPerTimeStep = true;
    globalSettings.ExportTargetDistanceFieldAsImage = true;
    globalSettings.ProcedureName = "ShrinkWrappingACirclePointCloud_NoInnerCurveNoRemeshing";
    globalSettings.OutputPath = dataOutPath + "core_tests\\";
    globalSettings.ExportResult = false;

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
        EXPECT_LT(norm(vPos), 2.0f);
        EXPECT_GT(norm(vPos), 1.0f);
    }
}

TEST(ManifoldEvolverTests_ManifoldCurveSuite, ShrinkWrappingAnIncompleteCirclePointCloud_NoInnerCircleNoRemeshing)
{
    // Arrange
    pmp::ManifoldCurve2D targetCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 0.75f, 16);
    auto targetPts = targetCurve.positions();
    targetPts.erase(targetPts.begin());
    targetPts.erase(targetPts.begin());
    targetPts.erase(targetPts.begin());
    targetPts.erase(targetPts.begin());

    ManifoldEvolutionSettings strategySettings;
    strategySettings.UseInnerManifolds = false;
    strategySettings.OuterManifoldEpsilon = STANDARD_EPSILON;
    strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
    {
        if (distance >= Geometry::DEFAULT_SCALAR_GRID_INIT_VAL)
            return 0.0;
        return 2.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
    };
    strategySettings.TimeStep = 0.01;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 40;
    globalSettings.DoRemeshing = false;
    globalSettings.ExportPerTimeStep = true;
    globalSettings.ExportTargetDistanceFieldAsImage = true;
    globalSettings.ProcedureName = "ShrinkWrappingAnIncompleteCirclePointCloud_NoInnerCircleNoRemeshing";
    globalSettings.OutputPath = dataOutPath + "core_tests\\";
    globalSettings.ExportResult = false;

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

    size_t nInvalidVertices = 0;
    for (const auto vPos : resultOuterCurve->positions())
    {
        if (norm(vPos) > 2.0f)
            nInvalidVertices++;
        else if (norm(vPos) < 1.0f)
            nInvalidVertices++;
    }
    EXPECT_LT(static_cast<double>(nInvalidVertices) / resultOuterCurve->n_vertices(), 0.25);
}

TEST(ManifoldEvolverTests_ManifoldCurveSuite, ShrinkWrappingAnIncompleteCirclePointCloud_NoRemeshing)
{
    // Arrange
    pmp::ManifoldCurve2D targetCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 0.75f, 16);
    auto targetPts = targetCurve.positions();
    targetPts.erase(targetPts.begin());
    targetPts.erase(targetPts.begin());
    targetPts.erase(targetPts.begin());
    targetPts.erase(targetPts.begin());

    ManifoldEvolutionSettings strategySettings;
    strategySettings.UseInnerManifolds = true;
    strategySettings.OuterManifoldEpsilon = STANDARD_EPSILON;
    strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
    {
        if (distance >= Geometry::DEFAULT_SCALAR_GRID_INIT_VAL)
            return 0.0;
        return 2.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
    };
    strategySettings.InnerManifoldEpsilon = strategySettings.OuterManifoldEpsilon;
    strategySettings.InnerManifoldEta = strategySettings.OuterManifoldEta;
    strategySettings.TimeStep = 0.01;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 40;
    globalSettings.DoRemeshing = false;
    globalSettings.ExportPerTimeStep = true;
    globalSettings.ExportTargetDistanceFieldAsImage = true;
    globalSettings.ProcedureName = "ShrinkWrappingAnIncompleteCirclePointCloud_NoRemeshing";
    globalSettings.OutputPath = dataOutPath + "core_tests\\";
    globalSettings.ExportResult = false;

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
    ASSERT_FALSE(resultInnerCurves.empty());

    size_t nInvalidVertices = 0;
    for (const auto vPos : resultOuterCurve->positions())
    {
        if (norm(vPos) > 2.0f)
            nInvalidVertices++;
        else if (norm(vPos) < 1.0f)
            nInvalidVertices++;
    }
    EXPECT_LT(static_cast<double>(nInvalidVertices) / resultOuterCurve->n_vertices(), 0.0);
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

TEST(ManifoldEvolverTests_ManifoldSurfaceSuite, ShrinkingSphere_BarycentricSemiImplicitNoRemeshing)
{
    // Arrange
    auto outerSurface = ConstructIcoSphere(pmp::Point(0.0f, 0.0f, 0.0f), 1.0f, 2);
    std::vector<pmp::SurfaceMesh> innerSurfaces;

    ManifoldEvolutionSettings strategySettings;
    strategySettings.TimeStep = 0.01;
    strategySettings.UseInnerManifolds = false;
    strategySettings.TangentialVelocityWeight = 0.0;
    strategySettings.UseStabilizationViaScaling = false;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 10;
    globalSettings.DoRemeshing = false;
    globalSettings.ExportPerTimeStep = true;
    globalSettings.ProcedureName = "ShrinkingSphere_BarycentricSemiImplicitNoRemeshing";
    globalSettings.OutputPath = dataOutPath + "core_tests\\";
    globalSettings.ExportResult = false;

    auto customStrategy = std::make_shared<CustomManifoldSurfaceEvolutionStrategy>(
        strategySettings, MeshLaplacian::Barycentric, outerSurface, innerSurfaces);
    ManifoldEvolver evolver(globalSettings, std::move(customStrategy));

    // Act
    evolver.Evolve();

    // Assert
    auto strategy = dynamic_cast<CustomManifoldSurfaceEvolutionStrategy*>(evolver.GetStrategy().get());
    auto resultOuterSurface = strategy->GetOuterSurfaceInOrigScale();
    auto resultInnerSurfaces = strategy->GetInnerSurfacesInOrigScale();

    ASSERT_TRUE(resultOuterSurface != nullptr);
    ASSERT_TRUE(resultInnerSurfaces.empty());

    for (int i = 0; i < 32; ++i)
    {
        EXPECT_NEAR(norm(resultOuterSurface->position(pmp::Vertex(i))), std::sqrt(1.0f - 4.0f * (globalSettings.NSteps * strategySettings.TimeStep)), 1e-2f);
    }
}

TEST(ManifoldEvolverTests_ManifoldSurfaceSuite, ShrinkingSphere_BarycentricExplicitNoRemeshing)
{
    // Arrange
    auto outerSurface = ConstructIcoSphere(pmp::Point(0.0f, 0.0f, 0.0f), 1.0f, 2);
    std::vector<pmp::SurfaceMesh> innerSurfaces;

    ManifoldEvolutionSettings strategySettings;
    strategySettings.TimeStep = 0.005; // smaller time step size due to instability
    strategySettings.UseSemiImplicit = false;
    strategySettings.UseInnerManifolds = false;
    strategySettings.TangentialVelocityWeight = 0.0;
    strategySettings.UseStabilizationViaScaling = false;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 6; // fewer time steps due to instability
    globalSettings.DoRemeshing = false;
    globalSettings.ExportPerTimeStep = true;
    globalSettings.ProcedureName = "ShrinkingSphere_BarycentricExplicitNoRemeshing";
    globalSettings.OutputPath = dataOutPath + "core_tests\\";
    globalSettings.ExportResult = false;

    auto customStrategy = std::make_shared<CustomManifoldSurfaceEvolutionStrategy>(
        strategySettings, MeshLaplacian::Barycentric, outerSurface, innerSurfaces);
    ManifoldEvolver evolver(globalSettings, std::move(customStrategy));

    // Act
    evolver.Evolve();

    // Assert
    auto strategy = dynamic_cast<CustomManifoldSurfaceEvolutionStrategy*>(evolver.GetStrategy().get());
    auto resultOuterSurface = strategy->GetOuterSurfaceInOrigScale();
    auto resultInnerSurfaces = strategy->GetInnerSurfacesInOrigScale();

    ASSERT_TRUE(resultOuterSurface != nullptr);
    ASSERT_TRUE(resultInnerSurfaces.empty());

    for (int i = 0; i < 32; ++i)
    {
        EXPECT_NEAR(norm(resultOuterSurface->position(pmp::Vertex(i))), std::sqrt(1.0f - 4.0f * (globalSettings.NSteps * strategySettings.TimeStep)), 1e-2f);
    }
}

TEST(ManifoldEvolverTests_ManifoldSurfaceSuite, ShrinkingSphere_VoronoiSemiImplicitNoRemeshing)
{
    // Arrange
    auto outerSurface = ConstructIcoSphere(pmp::Point(0.0f, 0.0f, 0.0f), 1.0f, 2);
    std::vector<pmp::SurfaceMesh> innerSurfaces;

    ManifoldEvolutionSettings strategySettings;
    strategySettings.TimeStep = 0.01;
    strategySettings.UseInnerManifolds = false;
    strategySettings.TangentialVelocityWeight = 0.0;
    strategySettings.UseStabilizationViaScaling = false;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 10;
    globalSettings.DoRemeshing = false;
    globalSettings.ExportPerTimeStep = true;
    globalSettings.ProcedureName = "ShrinkingSphere_VoronoiSemiImplicitNoRemeshing";
    globalSettings.OutputPath = dataOutPath + "core_tests\\";
    globalSettings.ExportResult = false;

    auto customStrategy = std::make_shared<CustomManifoldSurfaceEvolutionStrategy>(
        strategySettings, MeshLaplacian::Voronoi, outerSurface, innerSurfaces);
    ManifoldEvolver evolver(globalSettings, std::move(customStrategy));

    // Act
    evolver.Evolve();

    // Assert
    auto strategy = dynamic_cast<CustomManifoldSurfaceEvolutionStrategy*>(evolver.GetStrategy().get());
    auto resultOuterSurface = strategy->GetOuterSurfaceInOrigScale();
    auto resultInnerSurfaces = strategy->GetInnerSurfacesInOrigScale();

    ASSERT_TRUE(resultOuterSurface != nullptr);
    ASSERT_TRUE(resultInnerSurfaces.empty());

    for (int i = 0; i < 32; ++i)
    {
        EXPECT_NEAR(norm(resultOuterSurface->position(pmp::Vertex(i))), std::sqrt(1.0f - 4.0f * (globalSettings.NSteps * strategySettings.TimeStep)), 1e-2f);
    }
}

TEST(ManifoldEvolverTests_ManifoldSurfaceSuite, ShrinkingSphere_VoronoiExplicitNoRemeshing)
{
    // Arrange
    auto outerSurface = ConstructIcoSphere(pmp::Point(0.0f, 0.0f, 0.0f), 1.0f, 2);
    std::vector<pmp::SurfaceMesh> innerSurfaces;

    ManifoldEvolutionSettings strategySettings;
    strategySettings.TimeStep = 0.005; // smaller time step size due to instability
    strategySettings.UseSemiImplicit = false;
    strategySettings.UseInnerManifolds = false;
    strategySettings.TangentialVelocityWeight = 0.0;
    strategySettings.UseStabilizationViaScaling = false;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 6; // fewer time steps due to instability
    globalSettings.DoRemeshing = false;
    globalSettings.ExportPerTimeStep = true;
    globalSettings.ProcedureName = "ShrinkingSphere_VoronoiExplicitNoRemeshing";
    globalSettings.OutputPath = dataOutPath + "core_tests\\";
    globalSettings.ExportResult = false;

    auto customStrategy = std::make_shared<CustomManifoldSurfaceEvolutionStrategy>(
        strategySettings, MeshLaplacian::Barycentric, outerSurface, innerSurfaces);
    ManifoldEvolver evolver(globalSettings, std::move(customStrategy));

    // Act
    evolver.Evolve();

    // Assert
    auto strategy = dynamic_cast<CustomManifoldSurfaceEvolutionStrategy*>(evolver.GetStrategy().get());
    auto resultOuterSurface = strategy->GetOuterSurfaceInOrigScale();
    auto resultInnerSurfaces = strategy->GetInnerSurfacesInOrigScale();

    ASSERT_TRUE(resultOuterSurface != nullptr);
    ASSERT_TRUE(resultInnerSurfaces.empty());

    for (int i = 0; i < 32; ++i)
    {
        EXPECT_NEAR(norm(resultOuterSurface->position(pmp::Vertex(i))), std::sqrt(1.0f - 4.0f * (globalSettings.NSteps * strategySettings.TimeStep)), 1e-2f);
    }
}

TEST(ManifoldEvolverTests_ManifoldSurfaceSuite, ShrinkWrappingAnIncompleteSurfacePointCloud_NoInnerSurfaceNoRemeshing)
{
    // Arrange
    auto targetIcoSphere = ConstructIcoSphere(pmp::Point(0.0f, 0.0f, 0.0f), 0.75f, 2);
    auto targetPts = targetIcoSphere.positions();
    targetPts.erase(std::remove_if(targetPts.begin(), targetPts.end(),
        [](const pmp::Point& p) {
            // Filter points in the x > 0 && y > 0 && z > 0 octant
            return (p[0] > 0.0f && p[1] > 0.0f && p[2] > 0.0f);
        }),
    targetPts.end());    

    ManifoldEvolutionSettings strategySettings;
    strategySettings.UseInnerManifolds = false;
    strategySettings.OuterManifoldEpsilon = [](double distance) {
        return 0.1 * (1.0 - exp(-distance * distance / 0.0125));
    };
    strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
    {
        if (distance >= Geometry::DEFAULT_SCALAR_GRID_INIT_VAL)
            return 0.0;
        return 2.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
    };
    strategySettings.TimeStep = 0.01;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 40;
    globalSettings.DoRemeshing = false;
    globalSettings.ExportPerTimeStep = true;
    globalSettings.ExportTargetDistanceFieldAsImage = true;
    globalSettings.ProcedureName = "ShrinkWrappingAnIncompleteSurfacePointCloud_NoInnerSurfaceNoRemeshing";
    globalSettings.OutputPath = dataOutPath + "core_tests\\";
    globalSettings.ExportResult = false;

    // Set up the evolution strategy
    auto surfaceStrategy = std::make_shared<ManifoldSurfaceEvolutionStrategy>(
        strategySettings, MeshLaplacian::Barycentric,
        std::make_shared<std::vector<pmp::Point>>(targetPts));
    ManifoldEvolver evolver(globalSettings, std::move(surfaceStrategy));

    // Act
    evolver.Evolve();

    // Assert
    auto strategy = dynamic_cast<ManifoldSurfaceEvolutionStrategy*>(evolver.GetStrategy().get());
    auto resultOuterSurface = strategy->GetOuterSurfaceInOrigScale();
    auto resultInnerSurfaces = strategy->GetInnerSurfacesInOrigScale();
    ASSERT_TRUE(resultOuterSurface != nullptr);
    ASSERT_TRUE(resultInnerSurfaces.empty());
    size_t countInSphere = 0;
    for (const auto& v : resultOuterSurface->vertices())
    {
        const auto& pos = resultOuterSurface->position(v);
        if (norm(pos) > 0.75f)
            continue;
        countInSphere++;
    }
    EXPECT_LT(static_cast<double>(countInSphere) / resultOuterSurface->n_vertices(), 0.25);
}

TEST(ManifoldEvolverTests_ManifoldSurfaceSuite, ShrinkWrappingAnIncompleteSurfacePointCloud_NoRemeshing)
{
    // Arrange
    auto targetIcoSphere = ConstructIcoSphere(pmp::Point(0.0f, 0.0f, 0.0f), 0.75f, 2);
    auto targetPts = targetIcoSphere.positions();
    targetPts.erase(std::remove_if(targetPts.begin(), targetPts.end(),
        [](const pmp::Point& p) {
            // Filter points in the x > 0 && y > 0 && z > 0 octant
            return (p[0] > 0.0f && p[1] > 0.0f && p[2] > 0.0f);
        }),
        targetPts.end());

    ManifoldEvolutionSettings strategySettings;
    strategySettings.UseInnerManifolds = true;
    strategySettings.OuterManifoldEpsilon = [](double distance) {
        return 0.1 * (1.0 - exp(-distance * distance / 0.0125));
    };
    strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
    {
        if (distance >= Geometry::DEFAULT_SCALAR_GRID_INIT_VAL)
            return 0.0;
        return 2.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
    };
    strategySettings.InnerManifoldEpsilon = strategySettings.OuterManifoldEpsilon;
    strategySettings.InnerManifoldEta = strategySettings.OuterManifoldEta;
    strategySettings.TimeStep = 0.01;
    GlobalManifoldEvolutionSettings globalSettings;
    globalSettings.NSteps = 40;
    globalSettings.DoRemeshing = false;
    globalSettings.ExportPerTimeStep = true;
    globalSettings.ExportTargetDistanceFieldAsImage = true;
    globalSettings.ProcedureName = "ShrinkWrappingAnIncompleteSurfacePointCloud_NoRemeshing";
    globalSettings.OutputPath = dataOutPath + "core_tests\\";
    globalSettings.ExportResult = false;

    // Set up the evolution strategy
    auto surfaceStrategy = std::make_shared<ManifoldSurfaceEvolutionStrategy>(
        strategySettings, MeshLaplacian::Barycentric,
        std::make_shared<std::vector<pmp::Point>>(targetPts));
    ManifoldEvolver evolver(globalSettings, std::move(surfaceStrategy));

    // Act
    evolver.Evolve();

    // Assert
    auto strategy = dynamic_cast<ManifoldSurfaceEvolutionStrategy*>(evolver.GetStrategy().get());
    auto resultOuterSurface = strategy->GetOuterSurfaceInOrigScale();
    auto resultInnerSurfaces = strategy->GetInnerSurfacesInOrigScale();
    ASSERT_TRUE(resultOuterSurface != nullptr);
    ASSERT_FALSE(resultInnerSurfaces.empty());
    size_t countInSphere = 0;
    for (const auto& v : resultOuterSurface->vertices())
    {
        const auto& pos = resultOuterSurface->position(v);
        if (norm(pos) > 0.75f)
            continue;
        countInSphere++;
    }
    EXPECT_LT(static_cast<double>(countInSphere) / resultOuterSurface->n_vertices(), 0.0);
}