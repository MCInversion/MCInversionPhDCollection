
#include "pmp/algorithms/Normals.h"

#include "geometry/IcoSphereBuilder.h"

#include "sdf/SDF.h"

#include "IMB_ShrinkWrapper.h"


bool IMB_ShrinkWrapper::ComputeAmbientFields(pmp::Scalar& minTargetSize, pmp::Scalar& maxTargetSize, pmp::Point& targetBoundsCenter)
{
    if (m_Points.empty())
        return false;

    const pmp::BoundingBox ptCloudBBox(m_Points);
    targetBoundsCenter = ptCloudBBox.center();
    const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
    minTargetSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
    maxTargetSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
    const auto cellSize = minTargetSize / static_cast<pmp::Scalar>(m_Settings.FieldSettings.NVoxelsPerMinDimension);
    
    if (m_Settings.FieldSettings.FieldIsoLevel > 0.0)
    {
        // adjust iso-level according to point cloud distance field cell size.
        m_Settings.FieldSettings.FieldIsoLevel *= sqrt(3.0) / 2.0 * static_cast<double>(cellSize);
    }

    const SDF::PointCloudDistanceFieldSettings dfSettings{
        cellSize,
        m_Settings.FieldSettings.FieldExpansionFactor,
        DBL_MAX
    };
    m_DistanceField = std::make_shared<Geometry::ScalarGrid>(
        SDF::PointCloudDistanceFieldGenerator::Generate(m_Points, dfSettings));
    m_DFNegNormalizedGradient = std::make_shared<Geometry::VectorGrid>(ComputeNormalizedNegativeGradient(*m_DistanceField));

    return true;
}

void IMB_ShrinkWrapper::ConstructInitialSurface(const pmp::Scalar& minTargetSize, const pmp::Scalar& maxTargetSize, const pmp::Point& targetBoundsCenter)
{
    const pmp::Scalar outerSphereRadius = 0.5 * SPHERE_RADIUS_FACTOR *
        (minTargetSize + (0.5 + m_Settings.FieldSettings.FieldExpansionFactor) * maxTargetSize);

    m_InitialSphere = Geometry::Sphere3D{ targetBoundsCenter, outerSphereRadius };
    m_Surface = std::make_shared<pmp::SurfaceMesh>(ConstructIcoSphere(m_InitialSphere, m_Settings.LevelOfDetail));
}

void IMB_ShrinkWrapper::AssignRemeshingSettings()
{
    m_RemeshingSettingsCustom = CollectRemeshingSettingsFromIcoSphere_OLD(
        m_Settings.LevelOfDetail,
        m_InitialSphere.Radius * m_ScalingFactor,
        m_Settings.RemeshingSettings);
}

//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

bool IMB_ShrinkWrapper::Preprocess()
{
    pmp::Scalar minTargetSize{ FLT_MAX };
    pmp::Scalar maxTargetSize{ -FLT_MAX };
    pmp::Point targetCenter{};

    if (!ComputeAmbientFields(minTargetSize, maxTargetSize, targetCenter))
        return false;

    ConstructInitialSurface(minTargetSize, maxTargetSize, targetCenter);
    StabilizeGeometries(minTargetSize);
    AssignRemeshingSettings();
    return true;
}

/// \brief The power of the stabilizing scale factor.
constexpr pmp::Scalar SCALE_FACTOR_POWER_2D = 1.0 / 2.0;
/// \brief the reciprocal value of how many times the surface area element shrinks during evolution.
constexpr pmp::Scalar INV_SHRINK_FACTOR_2D = 5.0;

void IMB_ShrinkWrapper::StabilizeGeometries(const pmp::Scalar& minTargetSize, const pmp::Scalar& stabilizationFactor)
{
    const pmp::Scalar radius = stabilizationFactor * minTargetSize + (1.0 - stabilizationFactor) * m_InitialSphere.Radius;

    const unsigned int expectedVertexCount = (N_ICO_EDGES_0 * static_cast<unsigned int>(pow(4, m_Settings.LevelOfDetail) - 1) + 3 * N_ICO_VERTS_0) / 3;
    const pmp::Scalar expectedMeanCoVolArea = (4.0 * static_cast<pmp::Scalar>(M_PI) * radius * radius / static_cast<pmp::Scalar>(expectedVertexCount));
    const auto scalingFactor = pow(static_cast<pmp::Scalar>(m_Settings.TimeStep) / expectedMeanCoVolArea * INV_SHRINK_FACTOR_2D, SCALE_FACTOR_POWER_2D);
    m_ScalingFactor = scalingFactor;

    // -----------------------------------------------------------------------------------------------
    // All geometric quantities need to be scaled by the scalingFactor to ensure numerical stability.
    // For spatial conveniences, geometries must also be centered, i.e.: translated by -origin.
    // Geometries that are already centered at (0, 0, 0) are not translated.
    // -----------------------------------------------------------------------------------------------
    m_Settings.FieldSettings.FieldIsoLevel *= scalingFactor;
    m_Settings.PointActivationRadius *= scalingFactor;

    const pmp::mat4 transfMatrixGeomScale{
        scalingFactor, 0.0, 0.0, 0.0,
        0.0, scalingFactor, 0.0, 0.0,
        0.0, 0.0, scalingFactor, 0.0,
        0.0, 0.0, 0.0, 1.0
    };
    const pmp::Point& origin = m_InitialSphere.Center;
    const pmp::mat4 transfMatrixGeomMove{
        1.0, 0.0, 0.0, -origin[0],
        0.0, 1.0, 0.0, -origin[1],
        0.0, 0.0, 1.0, -origin[2],
        0.0, 0.0, 0.0, 1.0
    };
    const auto transfMatrixFull = transfMatrixGeomScale * transfMatrixGeomMove;
    m_TransformToOriginal = inverse(transfMatrixFull);

    // transform geometries
    (*m_Surface) *= transfMatrixGeomScale;

    // test box for geometry validation
    const pmp::Scalar evolBoxFactor = 5.0 * scalingFactor;
    m_EvolBox = pmp::BoundingBox(
        pmp::Point{ -radius, -radius, -radius } *evolBoxFactor,
        pmp::Point{ radius, radius, radius } *evolBoxFactor);

    (*m_DistanceField) *= transfMatrixFull;
    (*m_DistanceField) *= static_cast<double>(scalingFactor); // scale also the distance values.
    (*m_DFNegNormalizedGradient) *= transfMatrixFull;
    // values are supposed to be unit vectors regardless of scaling

    // Note: m_InitialSphere will remain in the original scale and position.
}

void IMB_ShrinkWrapper::Remesh()
{
    if (!m_ShouldRemesh)
        return; // nothing to do.

    pmp::Remeshing remesher(*m_Surface);
    remesher.adaptive_remeshing(m_RemeshingSettingsCustom);
    m_ShouldRemesh = false; // reset flag
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


bool IMB_ShrinkWrapper::PerformEvolutionStep(unsigned int stepId)
{
    const auto NVertices = static_cast<unsigned int>(m_Surface->n_vertices());
    SparseMatrix sysMat(NVertices, NVertices);
    Eigen::MatrixXd sysRhs(NVertices, 3);

    const auto tStep = m_Settings.TimeStep;

    pmp::Normals::compute_vertex_normals(*m_Surface);
    auto vNormalsProp = m_Surface->get_vertex_property<pmp::vec3>("v:normal");

    // prepare matrix & rhs for m_Surface:
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(static_cast<size_t>(NVertices) * 6);  // Assuming an average of 6 entries per vertex

    for (const auto v : m_Surface->vertices())
    {
        const auto& vPosToUpdate = m_Surface->position(v);

        if (m_Surface->is_boundary(v))
        {
            // freeze boundary/feature vertices
            const Eigen::Vector3d vertexRhs = vPosToUpdate;
            sysRhs.row(v.idx()) = vertexRhs;
            tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0));
            continue;
        }

        double vDistanceToTarget = m_DistanceField ? m_ScalarInterpolate(vPosToUpdate, *m_DistanceField) : DBL_MAX;
        vDistanceToTarget -= m_Settings.FieldSettings.FieldIsoLevel;
        const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec3(0, 0, 0);

        const double epsilonCtrlWeight = m_Settings.Epsilon(vDistanceToTarget);
        const auto& vNormal = static_cast<const pmp::vec3&>(vNormalsProp[v]); // vertex unit normal

        const auto negGradDotNormal = std::clamp(pmp::ddot(vNegGradDistanceToTarget, vNormal), -1.0, 1.0);
        const double advectionDistance = vDistanceToTarget;
        const double etaCtrlWeight = m_Settings.Eta(advectionDistance, negGradDotNormal);

        const Eigen::Vector3d vertexRhs = vPosToUpdate + tStep * etaCtrlWeight * vNormal;
        sysRhs.row(v.idx()) = vertexRhs;
        const auto tanRedistWeight = static_cast<double>(m_Settings.TangentialVelocityWeight) * std::abs(epsilonCtrlWeight);
        if (tanRedistWeight > 0.0)
        {
            // compute tangential velocity
            const auto vTanVelocity = ComputeTangentialUpdateVelocityAtVertex(*m_Surface, v, vNormal, tanRedistWeight);
            sysRhs.row(v.idx()) += tStep * Eigen::Vector3d(vTanVelocity);
        }

        const auto laplaceWeightInfo = m_ImplicitLaplacianFunction(*m_Surface, v); // Laplacian weights
        tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0 + tStep * epsilonCtrlWeight * static_cast<double>(laplaceWeightInfo.weightSum)));

        for (const auto& [w, weight] : laplaceWeightInfo.vertexWeights)
        {
            tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), w.idx(), -1.0 * tStep * epsilonCtrlWeight * static_cast<double>(weight)));
        }
    }

    // After the loop
    sysMat.setFromTriplets(tripletList.begin(), tripletList.end());
    m_ShouldRemesh = IsRemeshingNecessary(*m_Surface, m_Settings.QualitySettings.FaceQualityFunc, m_Settings.QualitySettings.Range);

    // solve
    Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double>> solver(sysMat);
    Eigen::MatrixXd x = solver.solve(sysRhs);
    if (solver.info() != Eigen::Success)
    {
        const std::string msg = "\nIMB_ShrinkWrapper::PerformEvolutionStep: solver.info() != Eigen::Success for time step id: "
            + std::to_string(stepId) + ", Error code: " + InterpretSolverErrorCode(solver.info()) + "\n";
        std::cerr << msg;
        return false;
    }

    // update vertex positions & verify mesh within bounds
    for (unsigned int i = 0; i < NVertices; i++)
    {
        const auto newPos = x.row(i);
        if (!m_EvolBox.Contains(newPos))
        {
            const std::string msg = "\nIMB_ShrinkWrapper::PerformEvolutionStep: vertex " + std::to_string(i)
                + ": (" + std::to_string(newPos[0]) + ", " + std::to_string(newPos[1]) + ", " + std::to_string(newPos[2]) + ") outside m_EvolBox: {"
                + "min_=(" + std::to_string(m_EvolBox.min()[0]) + ", " + std::to_string(m_EvolBox.min()[1]) + ", " + std::to_string(m_EvolBox.min()[2])
                + "), max_=(" + std::to_string(m_EvolBox.max()[0]) + ", " + std::to_string(m_EvolBox.max()[1]) + ", " + std::to_string(m_EvolBox.max()[2])
                + ")} for time step id: " + std::to_string(stepId) + "!\n";
            std::cerr << msg;
            return false;
        }
        m_Surface->position(pmp::Vertex(i)) = newPos;
    }

    return true;
}

//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

double IMB_ShrinkWrapper::GetVertexCoveragePercentage() const
{
    if (!m_Surface || !m_DistanceField)
        return -1.0;

    unsigned int count = 0;
    for (const auto& vPos : m_Surface->positions())
    {
        if (m_ScalarInterpolate(vPos, *m_DistanceField) > m_Settings.ActivatedPointPercentageThreshold)
            continue;

        count++;
    }

    return static_cast<double>(count) / static_cast<double>(m_Surface->n_vertices());
}

//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

std::optional<pmp::SurfaceMesh> IMB_ShrinkWrapper::GetSurface() const
{
    if (!m_Surface)
        return {};

    return *m_Surface;
}

//
// ========================================================
//

IMB_ShrinkWrapper::IMB_ShrinkWrapper(const IMB_ShrinkWrapperSettings& settings, const std::vector<pmp::Point>& points)
    : m_Settings(settings), m_Points(points)
{
    m_LaplacianAreaFunction = (m_Settings.LaplacianType == MeshLaplacian::Barycentric ?
        pmp::barycentric_area : pmp::voronoi_area);
    m_ImplicitLaplacianFunction = (m_Settings.LaplacianType == MeshLaplacian::Barycentric ?
        pmp::laplace_implicit_barycentric : pmp::laplace_implicit_voronoi);

    if (m_Settings.UseLinearGridInterpolation)
    {
        m_ScalarInterpolate = Geometry::TrilinearInterpolateScalarValue;
        m_VectorInterpolate = Geometry::TrilinearInterpolateVectorValue;
    }
    else
    {
        m_ScalarInterpolate = Geometry::GetNearestNeighborScalarValue;
        m_VectorInterpolate = Geometry::GetNearestNeighborVectorValue;
    }
}

// ---------------------------------------------------------------------------------

std::optional<pmp::SurfaceMesh> IMB_ShrinkWrapper::Perform()
{
    if (!Preprocess())
        return {};

    unsigned int step = 1;
    while (GetVertexCoveragePercentage() < m_Settings.ActivatedPointPercentageThreshold ||
        step < m_Settings.MaxSteps)
    {
        if (!PerformEvolutionStep(step))
            return GetSurface(); // integration error

        Remesh();
        step++;
    }

    return GetSurface();
}