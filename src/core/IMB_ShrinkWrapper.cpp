
#include "utils/TimingUtils.h"

#include "IMB_ShrinkWrapper.h"


bool IMB_ShrinkWrapper::Preprocess()
{
	return false;
}

void IMB_ShrinkWrapper::StabilizeGeometries(pmp::Scalar stabilizationFactor)
{
}

bool IMB_ShrinkWrapper::Remesh()
{
    if (!m_ShouldRemesh)
        return true; // nothing to do.

    if (!m_Surface)
        return false;

    // TODO: play with this
    //return Utils::DeadlinedScope::Run([this]() {
    //    pmp::Remeshing remesher(*m_Surface);
    //    remesher.adaptive_remeshing(m_RemeshingSettings);
    //    },
    //    std::chrono::milliseconds(500) /* TODO: use customized time limit */);

    pmp::Remeshing remesher(*m_Surface);
    remesher.adaptive_remeshing(m_RemeshingSettings);
    return true;
}

bool IMB_ShrinkWrapper::PerformEvolutionStep(unsigned int stepId)
{
    //if (!IsRemeshingNecessary(*m_Surface, m_Settings.QualitySettings.FaceQualityFunc, m_Settings.QualitySettings.Range))
    //    return true;

    return false;
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

std::optional<pmp::SurfaceMesh> IMB_ShrinkWrapper::Perform()
{
    if (!Preprocess())
        return {};

    unsigned int step = 1;
    while (GetVertexCoveragePercentage() < m_Settings.ActivatedPointPercentageThreshold ||
        step > m_Settings.MaxSteps)
    {
        if (!PerformEvolutionStep(step))
            return GetSurface(); // integration error

        if (!Remesh())
            return GetSurface(); // terminate if it takes too long

        step++;
    }

    return GetSurface();
}
