#include "ManifoldEvolver.h"

void ManifoldEvolver::PerformEvolutionStep(unsigned int step)
{
    // Implement the evolution logic for a single step
    // Example:
    // m_GeometryAdapter.AdvanceVertices(m_Settings.TimeStep);
}

void ManifoldEvolver::AdaptiveRemeshing()
{
    // Implement adaptive remeshing logic
    // Example:
    // m_GeometryAdapter.AdaptiveRemesh();
}

void ManifoldEvolver::ExportCurrentState(unsigned int step)
{
    /**
     * \brief Writes m_GeometryAdapter using the configured output path.
     * \param step The current time step index.
     */
    std::string outputPath = m_Settings.OutputPath + "/state_" + std::to_string(step) + m_OutputExtension;
    //m_GeometryAdapter.Export(outputPath);
}

void ManifoldEvolver::ExportFinalResult()
{
    /**
     * \brief Writes the final result after all time steps are completed.
     */
    //m_GeometryAdapter.Export(m_Settings.OutputPath + "/final_result.obj");
}

// TODO: change for something that works
[[nodiscard]] bool ShouldPerformRemeshing(const Geometry::ManifoldGeometryAdapter& outerGeom, const std::vector<std::shared_ptr<Geometry::ManifoldGeometryAdapter>>& innerGeoms)
{
    return false;
}

void ManifoldEvolver::Evolve()
{


    for (unsigned int step = 0; step < m_Settings.NSteps; ++step) 
    {
        PerformEvolutionStep(step);

        if (m_Settings.ExportPerTimeStep)
        {
            ExportCurrentState(step);
        }

        if (m_Settings.DoRemeshing && ShouldPerformRemeshing(*m_OuterManifoldAdapter, m_InnerManifoldAdapters))
        {
            AdaptiveRemeshing();
        }
    }

    if (m_Settings.ExportResult) 
    {
        ExportFinalResult();
    }
}
