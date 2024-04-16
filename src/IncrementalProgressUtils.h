#pragma once

#include "VertexSamplingStrategies.h"
#include "PointCloudMeshingStrategies.h"

#include <atomic>
#include <mutex>
#include <optional>
#include <vector>

namespace pmp
{
    class Point;
}


namespace IMB
{
    class IncrementalProgressTracker
	{
    public:
        IncrementalProgressTracker(const size_t& nTotal, const double& frequency)
            : m_nTotalExpectedVertices(nTotal), m_CompletionFrequency(frequency) {}

        void Update(const size_t& processedPerChunk);

        [[nodiscard]] bool ShouldTriggerUpdate();

    private:

        size_t m_nTotalExpectedVertices; //>! the total amount of expected mesh vertices.
        std::atomic<size_t> m_ProcessedVertices{ 0 }; //>! a counter for the total amount of processed vertices in all threads.
        double m_CompletionFrequency; //>! the frequency under which updates are triggered.
        std::mutex m_Mutex; //>! mutex for ensuring thread safety of this object.
    };

	class IncrementalMeshBuilderDispatcher
	{
	public:
        IncrementalMeshBuilderDispatcher(const size_t& totalExpectedVertices, const double& frequency,
            std::unique_ptr<VertexSamplingStrategy> vertStrategy, std::unique_ptr<PointCloudMeshingStrategy> meshStrategy)
            : m_ProgressTracker(totalExpectedVertices, frequency),
            m_VertexSamplingStrategy(std::move(vertStrategy)),
            m_MeshingStrategy(std::move(meshStrategy))
        {
        }

        void ProcessChunk(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed);

	private:

        IncrementalProgressTracker m_ProgressTracker;
        std::unique_ptr<VertexSamplingStrategy> m_VertexSamplingStrategy{ nullptr };
        std::unique_ptr<PointCloudMeshingStrategy> m_MeshingStrategy{ nullptr };
	};
	
} // namespace IMB
