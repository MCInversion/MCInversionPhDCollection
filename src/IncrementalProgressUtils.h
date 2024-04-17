#pragma once

#include "VertexSamplingStrategies.h"
#include "PointCloudMeshingStrategies.h"

#include <atomic>
#include <functional>
#include <mutex>
#include <optional>
#include <vector>

namespace pmp
{
    class Point;
}

using MeshUpdateCallback = std::function<void(const std::vector<pmp::Point>& , const std::vector<std::vector<unsigned int>>& )>;

namespace IMB
{
    class IncrementalProgressTracker
	{
    public:
        IncrementalProgressTracker(const size_t& nTotal, const double& frequency, const std::function<void()>& callback)
            : m_nTotalExpectedVertices(nTotal), m_CompletionFrequency(frequency)
        {
            SetDispatcherCallback(callback);
        }

        void Update(const size_t& nLocalVerts);

        [[nodiscard]] bool ShouldTriggerUpdate() const;

        void SetDispatcherCallback(const std::function<void()>& callback)
        {
            m_DispatcherCallback = callback;
        }

    private:

        std::function<void()> m_DispatcherCallback;

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
            : m_ProgressTracker(totalExpectedVertices, frequency, [this]() { ProcessMeshUpdate(); }),
            m_VertexSamplingStrategy(std::move(vertStrategy)),
            m_MeshingStrategy(std::move(meshStrategy))
        {
        }

        void ProcessChunk(const char* start, const char* end, const std::optional<unsigned int>& seed);

        void SetProgressCallback(const MeshUpdateCallback& callback)
        {
            m_ProgressCallback = callback;
        }

        void ProcessMeshUpdate();

	private:

        MeshUpdateCallback m_ProgressCallback;
        std::mutex m_UpdateMutex;
        std::vector<pmp::Point> m_ThreadResult;

        IncrementalProgressTracker m_ProgressTracker;
        std::unique_ptr<VertexSamplingStrategy> m_VertexSamplingStrategy{ nullptr };
        std::unique_ptr<PointCloudMeshingStrategy> m_MeshingStrategy{ nullptr };
	};
	
} // namespace IMB
