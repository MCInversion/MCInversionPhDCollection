#pragma once

#include "VertexSamplingStrategies.h"
#include "PointCloudMeshingStrategies.h"

#include <atomic>
#include <functional>
#include <mutex>
#include <optional>
#include <vector>
#include <queue>

namespace pmp
{
    class Point;
}

using MeshUpdateCallback = std::function<void(const std::vector<pmp::Point>&)>;

namespace IMB
{
    /// =======================================================================================
    /// \brief Monitors the progress of vertex processing across all threads and triggers mesh updates when necessary.
    ///
    /// \class IncrementalProgressTracker
    /// 
    /// =======================================================================================
    class IncrementalProgressTracker
	{
    public:
        IncrementalProgressTracker(const size_t& nTotal, const double& frequency, const std::function<void()>& callback)
            : m_nTotalExpectedVertices(nTotal), m_CompletionFrequency(frequency)
        {
            SetDispatcherCallback(callback);
        }

        void Update(const size_t& nLocalVerts);

        [[nodiscard]] bool ShouldTriggerUpdate(const size_t& currentCount) const;

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

    /// =======================================================================================
    /// \brief A simple queue for managing the dispatching of mesh update tasks to different threads.
    ///         
    /// \class MeshUpdateQueue
    /// 
    /// =======================================================================================
    class MeshUpdateQueue 
    {
    public:
        void Enqueue(const std::function<void()>& task);

        void ProcessTasks();

        void Shutdown();
    private:
        std::queue<std::function<void()>> m_Tasks;
        std::mutex m_QueueMutex;
        std::condition_variable m_Condition;
        bool m_ShutDown{ false };
    };

    /// =======================================================================================
    /// \brief Manages the dispatching of vertex processing tasks to different threads and coordinates the collection of results for mesh updates.
    /// 
    /// \class IncrementalMeshBuilderDispatcher
    /// 
    /// =======================================================================================
	class IncrementalMeshBuilderDispatcher
	{
	public:
        IncrementalMeshBuilderDispatcher(const size_t& totalExpectedVertices, const double& frequency,
            std::unique_ptr<VertexSamplingStrategy> vertStrategy)
            : m_ProgressTracker(totalExpectedVertices, frequency, [this]() { this->EnqueueMeshUpdate(); }),
            m_VertexSamplingStrategy(std::move(vertStrategy))
        {
            m_UpdateThread = std::thread([this] { m_UpdateQueue.ProcessTasks(); });
        }

        ~IncrementalMeshBuilderDispatcher();

        void ProcessChunk(const char* start, const char* end, const std::optional<unsigned int>& seed);

        void SetMeshUpdateCallback(const MeshUpdateCallback& callback)
        {
            m_MeshUpdateCallback = callback;
        }

        void ProcessMeshUpdate(const std::vector<pmp::Point>& data);

        void EnqueueMeshUpdate();

	private:

        void ProcessQueue()
        {
            m_UpdateQueue.ProcessTasks();
        }

        MeshUpdateCallback m_MeshUpdateCallback;
        MeshUpdateQueue m_UpdateQueue;
        std::thread m_UpdateThread;
        std::vector<pmp::Point> m_ThreadResult;

        IncrementalProgressTracker m_ProgressTracker;
        std::unique_ptr<VertexSamplingStrategy> m_VertexSamplingStrategy{ nullptr };
	};
	
} // namespace IMB
