#pragma once

#include "VertexSamplingStrategies.h"

#include <unordered_map>
#include <atomic>
#include <functional>
#include <mutex>
#include <optional>
#include <vector>
#include <queue>

#include "pmp/Types.h"

/// \brief A function to call when enough points are counted.
using MeshUpdateCallback = std::function<void(const std::vector<pmp::Point>&)>;

/// \brief A function to call when enough points are counted. Moves the result data to its inner scope.
/// TODO: Utilize this for dynamic update
using MeshUpdateMoveCallback = std::function<void(std::vector<pmp::Point>&&)>;

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
        IncrementalProgressTracker(const size_t& nTotal, const unsigned int& frequency, const std::function<void()>& callback)
            : m_nTotalExpectedVertices(nTotal), m_CompletionFrequency(frequency)
        {
            SetDispatcherCallback(callback);
        }

        void Update(const size_t& nLocalVerts, const bool& forceUpdate = false);

        [[nodiscard]] bool ShouldTriggerUpdate(const size_t& newCount) const;

        void SetDispatcherCallback(const std::function<void()>& callback)
        {
            m_DispatcherCallback = callback;
        }

    private:

        std::function<void()> m_DispatcherCallback;

        size_t m_nTotalExpectedVertices; //>! the total amount of expected mesh vertices.
        std::atomic<size_t> m_ProcessedVertices{ 0 }; //>! a counter for the total amount of processed vertices in all threads.
        unsigned int m_CompletionFrequency; //>! the frequency under which updates are triggered.
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

        void ShutDown();

        [[nodiscard]] size_t Size() const
        {
            return m_Tasks.size();
        }
    private:
        std::queue<std::function<void()>> m_Tasks;
        mutable std::mutex m_QueueMutex;
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
        IncrementalMeshBuilderDispatcher(const unsigned int& frequency, const VertexSelectionType& selectionType, const std::shared_ptr<IncrementalMeshFileHandler>& handler)
            : m_UpdateFrequency(frequency)
        {
            try
            {
				m_VertexSamplingStrategy = GetVertexSelectionStrategy(selectionType, frequency, handler);	            
	            m_ProgressTracker = std::make_unique<IncrementalProgressTracker>(m_VertexSamplingStrategy->GetVertexCountEstimate(), frequency, [this] { EnqueueMeshUpdate(); });
	            m_UpdateThread = std::thread([this] { ProcessQueue(); });
            }
            catch (const std::exception& e) {
                std::cerr << "IncrementalMeshBuilderDispatcher: Encountered an internal error: \n" << e.what() << '\n';
            }
        }

        [[nodiscard]] bool IsValid() const;

        ~IncrementalMeshBuilderDispatcher();

        void ProcessChunk(const char* start, const char* end, const std::optional<unsigned int>& seed);

        void SetMeshUpdateCallback(const MeshUpdateCallback& callback)
        {
            m_MeshUpdateCallback = callback;
        }

        void ProcessMeshUpdate(const std::vector<pmp::Point>& data) const;

        void EnqueueMeshUpdate();

	private:

        void ProcessQueue()
        {
            m_UpdateQueue.ProcessTasks();
        }

        MeshUpdateCallback m_MeshUpdateCallback;
        MeshUpdateQueue m_UpdateQueue;
        std::thread m_UpdateThread;
        std::unordered_map<std::thread::id, std::vector<pmp::Point>> m_ThreadResults;
        std::mutex m_ThreadResultMutex;

        const unsigned int& m_UpdateFrequency;
        std::atomic<unsigned int> m_UpdateCounter{ 0 };

        std::unique_ptr<IncrementalProgressTracker> m_ProgressTracker{ nullptr };
        std::unique_ptr<VertexSamplingStrategy> m_VertexSamplingStrategy{ nullptr };
	};
	
} // namespace IMB
