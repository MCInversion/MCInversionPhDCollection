#pragma once

#include "VertexSamplingStrategies.h"

#include "Utils/IncrementalUtils.h"

#include <unordered_map>
#include <atomic>
#include <functional>
#include <mutex>
#include <optional>
#include <vector>
#include <queue>
#include <thread>

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
        IncrementalProgressTracker(
            const size_t& nTotal, const unsigned int& frequency, 
            const size_t& maxVertexCount, const size_t& minVertexCount,
            const std::function<void()>& addJobCallback, const std::function<void()>& terminationCallback)
            : m_MinVertexCount(minVertexCount),
    	      m_nTotalExpectedVertices(std::min(nTotal, maxVertexCount)),
    	      m_CompletionFrequency(frequency)
        {
            m_GrowthRate = (log(m_nTotalExpectedVertices) - log(minVertexCount)) / static_cast<double>(m_nTotalExpectedVertices);
#if DEBUG_PRINT
            DBG_OUT << "IncrementalProgressTracker::IncrementalProgressTracker: m_GrowthRate  = " << m_GrowthRate << "\n";
            DBG_OUT << "IncrementalProgressTracker::IncrementalProgressTracker: m_nTotalExpectedVertices  = " << m_nTotalExpectedVertices << "\n";
#endif
            SetDispatcherAddJobCallback(addJobCallback);
            SetDispatcherTerminationCallback(terminationCallback);
        }

        void Update(const size_t& nLocalVerts, const bool& forceUpdate = false);

        void SetDispatcherAddJobCallback(const std::function<void()>& callback)
        {
            m_DispatcherAddJobCallback = callback;
        }


        void SetDispatcherTerminationCallback(const std::function<void()>& callback)
        {
            m_DispatcherTerminationCallback = callback;
        }

    private:
        [[nodiscard]] bool ShouldTriggerUpdate(const size_t& newCount) const;

        std::function<void()> m_DispatcherAddJobCallback;
        std::function<void()> m_DispatcherTerminationCallback;

        std::atomic<size_t> m_ProcessedVertices{ 0 }; //>! a counter for the total amount of processed vertices in all threads.
        std::mutex m_Mutex; //>! mutex for ensuring thread safety of this object.

        std::atomic<unsigned int> m_UpdateCount{ 0 };
        double m_GrowthRate{ 1.0 }; //>! a helper parameter for non-linear update rate
        const size_t m_MinVertexCount; //>! the minimum amount of vertices to be rendered
        const size_t m_nTotalExpectedVertices; //>! the total amount of expected mesh vertices.
        const unsigned int m_CompletionFrequency; //>! the frequency under which updates are triggered.
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

        //void ForcedTerminate();
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
        IncrementalMeshBuilderDispatcher(const unsigned int& frequency, 
            const size_t& maxVertexCount, const VertexSelectionType& selectionType, 
            const std::shared_ptr<IncrementalMeshFileHandler>& handler)
            : m_UpdateFrequency(frequency)
        {
            try
            {
				m_VertexSamplingStrategy = GetVertexSelectionStrategy(selectionType, frequency, maxVertexCount, handler);

                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                // TODO: remove
                if (selectionType == IMB::VertexSelectionType::SoftMaxUniform)
                {
                    IMB::GeometricSamplingParams params;
                    params.ErrorMetricGradientMultiplier = 1.0;
                    params.DistributionMetricGradientMultiplier = 10'000.0; //0.000'1;
                    params.TargetVertexDensity = -1.0;
                    params.AvgNeighborhoodDisplacementMultiplier = 1.0;
                    params.SelectionSharpness = 0.001;
                    params.ConfidenceGM1Imax = 0.05;
                    params.ConfidenceNormal = 0.8;

                    dynamic_cast<IMB::SoftmaxUniformVertexSamplingStrategy*>(m_VertexSamplingStrategy.get())->Params = params;
                }
                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	            m_ProgressTracker = std::make_unique<IncrementalProgressTracker>(
                    m_VertexSamplingStrategy->GetVertexCountEstimate(), frequency, 
                    m_VertexSamplingStrategy->GetVertexCap(), m_VertexSamplingStrategy->GetMinVertexCount(),
                    [this] { EnqueueMeshUpdate(); },
                    [this] { ShutDownQueue(); });
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

        void ShutDownQueue();

        MeshUpdateCallback m_MeshUpdateCallback;
        MeshUpdateQueue m_UpdateQueue;
        std::thread m_UpdateThread;
        std::unordered_map<std::thread::id, std::vector<pmp::Point>> m_ThreadResults;
        std::mutex m_ThreadResultMutex;

        const unsigned int& m_UpdateFrequency;
        std::atomic<unsigned int> m_UpdateCounter{ 0 };

        std::unique_ptr<IncrementalProgressTracker> m_ProgressTracker{ nullptr };
        std::unique_ptr<VertexSamplingStrategy> m_VertexSamplingStrategy{ nullptr };

        bool m_UpdateThreadTerminated{ false };
	};
	
} // namespace IMB
