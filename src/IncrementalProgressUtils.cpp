#include "IncrementalProgressUtils.h"

#include "utils/IncrementalUtils.h"

namespace IMB
{
	IncrementalMeshBuilderDispatcher::~IncrementalMeshBuilderDispatcher()
	{
		m_UpdateQueue.Shutdown();
		if (m_UpdateThread.joinable())
		{
			m_UpdateThread.join();
		}
	}

	void IncrementalMeshBuilderDispatcher::ProcessChunk(const char* start, const char* end, const std::optional<unsigned int>& seed)
	{
#if DEBUG_PRINT
		DBG_OUT << "IncrementalMeshBuilderDispatcher::ProcessChunk: [" << FormatAddresses(start, end) << "] ...\n";
#endif
		//std::vector<pmp::Point> localResults;
		m_VertexSamplingStrategy->Sample(start, end, m_ThreadResult, seed, m_ProgressTracker);
	}

	void IncrementalMeshBuilderDispatcher::ProcessMeshUpdate(std::vector<pmp::Point>&& data) const
	{
		// Process mesh update with isolated data
		m_MeshUpdateCallback(std::move(data));
	}

	void IncrementalMeshBuilderDispatcher::EnqueueMeshUpdate()
	{
#if DEBUG_PRINT
		DBG_OUT << "IncrementalMeshBuilderDispatcher::EnqueueMeshUpdate: with "<< m_ThreadResult.size() << " points ...\n";
#endif
		m_UpdateQueue.Enqueue([this](){ this->ProcessMeshUpdate(std::move(m_ThreadResult)); });
#if DEBUG_PRINT
		DBG_OUT << "IncrementalMeshBuilderDispatcher::EnqueueMeshUpdate: Job Queue contains " << m_UpdateQueue.Size() << " jobs.\n";
#endif
	}

	void IncrementalProgressTracker::Update(const size_t& nLocalVerts, const bool& forceUpdate)
	{
		if (forceUpdate)
		{
#if DEBUG_PRINT
			DBG_OUT << "IncrementalProgressTracker::Update: Forced update.\n";
#endif
			m_DispatcherCallback();
#if DEBUG_PRINT
			DBG_OUT << "IncrementalProgressTracker::Update: invoking m_DispatcherCallback.\n";
#endif
			return;
		}

		// increment a shared value for the amount of processed vertices
		const auto newCount = m_ProcessedVertices.fetch_add(nLocalVerts, std::memory_order_relaxed);
#if DEBUG_PRINT
		DBG_OUT << "IncrementalProgressTracker::Update: newCount  = " << newCount << "\n";
#endif
		if (!ShouldTriggerUpdate(newCount))
			return;
#if DEBUG_PRINT
		DBG_OUT << "IncrementalProgressTracker::Update: invoking m_DispatcherCallback.\n";
#endif
		m_DispatcherCallback();
	}

	constexpr unsigned int DEFAULT_MAX_VERTEX_CAPACITY = 1'000'000;

	bool IncrementalProgressTracker::ShouldTriggerUpdate(const size_t& newCount) const
	{
		// Calculate the threshold for updates based on the total expected vertices and completion frequency
		const double progressTrackerUpdateThreshold = static_cast<double>(m_nTotalExpectedVertices) / m_CompletionFrequency;
		return static_cast<double>(newCount) >= progressTrackerUpdateThreshold || m_ProcessedVertices + newCount >= DEFAULT_MAX_VERTEX_CAPACITY;
	}

	void MeshUpdateQueue::Enqueue(const std::function<void()>& task)
	{
		{ // ensure that the lock is held only for the duration of the operations that need synchronization
			std::lock_guard lock(m_QueueMutex);
			m_Tasks.push(task);
		}
		m_Condition.notify_one();
	}

	void MeshUpdateQueue::ProcessTasks()
	{
		while (!m_ShutDown) 
		{
			std::function<void()> task;
			{ // ensure that the lock is held only for the duration of the operations that need synchronization
				std::unique_lock lock(m_QueueMutex);
				m_Condition.wait(lock, [this] { return !m_Tasks.empty() || m_ShutDown; });
				if (m_ShutDown) break;
				task = std::move(m_Tasks.front());
				m_Tasks.pop();
			}
			task();
		}
	}

	void MeshUpdateQueue::Shutdown()
	{
		{ // ensure that the lock is held only for the duration of the operations that need synchronization
			std::lock_guard lock(m_QueueMutex);
			m_ShutDown = true;
		}
		m_Condition.notify_all();
	}
}