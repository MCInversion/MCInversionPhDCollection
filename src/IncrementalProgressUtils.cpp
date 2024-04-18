#include "IncrementalProgressUtils.h"

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
		std::vector<pmp::Point> localResults;
		m_VertexSamplingStrategy->Sample(start, end, localResults, seed, m_ProgressTracker);
	}

	void IncrementalMeshBuilderDispatcher::ProcessMeshUpdate(const std::vector<pmp::Point>& data)
	{
		// Process mesh update with isolated data
		m_MeshUpdateCallback(data);
		m_ThreadResult.clear();
	}

	void IncrementalMeshBuilderDispatcher::EnqueueMeshUpdate()
	{
		auto dataCopy = m_ThreadResult;  // Copy data for the task
		m_UpdateQueue.Enqueue([this, &dataCopy]() { this->ProcessMeshUpdate(dataCopy); });
	}

	void IncrementalProgressTracker::Update(const size_t& nLocalVerts)
	{
		// increment a shared value for the amount of processed vertices
		auto currentCount = m_ProcessedVertices.fetch_add(nLocalVerts, std::memory_order_relaxed);
		if (!ShouldTriggerUpdate(currentCount))
			return;
		m_DispatcherCallback();
	}

	constexpr unsigned int DEFAULT_MAX_VERTEX_CAPACITY = 1'000'000;

	bool IncrementalProgressTracker::ShouldTriggerUpdate(const size_t& currentCount) const
	{
		return (currentCount >= m_nTotalExpectedVertices * m_CompletionFrequency) ||
			(currentCount >= DEFAULT_MAX_VERTEX_CAPACITY);
	}

	void MeshUpdateQueue::Enqueue(const std::function<void()>& task)
	{
		{ // ensure that the lock is held only for the duration of the operations that need synchronization
			std::lock_guard<std::mutex> lock(m_QueueMutex);
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
				std::unique_lock<std::mutex> lock(m_QueueMutex);
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
			std::lock_guard<std::mutex> lock(m_QueueMutex);
			m_ShutDown = true;
		}
		m_Condition.notify_all();
	}
}