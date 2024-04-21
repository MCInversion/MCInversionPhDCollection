#include "IncrementalProgressUtils.h"

#include "utils/IncrementalUtils.h"

#include <ranges>

namespace IMB
{
	bool IncrementalMeshBuilderDispatcher::IsValid() const
	{
		return (m_ProgressTracker && m_VertexSamplingStrategy);
	}

	IncrementalMeshBuilderDispatcher::~IncrementalMeshBuilderDispatcher()
	{
		m_UpdateQueue.ShutDown();
		if (m_UpdateThread.joinable() && !m_UpdateThreadTerminated)
		{
			m_UpdateThread.join();
		}
	}

	void IncrementalMeshBuilderDispatcher::ProcessChunk(const char* start, const char* end, const std::optional<unsigned int>& seed)
	{
#if DEBUG_PRINT
		DBG_OUT << "IncrementalMeshBuilderDispatcher::ProcessChunk: [" << FormatAddresses(start, end) << "] ...\n";
#endif
		auto& localResults = m_ThreadResults[std::this_thread::get_id()];
		m_VertexSamplingStrategy->Sample(start, end, localResults, seed, *m_ProgressTracker);
	}

	void IncrementalMeshBuilderDispatcher::ProcessMeshUpdate(const std::vector<pmp::Point>& data) const
	{
		// Process mesh update with isolated data
		m_MeshUpdateCallback(data);
	}

	void IncrementalMeshBuilderDispatcher::EnqueueMeshUpdate()
	{
		// Gather point data from all worker threads.
		std::vector<pmp::Point> aggregatedData;
		{
			std::lock_guard lock(m_ThreadResultMutex);
			for (auto& results : m_ThreadResults | std::views::values) 
			{
				aggregatedData.insert(aggregatedData.end(), results.begin(), results.end());
				results.clear();  // Clear the thread-specific results after aggregating
			}
		}
		// ............................................

		if (aggregatedData.empty()) 
			return;  // Nothing to update

		//if (m_UpdateCounter.load() > m_UpdateFrequency)
		//{
#if DEBUG_PRINT
		//	DBG_OUT << "IncrementalMeshBuilderDispatcher::EnqueueMeshUpdate: Reached " << m_UpdateCounter.load() << " updates with " << aggregatedData.size() << " remaining points! Terminating.\n";
#endif
		//	m_UpdateQueue.Enqueue([this, dataCopy = std::move(aggregatedData)]() mutable {
		//		this->ProcessMeshUpdate(dataCopy);
		//	});
		//	return;
		//}
#if DEBUG_PRINT
		DBG_OUT << "IncrementalMeshBuilderDispatcher::EnqueueMeshUpdate: with "<< aggregatedData.size() << " points ...\n";
#endif

		m_UpdateQueue.Enqueue([this, dataCopy = std::move(aggregatedData)]() mutable {
			this->ProcessMeshUpdate(dataCopy);
		});
		m_UpdateCounter.fetch_add(1);

#if DEBUG_PRINT
		DBG_OUT << "IncrementalMeshBuilderDispatcher::EnqueueMeshUpdate: Job Queue contains " << m_UpdateQueue.Size() << " jobs.\n";
#endif
	}

	void IncrementalMeshBuilderDispatcher::ShutDownQueue()
	{
		{ // ensure that the lock is held only for the duration of the operations that need synchronization
			std::lock_guard lock(m_ThreadResultMutex);
			m_UpdateQueue.ShutDown();
			if (m_UpdateThread.joinable() && !m_UpdateThreadTerminated)
			{
				m_UpdateThread.join();
				m_UpdateThreadTerminated = true;
			}
		}
	}

	void IncrementalProgressTracker::Update(const size_t& nLocalVerts, const bool& forceUpdate)
	{
		// increment a shared value for the amount of processed vertices
		const auto newCount = m_ProcessedVertices.fetch_add(nLocalVerts, std::memory_order_relaxed);

		if (forceUpdate && newCount < m_nTotalExpectedVertices)
		{
#if DEBUG_PRINT
			DBG_OUT << "IncrementalProgressTracker::Update: Forced update with newCount  = " << newCount << "\n";
			DBG_OUT << "IncrementalProgressTracker::Update: invoking m_DispatcherCallback ... \n";
#endif
			m_DispatcherAddJobCallback();
			return;
		}

#if DEBUG_PRINT
		DBG_OUT << "IncrementalProgressTracker::Update: newCount  = " << newCount << "\n";
#endif
		if (!ShouldTriggerUpdate(newCount))
			return;
#if DEBUG_PRINT
		DBG_OUT << "IncrementalProgressTracker::Update: invoking m_DispatcherCallback.\n";
#endif
		m_DispatcherAddJobCallback();
	}

	bool IncrementalProgressTracker::ShouldTriggerUpdate(const size_t& newCount) const
	{
		if (newCount >= m_nTotalExpectedVertices)
		{
#if DEBUG_PRINT
			DBG_OUT << "IncrementalProgressTracker::ShouldTriggerUpdate: MAX CAPACITY REACHED!!!\n";
#endif
			m_DispatcherTerminationCallback();
			return false;
		}

		// Calculate the threshold for updates based on the total expected vertices and completion frequency
		const double progressTrackerUpdateThreshold = static_cast<double>(m_nTotalExpectedVertices) / m_CompletionFrequency;
		return static_cast<double>(newCount) >= progressTrackerUpdateThreshold;
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
#if DEBUG_PRINT
		DBG_OUT << "MeshUpdateQueue::ProcessTasks: ... \n";
#endif
		std::unique_lock lock(m_QueueMutex);
		while (!m_ShutDown || !m_Tasks.empty()) {
			m_Condition.wait(lock, [this] { return m_ShutDown || !m_Tasks.empty(); });
			if (!m_Tasks.empty()) {
				auto task = std::move(m_Tasks.front());
				m_Tasks.pop();
				lock.unlock();
				task();
				lock.lock();
			}
		}
#if DEBUG_PRINT
		DBG_OUT << "MeshUpdateQueue::ProcessTasks: ... done.\n";
		DBG_OUT << "MeshUpdateQueue::ProcessTasks: " << m_Tasks.size() << " tasks remaining!\n";
#endif
	}

	void MeshUpdateQueue::ShutDown()
	{
		{ // ensure that the lock is held only for the duration of the operations that need synchronization
			std::lock_guard lock(m_QueueMutex);
			m_ShutDown = true;
		}
		m_Condition.notify_all();
	}

	//void MeshUpdateQueue::ForcedTerminate()
	//{
	//	{// ensure that the lock is held only for the duration of the operations that need synchronization
	//		std::lock_guard lock(m_QueueMutex);
	//		while (!m_Tasks.empty())
	//			m_Tasks.pop();
	//	}
	//	m_Condition.notify_all();
	//}
}