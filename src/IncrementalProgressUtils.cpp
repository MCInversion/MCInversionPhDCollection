#include "IncrementalProgressUtils.h"

namespace IMB
{
	void IncrementalMeshBuilderDispatcher::ProcessChunk(const char* start, const char* end, const std::optional<unsigned int>& seed)
	{
		std::vector<pmp::Point> localResults;
		m_VertexSamplingStrategy->Sample(start, end, localResults, seed, m_ProgressTracker);

		// Locking to safely update shared resources
		{
			std::lock_guard lock(m_UpdateMutex);
			m_ThreadResult.insert(m_ThreadResult.end(), localResults.begin(), localResults.end());
		}

		// Check if it's time to update the mesh
		if (m_ProgressTracker.ShouldTriggerUpdate()) 
		{
			ProcessMeshUpdate();
		}
	}

	void IncrementalMeshBuilderDispatcher::ProcessMeshUpdate()
	{
		std::vector<std::vector<unsigned int>> resultVertexIds;
		m_MeshingStrategy->Process(m_ThreadResult, resultVertexIds);

		m_ProgressCallback(m_ThreadResult, resultVertexIds);
	}

	void IncrementalProgressTracker::Update(const size_t& nLocalVerts)
	{
		// increment a shared value for the amount of processed vertices
		auto currentCount = m_ProcessedVertices.fetch_add(nLocalVerts, std::memory_order_relaxed);
		if (!ShouldTriggerUpdate())
			return;

		m_DispatcherCallback();
	}

	constexpr unsigned int DEFAULT_MAX_VERTEX_CAPACITY = 1'000'000;

	bool IncrementalProgressTracker::ShouldTriggerUpdate() const
	{
		const size_t count = m_ProcessedVertices.load(std::memory_order_relaxed);

		if (count >= DEFAULT_MAX_VERTEX_CAPACITY)
			return true; // exceeding max mesh memory capacity

		if (static_cast<double>(count) < (static_cast<double>(m_nTotalExpectedVertices) / 100.0) * m_CompletionFrequency)
			return false;

		return true;
	}
}