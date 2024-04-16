#include "IncrementalProgressUtils.h"

namespace IMB
{
	void IncrementalMeshBuilderDispatcher::ProcessChunk(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed)
	{
		//m_VertexSamplingStrategy->Sample(start, end, result, seed, m_ProgressTracker);
		if (m_ProgressTracker.ShouldTriggerUpdate())
			return;

		// TODO: Trigger some action
	}

	void IncrementalProgressTracker::Update(const size_t& processedPerChunk)
	{
		m_ProcessedVertices += processedPerChunk;
	}

	bool IncrementalProgressTracker::ShouldTriggerUpdate()
	{
		size_t count = m_ProcessedVertices.load(std::memory_order_relaxed);

		if (count < (m_nTotalExpectedVertices / 100.0) * m_CompletionFrequency)
			return false;

		//std::lock_guard lock(m_Mutex);
		// Double-check locking pattern
		//if (m_ProcessedVertices >= (m_nTotalExpectedVertices / 100.0) * m_CompletionFrequency)
		//{
		//	m_ProcessedVertices -= (m_nTotalExpectedVertices / 100.0) * m_CompletionFrequency;
		//	return true;
		//}

		return true;
	}
}