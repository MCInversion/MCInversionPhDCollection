#include "VertexSamplingStrategies.h"

#include <numeric>
#include <random>
//#include <unordered_set>

#include "IncrementalProgressUtils.h"
#include "IncrementalMeshFileHandler.h"

#include "utils/IncrementalUtils.h"


namespace
{
	/// \brief A verification utility for the uniqueness of randomly generated indices.
	//[[nodiscard]] size_t CountNonUniqueIndices(const std::vector<size_t>& indices)
	//{
	//	const std::unordered_set uniqueCheck(indices.begin(), indices.end());
	//	return indices.size() - uniqueCheck.size();
	//}
	// USE CASE:
	//#if DEBUG_PRINT
	//	if (const auto nonUnique = CountNonUniqueIndices(resultIndices); nonUnique == 0)
	//		DBG_OUT << "RandomSampleIndices: Finished generating " << expectedCount << " unique indices.\n";
	//	else
	//		DBG_OUT << "RandomSampleIndices: [WARNING]: Generated " << nonUnique << "/" << expectedCount << " non-unique indices!\n";
	//#endif

	/// \brief Generates expectedCount indices from 0 to expectedCount - 1 into resultIndices sequentially.
	void SampleIndicesSequentially(const size_t& expectedCount, std::vector<size_t>& resultIndices)
	{
		resultIndices.clear();
		resultIndices.resize(expectedCount);
		std::iota(resultIndices.begin(), resultIndices.end(), 0);
	}

	/// \brief Generates expectedCount random indices from 0 to expectedCount - 1 into resultIndices from the unifom distribution type.
	void RandomSampleIndices(const size_t& expectedCount, std::vector<size_t>& resultIndices, const std::optional<unsigned int>& seed)
	{
		if (expectedCount == 0)
		{
			std::cerr << "RandomSampleIndices: expectedCount is 0!\n";
			return;
		}

#if DEBUG_PRINT
		DBG_OUT << "RandomSampleIndices: Starting...\n";
#endif

		resultIndices.clear();
		resultIndices.reserve(expectedCount);
		std::mt19937 gen(seed ? *seed : std::random_device{}());
		std::uniform_int_distribution<size_t> distrib(0, expectedCount - 1);

		std::vector usedIndices(expectedCount, false); 

		while (resultIndices.size() < expectedCount)
		{
			size_t newIndex = distrib(gen);
			if (!usedIndices[newIndex]) // Check if the index has not been used yet
			{
				usedIndices[newIndex] = true; // Mark the index as used
				resultIndices.push_back(newIndex);
			}
		}

#if DEBUG_PRINT
		DBG_OUT << "RandomSampleIndices: Finished generating " << expectedCount << " unique indices.\n";
#endif
	}	
}

namespace IMB
{
	void SequentialVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		std::vector<size_t> indices;
		SampleIndicesSequentially(m_FileHandler->GetLocalVertexCountEstimate(start, end), indices);
		m_FileHandler->Sample(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	void UniformRandomVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		std::vector<size_t> indices;
		RandomSampleIndices(m_FileHandler->GetLocalVertexCountEstimate(start, end), indices, seed);
		m_FileHandler->Sample(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	void NormalRandomVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		std::vector<size_t> indices;
		// TODO: verify utility before implementing
		RandomSampleIndices(m_FileHandler->GetLocalVertexCountEstimate(start, end), indices, seed);
		m_FileHandler->Sample(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	void SoftmaxFeatureDetectingVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		throw std::runtime_error("SampleVerticesWithSoftmaxFeatureDectection Not implemented\n");
	}

	constexpr unsigned int FREQUENCY_UPDATE_MULTIPLIER = 3;

	VertexSamplingStrategy::VertexSamplingStrategy(const unsigned int& completionFrequency, const size_t& maxVertexCount, const std::shared_ptr<IncrementalMeshFileHandler>& handler)
	{
		m_FileHandler = handler;
		if (!m_FileHandler)
		{
			throw std::invalid_argument("VertexSamplingStrategy::VertexSamplingStrategy: m_FileHandler could not be created!\n");
		}
		const auto totalExpectedVertices = std::min(m_FileHandler->GetGlobalVertexCountEstimate(), maxVertexCount);
		m_UpdateThreshold = static_cast<size_t>(std::round(static_cast<double>(totalExpectedVertices) /
			static_cast<double>(completionFrequency * FREQUENCY_UPDATE_MULTIPLIER)));
#if DEBUG_PRINT
		DBG_OUT << "VertexSamplingStrategy::VertexSamplingStrategy: completionFrequency = " << completionFrequency << " jobs per file load.\n";
		DBG_OUT << "VertexSamplingStrategy::VertexSamplingStrategy: m_UpdateThreshold = " << m_UpdateThreshold << " vertices.\n";
#endif
	}

	size_t VertexSamplingStrategy::GetVertexCountEstimate() const
	{
		if (!m_FileHandler)
		{
			return 0;
		}
		return m_FileHandler->GetGlobalVertexCountEstimate();
	}
} // namespace IMB