#include "VertexSamplingStrategies.h"

#include "IncrementalProgressUtils.h"

#include "utils/IncrementalUtils.h"

#include <random>

namespace
{
	void SampleIndicesUniformly(const char* start, const char* end, std::vector<size_t>& resultIndices, const std::optional<unsigned int>& seed)
	{
#if DEBUG_PRINT
		DBG_OUT << "SampleIndicesUniformly: ... \n";
#endif
		resultIndices.clear();
		std::mt19937 gen;
		if (seed.has_value())
		{
			gen.seed(seed.value());
		}
		else
		{
			std::random_device rd;
			gen = std::mt19937(rd());
		}
		const size_t nLines = std::count(start, end, '\n');
		std::uniform_int_distribution<> distrib(0, nLines - 1);
		resultIndices.reserve(nLines);
		for (size_t i = 0; i < nLines; i++)
		{
			resultIndices.push_back(distrib(gen));
		}
#if DEBUG_PRINT
		DBG_OUT << "SampleIndicesUniformly: from " << nLines << " lines ... done.\n";
#endif
	}

	void SampleIndicesNormalRandomly(const char* start, const char* end, std::vector<size_t>& resultIndices, const std::optional<unsigned int>& seed)
	{
#if DEBUG_PRINT
		DBG_OUT << "SampleIndicesNormalRandomly: ... \n";
#endif
		resultIndices.clear();
		std::mt19937 gen;
		if (seed.has_value())
		{
			gen.seed(seed.value());
		}
		else
		{
			std::random_device rd;
			gen = std::mt19937(rd());
		}
		const size_t nLines = std::count(start, end, '\n');
		std::normal_distribution<> distrib(nLines / 2, nLines / 6);
		resultIndices.reserve(nLines);
		for (size_t i = 0; i < nLines; i++)
		{
			resultIndices.push_back(distrib(gen));
		}
#if DEBUG_PRINT
		DBG_OUT << "SampleIndicesNormalRandomly: from " << nLines << " lines ... done.\n";
#endif
	}

	constexpr size_t APPROX_BYTES_PER_VERTEX = 24;

	void SampleVerticesFromIndices(const char* start, const char* end, const std::vector<size_t>& indices, const size_t& updateThreshold, std::vector<pmp::Point>& result, IMB::IncrementalProgressTracker& tracker)
	{
#if DEBUG_PRINT
		DBG_OUT << "SampleVerticesFromIndices: ... \n";
#endif
		size_t localVertexCount = 0;

		for (const auto index : indices)
		{
			// Calculate an approximate position to jump to
			const char* cursor = start + index * APPROX_BYTES_PER_VERTEX;

			// Adjust cursor to the start of the next line if not already at a newline
			while (cursor < end && *cursor != '\n') cursor++;
			if (cursor < end) cursor++;  // Move past the newline to the start of the next line

			// Ensure the cursor is within valid range after adjustment
			if (cursor >= end) continue;


			// Check if the line starts with "v " indicating a vertex
			if (strncmp(cursor, "v ", 2) == 0)
			{
				cursor += 2; // skip "v "

				pmp::vec3 vec;
				char* tempCursor;
				vec[0] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				vec[1] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				vec[2] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				result.push_back(vec);
				localVertexCount++;

				// Check if it's time to update the tracker
				if (localVertexCount >= updateThreshold) 
				{
#if DEBUG_PRINT
					DBG_OUT << "SampleVerticesFromIndices: Time to update the tracker with " << localVertexCount << " collected vertices.\n";
#endif
					tracker.Update(localVertexCount);
					localVertexCount = 0; // Reset local count after update
				}
			}
			
			// Otherwise, skip to the next line
			while (*cursor != '\n' && cursor < end) cursor++;
		}

#if DEBUG_PRINT
		DBG_OUT << "SampleVerticesFromIndices: ... done.\n";
#endif
		// Ensure any remaining vertices are accounted for
		if (localVertexCount > 0)
		{
#if DEBUG_PRINT
			DBG_OUT << "SampleVerticesFromIndices: Time for a final tracker update with " << localVertexCount << " collected vertices.\n";
#endif
			tracker.Update(localVertexCount, true);
		}
	}
}

namespace IMB
{
	void SequentialVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
#if DEBUG_PRINT
		DBG_OUT << "SequentialVertexSamplingStrategy::Sample: ... \n";
#endif
		size_t localVertexCount = 0;
		const char* cursor = start;

		while (cursor < end)
		{
			// If the current line is incomplete, skip to the next line
			if (*cursor == '\n')
			{
				cursor++;
				continue;
			}

			// If it's a vertex, parse the three floats
			if (strncmp(cursor, "v ", 2) == 0)
			{
				cursor += 2; // skip "v "

				pmp::vec3 vec;
				char* tempCursor;
				vec[0] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				vec[1] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				vec[2] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				result.push_back(vec);
				localVertexCount++;

				// Check if it's time to update the tracker
				if (localVertexCount >= m_UpdateThreshold)
				{
#if DEBUG_PRINT
					DBG_OUT << "SequentialVertexSamplingStrategy::Sample: Time to update the tracker with " << localVertexCount << " collected vertices.\n";
#endif
					tracker.Update(localVertexCount);
					localVertexCount = 0; // Reset local count after update
				}
			}
			else
			{
				// Skip to the next line if the current line isn't recognized
				while (*cursor != '\n' && cursor < end)
					cursor++;
			}
		}
#if DEBUG_PRINT
		DBG_OUT << "SequentialVertexSamplingStrategy::Sample: ... done.\n";
#endif

		// Ensure any remaining vertices are accounted for
		if (localVertexCount > 0) 
		{
#if DEBUG_PRINT
			DBG_OUT << "SequentialVertexSamplingStrategy::Sample: Time for a final tracker update with " << localVertexCount << " collected vertices.\n";
#endif
			tracker.Update(localVertexCount, true);
		}
	}

	void UniformRandomVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		std::vector<size_t> indices;
		SampleIndicesUniformly(start, end, indices, seed);
		SampleVerticesFromIndices(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	void NormalRandomVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		std::vector<size_t> indices;
		SampleIndicesNormalRandomly(start, end, indices, seed);
		SampleVerticesFromIndices(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	void SoftmaxFeatureDetectingVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		throw std::runtime_error("SampleVerticesWithSoftmaxFeatureDectection Not implemented\n");
	}

	constexpr unsigned int FREQUENCY_UPDATE_MULTIPLIER = 1;

	VertexSamplingStrategy::VertexSamplingStrategy(const unsigned int& completionFrequency, const size_t& totalExpectedVertices)
		: m_UpdateThreshold(
			static_cast<size_t>(std::round(static_cast<double>(totalExpectedVertices) / 
				static_cast<double>(completionFrequency * FREQUENCY_UPDATE_MULTIPLIER))))
	{
#if DEBUG_PRINT
		DBG_OUT << "VertexSamplingStrategy::VertexSamplingStrategy: completionFrequency = " << completionFrequency << " jobs per file load.\n";
		DBG_OUT << "VertexSamplingStrategy::VertexSamplingStrategy: totalExpectedVertices = " << totalExpectedVertices << "\n";
		DBG_OUT << "VertexSamplingStrategy::VertexSamplingStrategy: m_UpdateThreshold = " << m_UpdateThreshold << " vertices.\n";
#endif
	}
} // namespace IMB