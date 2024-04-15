#include "pmp/Types.h"

#include "VertexSamplingStrategies.h"

#include <random>

#include "IncrementalProgressTracker.h"

namespace
{
	void SampleIndicesUniformly(const char* start, const char* end, std::vector<size_t>& resultIndices, const std::optional<unsigned int>& seed)
	{
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
	}

	void SampleIndicesNormalRandomly(const char* start, const char* end, std::vector<size_t>& resultIndices, const std::optional<unsigned int>& seed)
	{
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
	}

	void SampleVerticesFromIndices(const char* start, const char* end, const std::vector<size_t>& indices, std::vector<pmp::Point>& result, IMB::IncrementalProgressTracker& tracker)
	{
		const char* cursor = start + indices[0];
		for (size_t i = 0; i < indices.size(); i++)
		{
			// skip if indices not in range
			if (indices[i] < 0 || indices[i] >= end - start)
				continue;

			// If the current line is incomplete, skip to the next line
			if (*cursor == '\n')
			{
				cursor = start + indices[(i + 1) % indices.size()];
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
				tracker.Update(result.size());
			}
			else
			{
				// Skip to the next line if the current line isn't recognized
				while (*cursor != '\n' && cursor < end)
					cursor = start + indices[(i + 1) % indices.size()];
			}
		}
	}
}

namespace IMB
{
	void SequentialVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
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
				tracker.Update(result.size());
			}
			else
			{
				// Skip to the next line if the current line isn't recognized
				while (*cursor != '\n' && cursor < end)
					cursor++;
			}
		}
	}

	void UniformRandomVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		std::vector<size_t> indices;
		SampleIndicesUniformly(start, end, indices, seed);
		SampleVerticesFromIndices(start, end, indices, result, tracker);
	}

	void NormalRandomVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		std::vector<size_t> indices;
		SampleIndicesNormalRandomly(start, end, indices, seed);
		SampleVerticesFromIndices(start, end, indices, result, tracker);
	}

	void SoftmaxFeatureDetectingVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		throw std::runtime_error("SampleVerticesWithSoftmaxFeatureDectection Not implemented\n");
	}
} // namespace IMB