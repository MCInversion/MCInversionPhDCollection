#include "VertexSamplingStrategies.h"

#include <numeric>
#include <random>
//#include <unordered_set>

#include "IncrementalProgressUtils.h"
#include "IncrementalMeshFileHandler.h"

#include "utils/IncrementalUtils.h"

#include "geometry/GeometryConversionUtils.h"


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

	// \brief Fills resultIndices with expectedCount unique random indices from [0..totalCount-1]
	///       by performing a "partial" Fisher-Yates shuffle. Requires O(totalCount) memory.
	/// \param expectedCount  How many distinct indices to sample.
	/// \param totalCount     The range of indices is [0..totalCount-1].
	/// \param resultIndices  Output vector: will be resized to expectedCount and filled with the sample.
	/// \param alreadyDrawn   The set of already drawn indices.
	/// \param seed           Optional seed; if not provided, a random_device-seeded mt19937 is used.
	void FisherYatesSampleIndices(
		const size_t& expectedCount,
		const size_t& totalCount,
		std::vector<size_t>& resultIndices,
		const std::set<size_t>& alreadyDrawn,
		const std::optional<unsigned int>& seed)
	{
		// Count how many indices remain available
		size_t excludedCount = alreadyDrawn.size();
		if (excludedCount >= totalCount) {
			std::cerr << "FisherYatesSampleIndices: no available indices (all already drawn)\n";
			return;
		}
		size_t availableCount = totalCount - excludedCount;
		if (expectedCount > availableCount) {
			std::cerr << "FisherYatesSampleIndices: expectedCount > available indices\n";
			return;
		}

		resultIndices.clear();
		resultIndices.reserve(expectedCount);

		// Build a buffer of all indices not in alreadyDrawn
		std::vector<size_t> buffer;
		buffer.reserve(availableCount);
		for (size_t i = 0; i < totalCount; ++i) {
			if (!alreadyDrawn.contains(i)) {
				buffer.push_back(i);
			}
		}

		// random engine
		std::mt19937 gen(seed ? *seed : std::random_device{}());

		// Perform a "partial" Fisher–Yates on buffer[0..availableCount-1]
		// Swap only up to expectedCount
		for (size_t i = 0; i < expectedCount; ++i) {
			std::uniform_int_distribution<size_t> distrib(i, availableCount - 1);
			size_t j = distrib(gen);
			std::swap(buffer[i], buffer[j]);
		}

		// Copy out the first expectedCount entries
		resultIndices.assign(buffer.begin(), buffer.begin() + expectedCount);
	}

	/// \brief Fills resultIndices with expectedCount unique random indices from [0..totalCount-1]
	///        by using reservoir sampling. Requires only O(expectedCount) memory.
	/// \param expectedCount  How many distinct indices to sample.
	/// \param totalCount     The range of indices is [0..totalCount-1].
	/// \param resultIndices  Output vector: will be resized to expectedCount and filled with the sample.
	/// \param alreadyDrawn   The set of already drawn indices.
	/// \param seed           Optional seed; if not provided, a random_device-seeded mt19937 is used.
	void ReservoirSampleIndices(
		const size_t& expectedCount,
		const size_t& totalCount,
		std::vector<size_t>& resultIndices,
		const std::set<size_t>& alreadyDrawn,
		const std::optional<unsigned int>& seed)
	{
		// Count how many indices remain available
		size_t excludedCount = alreadyDrawn.size();
		if (excludedCount >= totalCount) {
			std::cerr << "ReservoirSampleIndices: no available indices (all already drawn)\n";
			return;
		}
		size_t availableCount = totalCount - excludedCount;
		if (expectedCount > availableCount) {
			std::cerr << "ReservoirSampleIndices: expectedCount > available indices\n";
			return;
		}

		resultIndices.clear();
		resultIndices.reserve(expectedCount);

		std::mt19937 gen(seed ? *seed : std::random_device{}());
		size_t seen = 0;

		// Iterate through the entire [0..totalCount-1], skipping alreadyDrawn
		for (size_t i = 0; i < totalCount; ++i) {
			if (alreadyDrawn.contains(i)) {
				continue;
			}
			if (seen < expectedCount) {
				// Fill initial reservoir
				resultIndices.push_back(i);
			}
			else {
				// Decide whether to replace one of the reservoir slots
				std::uniform_int_distribution<size_t> distrib(0, seen);
				size_t j = distrib(gen);
				if (j < expectedCount) {
					resultIndices[j] = i;
				}
			}
			++seen;
		}
		// At this point, `seen == availableCount >= expectedCount`
	}

	/// \brief Generates expectedCount random indices from 0 to totalCount - 1 into resultIndices from the uniform distribution type.
	void RandomSampleIndices(
		const size_t& expectedCount, 
		const size_t& totalCount, 
		std::vector<size_t>& resultIndices, 
		const std::set<size_t>& alreadyDrawn,
		const std::optional<unsigned int>& seed)
	{
		// Determine how many indices remain available
		size_t excludedCount = alreadyDrawn.size();
		if (excludedCount >= totalCount) 
		{
			std::cerr << "ReservoirSampleIndices: no available indices (all already drawn)\n";
			return;
		}
		size_t availableCount = totalCount - excludedCount;
		if (expectedCount > availableCount) 
		{
			std::cerr << "ReservoirSampleIndices: expectedCount > available indices\n";
			return;
		}

#if DEBUG_PRINT
		DBG_OUT << "RandomSampleIndices: Starting...\n";
#endif

		if (HasEnoughMemoryForFisherYates(totalCount))
		{
			FisherYatesSampleIndices(expectedCount, totalCount, resultIndices, alreadyDrawn, seed);
		}
		else {
			ReservoirSampleIndices(expectedCount, totalCount, resultIndices, alreadyDrawn, seed);
		}

#if DEBUG_PRINT
		DBG_OUT << "RandomSampleIndices: Finished generating " << expectedCount << " unique indices.\n";
#endif
	}

} // anonymous namespace

	// .........................................
	// TODO: move to anonymous namespace vvvvv
	// .........................................
	struct KeyIdx {
		double key;
		size_t idx;
	};

	void SoftmaxUniformGeometricSampleIndices(
		const size_t& expectedCount, const size_t& totalCount,
		std::vector<size_t>& resultIndices, 
		const std::optional<unsigned int>& seed, 
		const std::set<size_t>& alreadyDrawn,
		const std::vector<pmp::Point>& prevIterPts,
		const IMB::GeometricSamplingParams& params, 
		IMB::IncrementalMeshFileHandler* fileHandler,
		const char* start, const char* end)
	{
		std::vector<size_t> candidateIndices;
		RandomSampleIndices(expectedCount, totalCount, candidateIndices, alreadyDrawn, seed);

		if (prevIterPts.empty())
		{
			// initial sample must be random
			resultIndices = candidateIndices;
			return;
		}

		// Pick a random subset of vertices of size expectedCount
		std::vector<pmp::Point> candidatePts;
		candidatePts.reserve(expectedCount);
		for (const auto& candidateId : candidateIndices)
		{
			candidatePts.push_back(fileHandler->SampleSinglePoint(start, end, candidateId));
		}
		// Build a KD-tree over candidatePts
		auto fullIndex = Geometry::Get3DPointSearchIndex(candidatePts);
		if (!fullIndex) 
		{
			// kd-tree construction error.
			resultIndices = candidateIndices;
			return;
		}
		auto& fullTree = fullIndex->tree;

		// For each candidate point, compute Dv = sum of distances to its k nearest neighbors
		const size_t kNeighbors = 6;
		std::vector<double> Dv(expectedCount, 0.0);
		std::vector<uint32_t>    ret_idx(kNeighbors);
		std::vector<pmp::Scalar> out_dist2(kNeighbors);

		for (size_t i = 0; i < expectedCount; ++i) 
		{
			pmp::Scalar queryPt[3] = { candidatePts[i][0], candidatePts[i][1], candidatePts[i][2] };
			size_t found = fullTree.knnSearch(&queryPt[0], kNeighbors, ret_idx.data(), out_dist2.data());
			if (found > 1) 
			{
				double sumLen = 0.0;
				for (size_t j = 1; j < found; ++j) {
					sumLen += std::sqrt(static_cast<double>(out_dist2[j]));
				}
				if (found < kNeighbors) {
					sumLen *= (static_cast<double>(kNeighbors - 1) / static_cast<double>(found - 1));
				}
				Dv[i] = sumLen;
			}
			else 
			{
				Dv[i] = 0.0;
			}
		}

		// Determine D_target
		double D_target = params.TargetVertexDensity;
		if (D_target < 0.0) 
		{
			double accum = 0.0;
			for (double v : Dv) accum += v;
			D_target = accum / static_cast<double>(expectedCount);
		}

		// Compute softmax-like keys and keep (key, originalIndex) pairs
		std::mt19937_64 rng(seed ? *seed : std::random_device{}());
		std::uniform_real_distribution<double> uniformReal(std::numeric_limits<double>::min(), 1.0);

		std::vector<KeyIdx> keys;
		keys.reserve(expectedCount);

		for (size_t i = 0; i < expectedCount; ++i) 
		{
			double score = params.DistributionMetricGradientMultiplier * (Dv[i] - D_target);
			double w_i = std::exp(score);
			if (w_i < 1e-300) w_i = 1e-300;
			if (w_i > 1e+300) w_i = 1e+300;

			double u = uniformReal(rng);
			double lnU = std::log(u);
			double key = std::exp(lnU / w_i);

			keys.emplace_back(key, candidateIndices[i]);
		}

		// Sort all candidates by descending key
		std::ranges::sort(keys, [](const KeyIdx& a, const KeyIdx& b) { return a.key > b.key; });

		// Extract the top expectedCount indices into resultIndices
		resultIndices.resize(expectedCount);
		for (size_t j = 0; j < expectedCount; ++j) { resultIndices[j] = keys[j].idx; }

		if (HasEnoughMemoryForFisherYates(expectedCount))
		{
			// Full Fisher–Yates on the `expectedCount` slots:
			for (size_t i = expectedCount - 1; i > 0; --i)
			{
				std::uniform_int_distribution<size_t> distrib(0, i);
				size_t swap_j = distrib(rng);
				std::swap(resultIndices[i], resultIndices[swap_j]);
			}
		}
		else
		{
			// Fallback shuffle that never allocates any extra buffer.
			// std::shuffle is guaranteed to work in-place, touching only resultIndices itself.
			std::shuffle(resultIndices.begin(), resultIndices.end(), rng);
		}
	}
	// .........................................
	// TODO: move to anonymous namespace ^^^^^
	// .........................................

namespace IMB
{
	void SequentialVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		std::vector<size_t> indices;
		size_t total = m_FileHandler->GetLocalVertexCountEstimate(start, end);
		size_t already = result.size();
		size_t remaining = (total > already ? total - already : 0);
		size_t toSample = std::min(remaining, m_UpdateThreshold);
		if (toSample == 0)
			return;
		SampleIndicesSequentially(toSample, indices);
		m_AlreadyDrawnIndices.insert(indices.begin(), indices.end());
		m_FileHandler->Sample(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	void UniformRandomVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		std::vector<size_t> indices;
		size_t total = m_FileHandler->GetLocalVertexCountEstimate(start, end);
		size_t already = result.size();
		size_t remaining = (total > already ? total - already : 0);
		size_t toSample = std::min(remaining, m_UpdateThreshold);
		if (toSample == 0) 
			return;
		RandomSampleIndices(toSample, total, indices, m_AlreadyDrawnIndices, seed);
		m_AlreadyDrawnIndices.insert(indices.begin(), indices.end());
		m_FileHandler->Sample(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	void SoftmaxUniformVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		std::vector<size_t> indices;
		size_t total = m_FileHandler->GetLocalVertexCountEstimate(start, end);
		size_t already = result.size();
		size_t remaining = (total > already ? total - already : 0);
		size_t toSample = std::min(remaining, m_UpdateThreshold);
		if (toSample == 0) 
			return;
		SoftmaxUniformGeometricSampleIndices(toSample, total, indices, seed, m_AlreadyDrawnIndices, result, Params, m_FileHandler.get(), start, end);
		m_AlreadyDrawnIndices.insert(indices.begin(), indices.end());
		m_FileHandler->Sample(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	void SoftmaxFeatureDetectingVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		std::vector<size_t> indices;
		size_t total = m_FileHandler->GetLocalVertexCountEstimate(start, end);
		size_t already = result.size();
		size_t remaining = (total > already ? total - already : 0);
		size_t toSample = std::min(remaining, m_UpdateThreshold);
		if (toSample == 0)
			return;
		RandomSampleIndices(toSample, total, indices, m_AlreadyDrawnIndices, seed);
		m_AlreadyDrawnIndices.insert(indices.begin(), indices.end());
		m_FileHandler->Sample(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	constexpr unsigned int FREQUENCY_UPDATE_MULTIPLIER = 2;

	constexpr double MIN_VERTEX_FRACTION = 0.001;

	VertexSamplingStrategy::VertexSamplingStrategy(const unsigned int& completionFrequency, const size_t& maxVertexCount, const std::shared_ptr<IncrementalMeshFileHandler>& handler)
	{
		m_FileHandler = handler;
		if (!m_FileHandler)
		{
			throw std::invalid_argument("VertexSamplingStrategy::VertexSamplingStrategy: m_FileHandler could not be created!\n");
		}
		m_VertexCap = std::min(m_FileHandler->GetGlobalVertexCountEstimate(), maxVertexCount);
		m_MinVertexCount = static_cast<size_t>(m_VertexCap * MIN_VERTEX_FRACTION);
		m_UpdateThreshold = static_cast<size_t>(std::round(static_cast<double>(m_VertexCap) / static_cast<double>(completionFrequency * FREQUENCY_UPDATE_MULTIPLIER)));
#if DEBUG_PRINT
		DBG_OUT << "VertexSamplingStrategy::VertexSamplingStrategy: completionFrequency = " << completionFrequency << " jobs per file load.\n";
		DBG_OUT << "VertexSamplingStrategy::VertexSamplingStrategy: m_MinVertexCount = " << m_MinVertexCount << " vertices.\n";
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