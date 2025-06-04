#include "VertexSamplingStrategies.h"

#include <numeric>
#include <random>
#include <unordered_set>

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
		for (size_t i = 0; i < totalCount; ++i) 
		{
			if (alreadyDrawn.contains(i))
				continue;
			
			buffer.push_back(i);
		}

		// random engine
		std::mt19937 gen(seed ? *seed : std::random_device{}());

		// Perform a "partial" Fisher–Yates on buffer[0..availableCount-1]
		// Swap only up to expectedCount
		for (size_t i = 0; i < expectedCount; ++i) 
		{
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

	[[nodiscard]] double ComputeStdDev(const std::vector<double>& values)
	{
		if (values.empty())
		{
			return 0.0;
		}
		double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
		double accum = 0.0;
		for (double v : values)
		{
			double diff = v - mean;
			accum += diff * diff;
		}
		return std::sqrt(accum / values.size());
	}

} // anonymous namespace

	//bool IsSupersetOfAlreadyDrawnIndices(
	//	const std::vector<size_t>& newIndices,
	//	const std::set<size_t>& alreadyDrawn)
	//{
	//	if (alreadyDrawn.empty())
	//		return false;

	//	// Build a hash set of the freshly generated indices:
	//	std::unordered_set<size_t> newSet;
	//	newSet.reserve(newIndices.size());
	//	for (size_t x : newIndices) {
	//		newSet.insert(x);
	//	}

	//	// Now check that every element of alreadyDrawn appears in newSet:
	//	return std::ranges::all_of(alreadyDrawn, [&newSet](const auto& i) { return newSet.contains(i); });
	//}

	void PrintSortedIndices(const std::vector<size_t>& indices, const size_t& amount)
	{
		std::cout << "sortedIndices:\n";

		std::vector<size_t> sortedIndices{ indices };
		std::ranges::sort(sortedIndices);

		for (int i = 0; i < amount; ++i)
			std::cout << "  [" << i << "]:  " << sortedIndices[i] << "\n";
	}

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
		size_t excludedCount = alreadyDrawn.size();
		size_t availableCount = (excludedCount < totalCount ? totalCount - excludedCount : 0);
		if (availableCount == 0 || expectedCount == 0)
		{
			resultIndices.clear();
			return;
		}

		size_t poolSize = expectedCount;
		if (!prevIterPts.empty())
		{
			// compute Dv for previous points to estimate sigma
			std::vector<double> prevDv;
			prevDv.reserve(prevIterPts.size());
			const size_t kNeighbors = 6;
			std::vector<uint32_t> ret_idx(kNeighbors);
			std::vector<pmp::Scalar> out_dist2(kNeighbors);
			auto prevIndex = Geometry::Get3DPointSearchIndex(prevIterPts);
			if (prevIndex)
			{
				auto& prevTree = prevIndex->tree;
				for (const auto& pt : prevIterPts)
				{
					pmp::Scalar queryPt[3] = { pt[0], pt[1], pt[2] };
					size_t found = prevTree.knnSearch(&queryPt[0], kNeighbors, ret_idx.data(), out_dist2.data());
					double sumLen = 0.0;
					if (found > 1)
					{
						for (size_t j = 1; j < found; ++j)
						{
							sumLen += std::sqrt(static_cast<double>(out_dist2[j]));
						}
						if (found < kNeighbors)
						{
							sumLen *= (static_cast<double>(kNeighbors - 1) / static_cast<double>(found - 1));
						}
					}
					prevDv.push_back(sumLen);
				}
			}
			else
			{
				prevDv.assign(prevIterPts.size(), 0.0);
			}
			double sigma = ComputeStdDev(prevDv);
			double Nk = static_cast<double>(prevIterPts.size());
			double Z = params.ConfidenceNormal;
			double rho = params.ConfidenceGM1Imax;
			double E = (Nk > 0.0 ? (Z * sigma / std::sqrt(Nk)) : 0.0);
			double ratio = (E > 0.0 ? (Z * sigma / E) : 1.0);
			double rawPool = (rho > 0.0 ? ((1.0 / rho) * (ratio * ratio)) : Nk);
			poolSize = static_cast<size_t>(std::ceil(rawPool));
		}

		if (poolSize > availableCount)
		{
			poolSize = availableCount;
		}

		std::vector<size_t> candidateIndices;
		candidateIndices.reserve(poolSize);
		RandomSampleIndices(poolSize, totalCount, candidateIndices, alreadyDrawn, seed);

		if (prevIterPts.empty())
		{
			std::mt19937_64 rng(seed ? *seed : std::random_device{}());
			std::shuffle(candidateIndices.begin(), candidateIndices.end(), rng);
			size_t takeCount = std::min(expectedCount, candidateIndices.size());
			resultIndices.assign(candidateIndices.begin(), candidateIndices.begin() + takeCount);
			return;
		}

		std::vector<pmp::Point> candidatePts;
		candidatePts.reserve(candidateIndices.size());
		for (size_t id : candidateIndices)
		{
			candidatePts.push_back(fileHandler->SampleSinglePoint(start, end, id));
		}

		auto fullIndex = Geometry::Get3DPointSearchIndex(candidatePts);
		if (!fullIndex)
		{
			std::mt19937_64 rng(seed ? *seed : std::random_device{}());
			std::shuffle(candidateIndices.begin(), candidateIndices.end(), rng);
			size_t takeCount = std::min(expectedCount, candidateIndices.size());
			resultIndices.assign(candidateIndices.begin(), candidateIndices.begin() + takeCount);
			return;
		}
		auto& fullTree = fullIndex->tree;

		const size_t kNeighbors = 10;
		std::vector<double> Dv(candidatePts.size(), 0.0);
		std::vector<uint32_t> ret_idx(kNeighbors);
		std::vector<pmp::Scalar> out_dist2(kNeighbors);
		for (size_t i = 0; i < candidatePts.size(); ++i)
		{
			pmp::Scalar queryPt[3] = { candidatePts[i][0], candidatePts[i][1], candidatePts[i][2] };
			size_t found = fullTree.knnSearch(&queryPt[0], kNeighbors, ret_idx.data(), out_dist2.data());
			double sumLen = 0.0;
			if (found > 1)
			{
				for (size_t j = 1; j < found; ++j)
				{
					sumLen += std::sqrt(static_cast<double>(out_dist2[j]));
				}
				if (found < kNeighbors)
				{
					sumLen *= (static_cast<double>(kNeighbors - 1) / static_cast<double>(found - 1));
				}
			}
			Dv[i] = sumLen;
		}

		double D_target = params.TargetVertexDensity;
		if (D_target < 0.0)
		{
			double accum = 0.0;
			for (double v : Dv)
			{
				accum += v;
			}
			D_target = (Dv.empty() ? 0.0 : (accum / static_cast<double>(Dv.size())));
		}

		std::mt19937_64 rng(seed ? *seed : std::random_device{}());
		std::uniform_real_distribution<double> uniformReal(std::numeric_limits<double>::min(), 1.0);

		std::vector<KeyIdx> keys;
		keys.reserve(candidateIndices.size());
		for (size_t i = 0; i < candidateIndices.size(); ++i)
		{
			double score = params.DistributionMetricGradientMultiplier * (Dv[i] - D_target);
			double w_i = std::exp(score);
			if (w_i < 1e-300)
			{
				w_i = 1e-300;
			}
			else if (w_i > 1e300)
			{
				w_i = 1e300;
			}
			double u = uniformReal(rng);
			double lnU = std::log(u);
			double key = std::exp(lnU / (params.SelectionSharpness * w_i));
			keys.push_back({ key, candidateIndices[i] });
		}

		std::ranges::sort(keys, [](const KeyIdx& a, const KeyIdx& b) { return a.key > b.key; });

		size_t keepCount = std::min(expectedCount, keys.size());
		resultIndices.resize(keepCount);
		for (size_t j = 0; j < keepCount; ++j)
		{
			resultIndices[j] = keys[j].idx;
		}

		if (HasEnoughMemoryForFisherYates(keepCount))
		{
			for (size_t i = keepCount - 1; i > 0; --i)
			{
				std::uniform_int_distribution<size_t> distrib(0, i);
				size_t swap_j = distrib(rng);
				std::swap(resultIndices[i], resultIndices[swap_j]);
			}
		}
		else
		{
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
		size_t total = m_FileHandler->GetLocalVertexCountEstimate(start, end);
		size_t already = result.size();
		size_t remaining = (total > already ? total - already : 0);
		size_t toSample = already == 0 ? m_MinVertexCount : std::min(remaining, m_UpdateThreshold);
		if (toSample == 0)
			return;

		std::vector<size_t> indices;
		SampleIndicesSequentially(toSample, indices);

		m_AlreadyDrawnIndices.insert(indices.begin(), indices.end());
		m_FileHandler->Sample(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	void UniformRandomVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		size_t total = m_FileHandler->GetLocalVertexCountEstimate(start, end);
		size_t already = result.size();
		size_t remaining = (total > already ? total - already : 0);
		size_t toSample = already == 0 ? m_MinVertexCount : std::min(remaining, m_UpdateThreshold);
		if (toSample == 0) 
			return;

		std::vector<size_t> indices;
		RandomSampleIndices(toSample, total, indices, m_AlreadyDrawnIndices, seed);

		std::cout << "----------------------------------------------\n";
		std::cout << "UniformRandomVertexSamplingStrategy::Sample:\n";
		PrintSortedIndices(indices, 8);
		std::cout << "----------------------------------------------\n";

		m_AlreadyDrawnIndices.insert(indices.begin(), indices.end());
		m_FileHandler->Sample(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	void SoftmaxUniformVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		size_t total = m_FileHandler->GetLocalVertexCountEstimate(start, end);
		size_t already = result.size();
		size_t remaining = (total > already ? total - already : 0);
		size_t toSample = already == 0 ? m_MinVertexCount : std::min(remaining, m_UpdateThreshold);
		if (toSample == 0) 
			return;

		std::vector<size_t> indices;
		SoftmaxUniformGeometricSampleIndices(toSample, total, indices, seed, m_AlreadyDrawnIndices, result, Params, m_FileHandler.get(), start, end);
		
		std::cout << "----------------------------------------------\n";
		std::cout << "SoftmaxUniformVertexSamplingStrategy::Sample:\n";
		PrintSortedIndices(indices, 8);
		std::cout << "----------------------------------------------\n";
		
		m_AlreadyDrawnIndices.insert(indices.begin(), indices.end());
		m_FileHandler->Sample(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	void SoftmaxFeatureDetectingVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		size_t total = m_FileHandler->GetLocalVertexCountEstimate(start, end);
		size_t already = result.size();
		size_t remaining = (total > already ? total - already : 0);
		size_t toSample = already == 0 ? m_MinVertexCount : std::min(remaining, m_UpdateThreshold);
		if (toSample == 0)
			return;

		std::vector<size_t> indices;
		// TODO: implement SoftmaxAdaptiveGeometricSampleIndices
		RandomSampleIndices(toSample, total, indices, m_AlreadyDrawnIndices, seed);

		m_AlreadyDrawnIndices.insert(indices.begin(), indices.end());
		m_FileHandler->Sample(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	constexpr unsigned int FREQUENCY_UPDATE_MULTIPLIER = 1;
	constexpr double MIN_VERTEX_FRACTION = 0.01;

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