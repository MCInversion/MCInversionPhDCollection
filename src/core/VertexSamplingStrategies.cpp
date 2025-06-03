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

	/// \brief Generates expectedCount random indices from 0 to expectedCount - 1 into resultIndices from the uniform distribution type.
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
		resultIndices.resize(expectedCount);

		// Create a generator with the given seed or a random device
		std::mt19937 gen(seed ? *seed : std::random_device{}());

		// Fill resultIndices with 0, 1, ..., expectedCount - 1
		std::iota(resultIndices.begin(), resultIndices.end(), 0);

		// Perform Fisher-Yates shuffle on resultIndices
		for (size_t i = expectedCount - 1; i > 0; --i)
		{
			std::uniform_int_distribution<size_t> distrib(0, i);
			size_t j = distrib(gen);
			std::swap(resultIndices[i], resultIndices[j]);
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

	void SoftmaxUniformGeometricSampleIndices(const size_t& expectedCount, std::vector<size_t>& resultIndices, 
		const std::optional<unsigned int>& seed, const std::vector<pmp::Point>& prevIterPts,
		const IMB::GeometricSamplingParams& params, 
		IMB::IncrementalMeshFileHandler* fileHandler,
		const char* start, const char* end)
	{
		if (prevIterPts.empty())
		{
			// initial sample must be random
			RandomSampleIndices(expectedCount, resultIndices, seed);
			return;
		}

		// Build a kdtree over prevIterPts
		auto pt3DIndex = Geometry::Get3DPointSearchIndex(prevIterPts);
		if (!pt3DIndex)
		{
#if DEBUG_PRINT
			DBG_OUT << "SoftmaxUniformGeometricSampleIndices: failed to build kdtree on prevIterPts\n";
#endif
		}

		auto& kdTree = pt3DIndex->tree;

		// Evaluate vertex weights and draw indices
		std::mt19937_64 rng(seed ? *seed : std::random_device{}());
		std::uniform_real_distribution<double> uniformReal(std::numeric_limits<double>::min(), 1.0);
		std::vector<KeyIdx> keys;
		keys.reserve(expectedCount);
		for (size_t i = 0; i < expectedCount; ++i)
		{
			const auto v_i = fileHandler->SampleSinglePoint(start, end, i);
			std::vector<uint32_t>    nnIdx(1);
			std::vector<pmp::Scalar> nnDistSq(1);
			{
				// find nearest neighbor from existing pts
				pmp::Scalar queryPt[3] = { v_i[0], v_i[1], v_i[2] };
				kdTree.knnSearch(&queryPt[0], 1, nnIdx.data(), nnDistSq.data());
			}

			// score(v_i) = DistributionMetricGradientMultiplier * (d_i - TargetVertexDensity)
			const double d_i = std::sqrt(static_cast<double>(nnDistSq[0]));
			const double score = params.DistributionMetricGradientMultiplier * (d_i - params.TargetVertexDensity);

			// clamp weight w[i] = exp(score(v_i))
			double w_i = std::exp(score);
			if (w_i < 1e-300) w_i = 1e-300;
			if (w_i > 1e+300) w_i = 1e+300;

			// Draw U in (0,1), compute key = exp( ln(U) / w_i )
			const double u = uniformReal(rng);
			const double lnU = std::log(u);
			const double key = std::exp(lnU / w_i);
			keys.emplace_back(key, i);
		}

		// Extract sorted indices into resultIndices
		resultIndices.resize(expectedCount);
		std::ranges::sort(keys, [](const KeyIdx& a, const KeyIdx& b) { return a.key > b.key; });
		for (size_t j = 0; j < expectedCount; ++j)
		{
			resultIndices[j] = keys[j].idx;
		}

		// Fisher–Yates shuffle among all expectedCount to break any ties
		for (size_t i = expectedCount - 1; i > 0; --i)
		{
			std::uniform_int_distribution<size_t> distrib(0, i);
			size_t swap_j = distrib(rng);
			std::swap(resultIndices[i], resultIndices[swap_j]);
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
		SampleIndicesSequentially(m_FileHandler->GetLocalVertexCountEstimate(start, end), indices);
		m_FileHandler->Sample(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	void UniformRandomVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		std::vector<size_t> indices;
		RandomSampleIndices(m_FileHandler->GetLocalVertexCountEstimate(start, end), indices, seed);
		m_FileHandler->Sample(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	void SoftmaxUniformVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		std::vector<size_t> indices;
		const auto nVertsToSample = m_FileHandler->GetLocalVertexCountEstimate(start, end);
		SoftmaxUniformGeometricSampleIndices(nVertsToSample, indices, seed, result, Params, m_FileHandler.get(), start, end);
		m_FileHandler->Sample(start, end, indices, m_UpdateThreshold, result, tracker);
	}

	void SoftmaxFeatureDetectingVertexSamplingStrategy::Sample(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker)
	{
		std::vector<size_t> indices;
		// TODO: implement index sampling
		RandomSampleIndices(m_FileHandler->GetLocalVertexCountEstimate(start, end), indices, seed);
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
		m_UpdateThreshold = static_cast<size_t>(std::round(static_cast<double>(m_VertexCap) /
			static_cast<double>(completionFrequency * FREQUENCY_UPDATE_MULTIPLIER)));
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