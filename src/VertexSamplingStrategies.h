#pragma once
#include <optional>
#include <vector>
#include <memory>

#include "pmp/Types.h"

namespace IMB
{
	// forward declarations
	class IncrementalProgressTracker;

	/// \brief enumerator for mesh simplification function type.
	enum class [[nodiscard]] VertexSelectionType
	{
		Sequential = 0, //>! selects vertices sequentially.
		UniformRandom = 1, //>! selects vertices uniformly at random.
		NormalRandom = 2, //>! selects vertices with a normal distribution.
		SoftMaxFeatureDetecting = 3, //>! selects vertices using a softmax function with feature detection.
	};

	class VertexSamplingStrategy
	{
	public:
		VertexSamplingStrategy(const unsigned int& completionFrequency, const size_t& totalExpectedVertices);

		virtual ~VertexSamplingStrategy() = default;

		virtual void Sample(const char* start, const char* end, 
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) = 0;
	protected:

		size_t m_UpdateThreshold;
	};

	class SequentialVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		SequentialVertexSamplingStrategy(const unsigned int& completionFrequency, const size_t& totalExpectedVertices)
			: VertexSamplingStrategy(completionFrequency, totalExpectedVertices) {}

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	class UniformRandomVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		UniformRandomVertexSamplingStrategy(const unsigned int& completionFrequency, const size_t& totalExpectedVertices)
			: VertexSamplingStrategy(completionFrequency, totalExpectedVertices) {}

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	// TODO: verify usefulness
	class NormalRandomVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		NormalRandomVertexSamplingStrategy(const unsigned int& completionFrequency, const size_t& totalExpectedVertices)
			: VertexSamplingStrategy(completionFrequency, totalExpectedVertices) {}

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	class SoftmaxFeatureDetectingVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		SoftmaxFeatureDetectingVertexSamplingStrategy(const unsigned int& completionFrequency, const size_t& totalExpectedVertices)
			: VertexSamplingStrategy(completionFrequency, totalExpectedVertices) {}

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	inline [[nodiscard]] std::unique_ptr<VertexSamplingStrategy> GetVertexSelectionStrategy(const VertexSelectionType& vertSelType, const unsigned int& completionFrequency, const size_t& totalExpectedVertices)
	{
		if (vertSelType == VertexSelectionType::Sequential)
			return std::make_unique<SequentialVertexSamplingStrategy>(completionFrequency, totalExpectedVertices);
		if (vertSelType == VertexSelectionType::UniformRandom)
			return std::make_unique<UniformRandomVertexSamplingStrategy>(completionFrequency, totalExpectedVertices);
		if (vertSelType == VertexSelectionType::NormalRandom)
			return std::make_unique<NormalRandomVertexSamplingStrategy>(completionFrequency, totalExpectedVertices);
		return std::make_unique<SoftmaxFeatureDetectingVertexSamplingStrategy>(completionFrequency, totalExpectedVertices);
	}

	inline [[nodiscard]] std::string GetVertexSelectionStrategyName(const VertexSelectionType& vertSelType)
	{
		if (vertSelType == VertexSelectionType::Sequential)
			return "VertexSelectionType::Sequential";
		if (vertSelType == VertexSelectionType::UniformRandom)
			return "VertexSelectionType::UniformRandom";
		if (vertSelType == VertexSelectionType::NormalRandom)
			return "VertexSelectionType::NormalRandom";
		return "VertexSelectionType::SoftMaxFeatureDetecting";
	}

} // namespace IMB
