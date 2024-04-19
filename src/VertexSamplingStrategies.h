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
		VertexSamplingStrategy(const double& completionFrequency, const size_t& totalExpectedVertices)
			: m_UpdateThreshold(static_cast<size_t>(totalExpectedVertices * completionFrequency))
		{ }

		virtual ~VertexSamplingStrategy() = default;

		virtual void Sample(const char* start, const char* end, 
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) = 0;
	protected:

		size_t m_UpdateThreshold;
	};

	class SequentialVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		SequentialVertexSamplingStrategy(const double& completionFrequency, const size_t& totalExpectedVertices)
			: VertexSamplingStrategy(completionFrequency, totalExpectedVertices) {}

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	class UniformRandomVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		UniformRandomVertexSamplingStrategy(const double& completionFrequency, const size_t& totalExpectedVertices)
			: VertexSamplingStrategy(completionFrequency, totalExpectedVertices) {}

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	// TODO: verify usefulness
	class NormalRandomVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		NormalRandomVertexSamplingStrategy(const double& completionFrequency, const size_t& totalExpectedVertices)
			: VertexSamplingStrategy(completionFrequency, totalExpectedVertices) {}

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	class SoftmaxFeatureDetectingVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		SoftmaxFeatureDetectingVertexSamplingStrategy(const double& completionFrequency, const size_t& totalExpectedVertices)
			: VertexSamplingStrategy(completionFrequency, totalExpectedVertices) {}

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	inline [[nodiscard]] std::unique_ptr<VertexSamplingStrategy> GetVertexSelectionStrategy(const VertexSelectionType& vertSelType, const double& completionFrequency, const size_t& totalExpectedVertices)
	{
		if (vertSelType == VertexSelectionType::Sequential)
			return std::make_unique<SequentialVertexSamplingStrategy>(completionFrequency, totalExpectedVertices);
		if (vertSelType == VertexSelectionType::UniformRandom)
			return std::make_unique<UniformRandomVertexSamplingStrategy>(completionFrequency, totalExpectedVertices);
		if (vertSelType == VertexSelectionType::NormalRandom)
			return std::make_unique<NormalRandomVertexSamplingStrategy>(completionFrequency, totalExpectedVertices);
		return std::make_unique<SoftmaxFeatureDetectingVertexSamplingStrategy>(completionFrequency, totalExpectedVertices);
	}

} // namespace IMB
