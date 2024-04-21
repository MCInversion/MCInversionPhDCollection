#pragma once
#include <optional>
#include <vector>
#include <memory>

#include "pmp/Types.h"

namespace IMB
{
	// forward declarations
	class IncrementalProgressTracker;
	class IncrementalMeshFileHandler;

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
		explicit VertexSamplingStrategy(const unsigned int& completionFrequency, const std::shared_ptr<IncrementalMeshFileHandler>& handler);

		virtual ~VertexSamplingStrategy() = default;

		virtual void Sample(const char* start, const char* end, 
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) = 0;

		[[nodiscard]] size_t GetVertexCountEstimate() const;

	protected:
		size_t m_UpdateThreshold;

		std::shared_ptr<IncrementalMeshFileHandler> m_FileHandler{ nullptr };
	};

	class SequentialVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		using VertexSamplingStrategy::VertexSamplingStrategy;

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	class UniformRandomVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		using VertexSamplingStrategy::VertexSamplingStrategy;

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	// TODO: verify usefulness
	class NormalRandomVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		using VertexSamplingStrategy::VertexSamplingStrategy;

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	class SoftmaxFeatureDetectingVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		using VertexSamplingStrategy::VertexSamplingStrategy;

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	inline [[nodiscard]] std::unique_ptr<VertexSamplingStrategy> GetVertexSelectionStrategy(const VertexSelectionType& vertSelType, const unsigned int& completionFrequency, const std::shared_ptr<IncrementalMeshFileHandler>& handler)
	{
		if (vertSelType == VertexSelectionType::Sequential)
			return std::make_unique<SequentialVertexSamplingStrategy>(completionFrequency, handler);
		if (vertSelType == VertexSelectionType::UniformRandom)
			return std::make_unique<UniformRandomVertexSamplingStrategy>(completionFrequency, handler);
		if (vertSelType == VertexSelectionType::NormalRandom)
			return std::make_unique<NormalRandomVertexSamplingStrategy>(completionFrequency, handler);
		return std::make_unique<SoftmaxFeatureDetectingVertexSamplingStrategy>(completionFrequency, handler);
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
