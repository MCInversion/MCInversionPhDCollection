#pragma once
#include <optional>
#include <vector>

namespace IMB
{
	// forward declarations
	class IncrementalProgressTracker;

	class VertexSamplingStrategy
	{
	public:
		VertexSamplingStrategy(double completionFrequency, size_t totalExpectedVertices)
			: m_UpdateThreshold(static_cast<size_t>(totalExpectedVertices* completionFrequency))
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
		SequentialVertexSamplingStrategy(double completionFrequency, size_t totalExpectedVertices)
			: VertexSamplingStrategy(completionFrequency, totalExpectedVertices) {}

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	class UniformRandomVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		UniformRandomVertexSamplingStrategy(double completionFrequency, size_t totalExpectedVertices)
			: VertexSamplingStrategy(completionFrequency, totalExpectedVertices) {}

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	// TODO: verify usefulness
	class NormalRandomVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		NormalRandomVertexSamplingStrategy(double completionFrequency, size_t totalExpectedVertices)
			: VertexSamplingStrategy(completionFrequency, totalExpectedVertices) {}

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	class SoftmaxFeatureDetectingVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		SoftmaxFeatureDetectingVertexSamplingStrategy(double completionFrequency, size_t totalExpectedVertices)
			: VertexSamplingStrategy(completionFrequency, totalExpectedVertices) {}

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

} // namespace IMB
