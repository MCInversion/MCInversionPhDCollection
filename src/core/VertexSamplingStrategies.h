#pragma once
#include <optional>
#include <vector>
#include <memory>
#include <set>

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
		SoftMaxUniform = 2, //>! selects vertices using a softmax function with uniform distribution across the surface.
		SoftMaxFeatureDetecting = 3, //>! selects vertices using a softmax function with feature detection.
		PoissonDisc = 4, //>! selects vertices using Poisson disc sampling
	};

	class VertexSamplingStrategy
	{
	public:
		explicit VertexSamplingStrategy(const unsigned int& completionFrequency, const size_t& maxVertexCount, const std::shared_ptr<IncrementalMeshFileHandler>& handler);

		virtual ~VertexSamplingStrategy() = default;

		virtual void Sample(const char* start, const char* end, 
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) = 0;

		[[nodiscard]] size_t GetVertexCountEstimate() const;

		[[nodiscard]] size_t GetVertexCap() const
		{
			return m_VertexCap;
		}

		[[nodiscard]] size_t GetMinVertexCount() const
		{
			return m_MinVertexCount;
		}

	protected:
		size_t m_UpdateThreshold;
		size_t m_VertexCap;
		size_t m_MinVertexCount;

		std::shared_ptr<IncrementalMeshFileHandler> m_FileHandler{ nullptr };
		std::set<size_t> m_AlreadyDrawnIndices;
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

	// -------------------------------------------------------------
	/// \brief A wrapper for geometric sampling parameters
	// -------------------------------------------------------------
	struct GeometricSamplingParams
	{
		double ErrorMetricGradientMultiplier{ 1.0 };   // >! alpha: error metric gradient multiplier
		double DistributionMetricGradientMultiplier{ 1.0 };   // >! beta: vertex distribution metric gradient multiplier
		double TargetVertexDensity{ 1.0 };             // >! Dtarget: target vertex density [verts/unit^2], if negative, compute as mean
		double AvgNeighborhoodDisplacementMultiplier{ 1.0 }; // >! lambda: average 1-ring neighborhood displacement multiplier
		double SelectionSharpness{ 1.0 };               // >! gamma: softmax selection sharpness
		double ConfidenceGM1Imax{ 0.05 };               // >! rho: confidence level of GM{1, Imax}
		double ConfidenceNormal{ 0.8 };                // >! Z: confidence level of N(0,1)
	};

	class SoftmaxUniformVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		using VertexSamplingStrategy::VertexSamplingStrategy;

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	
		GeometricSamplingParams Params;
	};

	class SoftmaxFeatureDetectingVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		using VertexSamplingStrategy::VertexSamplingStrategy;

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;
	};

	// -------------------------------------------------------------
	/// \brief A wrapper for geometric sampling parameters
	// -------------------------------------------------------------
	struct PoissonDiscSamplingParams
	{
		double MinDiscRadius{ 1.0 }; //>! Minimum geometric spacing (in mesh units) between any two sampled vertices.
		double MaxDiscRadius{ 1.5 }; //>! Maximum geometric spacing (in mesh units) between any two sampled vertices.
	
		double ExpectedPtsMultiplier{ 10.0 }; // >! multiplies the extpectedCount to randomly sample available points.
	};

	class PoissonDiscVertexSamplingStrategy : public VertexSamplingStrategy
	{
	public:
		using VertexSamplingStrategy::VertexSamplingStrategy;

		void Sample(const char* start, const char* end,
			std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed, IncrementalProgressTracker& tracker) override;

		PoissonDiscSamplingParams Params;
	};

	inline [[nodiscard]] std::unique_ptr<VertexSamplingStrategy> GetVertexSelectionStrategy(const VertexSelectionType& vertSelType, const unsigned int& completionFrequency, const size_t& maxVertexCount, const std::shared_ptr<IncrementalMeshFileHandler>& handler)
	{
		if (vertSelType == VertexSelectionType::Sequential)
			return std::make_unique<SequentialVertexSamplingStrategy>(completionFrequency, maxVertexCount, handler);
		if (vertSelType == VertexSelectionType::UniformRandom)
			return std::make_unique<UniformRandomVertexSamplingStrategy>(completionFrequency, maxVertexCount, handler);
		if (vertSelType == VertexSelectionType::SoftMaxUniform)
			return std::make_unique<SoftmaxUniformVertexSamplingStrategy>(completionFrequency, maxVertexCount, handler);
		if (vertSelType == VertexSelectionType::SoftMaxFeatureDetecting)
			return std::make_unique<SoftmaxFeatureDetectingVertexSamplingStrategy>(completionFrequency, maxVertexCount, handler);
		if (vertSelType == VertexSelectionType::PoissonDisc)
			return std::make_unique<PoissonDiscVertexSamplingStrategy>(completionFrequency, maxVertexCount, handler);
		return nullptr; // not supported
	}

	inline [[nodiscard]] std::string GetVertexSelectionStrategyName(const VertexSelectionType& vertSelType)
	{
		if (vertSelType == VertexSelectionType::Sequential)
			return "VertexSelectionType::Sequential";
		if (vertSelType == VertexSelectionType::UniformRandom)
			return "VertexSelectionType::UniformRandom";
		if (vertSelType == VertexSelectionType::SoftMaxUniform)
			return "VertexSelectionType::SoftMaxUniform";
		if (vertSelType == VertexSelectionType::SoftMaxFeatureDetecting)
			return "VertexSelectionType::SoftMaxFeatureDetecting";
		if (vertSelType == VertexSelectionType::PoissonDisc)
			return "VertexSelectionType::PoissonDisc";
		return "UNSUPPORTED!";
	}

} // namespace IMB
