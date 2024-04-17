#include "IncrementalMeshBuilder.h"

// #include <nanoflann.hpp> // TODO: use

#include "utils/FileMappingWrapper.h"

#include <iostream>

constexpr unsigned int DEFAULT_MAX_VERTEX_CAPACITY = 1'000'000;

// --------------------------------------------------------------------------------------------------------

[[nodiscard]] std::unique_ptr<IMB::PointCloudMeshingStrategy> GetReconstructionStrategy(const IMB::ReconstructionFunctionType& reconstructType)
{
	if (reconstructType == IMB::ReconstructionFunctionType::BallPivoting)
		return std::make_unique<IMB::BallPivotingMeshingStrategy>();
	if (reconstructType == IMB::ReconstructionFunctionType::Poisson)
		return std::make_unique<IMB::PoissonMeshingStrategy>();
	if (reconstructType == IMB::ReconstructionFunctionType::MarchingCubes)
		return std::make_unique<IMB::MarchingCubesMeshingStrategy>();
	return std::make_unique<IMB::LagrangianShrinkWrappingMeshingStrategy>();
}

[[nodiscard]] std::unique_ptr<IMB::VertexSamplingStrategy> GetVertexSelectionStrategy(const IMB::VertexSelectionType& vertSelType)
{
	if (vertSelType == IMB::VertexSelectionType::Sequential)
		return std::make_unique<IMB::SequentialVertexSamplingStrategy>();
	if (vertSelType == IMB::VertexSelectionType::UniformRandom)
		return std::make_unique<IMB::UniformRandomVertexSamplingStrategy>();
	if (vertSelType == IMB::VertexSelectionType::NormalRandom)
		return std::make_unique<IMB::NormalRandomVertexSamplingStrategy>();
	return std::make_unique<IMB::SoftmaxFeatureDetectingVertexSamplingStrategy>();
}

namespace IMB
{
	//
	// --------------------------------------------------------------------------------------------------------

	void IncrementalMeshBuilder::Init(
		const std::string& fileName, const unsigned int& completionFrequency, 
		const ReconstructionFunctionType& reconstructType, 
		const VertexSelectionType& vertSelType)
	{
		auto ptCloudMeshingStrategy = std::move(GetReconstructionStrategy(reconstructType));
		auto vertexSamplingStrategy = std::move(GetVertexSelectionStrategy(vertSelType));

		m_FileMapping = std::make_unique<Utils::FileMappingWrapper>(fileName);
		const auto fileSize = m_FileMapping->GetFileSize();
		const auto nExpectedVertices = fileSize / 3;

		m_Dispatcher = std::make_unique<IncrementalMeshBuilderDispatcher>(
			nExpectedVertices, completionFrequency, std::move(vertexSamplingStrategy), std::move(ptCloudMeshingStrategy));
	}

	void IncrementalMeshBuilder::ProcessVertices(const std::optional<unsigned int>& seed, const unsigned int& nThreads)
	{
		
	}
} // namespace IMB
