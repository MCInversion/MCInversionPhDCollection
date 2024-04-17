#include "IncrementalMeshBuilder.h"

// #include <nanoflann.hpp> // TODO: use

#include "utils/FileMappingWrapper.h"

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
		// m_MeshData (Geometry::BaseMeshGeometryData) will be filled
		const auto invokeDispatcherChunkProcess = [this](const char* start, const char* end, const std::optional<unsigned int>& seed)
		{ m_Dispatcher->ProcessChunk(start, end, seed);	};
		m_Dispatcher->SetProgressCallback([this](const std::vector<pmp::Point>& vertices, const std::vector<std::vector<unsigned int>>& polyIndices)
			{ UpdateMesh(vertices, polyIndices); });

		// initiate threads
		const size_t nAvailableThreads = std::thread::hardware_concurrency() - 1;
		const size_t threadCount = nThreads == 0 ? 1 : (nThreads >= nAvailableThreads ? nAvailableThreads : nThreads);
		const size_t chunkSize = m_FileMapping->GetFileSize() / threadCount;
		std::vector<std::thread> threads(threadCount);
		std::vector<std::vector<pmp::Point>> threadResults(threadCount);

		char* fileStart = m_FileMapping->GetFileMemory();
		char* fileEnd = fileStart + m_FileMapping->GetFileSize();
		for (size_t i = 0; i < threadCount; ++i)
		{
			char* chunkStart = fileStart + (i * chunkSize);
			char* chunkEnd = (i == threadCount - 1) ? fileEnd : chunkStart + chunkSize;

			// Adjust chunk_end to point to the end of a line
			while (*chunkEnd != '\n' && chunkEnd < fileEnd) {
				chunkEnd++;
			}
			if (chunkEnd != fileEnd) {
				chunkEnd++;  // move past the newline character
			}

			// Start a thread to process this chunk
			threads[i] = std::thread(invokeDispatcherChunkProcess, chunkStart, chunkEnd, std::ref(threadResults[i]));
		}

		// Wait for all threads to finish
		for (auto& t : threads) {
			t.join();
		}

		// Combine results from all threads
		std::vector<pmp::Point> combinedVertices;
		std::vector<std::vector<unsigned int>> combinedPolyIndices;
		for (const auto& result : threadResults) {
			combinedVertices.insert(combinedVertices.end(), result.begin(), result.end());
		}

		// Locking not necessary here if this section is single-threaded after joining
		UpdateMesh(combinedVertices, combinedPolyIndices);  // This now calls the mesh strategy internally
	}

	void IncrementalMeshBuilder::UpdateMesh(const std::vector<pmp::Point>& vertices, const std::vector<std::vector<unsigned int>>& polyIndices)
	{
		std::lock_guard lock(m_MeshDataMutex);
		m_MeshData.Vertices.insert(m_MeshData.Vertices.end(), vertices.begin(), vertices.end());
		m_MeshData.PolyIndices.insert(m_MeshData.PolyIndices.end(), polyIndices.begin(), polyIndices.end());

		// Now invoke meshing strategy
		m_Dispatcher->ProcessMeshUpdate();

		m_RenderCallback(m_MeshData);
	}
} // namespace IMB
