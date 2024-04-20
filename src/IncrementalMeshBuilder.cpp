#include "IncrementalMeshBuilder.h"

// #include <nanoflann.hpp> // TODO: use

#include "utils/FileMappingWrapper.h"
#include "utils/IncrementalUtils.h"

#include <thread>

[[nodiscard]] bool InitParamsAreCorrect(const std::string& fileName, const unsigned int& completionFrequency)
{
	if (fileName.empty())
	{
		std::cerr << "InitParamsAreCorrect: fileName.empty()!\n";
		return false;
	}

	if (completionFrequency == 0)
	{
		std::cerr << "InitParamsAreCorrect: completionFrequency == 0!\n";
		return false;
	}

	return true;
}

/// \brief Well, this is merely an approximation. Ascii OBJ files may consist of larger vertex lines.
constexpr size_t APPROX_BYTES_PER_VERTEX = 24;

namespace IMB
{
	//
	// --------------------------------------------------------------------------------------------------------

	void IncrementalMeshBuilder::Init(
		const std::string& fileName, const unsigned int& completionFrequency, 
		const ReconstructionFunctionType& reconstructType, 
		const VertexSelectionType& vertSelType,
		const MeshRenderFunction& renderCallback)
	{
		if (!InitParamsAreCorrect(fileName, completionFrequency))
		{
			std::cerr << "IncrementalMeshBuilder::Init: Invalid parameters! Terminating.\n";
			return;
		}
#if DEBUG_PRINT
		DBG_OUT << "-------------------------------------------------------------------------------------\n";
		DBG_OUT << "IncrementalMeshBuilder::Init: of file: \"" << fileName << "\"\n";
		DBG_OUT << "completionFrequency: " << completionFrequency << " jobs / file load.\n";
		DBG_OUT << "surface reconstruction (meshing) type: " << GetReconstructionStrategyName(reconstructType) << "\n";
		DBG_OUT << "vertex sampling type: " << GetVertexSelectionStrategyName(vertSelType) << "\n";
		DBG_OUT << "-------------------------------------------------------------------------------------\n";
#endif
		m_RenderCallback = renderCallback;
		m_MeshingStrategy = GetReconstructionStrategy(reconstructType);
		m_FileMapping = std::make_unique<Utils::FileMappingWrapper>(fileName);

		// Estimate number of vertices based on [Botsch et al., 2010] ratio N_F = 2 * N_V.
		const auto nExpectedVertices = m_FileMapping->GetLineCount() / 3;

		m_Dispatcher = std::make_unique<IncrementalMeshBuilderDispatcher>(nExpectedVertices, completionFrequency, vertSelType);
		m_IsInitialized = true;
	}

	void IncrementalMeshBuilder::DispatchAndSyncWorkers(const std::optional<unsigned int>& seed, const unsigned int& nThreads)
	{
#if DEBUG_PRINT
		DBG_OUT << "-------------------------------------------------------------------------------------\n";
		DBG_OUT << "IncrementalMeshBuilder::DispatchAndSyncWorkers: \n";
		DBG_OUT << "nThreads : " << nThreads << " " << (nThreads == 0 ? "(will default to 1 worker thread)" : "") << "\n";
		DBG_OUT << "seed : " << (seed.has_value() ? std::to_string(seed.value()) : "std::nullopt") << "\n";
		DBG_OUT << ".....................................................................................\n";
#endif
		if (!m_IsInitialized)
		{
			std::cerr << "IncrementalMeshBuilder::DispatchAndSyncWorkers: This service needs initialization! Use IncrementalMeshBuilder::Init to initialize! Terminating.\n";
			return;
		}
		if (m_IsWorking.exchange(true)) 
		{
			std::cerr << "IncrementalMeshBuilder::DispatchAndSyncWorkers: Processing is already underway.\n";
			return;
		}

		m_IsWorking = true;

		// m_MeshData (Geometry::BaseMeshGeometryData) will be filled
		const auto invokeDispatcherChunkProcess = [this, &seed](const char* start, const char* end)
		{
			try {
				m_Dispatcher->ProcessChunk(start, end, seed);
			}
			catch (const std::exception& e) {
				std::cerr << "IncrementalMeshBuilder::DispatchAndSyncWorkers: Error processing chunk: [" << FormatAddresses(start, end) << "]: " << e.what() << '\n';
				throw;  // Rethrow to propagate the exception
			}
		};
		m_Dispatcher->SetMeshUpdateCallback([this](const std::vector<pmp::Point>& vertices)
			{ UpdateMesh(vertices); });

		try {
			// initiate worker threads
			const size_t nAvailableWorkerThreads = std::thread::hardware_concurrency() - 2;
			const size_t threadCount = nThreads == 0 ? 1 : (nThreads >= nAvailableWorkerThreads ? nAvailableWorkerThreads : nThreads);
			const size_t chunkSize = m_FileMapping->GetFileSize() / threadCount;
#if DEBUG_PRINT
			DBG_OUT << "IncrementalMeshBuilder::DispatchAndSyncWorkers: nAvailableWorkerThreads = " << nAvailableWorkerThreads << "\n";
			DBG_OUT << "IncrementalMeshBuilder::DispatchAndSyncWorkers: threadCount = " << threadCount << "\n";
			DBG_OUT << "IncrementalMeshBuilder::DispatchAndSyncWorkers: m_FileMapping->GetFileSize() = " << m_FileMapping->GetFileSize() << " bytes\n";
			DBG_OUT << "IncrementalMeshBuilder::DispatchAndSyncWorkers: chunkSize = " << chunkSize << " bytes\n";
#endif
			std::vector<std::thread> threads(threadCount);

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
				threads[i] = std::thread(invokeDispatcherChunkProcess, chunkStart, chunkEnd);
			}

			// Wait for all threads to finish
			for (auto& t : threads) {
				t.join();
			}

		}
		catch (...) {
			std::cerr << "IncrementalMeshBuilder::DispatchAndSyncWorkers: A thread encountered a severe error. Terminating all operations.\n";
		}

		m_IsWorking = false;
		// Terminate all owned objects
		Terminate();
	}

	void IncrementalMeshBuilder::UpdateMesh(const std::vector<pmp::Point>& vertices)
	{
#if DEBUG_PRINT
		DBG_OUT << "IncrementalMeshBuilder::UpdateMesh: About to triangulate m_MeshData with " << m_MeshData.Vertices.size() + vertices.size() << " vertices ... \n";
#endif
		std::lock_guard lock(m_MeshDataMutex);
		m_MeshData.Vertices.insert(m_MeshData.Vertices.end(), vertices.begin(), vertices.end());

		if (m_MeshingStrategy)
		{
			// computing m_MeshData.PolyIndices, and for some strategies also modifying m_MeshData.Vertices!
			m_MeshingStrategy->Process(m_MeshData.Vertices, m_MeshData.PolyIndices);
		}
#if DEBUG_PRINT
		DBG_OUT << "IncrementalMeshBuilder::UpdateMesh: done.\n";
		DBG_OUT << "IncrementalMeshBuilder::UpdateMesh: Rendering ... \n";
#endif
		m_RenderCallback(m_MeshData);
#if DEBUG_PRINT
		DBG_OUT << "IncrementalMeshBuilder::UpdateMesh: done.\n";
#endif
	}

	void IncrementalMeshBuilder::Terminate()
	{
		// Invoke destructor of all owned objects
		m_Dispatcher.reset();
		m_MeshingStrategy.reset();
		m_FileMapping.reset();
		m_IsInitialized = false;
	}

} // namespace IMB
