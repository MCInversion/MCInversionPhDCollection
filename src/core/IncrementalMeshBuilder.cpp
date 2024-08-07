#include "IncrementalMeshBuilder.h"

// #include <nanoflann.hpp> // TODO: use

#include "utils/FileMappingWrapper.h"
#include "utils/IncrementalUtils.h"

#include "IncrementalMeshFileHandler.h"

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

namespace IMB
{
	//
	// --------------------------------------------------------------------------------------------------------

	void IncrementalMeshBuilder::Init(
		const std::string& fileName, const unsigned int& completionFrequency, 
		const ReconstructionFunctionType& reconstructType, 
		const VertexSelectionType& vertSelType,
		const MeshRenderFunction& renderCallback,
		const size_t& maxVertexCount)
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
		if (!m_FileMapping->IsValid())
		{
			std::cerr << "IncrementalMeshBuilder::Init: Error during initialization of FileMappingWrapper!\n";
			return;
		}
		const char* fileStart = m_FileMapping->GetFileMemory();
		const char* fileEnd = fileStart + m_FileMapping->GetFileSize();
		m_FileHandler = CreateMeshFileHandler(FileHandlerParams{ fileName, fileStart, fileEnd });
		if (!m_FileHandler)
		{
			std::cerr << "IncrementalMeshBuilder::Init: Error during initialization of IncrementalMeshFileHandler!\n";
			return;
		}
		m_Dispatcher = std::make_unique<IncrementalMeshBuilderDispatcher>(completionFrequency, maxVertexCount, vertSelType, m_FileHandler);
		if (!m_Dispatcher->IsValid())
		{
			std::cerr << "IncrementalMeshBuilder::Init: Error during initialization of IncrementalMeshBuilderDispatcher!\n";
			return;
		}
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
		m_Dispatcher->SetMeshUpdateCallback([this](const std::vector<pmp::Point>& newVertices)
			{ UpdateMesh(newVertices); });

		try {
			// initiate worker threads
			const size_t nAvailableWorkerThreads = std::thread::hardware_concurrency() - 2;
			const size_t threadCount = nThreads == 0 ? 1 : (nThreads >= nAvailableWorkerThreads ? nAvailableWorkerThreads : nThreads);
#if DEBUG_PRINT
			DBG_OUT << "IncrementalMeshBuilder::DispatchAndSyncWorkers: nAvailableWorkerThreads = " << nAvailableWorkerThreads << "\n";
			DBG_OUT << "IncrementalMeshBuilder::DispatchAndSyncWorkers: threadCount = " << threadCount << "\n";
			DBG_OUT << "IncrementalMeshBuilder::DispatchAndSyncWorkers: m_FileMapping->GetFileSize() = " << m_FileMapping->GetFileSize() << " bytes\n";
#endif
			std::vector<std::thread> threads(threadCount);

#if DEBUG_PRINT
			const char* fileStart = m_FileHandler->GetMemoryStart();
			const char* fileEnd = m_FileHandler->GetMemoryEnd();
			const size_t dataSize = fileEnd - fileStart;
			const size_t chunkSize = dataSize / threadCount;
			DBG_OUT << "IncrementalMeshBuilder::DispatchAndSyncWorkers: (vertex)dataSize = " << dataSize << " bytes\n";
			DBG_OUT << "IncrementalMeshBuilder::DispatchAndSyncWorkers: chunkSize = " << chunkSize << " bytes\n";
#endif

			for (size_t i = 0; i < threadCount; ++i)
			{
				const auto [chunkStart, chunkEnd] = m_FileHandler->GetChunkBounds(i, threadCount);
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

	void IncrementalMeshBuilder::UpdateMesh(const std::vector<pmp::Point>& newVertices)
	{
		if (newVertices.empty())
		{
#if DEBUG_PRINT
			DBG_OUT << "IncrementalMeshBuilder::UpdateMesh: Nothing to triangulate. Terminating.\n";
#endif
			return;
		}

		{ // ensure that the lock is held only for the duration of the operations that need synchronization
			std::lock_guard lock(m_MeshDataMutex);
			m_MeshData.Vertices.insert(m_MeshData.Vertices.end(), newVertices.begin(), newVertices.end());
#if DEBUG_PRINT
		DBG_OUT << "IncrementalMeshBuilder::UpdateMesh: About to triangulate m_MeshData with " << m_MeshData.Vertices.size() << " vertices ... \n";
#endif

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
		}
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
