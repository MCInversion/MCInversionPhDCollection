#include "IncrementalMeshBuilder.h"

// #include <nanoflann.hpp> // TODO: use

#include "utils/FileMappingWrapper.h"

#include <thread>

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
		m_RenderCallback = renderCallback;
		m_MeshingStrategy = GetReconstructionStrategy(reconstructType);
		m_FileMapping = std::make_unique<Utils::FileMappingWrapper>(fileName);
		const auto fileSize = m_FileMapping->GetFileSize();
		const auto nExpectedVertices = fileSize / 3;
		m_Dispatcher = std::make_unique<IncrementalMeshBuilderDispatcher>(nExpectedVertices, completionFrequency, vertSelType);
	}

	void IncrementalMeshBuilder::DispatchAndSyncWorkers(const std::optional<unsigned int>& seed, const unsigned int& nThreads)
	{
		if (m_IsWorking)
			return;

		m_IsWorking = true;

		// m_MeshData (Geometry::BaseMeshGeometryData) will be filled
		const auto invokeDispatcherChunkProcess = [this, &seed](const char* start, const char* end)
		{
			try {
				m_Dispatcher->ProcessChunk(start, end, seed);
			}
			catch (const std::exception& e) {
				std::cerr << "Error processing chunk: [" << std::stoi(start) << ", " << std::stoi(end) << "]: " << e.what() << '\n';
				throw;  // Rethrow to propagate the exception
			}
		};
		m_Dispatcher->SetMeshUpdateCallback([this](const std::vector<pmp::Point>& vertices)
			{ UpdateMesh(vertices); });

		try {
			// initiate threads
			const size_t nAvailableThreads = std::thread::hardware_concurrency() - 1;
			const size_t threadCount = nThreads == 0 ? 1 : (nThreads >= nAvailableThreads ? nAvailableThreads : nThreads);
			const size_t chunkSize = m_FileMapping->GetFileSize() / threadCount;
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

			//// Combine results from all threads
			//std::vector<pmp::Point> combinedVertices;
			//std::vector<std::vector<unsigned int>> combinedPolyIndices;
			//for (const auto& result : threadResults) {
			//	combinedVertices.insert(combinedVertices.end(), result.begin(), result.end());
			//}

			//// Locking not necessary here if this section is single-threaded after joining
			//UpdateMesh(combinedVertices);  // This now calls the mesh strategy internally

		}
		catch (...) {
			std::cerr << "A thread encountered a severe error. Terminating all operations.\n";
		}

		m_IsWorking = false;
		// Terminate all owned objects
		Terminate();
	}

	void IncrementalMeshBuilder::UpdateMesh(const std::vector<pmp::Point>& vertices)
	{
		std::lock_guard lock(m_MeshDataMutex);
		m_MeshData.Vertices.insert(m_MeshData.Vertices.end(), vertices.begin(), vertices.end());

		if (m_MeshingStrategy) 
		{
			// computing m_MeshData.PolyIndices, and for some strategies also modifying m_MeshData.Vertices!
			m_MeshingStrategy->Process(m_MeshData.Vertices, m_MeshData.PolyIndices);
		}
		if (m_RenderCallback) 
		{
			m_RenderCallback(m_MeshData);
		}
	}

	void IncrementalMeshBuilder::Terminate()
	{
		// Invoke destructor of all owned objects
		m_Dispatcher.reset();
		m_MeshingStrategy.reset();
		m_FileMapping.reset();
	}

} // namespace IMB
