#pragma once

#include "pmp/Types.h"

#include <string>
#include <vector>
#include <memory>

namespace IMB
{
	// forward declarations
	class IncrementalProgressTracker;

	/// \brief A wrapper for file handler input params.
	struct FileHandlerParams
	{
		std::string FilePath{};
		const char* FileStart{ nullptr };
		const char* FileEnd{ nullptr };
	};
	
	class IncrementalMeshFileHandler
	{
	public:
		virtual void Sample(const char* start, const char* end, const std::vector<size_t>& indices, size_t updateThreshold, std::vector<pmp::Point>& result, IncrementalProgressTracker& tracker) = 0;

		[[nodiscard]] size_t GetGlobalVertexCountEstimate() const { return m_GlobalVertexCountEstimate; }

		virtual [[nodiscard]] size_t GetLocalVertexCountEstimate(const char* start, const char* end) const = 0;

		virtual [[nodiscard]] std::pair<const char*, const char*> GetChunkBounds(size_t chunkIndex, size_t totalChunks) const = 0;

		virtual ~IncrementalMeshFileHandler() = default;

		virtual void EstimateGlobalVertexCount(const char* start, const char* end) = 0;

		virtual void InitializeMemoryBounds(const char* fileStart, const char* fileEnd) = 0;

		[[nodiscard]] const char* GetMemoryStart() const
		{
			return m_VertexDataStart;
		}

		[[nodiscard]] const char* GetMemoryEnd() const
		{
			return m_VertexDataEnd;
		}

	protected:
		const char* m_VertexDataStart{ nullptr };
		const char* m_VertexDataEnd{ nullptr };
		size_t m_GlobalVertexCountEstimate{0};
	};

	class IncrementalASCIIOBJFileHandler : public IncrementalMeshFileHandler
	{
	public:
		using IncrementalMeshFileHandler::IncrementalMeshFileHandler; 

		void Sample(const char* start, const char* end, const std::vector<size_t>& indices, size_t updateThreshold, std::vector<pmp::Point>& result, IncrementalProgressTracker& tracker) override;

		void EstimateGlobalVertexCount(const char* start, const char* end) override;

		[[nodiscard]] std::pair<const char*, const char*> GetChunkBounds(size_t chunkIndex, size_t totalChunks) const override;

		void InitializeMemoryBounds(const char* fileStart, const char* fileEnd) override;

		[[nodiscard]] size_t GetLocalVertexCountEstimate(const char* start, const char* end) const override;
	};

	class IncrementalBinaryPLYFileHandler : public IncrementalMeshFileHandler
	{
	public:
		using IncrementalMeshFileHandler::IncrementalMeshFileHandler;

		void Sample(const char* start, const char* end, const std::vector<size_t>& indices, size_t updateThreshold, std::vector<pmp::Point>& result, IncrementalProgressTracker& tracker) override;

	    void EstimateGlobalVertexCount(const char* start, const char* end) override;

		[[nodiscard]] std::pair<const char*, const char*> GetChunkBounds(size_t chunkIndex, size_t totalChunks) const override;

		void InitializeMemoryBounds(const char* fileStart, const char* fileEnd) override;

		[[nodiscard]] size_t GetLocalVertexCountEstimate(const char* start, const char* end) const override;
	private:

		void ParseHeader(const char* fileStart);
		
		size_t m_VertexItemSize{ 3 }; //>! vertex data could contain additional information besides coordinates, like color.
	};

	/**
	 * \brief Creates a mesh file handler based on the file's extension.
	 *
	 * \param handlerParams input params for file handler.
	 *
	 * \return A unique_ptr to the created file handler, or nullptr on failure (invalid input or unsupported format).
	 */
	[[nodiscard]] std::shared_ptr<IncrementalMeshFileHandler> CreateMeshFileHandler(const FileHandlerParams& handlerParams);
	
} // namespace IMB
