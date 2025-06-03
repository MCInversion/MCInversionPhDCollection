#pragma once

#include "pmp/Types.h"

#include <string>
#include <vector>
#include <memory>

// ================================================================
//                        PLY Types

// -----------------------------------------------------------
enum class PlyFormat
{
	ASCII,
	BINARY_LITTLE_ENDIAN,
	BINARY_BIG_ENDIAN,
	UNKNOWN
};

enum class PlyDataType : int
{
	INVALID = -1,
	INT8,    // char / int8_t
	UINT8,   // uchar / uint8_t
	INT16,   // short / int16_t
	UINT16,  // ushort / uint16_t
	INT32,   // int / int32_t
	UINT32,  // uint / uint32_t
	FLOAT32, // float
	FLOAT64  // double
};

// Helper: map from string -> PlyDataType
static PlyDataType TypeFromString(const std::string& token)
{
	// based on common PLY type names
	static const std::unordered_map<std::string, PlyDataType> s_map = {
		{ "char",    PlyDataType::INT8    },
		{ "int8",    PlyDataType::INT8    },
		{ "uchar",   PlyDataType::UINT8   },
		{ "uint8",   PlyDataType::UINT8   },
		{ "short",   PlyDataType::INT16   },
		{ "int16",   PlyDataType::INT16   },
		{ "ushort",  PlyDataType::UINT16  },
		{ "uint16",  PlyDataType::UINT16  },
		{ "int",     PlyDataType::INT32   },
		{ "int32",   PlyDataType::INT32   },
		{ "uint",    PlyDataType::UINT32  },
		{ "uint32",  PlyDataType::UINT32  },
		{ "float",   PlyDataType::FLOAT32 },
		{ "float32", PlyDataType::FLOAT32 },
		{ "double",  PlyDataType::FLOAT64 },
		{ "float64", PlyDataType::FLOAT64 }
	};

	auto it = s_map.find(token);
	return (it == s_map.end()) ? PlyDataType::INVALID : it->second;
}

// One property (either scalar or list)
struct PlyProperty
{
	std::string       name;
	bool              isList = false;
	PlyDataType       dataType = PlyDataType::INVALID; // for scalars
	PlyDataType       listCountType = PlyDataType::INVALID; // only used if isList == true
	PlyDataType       listDataType = PlyDataType::INVALID; // only used if isList == true

	// convenience:
	bool IsValid() const
	{
		if (!isList)
			return dataType != PlyDataType::INVALID;
		else
			return (listCountType != PlyDataType::INVALID) && (listDataType != PlyDataType::INVALID);
	}
};

// One element, with name, count, and a vector of properties
struct PlyElement
{
	std::string             name;
	size_t                  count = 0;
	std::vector<PlyProperty> properties;
};

/// ===============================================================

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
		virtual void Sample(const char* start, const char* end, const std::vector<size_t>& indices, size_t updateThreshold, std::vector<pmp::Point>& result, IncrementalProgressTracker& tracker) const = 0;

		virtual [[nodiscard]] pmp::Point SampleSinglePoint(const char* start, const char* end, const size_t& idx) const = 0;

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

		void Sample(const char* start, const char* end, const std::vector<size_t>& indices, size_t updateThreshold, std::vector<pmp::Point>& result, IncrementalProgressTracker& tracker) const override;

		[[nodiscard]] pmp::Point SampleSinglePoint(const char* start, const char* end, const size_t& idx) const override;

		void EstimateGlobalVertexCount(const char* start, const char* end) override;

		[[nodiscard]] std::pair<const char*, const char*> GetChunkBounds(size_t chunkIndex, size_t totalChunks) const override;

		void InitializeMemoryBounds(const char* fileStart, const char* fileEnd) override;

		[[nodiscard]] size_t GetLocalVertexCountEstimate(const char* start, const char* end) const override;
	};

	class IncrementalBinaryPLYFileHandler : public IncrementalMeshFileHandler
	{
	public:
		using IncrementalMeshFileHandler::IncrementalMeshFileHandler;

		void Sample(const char* start, const char* end, const std::vector<size_t>& indices, size_t updateThreshold, std::vector<pmp::Point>& result, IncrementalProgressTracker& tracker) const override;

		[[nodiscard]] pmp::Point SampleSinglePoint(const char* start, const char* end, const size_t& idx) const override;

	    void EstimateGlobalVertexCount(const char* start, const char* end) override;

		[[nodiscard]] std::pair<const char*, const char*> GetChunkBounds(size_t chunkIndex, size_t totalChunks) const override;

		void InitializeMemoryBounds(const char* fileStart, const char* fileEnd) override;

		[[nodiscard]] size_t GetLocalVertexCountEstimate(const char* start, const char* end) const override;
	private:

		void ParseHeader(const char* fileStart);

		PlyFormat                m_Format = PlyFormat::UNKNOWN;
		double                   m_Version = 0.0;
		std::vector<std::string> m_Comments;       // all "comment" or "obj_info" lines
		std::vector<PlyElement>  m_Elements;       // each element and its properties
		size_t                   m_HeaderSizeBytes = 0; // how many bytes from fileStart to end_header + EOL
		
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
