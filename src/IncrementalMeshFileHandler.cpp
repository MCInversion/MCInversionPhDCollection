
#include "IncrementalMeshFileHandler.h"

#include "utils/IncrementalUtils.h"
#include "utils/StringUtils.h"
#include "IncrementalProgressUtils.h"

namespace
{
	// Function to check if the system is little-endian
	bool IsSystemLittleEndian()
	{
		uint16_t num = 0x1;
		char* numPtr = reinterpret_cast<char*>(&num);
		return (numPtr[0] == 1);
	}

	// Function to convert float from little-endian to system endianness if needed
	float ConvertFloatFromLittleEndian(float input)
	{
		//if (!IsSystemLittleEndian())
		//{
			char* floatPtr = reinterpret_cast<char*>(&input);
			std::reverse(floatPtr, floatPtr + sizeof(float));  // Reverse the bytes of the float
		//}
		return input;
	}
}

namespace IMB
{
	constexpr size_t APPROX_BYTES_PER_OBJ_VERTEX = 24;

	void IncrementalASCIIOBJFileHandler::Sample(const char* start, const char* end, const std::vector<size_t>& indices, size_t updateThreshold, std::vector<pmp::Point>& result, IncrementalProgressTracker& tracker)
	{
#if DEBUG_PRINT
		DBG_OUT << "IncrementalASCIIOBJFileHandler::Sample: ... \n";
#endif
		size_t localVertexCount = 0;

		for (const auto index : indices)
		{
			// Calculate an approximate position to jump to
			const char* cursor = start + index * APPROX_BYTES_PER_OBJ_VERTEX;

			// Adjust cursor to the start of the next line if not already at a newline
			while (cursor < end && *cursor != '\n') cursor++;
			if (cursor < end) cursor++;  // Move past the newline to the start of the next line

			// Ensure the cursor is within valid range after adjustment
			if (cursor >= end) continue;


			// Check if the line starts with "v " indicating a vertex
			if (strncmp(cursor, "v ", 2) == 0)
			{
				cursor += 2; // skip "v "

				pmp::vec3 vec;
				char* tempCursor;
				vec[0] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				vec[1] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				vec[2] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				result.push_back(vec);
				localVertexCount++;

				// Check if it's time to update the tracker
				if (localVertexCount >= updateThreshold)
				{
#if DEBUG_PRINT
					DBG_OUT << "IncrementalASCIIOBJFileHandler::Sample: Time to update the tracker with " << localVertexCount << " collected vertices.\n";
#endif
					tracker.Update(localVertexCount);
					localVertexCount = 0; // Reset local count after update
				}
			}

			// Otherwise, skip to the next line
			while (*cursor != '\n' && cursor < end) cursor++;
		}

#if DEBUG_PRINT
		DBG_OUT << "IncrementalASCIIOBJFileHandler::Sample: ... done.\n";
#endif
		// Ensure any remaining vertices are accounted for
		if (localVertexCount > 0)
		{
#if DEBUG_PRINT
			DBG_OUT << "IncrementalASCIIOBJFileHandler::Sample: Time for a final tracker update with " << localVertexCount << " collected vertices.\n";
#endif
			tracker.Update(localVertexCount, true);
		}
	}

	void IncrementalASCIIOBJFileHandler::EstimateGlobalVertexCount(const char* start, const char* end)
	{
		// Rough estimation by counting the number of lines and assuming every third line might be a vertex
		const size_t nLines = std::count(start, end, '\n');
		m_GlobalVertexCountEstimate = nLines / 3;  // Estimate number of vertices based on [Botsch et al., 2010] ratio N_F = 2 * N_V.
	}

	std::pair<const char*, const char*> IncrementalASCIIOBJFileHandler::GetChunkBounds(size_t chunkIndex, size_t totalChunks) const
	{
		const size_t totalSize = m_VertexDataEnd - m_VertexDataStart;
		const size_t chunkSize = totalSize / totalChunks;

		const char* start = m_VertexDataStart + chunkIndex * chunkSize;
		const char* end = (chunkIndex + 1 == totalChunks) ? m_VertexDataEnd : start + chunkSize;

		// Adjust to nearest line ending
		while (end < m_VertexDataEnd && *end != '\n') ++end;
		if (end != m_VertexDataEnd) ++end; // include the newline character

		return { start, end };
	}

	void IncrementalASCIIOBJFileHandler::InitializeMemoryBounds(const char* fileStart, const char* fileEnd)
	{
		// ASCII data processing might start directly
		m_VertexDataStart = fileStart;
		m_VertexDataEnd = fileEnd;
	}

	size_t IncrementalASCIIOBJFileHandler::GetLocalVertexCountEstimate(const char* start, const char* end) const
	{
		const size_t nLines = std::count(start, end, '\n');
		return nLines / 3; 
	}

	//
	// ===================================================================================================
	//

	void IncrementalBinaryPLYFileHandler::Sample(const char* start, const char* end, const std::vector<size_t>& indices, size_t updateThreshold, std::vector<pmp::Point>& result, IncrementalProgressTracker& tracker)
	{
#if DEBUG_PRINT
		DBG_OUT << "IncrementalBinaryPLYFileHandler::Sample: Starting sample process...\n";
		//DBG_OUT << "IncrementalBinaryPLYFileHandler::Sample: System endianness: " << (IsSystemLittleEndian() ? "little endian" : "big endian") << "\n";
#endif
		size_t localVertexCount = 0;
		constexpr size_t vertexSize = sizeof(float) * 3; // Assuming each vertex is represented by three floats
		result.reserve(m_GlobalVertexCountEstimate);

		for (const auto index : indices)
		{
			const char* vertexData = start + index * vertexSize;
			if (vertexData + vertexSize > end) continue; // Ensure the vertex data is within bounds

			// Interpret the bytes as floats directly
			const auto* coords = reinterpret_cast<const float*>(vertexData);

			pmp::Point vec(coords[0], coords[1], coords[2]);
			result.push_back(vec);
			localVertexCount++;

			// Check if it's time to update the tracker
			if (localVertexCount >= updateThreshold)
			{
#if DEBUG_PRINT
				DBG_OUT << "IncrementalBinaryPLYFileHandler::Sample: Updating tracker with " << localVertexCount << " vertices.\n";
#endif
				tracker.Update(localVertexCount);
				localVertexCount = 0; // Reset local count after update
			}
		}

#if DEBUG_PRINT
		DBG_OUT << "IncrementalBinaryPLYFileHandler::Sample: Finished sampling " << result.size() << " vertices.\n";
#endif
		// Ensure any remaining vertices are accounted for
		if (localVertexCount > 0)
		{
#if DEBUG_PRINT
			DBG_OUT << "IncrementalBinaryPLYFileHandler::Sample: Final tracker update with " << localVertexCount << " vertices.\n";
#endif
			tracker.Update(localVertexCount, true);
		}
	}

	void IncrementalBinaryPLYFileHandler::EstimateGlobalVertexCount(const char* start, const char* end)
	{
		const char* cursor = start;
		std::string line;
		m_GlobalVertexCountEstimate = 0;  // Default if no vertex count is found

		while (cursor < end) {
			line.clear();
			// Read until newline or end of the header section
			while (cursor < end && *cursor != '\n') {
				line.append(1, *cursor++);
			}
			cursor++;  // Move past the newline character

			// Check if we've reached the end of the header
			if (line == "end_header") {
				break;
			}

			// Check for the vertex element line
			if (line.find("element vertex ") != std::string::npos) {
				std::istringstream iss(line);
				std::string element, vertex;
				size_t count;
				iss >> element >> vertex >> count;
				if (vertex == "vertex") {
					m_GlobalVertexCountEstimate = count;
					break;  // Exit once the vertex count is found
				}
			}
		}

		if (m_GlobalVertexCountEstimate == 0) {
			std::cerr << "No vertex count found in the PLY header." << std::endl;
		}
	}

	std::pair<const char*, const char*> IncrementalBinaryPLYFileHandler::GetChunkBounds(size_t chunkIndex, size_t totalChunks) const
	{
		const size_t totalSize = m_VertexDataEnd - m_VertexDataStart;
		size_t chunkSize = (totalSize / totalChunks);
		chunkSize -= chunkSize % (sizeof(float) * 3); // Ensure chunkSize is multiple of vertex size

		const char* start = m_VertexDataStart + chunkIndex * chunkSize;
		const char* end = (chunkIndex + 1 == totalChunks) ? m_VertexDataEnd : start + chunkSize;

		return { start, end };
	}

	void IncrementalBinaryPLYFileHandler::InitializeMemoryBounds(const char* fileStart, const char* fileEnd)
	{
		// Find "end_header" and move to the start of the vertex data
		const std::string headerEndTag = "end_header";
		const char* cursor = fileStart;
		while (cursor + headerEndTag.length() < fileEnd) {
			if (strncmp(cursor, headerEndTag.c_str(), headerEndTag.length()) == 0) {
				cursor += headerEndTag.length();
				while (cursor < fileEnd && (*cursor == '\n' || *cursor == '\r')) {
					cursor++; // Skip the newline character(s)
				}
				break;
			}
			cursor++;
		}

		m_VertexDataStart = cursor;
		constexpr size_t vertexSize = sizeof(float) * 3; // Adjust based on actual per-vertex data size and types
		m_VertexDataEnd = m_VertexDataStart + (m_GlobalVertexCountEstimate * vertexSize);
	}

	size_t IncrementalBinaryPLYFileHandler::GetLocalVertexCountEstimate(const char* start, const char* end) const
	{
		const auto dataSize = static_cast<size_t>(end - start);
		constexpr size_t vertexSize = sizeof(float) * 3; // Assuming each vertex consists of three floats
		const size_t localVertexCount = dataSize / vertexSize; // Calculate how many complete vertices are within this block
		return localVertexCount;
	}

	//
	// ===================================================================================================
	//

	std::shared_ptr<IncrementalMeshFileHandler> CreateMeshFileHandler(const FileHandlerParams& handlerParams)
	{
		if (handlerParams.FilePath.empty())
		{
			std::cerr << "IMB::CreateMeshFileHandler: filePath.empty()!\n";
			return nullptr;
		}

		if (!handlerParams.FileStart || !handlerParams.FileEnd || handlerParams.FileStart >= handlerParams.FileEnd)
		{
			std::cerr << "IMB::CreateMeshFileHandler: Invalid file pointers provided!\n";
		}

		const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(handlerParams.FilePath);
		std::shared_ptr<IncrementalMeshFileHandler> handler;

		// Determine the correct handler based on the file extension
		if (extension == "obj")
		{
			handler = std::make_shared<IncrementalASCIIOBJFileHandler>();
		}
		else if (extension == "ply")
		{
			handler = std::make_shared<IncrementalBinaryPLYFileHandler>();
		}
		else
		{
			std::cerr << "IMB::CreateMeshFileHandler: Unsupported file format: " << extension << "!\n";
			return nullptr;
		}
		handler->EstimateGlobalVertexCount(handlerParams.FileStart, handlerParams.FileEnd);
		handler->InitializeMemoryBounds(handlerParams.FileStart, handlerParams.FileEnd);
		return handler;
	}
} // namespace IMB