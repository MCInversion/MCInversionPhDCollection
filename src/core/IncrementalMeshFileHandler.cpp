#include "IncrementalMeshFileHandler.h"
#include "IncrementalMeshFileHandler.h"

#include "IncrementalMeshFileHandler.h"

#include "utils/IncrementalUtils.h"
#include "utils/StringUtils.h"

#include "IncrementalProgressUtils.h"

#define CHECK_LARGE_COORDS false

#if CHECK_LARGE_COORDS
constexpr pmp::Scalar MAX_ABS_COORD = 1e+5;
#endif

namespace
{
	//// Function to check if the system is little-endian
	//bool IsSystemLittleEndian()
	//{
	//	uint16_t num = 0x1;
	//	char* numPtr = reinterpret_cast<char*>(&num);
	//	return (numPtr[0] == 1);
	//}

	//// Function to convert float from little-endian to system endianness if needed
	//float ConvertFloatFromLittleEndian(float input)
	//{
	//	//if (!IsSystemLittleEndian())
	//	//{
	//		char* floatPtr = reinterpret_cast<char*>(&input);
	//		std::reverse(floatPtr, floatPtr + sizeof(float));  // Reverse the bytes of the float
	//	//}
	//	return input;
	//}

	//
	// ——————————————————————————————————————————————————————————————————
	// Step 0: Helper routines for endianness and type‐sizes
	// ——————————————————————————————————————————————————————————————————
	//

	// A. Determine at runtime if we need to swap 32‐bit values.
	//    We assume the host is little‐endian. If the file is BIG_ENDIAN, we swap.
	//    (If you ever need to detect a big‐endian host, you can add a check via std::endian.)
	static bool NeedSwap(PlyFormat fileFormat)
	{
		// If file is big‐endian and host is little, we must swap.
		return (fileFormat == PlyFormat::BINARY_BIG_ENDIAN);
	}

	// B. Swap a uint32_t (e.g. to go from big‐endian to little)
	static uint32_t SwapUInt32(uint32_t x)
	{
		return  ((x & 0xFF000000u) >> 24)
			| ((x & 0x00FF0000u) >> 8)
			| ((x & 0x0000FF00u) << 8)
			| ((x & 0x000000FFu) << 24);
	}

	// C. Map PlyDataType → size in bytes
	static size_t SizeOfDataType(PlyDataType dt)
	{
		switch (dt)
		{
		case PlyDataType::INT8:    return 1;
		case PlyDataType::UINT8:   return 1;
		case PlyDataType::INT16:   return 2;
		case PlyDataType::UINT16:  return 2;
		case PlyDataType::INT32:   return 4;
		case PlyDataType::UINT32:  return 4;
		case PlyDataType::FLOAT32: return 4;
		case PlyDataType::FLOAT64: return 8;
		default:
			throw std::runtime_error("SizeOfDataType: Unknown or unsupported PlyDataType");
		}
	}
}

namespace IMB
{
	constexpr size_t APPROX_BYTES_PER_OBJ_VERTEX = 24;

	void IncrementalASCIIOBJFileHandler::Sample(const char* start, const char* end, const std::vector<size_t>& indices, size_t updateThreshold, std::vector<pmp::Point>& result, IncrementalProgressTracker& tracker) const
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

	pmp::Point IncrementalASCIIOBJFileHandler::SampleSinglePoint(const char* start, const char* end, const size_t& idx) const
	{
#if DEBUG_PRINT
		DBG_OUT << "IncrementalASCIIOBJFileHandler::SampleSinglePoint: ... \n";
#endif
		// Calculate an approximate position to jump to
		const char* cursor = start + idx * APPROX_BYTES_PER_OBJ_VERTEX;

		// Adjust cursor to the start of the next line if not already at a newline
		while (cursor < end && *cursor != '\n') cursor++;
		if (cursor < end) cursor++;  // Move past the newline to the start of the next line

		// Ensure the cursor is within valid range after adjustment
		if (cursor >= end)
		{
			throw std::runtime_error("IncrementalASCIIOBJFileHandler::SampleSinglePoint: cursor >= end!\n");
		}

		// Check if the line starts with "v " indicating a vertex
		if (strncmp(cursor, "v ", 2) != 0)
		{
			throw std::runtime_error("IncrementalASCIIOBJFileHandler::SampleSinglePoint: no \"v\" element!\n");
		}

		cursor += 2; // skip "v "

		pmp::vec3 vec;
		char* tempCursor;
		vec[0] = std::strtof(cursor, &tempCursor);
		cursor = tempCursor;

		vec[1] = std::strtof(cursor, &tempCursor);
		cursor = tempCursor;

		vec[2] = std::strtof(cursor, &tempCursor);
		cursor = tempCursor;

		return vec;
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

	void IncrementalBinaryPLYFileHandler::Sample(const char* start, const char* end, const std::vector<size_t>& indices, size_t updateThreshold, std::vector<pmp::Point>& result, IncrementalProgressTracker& tracker) const
	{
#if DEBUG_PRINT
		DBG_OUT << "IncrementalBinaryPLYFileHandler::Sample: Starting sample process...\n";
#endif

		// Find the "vertex" element in m_Elements.
		const PlyElement* vertexElem = nullptr;
		for (auto const& elem : m_Elements)
		{
			if (elem.name == "vertex")
			{
				vertexElem = &elem;
				break;
			}
		}
		if (!vertexElem)
		{
			throw std::runtime_error("Sample: no \"vertex\" element in header.");
		}

		// Compute byte‐offsets of x, y, z within each vertex record, and compute total record size.
		size_t byteCursor = 0;
		size_t xOffset = 0, yOffset = 0, zOffset = 0;
		bool   foundX = false, foundY = false, foundZ = false;

		for (auto const& prop : vertexElem->properties)
		{
			if (!prop.isList)
			{
				// scalar property: size is known
				size_t sz = SizeOfDataType(prop.dataType);

				if (prop.name == "x")
				{
					xOffset = byteCursor;
					foundX = true;
				}
				else if (prop.name == "y")
				{
					yOffset = byteCursor;
					foundY = true;
				}
				else if (prop.name == "z")
				{
					zOffset = byteCursor;
					foundZ = true;
				}

				byteCursor += sz;
			}
			else
			{
				// list‐typed property: at header time we only know the size of the count field.
				// If any list appears before x/y/z, we cannot compute a fixed offset for x, etc.
				if (prop.name == "x" || prop.name == "y" || prop.name == "z")
				{
					throw std::runtime_error("Sample: found list-typed property named x, y, or z");
				}
				// add size of count field (e.g. uchar or uint32)
				size_t cntSz = SizeOfDataType(prop.listCountType);
				byteCursor += cntSz;
				// The actual data block for that list is variable per‐record; we assume
				// no list appears before x,y,z, and that lists come after all coordinate scalars.
			}
		}

		if (!foundX || !foundY || !foundZ)
		{
			throw std::runtime_error("Sample: could not locate all of x, y, z properties.");
		}

		// 'vertexRecordSize' is the total number of bytes occupied by all scalar props (and
		// count fields for lists) in one vertex record.  Any variable‐length data (list entries)
		// must be read at runtime if they precede coordinates.  Here we assume lists, if any,
		// come after x,y,z.
		size_t vertexRecordSize = byteCursor;

		// Determine if we need to swap bytes (file is big‐endian and host is little‐endian).
		bool swapBytes = (m_Format == PlyFormat::BINARY_BIG_ENDIAN);

		// Reserve space in result.
		result.reserve(result.size() + updateThreshold);

		size_t localVertexCount = 0;

		// Loop over each requested index, read x,y,z, and accumulate into 'result'.
		for (size_t idx : indices)
		{
			const char* vertexData = start + idx * vertexRecordSize;
			if (vertexData + vertexRecordSize > end)
			{
				// Out of bounds: skip
				continue;
			}

			// Copy raw 4 bytes for x, y, z
			uint32_t rawX = 0, rawY = 0, rawZ = 0;
			std::memcpy(&rawX, vertexData + xOffset, sizeof(uint32_t));
			std::memcpy(&rawY, vertexData + yOffset, sizeof(uint32_t));
			std::memcpy(&rawZ, vertexData + zOffset, sizeof(uint32_t));

			if (swapBytes)
			{
				rawX = SwapUInt32(rawX);
				rawY = SwapUInt32(rawY);
				rawZ = SwapUInt32(rawZ);
			}

			float fx = *reinterpret_cast<float*>(&rawX);
			float fy = *reinterpret_cast<float*>(&rawY);
			float fz = *reinterpret_cast<float*>(&rawZ);

#if CHECK_LARGE_COORDS
			if (std::isnan(fx) || std::isnan(fy) || std::isnan(fz) ||
				std::abs(fx) > MAX_ABS_COORD ||
				std::abs(fy) > MAX_ABS_COORD ||
				std::abs(fz) > MAX_ABS_COORD)
			{
				continue;
			}
#endif

			result.emplace_back(fx, fy, fz);
			localVertexCount++;

			if (localVertexCount >= updateThreshold)
			{
#if DEBUG_PRINT
				DBG_OUT << "IncrementalBinaryPLYFileHandler::Sample: Updating tracker with " << localVertexCount << " vertices.\n";
#endif
				tracker.Update(localVertexCount);
				localVertexCount = 0;
			}
		}

		// Final tracker update (if any vertices remain).
		if (localVertexCount > 0)
		{
#if DEBUG_PRINT
			DBG_OUT << "IncrementalBinaryPLYFileHandler::Sample: Final tracker update with " << localVertexCount << " vertices.\n";
#endif
			tracker.Update(localVertexCount, /*finalChunk=*/true);
		}

#if DEBUG_PRINT
		DBG_OUT << "IncrementalBinaryPLYFileHandler::Sample: Finished sampling " << updateThreshold << " vertices.\n";
#endif
	}

	pmp::Point IncrementalBinaryPLYFileHandler::SampleSinglePoint(const char* start, const char* end, const size_t& idx) const
	{
		// 1) Find the "vertex" element in m_Elements
		const PlyElement* vertexElem = nullptr;
		for (auto const& elem : m_Elements)
		{
			if (elem.name == "vertex")
			{
				vertexElem = &elem;
				break;
			}
		}
		if (!vertexElem)
		{
			throw std::runtime_error("ReadVertexPositionByIndex: no \"vertex\" element in header");
		}

		// 2) Compute byte‐offsets of x,y,z inside one vertex record
		size_t byteCursor = 0;
		size_t xOffset = 0, yOffset = 0, zOffset = 0;
		bool   foundX = false, foundY = false, foundZ = false;

		for (auto const& prop : vertexElem->properties)
		{
			if (!prop.isList)
			{
				size_t sz = SizeOfDataType(prop.dataType);

				if (prop.name == "x")
				{
					xOffset = byteCursor;
					foundX = true;
				}
				else if (prop.name == "y")
				{
					yOffset = byteCursor;
					foundY = true;
				}
				else if (prop.name == "z")
				{
					zOffset = byteCursor;
					foundZ = true;
				}

				byteCursor += sz;
			}
			else
			{
				// If a list precedes x/y/z, we can’t compute a fixed offset
				if (prop.name == "x" || prop.name == "y" || prop.name == "z")
				{
					throw std::runtime_error("ReadVertexPositionByIndex: found list-typed property named x, y, or z");
				}
				byteCursor += SizeOfDataType(prop.listCountType);
				// We assume lists come after x,y,z
			}
		}

		if (!foundX || !foundY || !foundZ)
		{
			throw std::runtime_error(
				"ReadVertexPositionByIndex: could not locate all of x,y,z properties");
		}

		size_t vertexRecordSize = byteCursor;

		// 3) Check bounds: idx must be < m_GlobalVertexCountEstimate
		if (idx >= vertexElem->count)
		{
			throw std::runtime_error("ReadVertexPositionByIndex: idx out of range");
		}

		// 4) Compute the pointer to this vertex’s data
		const char* vertexData = start + idx * vertexRecordSize;
		if (vertexData + vertexRecordSize > end)
		{
			throw std::runtime_error("ReadVertexPositionByIndex: vertex data out of bounds");
		}

		// 5) Copy raw bytes for x,y,z, handle endianness
		bool swapBytes = (m_Format == PlyFormat::BINARY_BIG_ENDIAN);

		uint32_t rawX = 0u, rawY = 0u, rawZ = 0u;
		std::memcpy(&rawX, vertexData + xOffset, sizeof(uint32_t));
		std::memcpy(&rawY, vertexData + yOffset, sizeof(uint32_t));
		std::memcpy(&rawZ, vertexData + zOffset, sizeof(uint32_t));

		if (swapBytes)
		{
			rawX = SwapUInt32(rawX);
			rawY = SwapUInt32(rawY);
			rawZ = SwapUInt32(rawZ);
		}

		float fx = *reinterpret_cast<float*>(&rawX);
		float fy = *reinterpret_cast<float*>(&rawY);
		float fz = *reinterpret_cast<float*>(&rawZ);

		return pmp::Point{ fx, fy, fz };
	}

	void IncrementalBinaryPLYFileHandler::EstimateGlobalVertexCount(const char* start, const char* end)
	{
		ParseHeader(start);
		if (m_GlobalVertexCountEstimate == 0) 
		{
			std::cerr << "No vertex count found in the PLY header." << std::endl;
		}
	}

	std::pair<const char*, const char*> IncrementalBinaryPLYFileHandler::GetChunkBounds(size_t chunkIndex, size_t totalChunks) const
	{
		const size_t totalSize = m_VertexDataEnd - m_VertexDataStart;
		size_t chunkSize = totalSize / totalChunks;
		chunkSize -= chunkSize % m_VertexItemSize; // Adjust chunk size based on the actual vertex size

		const char* start = m_VertexDataStart + chunkIndex * chunkSize;
		const char* end = (chunkIndex + 1 == totalChunks) ? m_VertexDataEnd : start + chunkSize;

		return { start, end };
	}

	void IncrementalBinaryPLYFileHandler::InitializeMemoryBounds(const char* fileStart, const char* fileEnd)
	{
		// Find "end_header" and move to the start of the vertex data
		const std::string headerEndTag = "end_header";
		const char* cursor = fileStart;
		while (cursor + headerEndTag.length() < fileEnd) 
		{
			if (strncmp(cursor, headerEndTag.c_str(), headerEndTag.length()) == 0) 
			{
				cursor += headerEndTag.length();
				while (cursor < fileEnd && (*cursor == '\n' || *cursor == '\r'))
				{
					cursor++; // Skip the newline character(s)
				}
				break;
			}
			cursor++;
		}

		m_VertexDataStart = cursor;
		m_VertexDataEnd = m_VertexDataStart + (m_GlobalVertexCountEstimate * m_VertexItemSize);
	}

	size_t IncrementalBinaryPLYFileHandler::GetLocalVertexCountEstimate(const char* start, const char* end) const
	{
		const auto dataSize = static_cast<size_t>(end - start);
		constexpr size_t vertexSize = sizeof(float) * 3; // Assuming each vertex consists of three floats
		const size_t localVertexCount = dataSize / vertexSize; // Calculate how many complete vertices are within this block
		return localVertexCount;
	}

	void IncrementalBinaryPLYFileHandler::ParseHeader(const char* fileStart)
	{
		if (!fileStart)
		{
			throw std::runtime_error("ParseHeader: null fileStart pointer");
		}

		const char* cursor = fileStart;
		const char* lineStart = cursor;

		// 1) Read the very first line: must be "ply"
		{
			const char* eol = cursor;
			while (*eol && *eol != '\r' && *eol != '\n')
			{
				++eol;
			}
			std::string firstLine(cursor, eol);
			if (firstLine != "ply")
			{
				throw std::runtime_error("ParseHeader: missing 'ply' as first token");
			}

			// skip over the EOL markers
			cursor = eol;
			if (*cursor == '\r') ++cursor;
			if (*cursor == '\n') ++cursor;
		}

		// 2) Next lines: keep reading until we see "end_header"
		PlyElement* currentElement = nullptr;

		while (true)
		{
			if (*cursor == '\0')
			{
				throw std::runtime_error("ParseHeader: reached EOF before end_header");
			}

			// 2.a) Extract one line
			lineStart = cursor;
			const char* eol = cursor;
			while (*eol && *eol != '\r' && *eol != '\n')
			{
				++eol;
			}
			std::string line(lineStart, eol);

			// Move cursor past the EOL
			cursor = eol;
			if (*cursor == '\r') ++cursor;
			if (*cursor == '\n') ++cursor;

			// 2.b) Tokenize by whitespace
			std::istringstream ss(line);
			std::string token;
			if (!(ss >> token))
			{
				// empty line? skip it
				continue;
			}

			if (token == "comment" || token == "obj_info")
			{
				std::string rest;
				std::getline(ss, rest);
				if (!rest.empty() && (rest.front() == ' ' || rest.front() == '\t'))
				{
					rest.erase(rest.begin());
				}
				m_Comments.push_back(rest);
				continue;
			}

			if (token == "format")
			{
				std::string fmtstr;
				if (!(ss >> fmtstr))
				{
					throw std::runtime_error("ParseHeader: 'format' line missing format token");
				}
				if (fmtstr == "ascii")
				{
					m_Format = PlyFormat::ASCII;
					throw std::runtime_error("ParseHeader: ASCII not supported!");
				}
				else if (fmtstr == "binary_little_endian")
				{
					m_Format = PlyFormat::BINARY_LITTLE_ENDIAN;
				}
				else if (fmtstr == "binary_big_endian")
				{
					m_Format = PlyFormat::BINARY_BIG_ENDIAN;
				}
				else
				{
					throw std::runtime_error("ParseHeader: unrecognized PLY format '" + fmtstr + "'");
				}

				double v = 0.0;
				if (!(ss >> v))
				{
					throw std::runtime_error("ParseHeader: 'format' line missing version number");
				}
				m_Version = v;
				continue;
			}

			if (token == "element")
			{
				std::string elemName;
				size_t elemCount;
				if (!(ss >> elemName))
				{
					throw std::runtime_error("ParseHeader: 'element' line missing element name");
				}
				if (!(ss >> elemCount))
				{
					throw std::runtime_error("ParseHeader: 'element' line missing element count");
				}

				m_Elements.push_back(PlyElement{ elemName, elemCount, {} });
				currentElement = &m_Elements.back();

				if (elemName == "vertex")
				{
					m_GlobalVertexCountEstimate = elemCount;
				}
				continue;
			}

			if (token == "property")
			{
				if (!currentElement)
				{
					throw std::runtime_error("ParseHeader: 'property' outside any element");
				}

				std::string maybeType;
				if (!(ss >> maybeType))
				{
					throw std::runtime_error("ParseHeader: 'property' missing type token");
				}

				PlyProperty prop;
				if (maybeType == "list")
				{
					// format: property list <count_type> <data_type> <property_name>
					std::string countTypeStr, dataTypeStr, propName;
					if (!(ss >> countTypeStr >> dataTypeStr >> propName))
					{
						throw std::runtime_error("ParseHeader: malformed 'property list' line");
					}
					prop.isList = true;
					prop.listCountType = TypeFromString(countTypeStr);
					prop.listDataType = TypeFromString(dataTypeStr);
					prop.name = propName;

					if (!prop.IsValid())
					{
						throw std::runtime_error("ParseHeader: invalid PLY data types in 'property list'");
					}
				}
				else
				{
					// scalar property: <type> <name>
					std::string propName;
					if (!(ss >> propName))
					{
						throw std::runtime_error("ParseHeader: malformed 'property' line");
					}
					prop.isList = false;
					prop.dataType = TypeFromString(maybeType);
					prop.name = propName;

					if (!prop.IsValid())
					{
						throw std::runtime_error("ParseHeader: invalid PLY data type '" + maybeType + "'");
					}
				}

				currentElement->properties.push_back(std::move(prop));
				continue;
			}

			if (token == "end_header")
			{
				m_HeaderSizeBytes = static_cast<size_t>(cursor - fileStart);
				break;
			}

			throw std::runtime_error("ParseHeader: unrecognized header token '" + token + "'");
		}

		// 3) After parsing header, estimate m_VertexItemSize from the 'vertex' element
		size_t vertexSize = 0;
		for (auto const& elem : m_Elements)
		{
			if (elem.name == "vertex")
			{
				for (auto const& prop : elem.properties)
				{
					if (!prop.isList)
					{
						// scalar property
						vertexSize += SizeOfDataType(prop.dataType);
					}
					else
					{
						// list: only count field size is fixed
						vertexSize += SizeOfDataType(prop.listCountType);
					}
				}
				break;
			}
		}
		// If no vertex element was found, vertexSize remains 0; that is acceptable fallback
		m_VertexItemSize = vertexSize;
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