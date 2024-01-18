#include "GeometryConversionUtils.h"

#include "utils/StringUtils.h"

#include <set>
#include <fstream>
#include <thread>

#ifdef _WINDOWS
// Windows-specific headers
#include <windows.h>
#include <fcntl.h>
#else
// Unsupported platform
#error "Unsupported platform"
#endif

// TODO: different representations, e.g.: ProgressiveMeshData which would then be exported into a disk file

namespace
{
	struct ChunkData
	{
		std::vector<pmp::vec3> Vertices;
		std::vector<pmp::vec3> VertexNormals;
		std::vector<std::vector<unsigned int>> PolyIndices;
	};

	void ParseChunk(const char* start, const char* end, ChunkData& data)
	{
		const char* cursor = start;

		while (cursor < end)
		{
			// If the current line is incomplete, skip to the next line
			if (*cursor == '\n')
			{
				cursor++;
				continue;
			}

			// If it's a vertex or normal, parse the three floats
			if (strncmp(cursor, "v ", 2) == 0 || strncmp(cursor, "vn ", 3) == 0)
			{
				pmp::vec3 vec;
				cursor += (cursor[1] == ' ') ? 2 : 3; // skip "v " or "vn "

				char* tempCursor;
				vec[0] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				vec[1] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				vec[2] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				if (cursor[-1] == 'n') // distinguish between vertex and normal
					data.VertexNormals.push_back(vec);
				else
					data.Vertices.push_back(vec);
			}
			// If it's a face, parse the vertex indices
			else if (strncmp(cursor, "f ", 2) == 0)
			{
				cursor += 2; // skip "f "
				std::vector<unsigned int> faceIndices;

				while (*cursor != '\n' && cursor < end)
				{
					char* tempCursor;
					const unsigned int index = std::strtoul(cursor, &tempCursor, 10);
					cursor = tempCursor;
					faceIndices.push_back(index - 1);

					while (*cursor == ' ' || *cursor == '/') // skip to next index or newline
						cursor++;
				}

				data.PolyIndices.push_back(faceIndices);
			}
			else
			{
				// Skip to the next line if the current line isn't recognized
				while (*cursor != '\n' && cursor < end)
					cursor++;
			}
		}
	}
} // anonymous namespace

namespace Geometry
{
	pmp::SurfaceMesh ConvertBufferGeomToPMPSurfaceMesh(const BaseMeshGeometryData& geomData)
	{
		pmp::SurfaceMesh result;

		// count edges
		std::set<std::pair<unsigned int, unsigned int>> edgeIdsSet;
		for (const auto& indexTuple : geomData.PolyIndices)
		{
			for (unsigned int i = 0; i < indexTuple.size(); i++)
			{
				unsigned int vertId0 = indexTuple[i];
				unsigned int vertId1 = indexTuple[(static_cast<size_t>(i) + 1) % indexTuple.size()];

				if (vertId0 > vertId1) std::swap(vertId0, vertId1);

				edgeIdsSet.insert({ vertId0, vertId1 });
			}
		}

		result.reserve(geomData.Vertices.size(), edgeIdsSet.size(), geomData.PolyIndices.size());
		for (const auto& v : geomData.Vertices)
			result.add_vertex(pmp::Point(v[0], v[1], v[2]));

		if (!geomData.VertexNormals.empty())
		{
			auto vNormal = result.vertex_property<pmp::Normal>("v:normal");
			for (auto v : result.vertices())
				vNormal[v] = geomData.VertexNormals[v.idx()];
		}

        for (const auto& indexTuple : geomData.PolyIndices)
        {
			std::vector<pmp::Vertex> vertices;
			vertices.reserve(indexTuple.size());
			for (const auto& vId : indexTuple)
				vertices.emplace_back(pmp::Vertex(vId));

			result.add_face(vertices);	        
        }

		return result;
	}

	pmp::SurfaceMesh ConvertMCMeshToPMPSurfaceMesh(const MC_Mesh& mcMesh)
	{
		pmp::SurfaceMesh result;

		// count edges
		std::set<std::pair<unsigned int, unsigned int>> edgeIdsSet;
		for (unsigned int i = 0; i < mcMesh.faceCount * 3; i += 3)
		{
			for (unsigned int j = 0; j < 3; j++)
			{
				unsigned int vertId0 = mcMesh.faces[i + j];
				unsigned int vertId1 = mcMesh.faces[i + (j + 1) % 3];

				if (vertId0 > vertId1) std::swap(vertId0, vertId1);

				edgeIdsSet.insert({ vertId0, vertId1 });
			}
		}

		result.reserve(mcMesh.vertexCount, edgeIdsSet.size(), mcMesh.faceCount);
		// MC produces normals by default
		auto vNormal = result.vertex_property<pmp::Normal>("v:normal");

		for (unsigned int i = 0; i < mcMesh.vertexCount; i++)
		{
			result.add_vertex(
			   pmp::Point(
				   mcMesh.vertices[i][0],
				   mcMesh.vertices[i][1],
				   mcMesh.vertices[i][2]
			));
			vNormal[pmp::Vertex(i)] = pmp::Normal{
				mcMesh.normals[i][0],
				mcMesh.normals[i][1],
				mcMesh.normals[i][2]
			};
		}

		for (unsigned int i = 0; i < mcMesh.faceCount * 3; i += 3)
		{
			std::vector<pmp::Vertex> vertices;
			vertices.reserve(3);
			for (unsigned int j = 0; j < 3; j++)
			{
				vertices.emplace_back(pmp::Vertex(mcMesh.faces[i + j]));
			}

			result.add_face(vertices);
		}

		return result;
	}

	bool ExportBaseMeshGeometryDataToOBJ(const BaseMeshGeometryData& geomData, const std::string& absFileName)
	{
		std::ofstream file(absFileName);
		if (!file.is_open())
		{
			std::cerr << "Failed to open file for writing: " << absFileName << std::endl;
			return false;
		}

		// Write vertices
		for (const auto& vertex : geomData.Vertices)
		{
			file << "v " << vertex[0] << ' ' << vertex[1] << ' ' << vertex[2] << '\n';
		}

		// Optionally, write vertex normals
		if (!geomData.VertexNormals.empty())
		{
			for (const auto& normal : geomData.VertexNormals)
			{
				file << "vn " << normal[0] << ' ' << normal[1] << ' ' << normal[2] << '\n';
			}
		}

		// Write faces
		for (const auto& indices : geomData.PolyIndices)
		{
			file << "f";
			for (unsigned int index : indices)
			{
				// OBJ indices start from 1, not 0
				file << ' ' << (index + 1);
			}
			file << '\n';
		}

		file.close();
		return true;
	}

	std::optional<BaseMeshGeometryData> ImportOBJMeshGeometryData(const std::string& absFileName, const bool& importInParallel)
	{
		const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(absFileName);
		if (extension != "obj")
			return {};

		const char* file_path = absFileName.c_str();

		const HANDLE file_handle = CreateFile(file_path, GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);

		if (file_handle == INVALID_HANDLE_VALUE) 
		{
			std::cerr << "ImportOBJMeshGeometryData [ERROR]: Failed to open the file.\n";
			return {};
		}

		// Get the file size
		const DWORD file_size = GetFileSize(file_handle, nullptr);

		// Create a file mapping object
		const HANDLE file_mapping = CreateFileMapping(file_handle, nullptr, PAGE_READONLY, 0, 0, nullptr);
		if (file_mapping == nullptr) 
		{
			std::cerr << "ImportOBJMeshGeometryData [ERROR]: Failed to create file mapping.\n";
			CloseHandle(file_handle);
			return {};
		}

		// Map the file into memory
		const LPVOID file_memory = MapViewOfFile(file_mapping, FILE_MAP_READ, 0, 0, 0);
		if (file_memory == nullptr) 
		{
			std::cerr << "ImportOBJMeshGeometryData [ERROR]: Failed to map the file.\n";
			CloseHandle(file_mapping);
			CloseHandle(file_handle);
			return {};
		}

		BaseMeshGeometryData resultData;

		// Determine the number of threads
		const size_t thread_count = importInParallel ? std::thread::hardware_concurrency() : 1;
		const size_t chunk_size = file_size / thread_count;
		std::vector<std::thread> threads(thread_count);
		std::vector<ChunkData> threadResults(thread_count);

		char* file_start = static_cast<char*>(file_memory);
		char* file_end = file_start + file_size;

		for (size_t i = 0; i < thread_count; ++i)
		{
			char* chunk_start = file_start + (i * chunk_size);
			char* chunk_end = (i == thread_count - 1) ? file_end : chunk_start + chunk_size;

			// Adjust chunk_end to point to the end of a line
			while (*chunk_end != '\n' && chunk_end < file_end) {
				chunk_end++;
			}
			chunk_end++;  // move past the newline character

			// Start a thread to process this chunk
			threads[i] = std::thread(ParseChunk, chunk_start, chunk_end, std::ref(threadResults[i]));
		}

		// Wait for all threads to finish
		for (auto& t : threads) {
			t.join();
		}

		for (const auto& result : threadResults)
		{
			resultData.Vertices.insert(resultData.Vertices.end(), result.Vertices.begin(), result.Vertices.end());
			resultData.PolyIndices.insert(resultData.PolyIndices.end(), result.PolyIndices.begin(), result.PolyIndices.end());
			resultData.VertexNormals.insert(resultData.VertexNormals.end(), result.VertexNormals.begin(), result.VertexNormals.end());
		}

		// Clean up
		UnmapViewOfFile(file_memory);
		CloseHandle(file_mapping);
		CloseHandle(file_handle);

		return std::move(resultData);
	}

} // namespace Geometry