#include "GeometryConversionUtils.h"

#include "utils/StringUtils.h"

#include <set>
#include <fstream>
#include <random>
#include <thread>
#include <unordered_set>
#include <cassert>
#include <memory>
#include <numeric>

#include "pmp/algorithms/Normals.h"

#include <tetgen.h>
#include <triangle.h>

#include "VCGAdapter.h"
#include "PoissonAdapter.h"

#include "sawhney_mat/MedialAxisTransform.h"
#include "sawhney_mat/BoundaryElement.h"
#include "sawhney_mat/BoundaryGenerator.h"
#include "sawhney_mat/Path.h"

#ifdef _WINDOWS
// Windows-specific headers
#include <windows.h>
#include <fcntl.h>
#else
// Unsupported platform
#error "Unsupported platform"
#endif


namespace
{
	void FillVCGMeshWithPoints(VCG_Mesh& vcgMesh, 
		const std::vector<pmp::Point>& points, const std::optional<std::vector<pmp::Normal>>& normals = std::nullopt)
	{
		if (normals && normals->size() != points.size())
		{
			std::cerr << "FillVCGMeshWithPoints: normals && normals->size() != points.size()!\n";
			return;
		}

		vcgMesh.Clear();  // Clear existing mesh data
		vcgMesh.vert.reserve(points.size());

		// Add vertices to the VCG_Mesh
		for (const auto& point : points) {
			VCG_Vertex v;
			v.P() = vcg::Point3f(point[0], point[1], point[2]);
			vcgMesh.vert.push_back(v);
		}
		vcgMesh.vn = vcgMesh.vert.size();  // Update the vertex count
		vcg::tri::UpdateBounding<VCG_Mesh>::Box(vcgMesh);

		if (normals.has_value())
		{
			// Add vertex normals to each vertex
			unsigned int i = 0;
			for (VCG_Mesh::VertexIterator vi = vcgMesh.vert.begin(); vi != vcgMesh.vert.end(); ++vi)
			{
				vi->N() = vcg::Point3f((*normals)[i][0], (*normals)[i][1], (*normals)[i][2]);
				++i;
			}
		}
	}

	[[nodiscard]] std::vector<std::vector<unsigned int>> ExtractVertexIndicesFromVCGMesh(const VCG_Mesh& vcgMesh)
	{
		std::vector<std::vector<unsigned int>> indices;
		indices.reserve(vcgMesh.face.size());
		for (const auto& f : vcgMesh.face) 
		{
			std::vector<unsigned int> faceIndices;
			faceIndices.reserve(f.VN()); // Assuming VN() is the number of vertices per face, typically 3 for a triangular mesh

			for (int i = 0; i < f.VN(); ++i) 
			{
				if (f.cV(i)) 
				{ // Check if the vertex pointer is not null
					faceIndices.push_back(vcg::tri::Index(vcgMesh, f.cV(i))); // Get the index of the vertex in the mesh
				}
			}
			indices.push_back(faceIndices);
		}
		return indices;
	}

	[[nodiscard]] std::vector<pmp::Normal> ExtractVertexNormalsFromVCGMesh(const VCG_Mesh& vcgMesh, const pmp::Scalar& orient = 1.0)
	{
		std::vector<pmp::Normal> normals;
		std::ranges::transform(vcgMesh.vert, std::back_inserter(normals), [&orient](const auto& vcgVert) {
			return pmp::Point{ orient * vcgVert.cN()[0], orient *vcgVert.cN()[1], orient * vcgVert.cN()[2] }; });
		return normals;
	}

	// TODO: different representations, e.g.: ProgressiveMeshData which would then be exported into a disk file
	/**
	 * \brief a thread-specific wrapper for mesh data.
	 */
	struct ChunkData
	{
		std::vector<pmp::vec3> Vertices;
		std::vector<pmp::vec3> VertexNormals;
		std::vector<std::vector<unsigned int>> PolyIndices;
	};

	/**
	 * \brief Parses a chunk of OBJ data. This function is run for each thread.
	 * \param start     start address for this thread.
	 * \param end       end address for this thread.
	 * \param data      preallocated chunk data.
	 */
	void ParseChunk(const char* start, const char* end, ChunkData& data)
	{
		const char* cursor = start;

		while (cursor < end)
		{
			// If the current line is incomplete, skip to the next line
			if (*cursor == '\n' || *cursor == '\r')
			{
				cursor++;
				continue;
			}

			// If it's a vertex or normal, parse the three floats
			if (strncmp(cursor, "v ", 2) == 0 || strncmp(cursor, "vn ", 3) == 0)
			{
				const bool isNormal = strncmp(cursor, "vn ", 3) == 0;
				cursor += isNormal ? 3 : 2; // skip "v " or "vn "

				pmp::vec3 vec;
				char* tempCursor;
				vec[0] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				vec[1] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				vec[2] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				if (isNormal)
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
					// Parse the vertex index
					const unsigned int vertexIndex = std::strtoul(cursor, &tempCursor, 10);
					if (cursor == tempCursor) // No progress in parsing, break to avoid infinite loop
						break;

					cursor = tempCursor;
					if (vertexIndex == 0) {
						// Underflow occurred, meaning strtoul failed to parse a valid number
						break;
					}
					faceIndices.push_back(vertexIndex - 1);

					if (*cursor == '/') // Check for additional indices
					{
						cursor++; // Skip the first '/'
						if (*cursor != '/') // Texture coordinate index is present
						{
							std::strtoul(cursor, &tempCursor, 10); // Parse and discard the texture index
							cursor = tempCursor;
						}

						if (*cursor == '/') // Normal index is present
						{
							cursor++; // Skip the second '/'
							std::strtoul(cursor, &tempCursor, 10); // Parse and discard the normal index
							cursor = tempCursor;
						}
					}

					// Skip to next index, newline, or end
					while (*cursor != ' ' && *cursor != '\n' && cursor < end)
						cursor++;

					if (*cursor == ' ')
						cursor++; // Skip space to start of next index
				}

				if (!faceIndices.empty())
					data.PolyIndices.push_back(faceIndices);

				// Move to the next line
				while (*cursor != '\n' && cursor < end)
					cursor++;
				if (cursor < end)
					cursor++; // Move past the newline character
			}
			else
			{
				// Skip to the next line if the current line isn't recognized
				while (*cursor != '\n' && cursor < end)
					cursor++;
			}
		}
	}

	// Assume average bytes per vertex and face line (these are just example values)
	//constexpr size_t avgBytesPerVertex = 13; // Example average size of a vertex line in bytes
	//constexpr size_t avgBytesPerFace = 7; // Example average size of a face line in bytes

	/**
	 * \brief A simple estimation utility for large OBJ files
	 * \param fileSize        file size in bytes.
	 * \param expectNormals   if true we expect the OBJ file to also contain vertex normals.
	 * \return pair { # of estimated vertices, # of estimated faces }.
	 */
	//std::pair<size_t, size_t> EstimateVertexAndFaceCapacitiesFromOBJFileSize(const size_t& fileSize, const bool& expectNormals = false)
	//{

	//	// Estimate total number of vertex and face lines
	//	const size_t estimatedTotalLines = fileSize / (((expectNormals ? 2 : 1) * avgBytesPerVertex + avgBytesPerFace) / 2);

	//	// Estimate number of vertices based on Botsch 2010 ratio N_F = 2 * N_V
	//	size_t estimatedVertices = estimatedTotalLines / (expectNormals ? 4 : 3); // Since N_F + N_V = 3 * N_V and N_V = N_N
	//	size_t estimatedFaces = 2 * estimatedVertices;

	//	return std::make_pair(estimatedVertices, estimatedFaces);
	//}

	/**
	 * \brief Parses a chunk of PLY point data. This function is run for each thread.
	 * \param start     start address for this thread.
	 * \param end       end address for this thread.
	 * \param data      preallocated chunk data.
	 */
	void ParsePointCloudChunk(const char* start, const char* end, std::vector<pmp::vec3>& data) {
		const char* cursor = start;

		while (cursor < end) 
		{
			// Skip any leading whitespace
			while ((*cursor == ' ' || *cursor == '\n') && cursor < end)
			{
				cursor++;
			}

			if (cursor >= end)
			{
				break; // Reached the end of the chunk
			}

			// Start of a new line
			const char* lineStart = cursor;

			// Find the end of the line
			while (*cursor != '\n' && cursor < end) 
			{
				cursor++;
			}

			// Check if the end of the line is within the chunk
			if (*cursor != '\n' && cursor == end) 
			{
				break; // The line is incomplete, so stop parsing
			}

			// Extract the line
			std::string line(lineStart, cursor);

			// Remove carriage return if present
			if (!line.empty() && line.back() == '\r') 
			{
				line.pop_back();
			}

			// Increment cursor to start the next line in the next iteration
			if (cursor < end) 
			{
				cursor++;
			}

			// Parse the line to extract vertex coordinates
			std::istringstream iss(line);
			pmp::vec3 vertex;
			if (!(iss >> vertex[0] >> vertex[1] >> vertex[2])) 
			{
				std::cerr << "ParsePointCloudChunk: Error parsing line: " << line << std::endl;
				continue;
			}

			data.push_back(vertex);
		}
	}

	/**
	 * \brief Reads the header of the PLY file for vertices.
	 * \param start     header start position in memory.
	 * \return pair { number of vertices read from the header, the memory position where the point data starts }
	 */
	[[nodiscard]] std::pair<size_t, char*> ReadPLYVertexHeader(const char* start)
	{
		const char* cursor = start;
		std::string line;
		size_t vertexCount = 0;

		while (*cursor != '\0') 
		{
			// Extract the line
			const char* lineStart = cursor;
			while (*cursor != '\n' && *cursor != '\0') 
			{
				cursor++;
			}
			line.assign(lineStart, cursor);

			// Trim trailing carriage return if present
			if (!line.empty() && line.back() == '\r') 
			{
				line.pop_back();
			}

			// Move to the start of the next line
			if (*cursor != '\0') 
			{
				cursor++;
			}

			// Check for the vertex count line
			if (line.rfind("element vertex", 0) == 0) 
			{
				std::istringstream iss(line);
				std::string element, vertex;
				if (!(iss >> element >> vertex >> vertexCount))
				{
					std::cerr << "ReadPLYVertexHeader [WARNING]: Failed to parse vertex count line: " << line << "\n";
					continue;
				}
				// Successfully parsed the vertex count line. Continue to look for the end of the header.
			}

			// Check for the end of the header
			if (line == "end_header") 
			{
				break;
			}
		}

		if (vertexCount == 0) 
		{
			std::cerr << "ReadPLYVertexHeader [ERROR]: Vertex count not found in the header.\n";
		}

		return { vertexCount, const_cast<char*>(cursor) }; // Return the vertex count and the cursor position
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

	BaseMeshGeometryData ConvertPMPSurfaceMeshToBaseMeshGeometryData(const pmp::SurfaceMesh& pmpMesh)
	{
		BaseMeshGeometryData geomData;

		// Extract vertices
		auto points = pmpMesh.get_vertex_property<pmp::Point>("v:point");
		for (const auto v : pmpMesh.vertices())
			geomData.Vertices.emplace_back(points[v][0], points[v][1], points[v][2]);

		// Extract vertex normals if available
		if (pmpMesh.has_vertex_property("v:normal")) 
		{
			auto normals = pmpMesh.get_vertex_property<pmp::Normal>("v:normal");
			for (const auto v : pmpMesh.vertices())
				geomData.VertexNormals.emplace_back(normals[v][0], normals[v][1], normals[v][2]);
		}

		// Extract face indices
		for (const auto f : pmpMesh.faces())
		{
			std::vector<unsigned int> faceIndices;
			for (auto v : pmpMesh.vertices(f))
			{
				faceIndices.push_back(v.idx());
			}
			geomData.PolyIndices.push_back(faceIndices);
		}

		return geomData;
	}

	pmp::ManifoldCurve2D ConvertBufferGeomToPMPManifoldCurve2D(const BaseCurveGeometryData& geomData)
	{
		pmp::ManifoldCurve2D result;

		for (const auto& v : geomData.Vertices)
			result.add_vertex(v);

		result.reserve(geomData.Vertices.size(), geomData.EdgeIndices.size());
		if (!geomData.VertexNormals.empty())
		{
			auto vNormal = result.vertex_property<pmp::Normal2>("v:normal");
			for (auto v : result.vertices())
				vNormal[v] = geomData.VertexNormals[v.idx()];
		}

		for (const auto& [vertId0, vertId1] : geomData.EdgeIndices)
		{
			result.add_edge(pmp::Vertex(vertId0), pmp::Vertex(vertId1));
		}

		return result;
	}

	BaseCurveGeometryData ConvertPMPManifoldCurve2DToBaseCurveGeometryData(const pmp::ManifoldCurve2D& pmpCurve)
	{
		BaseCurveGeometryData result;
		result.Vertices.reserve(pmpCurve.n_vertices());
		result.EdgeIndices.reserve(pmpCurve.n_vertices());
		for (const auto& vPos : pmpCurve.positions())
		{
			result.Vertices.push_back(vPos);
		}

		if (pmpCurve.has_vertex_property("v:normal"))
		{
			auto normals = pmpCurve.get_vertex_property<pmp::Normal2>("v:normal");
			for (const auto v : pmpCurve.vertices())
				result.VertexNormals.emplace_back(normals[v][0], normals[v][1]);
		}
		
		for (const auto e : pmpCurve.edges())
		{
			const auto [v0, v1] = pmpCurve.vertices(e);
			result.EdgeIndices.emplace_back(
				static_cast<unsigned int>(v0.idx()), 
				static_cast<unsigned int>(v1.idx())
			);
		}

		return result;
	}

	pmp::SurfaceMesh ConvertMCMeshToPMPSurfaceMesh(const IlatsikMC::MC_Mesh& mcMesh)
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

	bool ExportBaseMeshGeometryDataToVTK(const BaseMeshGeometryData& geomData, const std::string& absFileName)
	{
		const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(absFileName);
		if (extension != "vtk")
		{
			std::cerr << absFileName << "ExportBaseMeshGeometryDataToVTK: Invalid file extension!" << std::endl;
			return false;
		}

		std::ofstream file(absFileName);
		if (!file.is_open())
		{
			std::cerr << "ExportBaseMeshGeometryDataToVTK: Failed to open file for writing: " << absFileName << std::endl;
			return false;
		}

		// Header
		file << "# vtk DataFile Version 3.0\n";
		file << "VTK output from mesh data\n";
		file << "ASCII\n";
		file << "DATASET POLYDATA\n";

		// Write vertices
		file << "POINTS " << geomData.Vertices.size() << " float\n";
		for (const auto& vertex : geomData.Vertices)
		{
			file << vertex[0] << ' ' << vertex[1] << ' ' << vertex[2] << '\n';
		}

		// Write polygons (faces)
		size_t numIndices = 0;
		for (const auto& indices : geomData.PolyIndices)
		{
			numIndices += indices.size() + 1;  // +1 because we also write the size of the polygon
		}
		file << "POLYGONS " << geomData.PolyIndices.size() << ' ' << numIndices << '\n';
		for (const auto& indices : geomData.PolyIndices)
		{
			file << indices.size();  // Number of vertices in this polygon
			for (unsigned int index : indices)
			{
				file << ' ' << index;  // VTK indices start from 0
			}
			file << '\n';
		}

		// Optionally, write vertex normals
		if (!geomData.VertexNormals.empty())
		{
			file << "POINT_DATA " << geomData.VertexNormals.size() << '\n';
			file << "NORMALS normals float\n";
			for (const auto& normal : geomData.VertexNormals)
			{
				file << normal[0] << ' ' << normal[1] << ' ' << normal[2] << '\n';
			}
		}

		file.close();
		return true;
	}

	bool ExportBaseTetraMeshGeometryDataToVTK(const BaseTetraMeshGeometryData& geomData, const std::string& absFileName)
	{
		const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(absFileName);
		if (extension != "vtk")
		{
			std::cerr << absFileName << "ExportBaseTetraMeshGeometryDataToVTK: Invalid file extension!" << std::endl;
			return false;
		}

		std::ofstream file(absFileName);
		if (!file.is_open()) 
		{
			std::cerr << "ExportBaseTetraMeshGeometryDataToVTK: Failed to open file for writing: " << absFileName << std::endl;
			return false;
		}

		// Write VTK header for an unstructured grid.
		file << "# vtk DataFile Version 3.0\n";
		file << "VTK output from tet mesh data\n";
		file << "ASCII\n";
		file << "DATASET UNSTRUCTURED_GRID\n";

		// Write points.
		file << "POINTS " << geomData.Vertices.size() << " float\n";
		for (const auto& vertex : geomData.Vertices)
		{
			file << vertex[0] << " " << vertex[1] << " " << vertex[2] << "\n";
		}

		// Write cells.
		// Each tetrahedron is written as: 4 <pt0> <pt1> <pt2> <pt3>
		size_t numTets = geomData.TetrahedraIndices.size();
		size_t totalCellEntries = numTets * 5;
		file << "CELLS " << numTets << " " << totalCellEntries << "\n";
		for (const auto& tet : geomData.TetrahedraIndices)
		{
			file << "4 " << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << "\n";
		}

		// Write cell types: for tetrahedra, VTK cell type is 10.
		file << "CELL_TYPES " << numTets << "\n";
		for (size_t i = 0; i < numTets; ++i)
		{
			file << "10\n";
		}

		// Optionally, add a CELL_DATA section so that volume rendering shows something.
		// Here we simply output a dummy scalar equal to the cell's index.
		file << "CELL_DATA " << numTets << "\n";
		file << "SCALARS cellID int 1\n";
		file << "LOOKUP_TABLE default\n";
		for (size_t i = 0; i < numTets; ++i)
		{
			file << i << "\n";
		}

		file.close();
		return true;
	}

	namespace
	{
		// A hash functor for an array of three unsigned ints.
		struct FaceHash {
			std::size_t operator()(const std::array<unsigned int, 3>& face) const {
				std::size_t h1 = std::hash<unsigned int>{}(face[0]);
				std::size_t h2 = std::hash<unsigned int>{}(face[1]);
				std::size_t h3 = std::hash<unsigned int>{}(face[2]);
				return ((h1 * 31 + h2) * 31) + h3;
			}
		};

		// Equality functor for faces.
		struct FaceEqual {
			bool operator()(const std::array<unsigned int, 3>& a, const std::array<unsigned int, 3>& b) const {
				return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
			}
		};
	} // anonymous namespace

	bool ExportBaseTetraMeshGeometryDataToVTKPoly(const BaseTetraMeshGeometryData& geomData, const std::string& absFileName)
	{
		const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(absFileName);
		if (extension != "vtk")
		{
			std::cerr << absFileName << "ExportBaseTetraMeshGeometryDataToVTKPoly: Invalid file extension!" << std::endl;
			return false;
		}

		std::ofstream file(absFileName);
		if (!file.is_open())
		{
			std::cerr << "ExportBaseTetraMeshGeometryDataToVTKPoly: Failed to open file for writing: " << absFileName << std::endl;
			return false;
		}

		// Use an unordered_map to count faces.
	// Each face is stored as a sorted array of three vertex indices.
		std::unordered_map<std::array<unsigned int, 3>, int, FaceHash, FaceEqual> faceCount;

		// For each tetrahedron, extract its four faces.
		// A tetrahedron has vertices tet[0], tet[1], tet[2], tet[3]. Its faces are:
		//   {0,1,2}, {0,1,3}, {0,2,3}, and {1,2,3}.
		for (const auto& tet : geomData.TetrahedraIndices) {
			std::array<std::array<unsigned int, 3>, 4> faces = { {
				{ tet[0], tet[1], tet[2] },
				{ tet[0], tet[1], tet[3] },
				{ tet[0], tet[2], tet[3] },
				{ tet[1], tet[2], tet[3] }
			} };

			// Sort each face (so that identical faces get the same ordering)
			for (auto& face : faces) {
				std::sort(face.begin(), face.end());
				faceCount[face] += 1;
			}
		}

		// We now have a collection of unique faces (both interior and boundary).
		std::vector<std::array<unsigned int, 3>> uniqueFaces;
		uniqueFaces.reserve(faceCount.size());
		for (const auto& kv : faceCount) {
			uniqueFaces.push_back(kv.first);
		}

		// Write VTK header for POLYDATA.
		file << "# vtk DataFile Version 3.0\n";
		file << "VTK output from tetrahedral mesh (all faces)\n";
		file << "ASCII\n";
		file << "DATASET POLYDATA\n";

		// Write points.
		file << "POINTS " << geomData.Vertices.size() << " float\n";
		for (const auto& vertex : geomData.Vertices) {
			file << vertex[0] << " " << vertex[1] << " " << vertex[2] << "\n";
		}

		// Write polygons.
		// For each unique face (triangle), we write a line with 4 integers: "3 v0 v1 v2".
		size_t numFaces = uniqueFaces.size();
		size_t totalEntries = numFaces * 4; // each face: one count + 3 indices.
		file << "POLYGONS " << numFaces << " " << totalEntries << "\n";
		for (const auto& face : uniqueFaces) {
			file << "3 " << face[0] << " " << face[1] << " " << face[2] << "\n";
		}

		file.close();
		return true;
	}

	std::optional<BaseMeshGeometryData> ImportOBJMeshGeometryData(const std::string& absFileName, const bool& importInParallel, std::optional<std::vector<pmp::Scalar>*> chunkIdsVertexPropPtrOpt)
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

		for (size_t i = 0; i < thread_count; ++i) {
			char* chunk_start = file_start + (i * chunk_size);
			char* chunk_end = (i == thread_count - 1) ? file_end : chunk_start + chunk_size;

			// Adjust chunk_end to point to the end of a line
			while (*chunk_end != '\n' && chunk_end < file_end) {
				chunk_end++;
			}
			if (chunk_end != file_end) {
				chunk_end++;  // move past the newline character
			}

			// Start a thread to process this chunk
			threads[i] = std::thread(ParseChunk, chunk_start, chunk_end, std::ref(threadResults[i]));
		}

		// Wait for all threads to finish
		for (auto& t : threads) {
			t.join();
		}

		if (chunkIdsVertexPropPtrOpt.has_value() && !chunkIdsVertexPropPtrOpt.value()->empty())
		{
			chunkIdsVertexPropPtrOpt.value()->clear();
		}
		for (int threadId = 0; const auto& result : threadResults)
		{
			if (chunkIdsVertexPropPtrOpt.has_value())
			{
				const auto nVertsPerChunk = result.Vertices.size();
				chunkIdsVertexPropPtrOpt.value()->insert(chunkIdsVertexPropPtrOpt.value()->end(), nVertsPerChunk, static_cast<pmp::Scalar>(threadId));
			}
			resultData.Vertices.insert(resultData.Vertices.end(), result.Vertices.begin(), result.Vertices.end());
			resultData.PolyIndices.insert(resultData.PolyIndices.end(), result.PolyIndices.begin(), result.PolyIndices.end());
			resultData.VertexNormals.insert(resultData.VertexNormals.end(), result.VertexNormals.begin(), result.VertexNormals.end());

			++threadId;
		}

		// Clean up
		UnmapViewOfFile(file_memory);
		CloseHandle(file_mapping);
		CloseHandle(file_handle);

		return std::move(resultData);
	}

	std::optional<std::vector<pmp::vec3>> ImportPLYPointCloudData(const std::string& absFileName, const bool& importInParallel)
	{
		const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(absFileName);
		if (extension != "ply")
			return {};

		const char* file_path = absFileName.c_str();

		const HANDLE file_handle = CreateFile(file_path, GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);

		if (file_handle == INVALID_HANDLE_VALUE)
		{
			std::cerr << "ImportPLYPointCloudData [ERROR]: Failed to open the file.\n";
			return {};
		}

		// Get the file size
		const DWORD file_size = GetFileSize(file_handle, nullptr);

		// Create a file mapping object
		const HANDLE file_mapping = CreateFileMapping(file_handle, nullptr, PAGE_READONLY, 0, 0, nullptr);
		if (file_mapping == nullptr)
		{
			std::cerr << "ImportPLYPointCloudData [ERROR]: Failed to create file mapping.\n";
			CloseHandle(file_handle);
			return {};
		}

		// Map the file into memory
		const LPVOID file_memory = MapViewOfFile(file_mapping, FILE_MAP_READ, 0, 0, 0);
		if (file_memory == nullptr)
		{
			std::cerr << "ImportPLYPointCloudData [ERROR]: Failed to map the file.\n";
			CloseHandle(file_mapping);
			CloseHandle(file_handle);
			return {};
		}

		std::vector<pmp::vec3> resultData;

		char* file_start = static_cast<char*>(file_memory);
		char* file_end = file_start + file_size;

		// Read the PLY header to get the number of vertices and start position of vertex data
		const auto [vertexCount, vertexDataStart] = ReadPLYVertexHeader(file_start);
		if (vertexDataStart == 0)
		{
			std::cerr << "ImportPLYPointCloudData [ERROR]: Failed to read PLY header or no vertices found.\n";
			UnmapViewOfFile(file_memory);
			CloseHandle(file_mapping);
			CloseHandle(file_handle);
			return {};
		}

		// Adjust file_start to point to the beginning of vertex data
		file_start = vertexDataStart;

		// Ensure there is vertex data to process
		if (file_start >= file_end) 
		{
			std::cerr << "ImportPLYPointCloudData [ERROR]: No vertex data to process.\n";
			UnmapViewOfFile(file_memory);
			CloseHandle(file_mapping);
			CloseHandle(file_handle);
			return {};
		}

		// Adjust the file size for chunk calculation
		const size_t adjusted_file_size = file_end - file_start;

		// Determine the number of threads and chunk size
		const size_t thread_count = importInParallel ? std::thread::hardware_concurrency() : 1;
		const size_t chunk_size = std::max<size_t>(adjusted_file_size / thread_count, 1ul); // Ensure chunk size is at least 1
		std::vector<std::thread> threads(thread_count);
		std::vector<std::vector<pmp::vec3>> threadResults(thread_count);

		for (size_t i = 0; i < thread_count; ++i) 
		{
			char* chunk_start = file_start + (i * chunk_size);
			char* chunk_end = (i == thread_count - 1) ? file_end : chunk_start + chunk_size;

			// Adjust chunk_start to the beginning of a line (for all chunks except the first)
			if (i > 0) 
			{
				while (*chunk_start != '\n' && chunk_start < file_end) 
				{
					chunk_start++;
				}
				if (chunk_start < file_end)
				{
					chunk_start++;  // Move past the newline character
				}
			}

			// Adjust chunk_end to the end of a line
			while (*chunk_end != '\n' && chunk_end < file_end)
			{
				chunk_end++;
			}
			if (chunk_end < file_end)
			{
				chunk_end++;  // Move past the newline character
			}

			// Start a thread to process this chunk
			threads[i] = std::thread(ParsePointCloudChunk, chunk_start, chunk_end, std::ref(threadResults[i]));
		}

		// Wait for all threads to finish
		for (auto& t : threads) 
		{
			t.join();
		}

		for (const auto& result : threadResults)
		{
			resultData.insert(resultData.end(), result.begin(), result.end());
		}

		// Clean up
		UnmapViewOfFile(file_memory);
		CloseHandle(file_mapping);
		CloseHandle(file_handle);

		return std::move(resultData);
	}

	std::optional<std::vector<pmp::vec3>> ImportPLYPointCloudDataMainThread(const std::string& absFileName)
	{
		const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(absFileName);
		if (extension != "ply")
		{
			std::cerr << absFileName << " has invalid extension!" << std::endl;
			return {};
		}

		std::ifstream file(absFileName);
		if (!file.is_open()) 
		{
			std::cerr << "Failed to open the file." << std::endl;
			return {};
		}

		std::string line;
		size_t vertexCount = 0;
		bool headerEnded = false;

		// Read header to find the vertex count
		while (std::getline(file, line) && !headerEnded) 
		{
			std::istringstream iss(line);
			std::string token;
			iss >> token;

			if (token == "element") {
				iss >> token;
				if (token == "vertex") {
					iss >> vertexCount;
				}
			}
			else if (token == "end_header") {
				headerEnded = true;
			}
		}

		if (!headerEnded || vertexCount == 0) {
			std::cerr << "Invalid PLY header or no vertices found." << std::endl;
			return {};
		}

		std::vector<pmp::vec3> vertices;
		vertices.reserve(vertexCount);

		// Read vertex data
		while (std::getline(file, line)) {
			std::istringstream iss(line);
			pmp::vec3 vertex;
			if (!(iss >> vertex[0] >> vertex[1] >> vertex[2])) {
				std::cerr << "Error parsing vertex data: " << line << std::endl;
				continue;
			}

			vertices.push_back(vertex);
		}

		return vertices;
	}

	std::vector<pmp::Point> SamplePointsFromMeshData(const BaseMeshGeometryData& meshData, size_t nVerts, const std::optional<unsigned int>& seed)
	{
		const unsigned int nVertsTotal = meshData.Vertices.size();

		// Determine the actual number of vertices to sample
		const size_t numVertsToSample = std::min(nVerts, static_cast<size_t>(nVertsTotal));

		// Generate a list of all vertex indices
		std::vector<size_t> indices(nVertsTotal);
		std::iota(indices.begin(), indices.end(), 0);

		// Create a generator with the given seed or a random device
		std::mt19937 gen(seed ? *seed : std::random_device{}());

		// Perform a partial Fisher-Yates shuffle to get numVertsToSample indices
		for (size_t i = 0; i < numVertsToSample; ++i)
		{
			std::uniform_int_distribution<size_t> distrib(i, nVertsTotal - 1);
			size_t j = distrib(gen);
			std::swap(indices[i], indices[j]);
		}

		// Extract the sampled points based on the first numVertsToSample indices
		std::vector<pmp::Point> resultPts;
		resultPts.reserve(numVertsToSample);

		for (size_t i = 0; i < numVertsToSample; ++i)
		{
			resultPts.push_back(meshData.Vertices[indices[i]]);
		}

		return resultPts;
	}

	std::vector<std::pair<pmp::Point, pmp::vec3>> SamplePointsWithNormalsFromMeshData(const BaseMeshGeometryData& meshData, size_t nVerts, const std::optional<unsigned int>& seed)
	{
		const unsigned int nVertsTotal = meshData.Vertices.size();

		// Determine the actual number of vertices to sample
		const size_t numVertsToSample = std::min(nVerts, static_cast<size_t>(nVertsTotal));

		// Generate a list of all vertex indices
		std::vector<size_t> indices(nVertsTotal);
		std::iota(indices.begin(), indices.end(), 0);

		// Create a generator with the given seed or a random device
		std::mt19937 gen(seed ? *seed : std::random_device{}());

		// Perform a partial Fisher-Yates shuffle to get numVertsToSample indices
		for (size_t i = 0; i < numVertsToSample; ++i)
		{
			std::uniform_int_distribution<size_t> distrib(i, nVertsTotal - 1);
			size_t j = distrib(gen);
			std::swap(indices[i], indices[j]);
		}

		// Extract the sampled points based on the first numVertsToSample indices
		std::vector<std::pair<pmp::Point, pmp::vec3>> resultPts;
		resultPts.reserve(numVertsToSample);

		for (size_t i = 0; i < numVertsToSample; ++i)
		{
			resultPts.push_back({ meshData.Vertices[indices[i]], meshData.VertexNormals[indices[i]] });
		}

		return resultPts;
	
	}

	bool ExportSampledVerticesToPLY(const BaseMeshGeometryData& meshData, size_t nVerts, const std::string& absFileName, const std::optional<unsigned int>& seed)
	{
		const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(absFileName);
		if (extension != "ply")
		{
			std::cerr << "ExportSampledVerticesToPLY: " << absFileName << " has invalid extension!" << std::endl;
			return false;
		}

		std::ofstream file(absFileName);
		if (!file.is_open())
		{
			std::cerr << "ExportSampledVerticesToPLY: Failed to open file for writing: " << absFileName << std::endl;
			return false;
		}

		// Generate nVerts random indices
		std::vector<size_t> indices;
		std::mt19937 gen;
		if (seed.has_value())
		{
			gen.seed(seed.value());
		}
		else
		{
			std::random_device rd;
			gen = std::mt19937(rd());
		}
		std::uniform_int_distribution<> distrib(0, meshData.Vertices.size() - 1);

		for (size_t i = 0; i < nVerts; ++i)
		{
			indices.push_back(distrib(gen));
		}

		// Write to a .ply file
		std::ofstream outFile(absFileName);
		if (!outFile.is_open()) 
		{
			std::cerr << "ExportSampledVerticesToPLY: Unable to open file: " << absFileName << std::endl;
			return false;
		}

		// PLY header
		outFile << "ply\n";
		outFile << "format ascii 1.0\n";
		outFile << "element vertex " << nVerts << "\n";
		outFile << "property float x\n";
		outFile << "property float y\n";
		outFile << "property float z\n";
		outFile << "end_header\n";

		// Write sampled vertices
		for (auto idx : indices) 
		{
			const auto& vertex = meshData.Vertices[idx];
			outFile << vertex[0] << " " << vertex[1] << " " << vertex[2] << "\n";
		}

		outFile.close();
		return true;
	}

	bool ExportSampledVerticesWithNormalsToVTK(const BaseMeshGeometryData& meshData, size_t nVerts, const std::string& absFileName, const std::optional<unsigned int>& seed)
	{
		if (meshData.VertexNormals.empty())
		{
			std::cerr << "ExportSampledVerticesWithNormalsToVTK: Vertex normals are not available!" << std::endl;
			return false;
		}

		if (meshData.Vertices.empty())
		{
			std::cerr << "ExportSampledVerticesWithNormalsToVTK: Vertices are not available!" << std::endl;
			return false;
		}

		if (meshData.Vertices.size() != meshData.VertexNormals.size())
		{
			std::cerr << "ExportSampledVerticesWithNormalsToVTK: Vertices and VertexNormals are of different size!" << std::endl;
			return false;
		}

		const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(absFileName);
		if (extension != "vtk")
		{
			std::cerr << absFileName << "ExportSampledVerticesWithNormalsToVTK: Invalid file extension!" << std::endl;
			return false;
		}

		std::ofstream file(absFileName);
		if (!file.is_open())
		{
			std::cerr << "ExportSampledVerticesWithNormalsToVTK: Failed to open file for writing: " << absFileName << std::endl;
			return false;
		}

		const size_t nVertices = meshData.Vertices.size();

		// Generate nVerts random indices
		std::vector<size_t> indices;
		std::mt19937 gen;
		if (seed.has_value())
		{
			gen.seed(seed.value());
		}
		else
		{
			std::random_device rd;
			gen = std::mt19937(rd());
		}
		std::uniform_int_distribution<> distrib(0, nVertices - 1);

		for (size_t i = 0; i < nVerts; ++i)
		{
			indices.push_back(distrib(gen));
		}

		// Write to a .vtk file
		std::ofstream outFile(absFileName);
		if (!outFile.is_open())
		{
			std::cerr << "ExportSampledVerticesWithNormalsToVTK: Unable to open file: " << absFileName << std::endl;
			return false;
		}

		// VTK header
		outFile << "# vtk DataFile Version 3.0\n";
		outFile << "Sampled vertex normals\n";
		outFile << "ASCII\n";
		outFile << "DATASET UNSTRUCTURED_GRID\n";
		outFile << "POINTS " << nVerts << " float\n";

		// Write sampled vertices
		for (auto idx : indices)
		{
			const auto& vertex = meshData.Vertices[idx];
			outFile << vertex[0] << " " << vertex[1] << " " << vertex[2] << "\n";
		}

		outFile << "\nPOINT_DATA " << nVerts << "\n";
		outFile << "VECTORS normals float\n";

		// Write sampled vertex normals
		for (auto idx : indices)
		{
			const auto& normal = meshData.VertexNormals[idx];
			outFile << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
		}

		outFile.close();
		return true;
	}

	bool ExportSampledVerticesWithNormalsToPLY(const BaseMeshGeometryData& meshData, size_t nVerts, const std::string& absFileName, const std::optional<unsigned int>& seed)
	{
		if (meshData.VertexNormals.empty())
		{
			std::cerr << "ExportSampledVerticesWithNormalsToPLY: Vertex normals are not available!" << std::endl;
			return false;
		}

		if (meshData.Vertices.empty())
		{
			std::cerr << "ExportSampledVerticesWithNormalsToPLY: Vertices are not available!" << std::endl;
			return false;
		}

		if (meshData.Vertices.size() != meshData.VertexNormals.size())
		{
			std::cerr << "ExportSampledVerticesWithNormalsToPLY: Vertices and VertexNormals are of different size!" << std::endl;
			return false;
		}

		const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(absFileName);
		if (extension != "ply")
		{
			std::cerr << "ExportSampledVerticesWithNormalsToPLY: Invalid file extension!" << std::endl;
			return false;
		}

		std::ofstream file(absFileName);
		if (!file.is_open())
		{
			std::cerr << "ExportSampledVerticesWithNormalsToPLY: Failed to open file for writing: " << absFileName << std::endl;
			return false;
		}

		const size_t nVertices = meshData.Vertices.size();

		// Generate nVerts random indices
		std::vector<size_t> indices;
		std::mt19937 gen;
		if (seed.has_value())
		{
			gen.seed(seed.value());
		}
		else
		{
			std::random_device rd;
			gen = std::mt19937(rd());
		}
		std::uniform_int_distribution<> distrib(0, nVertices - 1);

		for (size_t i = 0; i < nVerts; ++i)
		{
			indices.push_back(distrib(gen));
		}

		// Write to a .ply file
		std::ofstream outFile(absFileName);
		if (!outFile.is_open())
		{
			std::cerr << "ExportSampledVerticesWithNormalsToPLY: Unable to open file: " << absFileName << std::endl;
			return false;
		}

		// PLY header
		outFile << "ply\n";
		outFile << "format ascii 1.0\n";
		outFile << "element vertex " << nVerts << "\n";
		outFile << "property float x\n";
		outFile << "property float y\n";
		outFile << "property float z\n";
		outFile << "property float nx\n";
		outFile << "property float ny\n";
		outFile << "property float nz\n";
		outFile << "end_header\n";

		// Write sampled vertices and normals
		for (auto idx : indices)
		{
			const auto& vertex = meshData.Vertices[idx];
			const auto& normal = meshData.VertexNormals[idx];
			outFile << vertex[0] << " " << vertex[1] << " " << vertex[2] << " ";
			outFile << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
		}

		outFile.close();
		return true;
	}

	bool ExportPointsToPLY(const BaseMeshGeometryData& meshData, const std::string& absFileName)
	{
		const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(absFileName);
		if (extension != "ply")
		{
			std::cerr << "Geometry::ExportPointsToPLY:" << absFileName << " has invalid extension!" << std::endl;
			return false;
		}

		if (meshData.Vertices.empty())
		{
			std::cerr << "Geometry::ExportPointsToPLY: No vertices to write!\n";
			return false;
		}

		// Write to a .ply file
		std::ofstream outFile(absFileName);
		if (!outFile.is_open())
		{
			std::cerr << "Geometry::ExportPointsToPLY: Unable to open file: " << absFileName << std::endl;
			return false;
		}

		const unsigned int nVerts = meshData.Vertices.size();
		const bool hasNormals = !meshData.VertexNormals.empty();

		// PLY header
		outFile << "ply\n";
		outFile << "format ascii 1.0\n";
		outFile << "element vertex " << nVerts << "\n";
		outFile << "property float x\n";
		outFile << "property float y\n";
		outFile << "property float z\n";
		if (hasNormals)
		{
			outFile << "property float nx\n";
			outFile << "property float ny\n";
			outFile << "property float nz\n";
		}
		outFile << "end_header\n";

		// Write vertices and optionally normals
		for (unsigned int i = 0; i < nVerts; ++i)
		{
			const auto& vertex = meshData.Vertices[i];
			outFile << vertex[0] << " " << vertex[1] << " " << vertex[2];
			if (hasNormals)
			{
				const auto& normal = meshData.VertexNormals[i];
				outFile << " " << normal[0] << " " << normal[1] << " " << normal[2];
			}
			outFile << "\n";
		}

		outFile.close();
		return true;
	}

	bool ExportPolylinesToOBJ(const std::vector<std::vector<pmp::vec3>>& polylines, const std::string& absFileName)
	{
		const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(absFileName);
		if (extension != "obj")
		{
			std::cerr << absFileName << " has invalid extension!" << std::endl;
			return false;
		}

		std::ofstream file(absFileName);
		if (!file.is_open()) 
		{
			std::cerr << "Failed to open file for writing: " << absFileName << std::endl;
			return false;
		}

		// Write vertices
		for (const auto& polyline : polylines)
		{
			for (const auto& vertex : polyline) 
			{
				file << "v " << vertex[0] << " " << vertex[1] << " " << vertex[2] << "\n";
			}
		}

		// Write polyline connections as lines
		size_t indexOffset = 1; // OBJ files are 1-indexed
		for (const auto& polyline : polylines)
		{
			if (polyline.size() < 2) continue; // Ensure there are at least two points to form a segment
			for (size_t i = 0; i < polyline.size() - 1; ++i)
			{
				file << "l " << (i + indexOffset) << " " << (i + indexOffset + 1) << "\n";
			}
			indexOffset += polyline.size();
		}

		file.close();
		return true;
	}

	bool ExportBoundaryEdgesToPLY(const pmp::SurfaceMesh& mesh, const std::string& absFileName)
	{
		const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(absFileName);
		if (extension != "ply")
		{
			std::cerr << "Geometry::ExportBoundaryEdgesToPLY: " << absFileName << " has invalid extension!" << std::endl;
			return false;
		}

		// Create a file stream
		std::ofstream outFile(absFileName);
		if (!outFile.is_open())
		{
			std::cerr << "Geometry::ExportBoundaryEdgesToPLY: Unable to open file: " << absFileName << std::endl;
			return false;
		}

		// Collect boundary vertices and their indices
		std::vector<pmp::Vertex> boundaryVertices;
		std::unordered_map<unsigned int, unsigned int> vertexIndexMap;
		unsigned int vertexCounter = 0;

		for (const auto v : mesh.vertices())
		{
			if (mesh.is_boundary(v))
			{
				boundaryVertices.push_back(v);
				vertexIndexMap[v.idx()] = vertexCounter++;
			}
		}

		const unsigned int nBoundaryVertices = boundaryVertices.size();

		if (nBoundaryVertices == 0)
		{
			std::cerr << "Geometry::ExportBoundaryEdgesToPLY: No boundary vertices found!\n";
			outFile.close();
			return false;
		}

		// PLY header
		outFile << "ply\n";
		outFile << "format ascii 1.0\n";
		outFile << "element vertex " << nBoundaryVertices << "\n";
		outFile << "property float x\n";
		outFile << "property float y\n";
		outFile << "property float z\n";

		// Now count the boundary edges
		unsigned int nBoundaryEdges = 0;
		for (const auto e : mesh.edges())
		{
			if (mesh.is_boundary(e))
				++nBoundaryEdges;
		}

		// Write edge count in PLY format (each edge is represented as a polyline of two vertices)
		outFile << "element edge " << nBoundaryEdges << "\n";
		outFile << "property int vertex1\n";
		outFile << "property int vertex2\n";
		outFile << "end_header\n";

		// Write boundary vertex positions
		for (const auto& v : boundaryVertices)
		{
			const auto& point = mesh.position(v);
			outFile << point[0] << " " << point[1] << " " << point[2] << "\n";
		}

		// Write boundary edges as polylines using vertex indices
		for (const auto e : mesh.edges())
		{
			if (mesh.is_boundary(e))
			{
				const auto v0 = mesh.vertex(e, 0); // First vertex of the edge
				const auto v1 = mesh.vertex(e, 1); // Second vertex of the edge
				outFile << vertexIndexMap[v0.idx()] << " " << vertexIndexMap[v1.idx()] << "\n";
			}
		}

		outFile.close();
		return true;
	}


	std::optional<BaseMeshGeometryData> ComputeConvexHullFromPoints(const std::vector<pmp::Point>& points)
	{
		if (points.size() < 4)
		{
			std::cerr << "Geometry::ComputeConvexHullFromPoints: Not enough points to form a convex hull!\n";
			return {};
		}

		BaseMeshGeometryData resultData;

		VCG_Mesh ptsMesh;
		FillVCGMeshWithPoints(ptsMesh, points);
		VCG_Mesh result;
		result.face.EnableFFAdjacency();
		vcg::tri::UpdateTopology<VCG_Mesh>::FaceFace(result);

		if (!vcg::tri::ConvexHull<VCG_Mesh, VCG_Mesh>::ComputeConvexHull(ptsMesh, result))
		{
			std::cerr << "Geometry::ComputeConvexHullFromPoints: internal error!\n";
			return {};
		}

		std::ranges::transform(result.vert, std::back_inserter(resultData.Vertices), [](const auto& vcgVert) {
			return pmp::Point{ vcgVert.P()[0], vcgVert.P()[1], vcgVert.P()[2] }; });
		resultData.PolyIndices = ExtractVertexIndicesFromVCGMesh(result);

		return resultData;
	}

	std::optional<BaseCurveGeometryData> ComputeConvexHullFrom2DPoints(const std::vector<pmp::Point2>& points)
	{
		if (points.size() < 3)
		{
			std::cerr << "Geometry::ComputeConvexHullFrom2DPoints: Not enough points to form a convex hull!\n";
			return {};
		}

		const auto delaunayMesh = ComputeDelaunayMeshFrom2DPoints(points);
		if (!delaunayMesh.has_value())
		{
			std::cerr << "Geometry::ComputeConvexHullFrom2DPoints: ComputeDelaunayMeshFrom2DPoints failed!\n";
			return {};
		}

		const auto pmpMesh = ConvertBufferGeomToPMPSurfaceMesh(*delaunayMesh);

		BaseCurveGeometryData resultCurve;

		std::map<pmp::Vertex, pmp::Vertex> origToNewBoundaryVertexMap;
		for (int newVId = 0; const auto v : pmpMesh.vertices())
		{
			if (!pmpMesh.is_boundary(v))
				continue;

			const auto& vPos = pmpMesh.position(v);
			resultCurve.Vertices.emplace_back(vPos[0], vPos[1]);
			origToNewBoundaryVertexMap[v] = pmp::Vertex(newVId); // remember new index for edges
			newVId++;
		}

		for (const auto e : pmpMesh.edges())
		{
			if (!pmpMesh.is_boundary(e))
				continue;

			const auto v0 = pmpMesh.vertex(e, 0);
			const auto v1 = pmpMesh.vertex(e, 1);

			resultCurve.EdgeIndices.emplace_back(
				static_cast<unsigned int>(origToNewBoundaryVertexMap[v0].idx()),
				static_cast<unsigned int>(origToNewBoundaryVertexMap[v1].idx())
			);
		}

		return resultCurve;
	}

	std::optional<pmp::ManifoldCurve2D> ComputePMPConvexHullFrom2DPoints(const std::vector<pmp::Point2>& points)
	{
		const auto baseCurveOpt = ComputeConvexHullFrom2DPoints(points);
		if (!baseCurveOpt.has_value())
			return {};

		try {
			const auto pmpCurve = ConvertBufferGeomToPMPManifoldCurve2D(*baseCurveOpt);
			return pmpCurve;
		}
		catch (...)
		{
			std::cerr << "Geometry::ComputePMPConvexHullFrom2DPoints: an internal error occured in ConvertBufferGeomToPMPManifoldCurve2D! Possibly non-manifold geometry!\n";
		}

		return {};
	}
	
	std::optional<pmp::SurfaceMesh> ComputePMPConvexHullFromPoints(const std::vector<pmp::Point>& points)
	{
		const auto baseMeshOpt = ComputeConvexHullFromPoints(points); // Uses VCG
		if (!baseMeshOpt.has_value())
			return {};

		try {
			const auto pmpMesh = ConvertBufferGeomToPMPSurfaceMesh(*baseMeshOpt);
			return pmpMesh;
		}
		catch (...)
		{
			std::cerr << "Geometry::ComputePMPConvexHullFromPoints: an internal error occured in ConvertBufferGeomToPMPSurfaceMesh! Possibly non-manifold geometry!\n";
		}

		return {};
	}

	std::optional<BaseMeshGeometryData> ComputeDelaunayMeshFrom2DPoints(const std::vector<pmp::Point2>& points)
	{
		if (points.empty())
		{
			std::cerr << "Geometry::ComputeDelaunayMeshFrom2DPoints: points.empty()!\n";
			return {};
		}

		if (points.size() < 3)
		{
			std::cerr << "Geometry::ComputeDelaunayMeshFrom2DPoints: points.size() < 3!\n";
			return {};
		}

		// Initialize Triangle's input and output structures.
		triangulateio in, out;
		// Zero out all fields.
		std::memset(&in, 0, sizeof(triangulateio));
		std::memset(&out, 0, sizeof(triangulateio));

		// Set up the input points.
		in.numberofpoints = static_cast<int>(points.size());
		in.numberofpointattributes = 0; // No extra attributes.
		in.pointlist = new TRI_REAL[in.numberofpoints * 2];
		for (size_t i = 0; i < points.size(); ++i)
		{
			// Each point occupies two REALs: x and y.
			in.pointlist[2 * i + 0] = static_cast<TRI_REAL>(points[i][0]);
			in.pointlist[2 * i + 1] = static_cast<TRI_REAL>(points[i][1]);
		}
		in.pointmarkerlist = nullptr; // No markers.

		// No segments, holes, or region constraints.
		in.numberofsegments = 0;
		in.numberofholes = 0;
		in.numberofregions = 0;

		// Use switches:
		//   z : zero-based indexing,
		//   Q : quiet operation.
		char switches[] = "zQ";

		// Call Triangle to compute the Delaunay triangulation.
		// We don't require Voronoi output (pass NULL for vorout).
		triangulate(switches, &in, &out, nullptr);

		// Free our allocated input pointlist.
		delete[] in.pointlist;

		if (out.numberoftriangles <= 0)
		{
			std::cerr << "Geometry::ComputeDelaunayMeshFrom2DPoints: No triangles generated!\n";
			// Free Triangle's output arrays if allocated.
			if (out.pointlist) trifree((int*)out.pointlist);
			if (out.trianglelist) trifree(out.trianglelist);
			return {};
		}
		if (out.numberofcorners != 3)
		{
			std::cerr << "Geometry::ComputeDelaunayMeshFrom2DPoints: Expected 3 corners per triangle, got "
				<< out.numberofcorners << "!\n";
			if (out.pointlist) trifree((int*)out.pointlist);
			if (out.trianglelist) trifree(out.trianglelist);
			return {};
		}

		// Build the result.
		BaseMeshGeometryData result;
		result.Vertices.reserve(points.size());
		for (const auto& pt2d : points)
		{
			// Embed each 2D point as a 3D point with z = 0.
			result.Vertices.push_back(pmp::Point{ pt2d[0], pt2d[1], static_cast<pmp::Scalar>(0.0) });
		}
		result.PolyIndices.reserve(out.numberoftriangles);
		// Triangle's output: out.trianglelist is an array of out.numberoftriangles * 3 ints.
		for (int i = 0; i < out.numberoftriangles; ++i)
		{
			int idx0 = out.trianglelist[3 * i + 0];
			int idx1 = out.trianglelist[3 * i + 1];
			int idx2 = out.trianglelist[3 * i + 2];
			result.PolyIndices.push_back({
				static_cast<unsigned int>(idx0),
				static_cast<unsigned int>(idx1),
				static_cast<unsigned int>(idx2)
				});
		}

		// Free Triangle's output arrays.
		if (out.pointlist) trifree((int*)out.pointlist);
		if (out.trianglelist) trifree(out.trianglelist);
		if (out.triangleattributelist) trifree((int*)out.triangleattributelist);
		if (out.neighborlist) trifree(out.neighborlist);
		if (out.segmentlist) trifree(out.segmentlist);
		if (out.segmentmarkerlist) trifree(out.segmentmarkerlist);
		if (out.holelist) trifree((int*)out.holelist);
		if (out.regionlist) trifree((int*)out.regionlist);
		if (out.edgelist) trifree(out.edgelist);
		if (out.edgemarkerlist) trifree(out.edgemarkerlist);
		if (out.normlist) trifree((int*)out.normlist);

		return result;
	}

	std::optional<BaseTetraMeshGeometryData> ComputeDelaunayTetrahedralMeshFromPoints(const std::vector<pmp::Point>& points)
	{
		if (points.empty())
		{
			std::cerr << "Geometry::ComputeDelaunayTetrahedralMeshFromPoints: points.empty()!\n";
			return {};
		}

		if (points.size() < 4)
		{
			std::cerr << "Geometry::ComputeDelaunayTetrahedralMeshFromPoints: points.size() < 4!\n";
			return {};
		}

		// Prepare TetGen input and output structures.
		tetgenio in, out;
		// The tetgenio constructor automatically calls initialize(), which sets all pointers to NULL.
		// We set the first number to 0 (zero-based indexing) and the mesh dimension to 3.
		in.firstnumber = 0;
		in.mesh_dim = 3;

		// Allocate the input point list: 3 REALs per point.
		int numPts = static_cast<int>(points.size());
		in.numberofpoints = numPts;
		in.pointlist = new REAL[numPts * 3];
		for (int i = 0; i < numPts; ++i) {
			in.pointlist[3 * i + 0] = static_cast<REAL>(points[i][0]);
			in.pointlist[3 * i + 1] = static_cast<REAL>(points[i][1]);
			in.pointlist[3 * i + 2] = static_cast<REAL>(points[i][2]);
		}
		// No additional point attributes or markers.
		in.numberofpointattributes = 0;
		in.pointattributelist = nullptr;
		in.pointmarkerlist = nullptr;

		// For a simple Delaunay tetrahedralization from a point set,
		// we do not supply facets, segments, holes, or regions.

		// Set up TetGen behavior.
		tetgenbehavior behavior;
		behavior.quiet = 1;      // Run quietly.
		behavior.zeroindex = 1;  // Use zero-based indexing.

		// Call tetrahedralize()
		try {
			tetrahedralize(&behavior, &in, &out, nullptr);
		}
		catch (...) {
			std::cerr << "Geometry::ComputeDelaunayTetrahedralMeshFromPoints: tetrahedralize() threw an exception!\n";
			return {};
		}

		// Check that TetGen produced tetrahedra.
		if (out.numberoftetrahedra <= 0 || out.tetrahedronlist == nullptr) {
			std::cerr << "Geometry::ComputeDelaunayTetrahedralMeshFromPoints: No tetrahedra produced!\n";
			return {};
		}
		if (out.numberofcorners != 4) {
			std::cerr << "Geometry::ComputeDelaunayTetrahedralMeshFromPoints: Expected tetrahedra with 4 corners, got "
				<< out.numberofcorners << "\n";
			return {};
		}

		BaseTetraMeshGeometryData result;

		// Copy output vertices. TetGen writes its output points into out.pointlist
		// (which may include Steiner points). Use out.numberofpoints.
		int outNumPts = out.numberofpoints;
		result.Vertices.resize(outNumPts);
		for (int i = 0; i < outNumPts; ++i) {
			result.Vertices[i] = pmp::Point(out.pointlist[3 * i + 0],
				out.pointlist[3 * i + 1],
				out.pointlist[3 * i + 2]);
		}

		// Extract tetrahedra connectivity.
		int numTets = out.numberoftetrahedra;
		result.TetrahedraIndices.reserve(numTets);
		for (int i = 0; i < numTets; ++i) {
			unsigned int a = static_cast<unsigned int>(out.tetrahedronlist[4 * i + 0]);
			unsigned int b = static_cast<unsigned int>(out.tetrahedronlist[4 * i + 1]);
			unsigned int c = static_cast<unsigned int>(out.tetrahedronlist[4 * i + 2]);
			unsigned int d = static_cast<unsigned int>(out.tetrahedronlist[4 * i + 3]);
			result.TetrahedraIndices.push_back({ a, b, c, d });
		}

		// If needed, additional error checking can be done here (e.g., verifying nonzero tetrahedron volumes).

		// When 'in' and 'out' go out of scope, their destructors call deinitialize(),
		// which frees all memory allocated by TetGen.
		return result;

	}

	std::pair<pmp::Point, pmp::Scalar> ComputeMeshBoundingSphere(const pmp::SurfaceMesh& mesh)
	{
		if (mesh.is_empty())
		{
			throw std::invalid_argument("Geometry::ComputeMeshBoundingSphere: Vertex positions are not available in the mesh.\n");
		}

		const auto vFirst = mesh.vertices().begin();
		auto vCurrent = vFirst;
		pmp::Point center = mesh.position(*vFirst);
		pmp::Scalar radius = 0.0;

		// Find the initial point farthest from the arbitrary start point
		for (const auto v : mesh.vertices()) 
		{
			if (norm(mesh.position(v) - center) < radius)
				continue;

			radius = norm(mesh.position(v)- center);
			vCurrent = v;
		}

		// Set sphere center to the midpoint of the initial and farthest point, radius as half the distance
		center = (center + mesh.position(*vCurrent)) * 0.5;
		radius = norm(mesh.position(*vCurrent) - center);

		// Ensure all points are within the sphere, adjust if necessary
		for (const auto v : mesh.vertices()) 
		{
			const pmp::Scalar dist = norm(mesh.position(v) - center);
			if (dist < radius)
				continue;

			const pmp::Scalar newRadius = (radius + dist) / 2;
			const pmp::Scalar moveBy = newRadius - radius;
			const pmp::Point direction = normalize(mesh.position(v) - center);
			center += direction * moveBy;
			radius = newRadius;
		}

		return { center, radius };
	}

	std::pair<pmp::Point, pmp::Scalar> ComputePointCloudBoundingSphere(const std::vector<pmp::Point>& points)
	{
		if (points.empty())
		{
			throw std::invalid_argument("Geometry::ComputePointCloudBoundingSphere: points.empty()!\n");
		}

		// Start with the first point as the center
		pmp::Point center = points[0];
		pmp::Scalar radius = 0.0;

		// First pass: find the farthest point from the initial point to set a rough sphere
		for (const auto& point : points)
		{
			const pmp::Scalar dist = norm(point - center);
			if (dist < radius)
				continue;
			radius = dist;
		}

		// Set initial sphere
		center = (center + points[0] + radius * normalize(points[0] - center)) / 2.0;
		radius /= 2.0;

		// Second pass: expand sphere to include all points
		for (const auto& point : points)
		{
			const pmp::Scalar dist = norm(point - center);
			if (dist < radius) // Point is inside the sphere
				continue;
			const pmp::Scalar newRadius = (radius + dist) / 2;
			const pmp::Scalar moveBy = newRadius - radius;
			const pmp::Point direction = normalize(point - center);
			center += direction * moveBy;
			radius = newRadius;
		}

		return { center, radius };
	}

	std::optional<BaseMeshGeometryData> ComputeBallPivotingMeshFromPoints(const std::vector<pmp::Point>& points, const pmp::Scalar& ballRadius, const pmp::Scalar& clusteringPercentageOfBallRadius, const pmp::Scalar& angleThreshold)
	{
		if (points.empty())
		{
			std::cerr << "Geometry::ComputeBallPivotingMeshFromPoints: No points to triangulate.\n";
			return {};
		}
		if (ballRadius < FLT_EPSILON)
		{
			std::cerr << "Geometry::ComputeBallPivotingMeshFromPoints: Invalid radius value: " << ballRadius << ".\n";
			return {};
		}

		// prepare data
		BaseMeshGeometryData resultData;
		resultData.Vertices = points;

		const auto clustering = clusteringPercentageOfBallRadius / 100.0;
		const auto angleRad = angleThreshold / 180.0 * M_PI;

		VCG_Mesh ptsMesh;
		FillVCGMeshWithPoints(ptsMesh, points);
		vcg::tri::BallPivoting bpa(ptsMesh, ballRadius, clustering, angleRad);
		bpa.BuildMesh();

		resultData.PolyIndices = ExtractVertexIndicesFromVCGMesh(ptsMesh);

		return resultData;
	}

	//
	// =======================================================================
	// ......... Nanoflann Kd-Tree Utils .....................................
	//

	template <typename T>
	struct PointCloud
	{
		struct Point
		{
			T x, y, z;
		};

		using coord_t = T;  //!< The type of each coordinate

		std::vector<Point> pts;

		// Must return the number of data points
		[[nodiscard]] size_t kdtree_get_point_count() const { return pts.size(); }

		// Returns the dim'th component of the idx'th point in the class:
		// Since this is inlined and the "dim" argument is typically an immediate
		// value, the
		//  "if/else's" are actually solved at compile time.
		[[nodiscard]] T kdtree_get_pt(const size_t idx, const size_t dim) const
		{
			if (dim == 0) return pts[idx].x;
			if (dim == 1) return pts[idx].y;
			return pts[idx].z;
		}

		// Optional bounding-box computation: return false to default to a standard
		// bbox computation loop.
		//   Return true if the BBOX was already computed by the class and returned
		//   in "bb" so it can be avoided to redo it again. Look at bb.size() to
		//   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
		template <class BBOX>
		[[nodiscard]] bool kdtree_get_bbox(BBOX& /* bb */) const
		{
			return false;
		}
	};

	// And this is the "dataset to kd-tree" adaptor class:
	template <typename Derived>
	struct PointCloudAdaptor
	{
		using coord_t = typename Derived::coord_t;

		const Derived& obj;  //!< A const ref to the data set origin

		/// The constructor that sets the data set source
		PointCloudAdaptor(const Derived& obj_) : obj(obj_) {}

		/// CRTP helper method
		[[nodiscard]] const Derived& derived() const { return obj; }

		// Must return the number of data points
		[[nodiscard]] size_t kdtree_get_point_count() const
		{
			return derived().pts.size();
		}

		// Returns the dim'th component of the idx'th point in the class:
		// Since this is inlined and the "dim" argument is typically an immediate
		// value, the
		//  "if/else's" are actually solved at compile time.
		[[nodiscard]] coord_t kdtree_get_pt(const size_t idx, const size_t dim) const
		{
			if (dim == 0) return derived().pts[idx].x;
			if (dim == 1) return derived().pts[idx].y;
			return derived().pts[idx].z;
		}

		// Optional bounding-box computation: return false to default to a standard
		// bbox computation loop.
		//   Return true if the BBOX was already computed by the class and returned
		//   in "bb" so it can be avoided to redo it again. Look at bb.size() to
		//   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
		template <class BBOX>
		bool kdtree_get_bbox(BBOX& /*bb*/) const
		{
			return false;
		}

	};  // end of PointCloudAdaptor

	std::optional<BaseMeshGeometryData> ComputePoissonMeshFromOrientedPoints(
		const std::vector<pmp::Point>& points, 
		const std::vector<pmp::Normal>& normals, 
		const PoissonReconstructionParams& params)
	{
		if (points.empty() || normals.empty())
		{
			std::cerr << "Geometry::ComputePoissonMeshFromOrientedPoints: points.empty() || normals.empty()!\n";
			return {};
		}

		if (points.size() != normals.size())
		{
			std::cerr << "Geometry::ComputePoissonMeshFromOrientedPoints: points.size() != normals.size()!\n";
			return {};
		}

		PoissonParam<Scalarm> pp;
		pp.MaxDepthVal = params.depth;
		pp.FullDepthVal = params.fullDepth;
		pp.CGDepthVal = params.cgDepth;
		pp.ScaleVal = params.scale;
		pp.SamplesPerNodeVal = params.samplesPerNode;
		pp.PointWeightVal = params.pointWeight;
		pp.ItersVal = params.iters;
		pp.ConfidenceFlag = params.confidence;
		pp.DensityFlag = true;
		pp.CleanFlag = params.preClean;
		pp.ThreadsVal = params.threads;

		VCG_Mesh ptsMesh;
		FillVCGMeshWithPoints(ptsMesh, points, normals);
		PoissonClean(ptsMesh, pp.ConfidenceFlag, pp.CleanFlag);
		const bool goodNormal = HasGoodNormal(ptsMesh);
		MeshModelPointStream<Scalarm> meshStream(ptsMesh);

		VCG_Mesh poissonMesh;
		_Execute<Scalarm, 2, BOUNDARY_NEUMANN, PlyColorAndValueVertex<Scalarm> >(&meshStream, ptsMesh.bbox, poissonMesh, pp, nullptr);

		BaseMeshGeometryData resultData;
		std::ranges::transform(poissonMesh.vert, std::back_inserter(resultData.Vertices), [](const auto& vcgVert) {
				return pmp::Point{ vcgVert.P()[0], vcgVert.P()[1], vcgVert.P()[2] };
		});
		resultData.PolyIndices = ExtractVertexIndicesFromVCGMesh(poissonMesh);

		return resultData;
	}

	pmp::Scalar ComputeMinInterVertexDistance(const std::vector<pmp::Point>& points)
	{
		if (points.empty())
		{
			std::cerr << "Geometry::ComputeMinInterVertexDistance: points.empty()!\n";
			return -1.0;
		}

		// convert points to a compatible dataset
		PointCloud<pmp::Scalar> nfPoints;
		nfPoints.pts.reserve(points.size());
		for (const auto& p : points)
			nfPoints.pts.push_back({ p[0], p[1], p[2] });
		using PointCloudAdapter = PointCloudAdaptor<PointCloud<pmp::Scalar>>;

		// construct a kd-tree index:
		using PointCloudKDTreeIndexAdapter = nanoflann::KDTreeSingleIndexAdaptor<
			nanoflann::L2_Simple_Adaptor<pmp::Scalar, PointCloudAdapter>, PointCloudAdapter, 3 /* dim */>;
		PointCloudKDTreeIndexAdapter indexAdapter(3 /*dim*/, nfPoints, { 10 /* max leaf */ });

		auto do_knn_search = [&indexAdapter](const pmp::Point& p) {
			// do a knn search
			const pmp::Scalar query_pt[3] = { p[0], p[1], p[2] };
			size_t num_results = 2;
			std::vector<uint32_t> ret_index(num_results);
			std::vector<pmp::Scalar> out_dist_sqr(num_results);
			num_results = indexAdapter.knnSearch(
				&query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);
			ret_index.resize(num_results);
			out_dist_sqr.resize(num_results);
			return out_dist_sqr[1];
		};

		pmp::Scalar minDistance = std::numeric_limits<pmp::Scalar>::max();
		for (const auto& point : points)
		{
			const auto currDistSq = do_knn_search(point);
			if (currDistSq < minDistance) minDistance = currDistSq;
		}
		return std::sqrt(minDistance);
	}

	pmp::Scalar ComputeNearestNeighborMeanInterVertexDistance(const std::vector<pmp::Point>& points, const size_t& nNeighbors)
	{
		if (points.empty()) 
		{
			std::cerr << "Geometry::ComputeNearestNeighborMeanInterVertexDistance: points.empty()!\n";
			return -1.0;
		}

		// convert points to a compatible dataset
		PointCloud<pmp::Scalar> nfPoints;
		nfPoints.pts.reserve(points.size());
		for (const auto& p : points)
			nfPoints.pts.push_back({ p[0], p[1], p[2] });
		using PointCloudAdapter = PointCloudAdaptor<PointCloud<pmp::Scalar>>;

		// construct a kd-tree index:
		using PointCloudKDTreeIndexAdapter = nanoflann::KDTreeSingleIndexAdaptor<
			nanoflann::L2_Simple_Adaptor<pmp::Scalar, PointCloudAdapter>, PointCloudAdapter, 3 /* dim */>;
		PointCloudKDTreeIndexAdapter indexAdapter(3 /*dim*/, nfPoints, { 10 /* max leaf */ });

		auto do_knn_search = [&indexAdapter, &nNeighbors](const pmp::Point& p) -> pmp::Scalar {
			// do a knn search
			const pmp::Scalar query_pt[3] = { p[0], p[1], p[2] };
			std::vector<uint32_t> ret_index(nNeighbors);
			std::vector<pmp::Scalar> out_dist_sqr(nNeighbors);
			const size_t num_results = indexAdapter.knnSearch(
				&query_pt[0], nNeighbors, &ret_index[0], &out_dist_sqr[0]);
			ret_index.resize(num_results);
			out_dist_sqr.resize(num_results);
			if (num_results == 0) return 0.0;
			pmp::Scalar totalDistSq{ 0.0 };
			for (size_t i = 0; i < num_results; ++i)
				totalDistSq += out_dist_sqr[i];
			return sqrt(totalDistSq / (static_cast<pmp::Scalar>(num_results - 1.0)));
		};

		pmp::Scalar totalDistance = 0.0;
		for (const auto& point : points)
		{
			totalDistance += do_knn_search(point);
		}

		return totalDistance / static_cast<pmp::Scalar>(points.size());
	}

	//
	// ===========================================================================================
	//

	pmp::Scalar ComputeMinInterVertexDistanceBruteForce(const std::vector<pmp::Point>& points)
	{
		std::cout << "Geometry::ComputeMinInterVertexDistanceBruteForce [WARNING]: This function is brute-force. Not recommended for large data!\n";
		if (points.size() < 2)
		{
			std::cerr << "Geometry::ComputeMinInterVertexDistanceBruteForce: points.size() < 2! No meaningful distance can be computed!\n";
			return -1.0;
		}

		pmp::Scalar minDistance = std::numeric_limits<pmp::Scalar>::max();
		for (size_t i = 0; i < points.size(); ++i)
		{
			for (size_t j = i + 1; j < points.size(); ++j) 
			{
				const pmp::Scalar dist = norm(points[i] - points[j]);
				if (dist < minDistance) 
				{
					minDistance = dist;
				}
			}
		}
		return minDistance;
	}

	pmp::Scalar ComputeMaxInterVertexDistanceBruteForce(const std::vector<pmp::Point>& points)
	{
		std::cout << "Geometry::ComputeMaxInterVertexDistanceBruteForce [WARNING]: This function is brute-force. Not recommended for large data!\n";
		if (points.size() < 2)
		{
			std::cerr << "Geometry::ComputeMaxInterVertexDistanceBruteForce: points.size() < 2! No meaningful distance can be computed!\n";
			return -1.0;
		}

		pmp::Scalar maxDistance = 0;
		for (size_t i = 0; i < points.size(); ++i) 
		{
			for (size_t j = i + 1; j < points.size(); ++j)
			{
				const pmp::Scalar dist = norm(points[i] - points[j]);
				if (dist > maxDistance)
				{
					maxDistance = dist;
				}
			}
		}
		return maxDistance;
	}

	pmp::Scalar ComputeMeanInterVertexDistanceBruteForce(const std::vector<pmp::Point>& points)
	{
		std::cout << "Geometry::ComputeMeanInterVertexDistanceBruteForce [WARNING]: This function is brute-force. Not recommended for large data!\n";
		if (points.size() < 2)
		{
			std::cerr << "Geometry::ComputeMeanInterVertexDistanceBruteForce: points.size() < 2! No meaningful distance can be computed!\n";
			return -1.0;
		}

		pmp::Scalar totalDistance = 0.0;
		size_t count = 0;
		for (size_t i = 0; i < points.size(); ++i) 
		{
			for (size_t j = i + 1; j < points.size(); ++j) 
			{
				const pmp::Scalar dist = norm(points[i] - points[j]);
				totalDistance += dist;
				++count;
			}
		}
		return count > 0 ? (totalDistance / count) : 0;
	}

	int GetClosestPointIndex2D(const std::vector<pmp::Point>& points, const pmp::Point& sampledPoint)
	{
		if (points.empty())
		{
			std::cerr << "ComputeClosestPointIndex2D: points.empty()!\n";
			return -1;
		}

		// Convert points to a compatible dataset
		PointCloud<pmp::Scalar> nfPoints;
		nfPoints.pts.reserve(points.size());
		for (const auto& p : points)
			nfPoints.pts.push_back({ p[0], p[1], 0.0 });  // Nullify the z-coordinate

		using PointCloudAdapter = PointCloudAdaptor<PointCloud<pmp::Scalar>>;

		std::cout << "Points converted for KD-Tree." << std::endl;
		std::cout << "Initializing KD-Tree with " << nfPoints.pts.size() << " points." << std::endl;

		// Construct a 2D KD-tree index
		using PointCloudKDTreeIndexAdapter = nanoflann::KDTreeSingleIndexAdaptor<
			nanoflann::L2_Simple_Adaptor<pmp::Scalar, PointCloudAdapter>, PointCloudAdapter, 2 /* dim */>;
		PointCloudKDTreeIndexAdapter indexAdapter(2 /*dim*/, nfPoints, { 10 /* max leaf */ });
		std::cout << "KD-Tree initialized." << std::endl;

		// Do a knn search for the closest point to the sampledPoint in 2D
		const pmp::Scalar query_pt[2] = { sampledPoint[0], sampledPoint[1] };
		const size_t num_results = 1;
		std::vector<uint32_t> ret_index(num_results);
		std::vector<pmp::Scalar> out_dist_sqr(num_results);

		// Ensure query point is valid
		assert(query_pt != nullptr);

		std::cout << "Performing knnSearch..." << std::endl;
		try {
			size_t results_found = indexAdapter.knnSearch(&query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);
			std::cout << "knnSearch found " << results_found << " results." << std::endl;

			if (results_found == 0) {
				std::cerr << "knnSearch found no results." << std::endl;
				return -1;
			}

			return ret_index[0];
		}
		catch (const std::exception& e) {
			std::cerr << "Exception during knnSearch: " << e.what() << std::endl;
			return -1;
		}
	}

	int GetClosestPointIndex(const std::vector<pmp::Point>& points, const pmp::Point& sampledPoint)
	{
		if (points.empty())
		{
			std::cerr << "ComputeClosestPointIndex: points.empty()!\n";
			return -1;
		}

		// Convert points to a compatible dataset
		PointCloud<pmp::Scalar> nfPoints;
		nfPoints.pts.reserve(points.size());
		for (const auto& p : points)
			nfPoints.pts.push_back({ p[0], p[1], p[2] });
		using PointCloudAdapter = PointCloudAdaptor<PointCloud<pmp::Scalar>>;

		// Construct a KD-tree index
		using PointCloudKDTreeIndexAdapter = nanoflann::KDTreeSingleIndexAdaptor<
			nanoflann::L2_Simple_Adaptor<pmp::Scalar, PointCloudAdapter>, PointCloudAdapter, 3 /* dim */>;
		const PointCloudKDTreeIndexAdapter indexAdapter(3 /*dim*/, nfPoints, { 10 /* max leaf */ });

		// Do a knn search for the closest point to the sampledPoint
		const pmp::Scalar query_pt[3] = { sampledPoint[0], sampledPoint[1], sampledPoint[2] };
		size_t num_results = 1;
		std::vector<uint32_t> ret_index(num_results);
		std::vector<pmp::Scalar> out_dist_sqr(num_results);
		num_results = indexAdapter.knnSearch(
			&query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);
		ret_index.resize(num_results);
		out_dist_sqr.resize(num_results);

		// Return the index of the closest point
		if (num_results > 0)
			return ret_index[0];
		return -1;  // Indicate failure
	}

	//
	// =======================================================================================================================
	//

	std::optional<pmp::Scalar> GetDistanceToClosestPoint2DSquared(const PointCloud2DTree& kdTree, const pmp::Point2& sampledPoint)
	{
		if (kdTree.size_ == 0)
		{
			std::cerr << "GetDistanceToClosestPoint2DSquared: kdTree contains no points!\n";
			return {};
		}

		// Find the nearest neighbor distance using KD-tree
		const pmp::Scalar query_pt[2] = { sampledPoint[0], sampledPoint[1] };
		size_t nearest_index;
		pmp::Scalar out_dist_sqr;
		nanoflann::KNNResultSet<pmp::Scalar> resultSet(1);
		resultSet.init(&nearest_index, &out_dist_sqr);
		kdTree.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(10));

		if (resultSet.size() == 0)
		{
			std::cerr << "GetDistanceToClosestPoint2DSquared: resultSet contains no points!\n";
			return {};
		}

		return out_dist_sqr;
	}

	std::optional<pmp::Scalar> GetDistanceToClosestPoint3DSquared(const PointCloud3DTree& kdTree, const pmp::Point& sampledPoint)
	{
		if (kdTree.size_ == 0)
		{
			std::cerr << "GetDistanceToClosestPoint3DSquared: kdTree contains no points!\n";
			return {};
		}

		// Find the nearest neighbor distance using KD-tree
		const pmp::Scalar query_pt[3] = { sampledPoint[0], sampledPoint[1], sampledPoint[2] };
		size_t nearest_index;
		pmp::Scalar out_dist_sqr;
		nanoflann::KNNResultSet<pmp::Scalar> resultSet(1);
		resultSet.init(&nearest_index, &out_dist_sqr);
		kdTree.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(10));

		if (resultSet.size() == 0)
		{
			std::cerr << "GetDistanceToClosestPoint3DSquared: resultSet contains no points!\n";
			return {};
		}

		return out_dist_sqr;
	}

	pmp::Scalar ComputeNearestNeighborMeanInterVertexDistance2D(const std::vector<pmp::Point2>& points, const size_t& nNeighbors)
	{
		if (points.empty())
		{
			std::cerr << "Geometry::ComputeNearestNeighborMeanInterVertexDistance2D: points.empty()!\n";
			return -1.0;
		}

		// convert points to a compatible dataset
		PointCloud<pmp::Scalar> nfPoints;
		nfPoints.pts.reserve(points.size());
		for (const auto& p : points)
			nfPoints.pts.push_back({ p[0], p[1] });
		using PointCloudAdapter = PointCloudAdaptor<PointCloud<pmp::Scalar>>;

		// construct a kd-tree index:
		using PointCloudKDTreeIndexAdapter = nanoflann::KDTreeSingleIndexAdaptor<
			nanoflann::L2_Simple_Adaptor<pmp::Scalar, PointCloudAdapter>, PointCloudAdapter, 2 /* dim */>;
		PointCloudKDTreeIndexAdapter indexAdapter(2 /*dim*/, nfPoints, { 10 /* max leaf */ });

		auto do_knn_search = [&indexAdapter, &nNeighbors](const pmp::Point2& p) -> pmp::Scalar {
			// do a knn search
			const pmp::Scalar query_pt[2] = { p[0], p[1] };
			std::vector<uint32_t> ret_index(nNeighbors);
			std::vector<pmp::Scalar> out_dist_sqr(nNeighbors);
			const size_t num_results = indexAdapter.knnSearch(
				&query_pt[0], nNeighbors, &ret_index[0], &out_dist_sqr[0]);
			ret_index.resize(num_results);
			out_dist_sqr.resize(num_results);
			if (num_results == 0) return 0.0;
			pmp::Scalar totalDistSq{ 0.0 };
			for (size_t i = 0; i < num_results; ++i)
				totalDistSq += out_dist_sqr[i];
			return sqrt(totalDistSq / (static_cast<pmp::Scalar>(num_results - 1.0)));
		};

		pmp::Scalar totalDistance = 0.0;
		for (const auto& point : points)
		{
			totalDistance += do_knn_search(point);
		}

		return totalDistance / static_cast<pmp::Scalar>(points.size());
	}

	//
	// =======================================================================================================================
	//


	std::vector<pmp::Point2> GetSliceOfThePointCloud(const std::vector<pmp::Point>& points, const pmp::Point& planePt, const pmp::vec3& planeNormal, const pmp::Scalar& distTolerance)
	{
		if (std::fabs(sqrnorm(planeNormal) - 1.0) > FLT_EPSILON)
		{
			std::cerr << "GetSliceOfThePointCloud: planeNormal needs to be normalized!\n";
			return {};
		}

		std::vector<pmp::Point2> slicedPoints;

		// Find two orthogonal vectors in the plane
		pmp::vec3 axis2;
		if (std::abs(planeNormal[0]) > std::abs(planeNormal[1]))
		{
			axis2 = normalize(pmp::vec3(-planeNormal[2], 0, planeNormal[0]));
		}
		else
		{
			axis2 = normalize(pmp::vec3(0, planeNormal[2], -planeNormal[1]));
		}
		const pmp::vec3 axis1 = normalize(cross(planeNormal, -axis2));

		for (const auto& point : points)
		{
			// Calculate the vector from the plane point to the current point
			pmp::vec3 vec = point - planePt;

			// Calculate the signed distance from the point to the plane
			const pmp::Scalar signedDistanceToPlane = dot(vec, planeNormal);

			// Check if the point is within the distance tolerance
			if (std::abs(signedDistanceToPlane) <= distTolerance)
			{
				// Project the point onto the plane
				pmp::Point projectedPoint = point - signedDistanceToPlane * planeNormal;

				// Convert the projected point to 2D by using the plane axes
				pmp::Point2 point2D(dot(projectedPoint, axis1), dot(projectedPoint, axis2));

				// Add the 2D point to the result
				slicedPoints.push_back(point2D);
			}
		}

		return slicedPoints;
	}

	namespace
	{
		// Function to sample points along a parabolic segment
		[[nodiscard]] std::vector<pmp::Point2> SampleParabola(const MAT::Parabola& parabola, const Vector2d& start, const Vector2d& end, size_t numSamples)
		{
			std::vector<pmp::Point2> sampledPoints;
			double step = (end.x() - start.x()) / static_cast<double>(numSamples);

			for (size_t i = 0; i <= numSamples; ++i)
			{
				double x = start.x() + i * step;
				auto ys = parabola.getY(x); // Get the two possible y-values for the parabola at x
				double y = (i == 0) ? ys.first : ys.second;  // Pick the appropriate branch of the parabola

				sampledPoints.emplace_back(pmp::Point2(x, y));
			}

			return sampledPoints;
		}
	} // anonymous namespace

	std::optional<BaseCurveGeometryData> CalculateApproxMedialAxisFromCurve(const pmp::ManifoldCurve2D& curve)
	{
		if (curve.is_empty())
		{
			std::cerr << "CalculateApproxMedialAxisFromCurve: curve is empty!\n";
			return {};
		}

		if (!curve.is_closed())
		{
			std::cerr << "CalculateApproxMedialAxisFromCurve: curve is not closed!\n";
			return {};
		}

		// Convert the ManifoldCurve2D vertices to a vector of Vector2d in reverse order (CW)
		std::vector<Vector2d> customShapePoints;
		for (auto it = curve.vertices().end() - 1; it != curve.vertices().begin(); --it)
		{
			const pmp::Point2& point = curve.position(*it);
			customShapePoints.push_back(Vector2d(point[0], point[1]));
		}

		// Initialize the medial axis data structure
		Geometry::BaseCurveGeometryData medialAxisData;

		// Generate boundary elements using the BoundaryGenerator
		MAT::BoundaryGenerator boundaryGen;
		std::vector<MAT::BoundaryElement> boundaryElements = boundaryGen.getBoundaryElementsFromCurve(customShapePoints);

		if (boundaryElements.empty())
		{
			std::cerr << "CalculateApproxMedialAxisFromCurve: Failed to generate boundary elements!\n";
			return {};
		}

		// force transitions for boundaryElements
		for (size_t i = 0; i < boundaryElements.size(); ++i)
		{
			size_t iNext = (i + 1) % boundaryElements.size();
			size_t iPrev = (i + boundaryElements.size() - 1) % boundaryElements.size();
			boundaryElements[i].transForward = iNext;
			boundaryElements[iNext].transBack = iPrev;
		}

		// Run the medial axis transform using the generated boundary elements
		MAT::MedialAxisTransform mat;
		mat.setBoundaryElements(boundaryElements);
		std::vector<MAT::Path> medialPaths = mat.run();

		// debug medialPaths
		for (const auto& path : medialPaths)
		{
			std::cout << "Path: { " << path.keyPoint1.x() << ", " << path.keyPoint1.y() << " } -> { " << path.keyPoint2.x() << ", " << path.keyPoint2.y() << " }\n";
		}

		// For each path, sample the parabolic segments and store them
		size_t numSamplesPerParabola = 20;  // Adjust the number of samples per parabolic segment as needed
		for (const auto& path : medialPaths)
		{
			if (path.parabola.set) 
			{
				// Sample points along the parabola
				auto sampledPoints = SampleParabola(path.parabola, path.keyPoint1, path.keyPoint2, numSamplesPerParabola);

				// Add the sampled points to the medial axis data
				for (const auto& point : sampledPoints)
				{
					medialAxisData.Vertices.push_back(point);
				}

				// Add edges between consecutive sampled points
				for (size_t i = 0; i < sampledPoints.size() - 1; ++i)
				{
					medialAxisData.EdgeIndices.push_back({ medialAxisData.Vertices.size() - sampledPoints.size() + i,
						medialAxisData.Vertices.size() - sampledPoints.size() + i + 1 });
				}
			}
			else 
			{
				// Handle linear segments by directly connecting the key points
				const auto p1 = pmp::Point2(path.keyPoint1.x(), path.keyPoint1.y());
				const auto p2 = pmp::Point2(path.keyPoint2.x(), path.keyPoint2.y());
				if (norm(p1 - p2) < 1e-6) 
					continue; // Skip zero-length segments

				medialAxisData.Vertices.push_back(p1);
				medialAxisData.Vertices.push_back(p2);
				medialAxisData.EdgeIndices.push_back({ medialAxisData.Vertices.size() - 2, medialAxisData.Vertices.size() - 1 });
			}
		}

		return medialAxisData;
	}

	std::optional<BaseCurveGeometryData> GetMedialAxisOfSawhneysStupidMATAlgorithm(unsigned char shape)
	{
		MAT::BoundaryGenerator boundaryGen;
		std::vector<MAT::BoundaryElement> boundaryElements = boundaryGen.getBoundaryElements(shape);

		MAT::MedialAxisTransform mat;
		mat.setBoundaryElements(boundaryElements);
		std::vector<MAT::Path> medialPaths = mat.run();

		// debug medialPaths
		for (const auto& path : medialPaths)
		{
			std::cout << "Path: { " << path.keyPoint1.x() << ", " << path.keyPoint1.y() << " } -> { " << path.keyPoint2.x() << ", " << path.keyPoint2.y() << " }\n";
		}

		// Initialize the medial axis data structure
		Geometry::BaseCurveGeometryData medialAxisData;

		// For each path, sample the parabolic segments and store them
		size_t numSamplesPerParabola = 20;  // Adjust the number of samples per parabolic segment as needed
		for (const auto& path : medialPaths)
		{
			if (path.parabola.set)
			{
				// Sample points along the parabola
				auto sampledPoints = SampleParabola(path.parabola, path.keyPoint1, path.keyPoint2, numSamplesPerParabola);

				// Add the sampled points to the medial axis data
				for (const auto& point : sampledPoints)
				{
					medialAxisData.Vertices.push_back(point);
				}

				// Add edges between consecutive sampled points
				for (size_t i = 0; i < sampledPoints.size() - 1; ++i)
				{
					medialAxisData.EdgeIndices.push_back({ medialAxisData.Vertices.size() - sampledPoints.size() + i,
						medialAxisData.Vertices.size() - sampledPoints.size() + i + 1 });
				}
			}
			else
			{
				// Handle linear segments by directly connecting the key points
				const auto p1 = pmp::Point2(path.keyPoint1.x(), path.keyPoint1.y());
				const auto p2 = pmp::Point2(path.keyPoint2.x(), path.keyPoint2.y());
				if (norm(p1 - p2) < 1e-6)
					continue; // Skip zero-length segments

				medialAxisData.Vertices.push_back(p1);
				medialAxisData.Vertices.push_back(p2);
				medialAxisData.EdgeIndices.push_back({ medialAxisData.Vertices.size() - 2, medialAxisData.Vertices.size() - 1 });
			}
		}

		return medialAxisData;
	}

	//
	// =================================================================================================
	//

	std::vector<std::pair<unsigned int, unsigned int>> GetBoundaryPointsOfPointCloudGaps2D(const std::vector<pmp::Point2>& points)
	{
		if (points.empty())
		{
			std::cerr << "Geometry::GetBoundaryPointsOfPointCloudGaps2D: points.empty()!\n";
			return {};
		}

		// Build the KD-tree
		PointCloud2D pointCloud;
		pointCloud.points = points;
		PointCloud2DTree kdTree(2, pointCloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
		kdTree.buildIndex();

		// evaluate min distance
		auto do_knn_search = [&kdTree](const pmp::Point2& p) {
			// do a knn search
			const pmp::Scalar query_pt[2] = { p[0], p[1] };
			size_t num_results = 2;
			std::vector<uint32_t> ret_index(num_results);
			std::vector<pmp::Scalar> out_dist_sqr(num_results);
			num_results = kdTree.knnSearch(
				&query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);
			ret_index.resize(num_results);
			out_dist_sqr.resize(num_results);
			return out_dist_sqr[1];
		};

		pmp::Scalar minDistance = std::numeric_limits<pmp::Scalar>::max();
		for (const auto& point : points)
		{
			const auto currDistSq = do_knn_search(point);
			if (currDistSq < minDistance) minDistance = currDistSq;
		}

		std::vector<std::pair<unsigned int, unsigned int>> result;



		return result;
	}

	std::vector<std::vector<pmp::Point>> GetPointClusters(const std::vector<pmp::Point>& points, const pmp::Scalar& criticalRadius)
	{
		if (points.empty())
		{
			std::cerr << "Geometry::GetPointClusters: points.empty()!\n";
			return { points };
		}

		const uint32_t N = points.size();
		std::vector<std::vector<pmp::Point>> clusters;
		if (N == 0)
			return clusters;

		// 1) Fill the PointCloud3D
		PointCloud3D cloud;
		cloud.points = points;  // copy

		// 2) Build the KD-tree index
		PointCloud3DTree index(/*dim=*/3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(/*max leaf=*/10));
		index.buildIndex();

		// 3) Prepare visit flags and search params
		std::vector<bool> visited(N, false);
		const pmp::Scalar radius2 = criticalRadius * criticalRadius;
		std::vector<nanoflann::ResultItem<uint32_t, pmp::Scalar>> results;
		nanoflann::SearchParameters searchParams;

		// 4) Flood‐fill over epsilon‐neighborhoods
		for (uint32_t seed = 0; seed < N; ++seed)
		{
			if (visited[seed])
				continue;

			std::vector<pmp::Point> cluster;
			std::queue<uint32_t> q;
			visited[seed] = true;
			q.push(seed);

			while (!q.empty()) 
			{
				uint32_t idx = q.front();
				q.pop();
				cluster.push_back(points[idx]);

				// query neighbors within criticalRadius
				pmp::Scalar query_pt[3] = {
					points[idx][0],
					points[idx][1],
					points[idx][2]
				};
				results.clear();
				const uint32_t nMatches = index.radiusSearch(&query_pt[0], radius2, results, searchParams);

				for (uint32_t i = 0; i < nMatches; ++i)
				{
					uint32_t nbr = results[i].first;
					if (visited[nbr])
						continue;

					visited[nbr] = true;
					q.push(nbr);
				}
			}

			if (cluster.empty())
				continue; // no point found

			clusters.push_back(std::move(cluster));
		}

		return clusters;
	}

	void Get3DPointSearchIndex(
		const std::vector<pmp::Point>& points,
		PointCloud3D& outCloud,
		std::unique_ptr<PointCloud3DTree>& outTree)
	{
		outCloud.points = points;   // copy data into the cloud
		// construct the tree on the heap so we can hold a unique_ptr to it
		outTree = std::make_unique<PointCloud3DTree>(
			/*dim=*/3, outCloud,
			nanoflann::KDTreeSingleIndexAdaptorParams(/*max leaf=*/10)
			);
		outTree->buildIndex();
	}

	std::unique_ptr<PointSearchIndex3D> Get3DPointSearchIndex(const std::vector<pmp::Point>& points)
	{
		return std::make_unique<PointSearchIndex3D>(points);
	}

	pmp::Scalar ComputeNearestNeighborMeanInterVertexDistance(const PointCloud3D& cloud, PointCloud3DTree& tree, const size_t& nNeighbors)
	{
		if (cloud.points.empty())
		{
			std::cerr << "Geometry::ComputeNearestNeighborMeanInterVertexDistance: cloud.points.empty()!\n";
			return -1.0;
		}

		std::vector<uint32_t> ret_idx(nNeighbors);
		std::vector<pmp::Scalar> out_dist2(nNeighbors);
		pmp::Scalar sumMeanDists{ 0 };

		const size_t N = cloud.points.size();
		for (size_t i = 0; i < N; ++i) {
			// form query
			pmp::Scalar q[3] = {
				cloud.points[i][0],
				cloud.points[i][1],
				cloud.points[i][2]
			};
			// kNN search (self + neighbors)
			const size_t found = tree.knnSearch(&q[0], nNeighbors,
				ret_idx.data(),
				out_dist2.data());
			if (found > 1) 
			{
				// skip the zero‐distance to itself
				pmp::Scalar accSq = 0;
				for (size_t k = 1; k < found; ++k)
					accSq += out_dist2[k];
				sumMeanDists += std::sqrt(accSq / (found - 1));
			}
		}
		return sumMeanDists / static_cast<pmp::Scalar>(N);
	}

	std::vector<std::vector<pmp::Point>> GetPointClusters(const PointCloud3D& cloud, PointCloud3DTree& tree, const pmp::Scalar& criticalRadius)
	{
		if (cloud.points.empty())
		{
			std::cerr << "Geometry::GetPointClusters: cloud.points.empty()!\n";
			return { cloud.points };
		}

		std::vector<std::vector<pmp::Point>> clusters;

		// Prepare visit flags and search params
		const uint32_t N = static_cast<uint32_t>(cloud.points.size());
		std::vector<bool> visited(N, false);
		const pmp::Scalar radius2 = criticalRadius * criticalRadius;
		std::vector<nanoflann::ResultItem<uint32_t, pmp::Scalar>> results;
		nanoflann::SearchParameters searchParams;

		// Flood‐fill over epsilon‐neighborhoods
		for (uint32_t seed = 0; seed < N; ++seed)
		{
			if (visited[seed])
				continue;

			std::vector<pmp::Point> cluster;
			std::queue<uint32_t> q;
			visited[seed] = true;
			q.push(seed);

			while (!q.empty())
			{
				uint32_t idx = q.front();
				q.pop();
				cluster.push_back(cloud.points[idx]);

				// query neighbors within criticalRadius
				pmp::Scalar query_pt[3] = {
					cloud.points[idx][0],
					cloud.points[idx][1],
					cloud.points[idx][2]
				};
				results.clear();
				const uint32_t nMatches = tree.radiusSearch(&query_pt[0], radius2, results, searchParams);

				for (uint32_t i = 0; i < nMatches; ++i)
				{
					uint32_t nbr = results[i].first;
					if (visited[nbr])
						continue;

					visited[nbr] = true;
					q.push(nbr);
				}
			}

			if (cluster.empty())
				continue; // no point found

			clusters.push_back(std::move(cluster));
		}

		return clusters;
	}

	std::vector<pmp::Normal> EstimatePointCloudNormalsVCG(const std::vector<pmp::Point>& points, const size_t& fittingAdjNum, const size_t& nSmoothingIters, const pmp::Point& viewPoint, const bool& useViewPoint)
	{
		if (points.empty())
		{
			std::cerr << "Geometry::EstimatePointCloudNormalsVCG: No points to triangulate.\n";
			return {};
		}

		// prepare data
		BaseMeshGeometryData resultData;
		resultData.Vertices = points;
		VCG_Mesh ptsMesh;
		FillVCGMeshWithPoints(ptsMesh, points);

		vcg::tri::PointCloudNormal<VCG_Mesh>::Param p;
		p.fittingAdjNum = fittingAdjNum;
		p.smoothingIterNum = nSmoothingIters;
		p.viewPoint = vcg::Point3f(viewPoint[0], viewPoint[1], viewPoint[2]);
		p.useViewPoint = useViewPoint;
		vcg::tri::PointCloudNormal<VCG_Mesh>::Compute(ptsMesh, p);

		resultData.VertexNormals = ExtractVertexNormalsFromVCGMesh(ptsMesh, -1.0 /* for some reason they're flipped */);
		return resultData.VertexNormals;
	}

} // namespace Geometry