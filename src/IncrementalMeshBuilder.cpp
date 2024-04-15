#include "IncrementalMeshBuilder.h"

// #include <nanoflann.hpp> // TODO: use

#include "utils/FileMappingWrapper.h"

#include <random>
#include <thread>
#include <iostream>

namespace
{
	void TriangulatePointsWithBallPivoting(std::vector<unsigned int>& resultTris, const std::vector<pmp::Point>& vertices, const std::vector<pmp::Normal>& /*normals*/)
	{
		std::cerr << "TriangulatePointsWithBallPivoting Not implemented\n";
	}

	void TriangulatePointsWithPoisson(std::vector<unsigned int>& resultTris, const std::vector<pmp::Point>& vertices, const std::vector<pmp::Normal>& normals)
	{
		std::cerr << "TriangulatePointsWithPoisson Not implemented\n";
	}

	void TriangulatePointsWithMarchingCubes(std::vector<unsigned int>& resultTris, const std::vector<pmp::Point>& vertices, const std::vector<pmp::Normal>& /*normals*/)
	{
		std::cerr << "TriangulatePointsWithMarchingCubes Not implemented\n";
	}

	void TriangulatePointsWithLagrangianShrinkWrapping(std::vector<unsigned int>& resultTris, const std::vector<pmp::Point>& vertices, const std::vector<pmp::Normal>& /*normals*/)
	{
		std::cerr << "TriangulatePointsWithLagrangianShrinkWrapping Not implemented\n";
	}

	// ========================================================================================================

	void SampleVerticesSequentially(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& /* seed */)
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
			if (strncmp(cursor, "v ", 2) == 0)
			{
				const bool isNormal = strncmp(cursor, "vn ", 3) == 0;
				cursor += isNormal ? 3 : 2; // skip "v "

				pmp::vec3 vec;
				char* tempCursor;
				vec[0] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				vec[1] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				vec[2] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				result.push_back(vec);
			}
			else
			{
				// Skip to the next line if the current line isn't recognized
				while (*cursor != '\n' && cursor < end)
					cursor++;
			}
		}
	}

	void SampleIndicesUniformly(const char* start, const char* end, std::vector<size_t>& resultIndices, const std::optional<unsigned int>& seed)
	{
		resultIndices.clear();
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
		size_t nLines = std::count(start, end, '\n');
		std::uniform_int_distribution<> distrib(0, nLines - 1);
		resultIndices.reserve(nLines);
		for (size_t i = 0; i < nLines; i++)
		{
			resultIndices.push_back(distrib(gen));
		}
	}

	void SampleIndicesNormalRandomly(const char* start, const char* end, std::vector<size_t>& resultIndices, const std::optional<unsigned int>& seed)
	{
		resultIndices.clear();
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
		size_t nLines = std::count(start, end, '\n');
		std::normal_distribution<> distrib(nLines / 2, nLines / 6);
		resultIndices.reserve(nLines);
		for (size_t i = 0; i < nLines; i++)
		{
			resultIndices.push_back(distrib(gen));
		}
	}

	void SampleVerticesFromIndices(const char* start, const char* end, const std::vector<size_t>& indices, std::vector<pmp::Point>& result)
	{
		const char* cursor = start + indices[0];
		for (size_t i = 0; i < indices.size(); i++)
		{
			// skip if indices not in range
			if (indices[i] < 0 || indices[i] >= end - start)
				continue;

			// If the current line is incomplete, skip to the next line
			if (*cursor == '\n')
			{
				cursor = start + indices[(i + 1) % indices.size()];
				continue;
			}

			// If it's a vertex, parse the three floats
			if (strncmp(cursor, "v ", 2) == 0)
			{
				const bool isNormal = strncmp(cursor, "vn ", 3) == 0;
				cursor += isNormal ? 3 : 2; // skip "v "

				pmp::vec3 vec;
				char* tempCursor;
				vec[0] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				vec[1] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				vec[2] = std::strtof(cursor, &tempCursor);
				cursor = tempCursor;

				result.push_back(vec);
			}
			else
			{
				// Skip to the next line if the current line isn't recognized
				while (*cursor != '\n' && cursor < end)
					cursor = start + indices[(i + 1) % indices.size()];
			}
		}
	}	

	void SampleVerticesUniformRandomly(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed)
	{
		std::vector<size_t> indices;
		SampleIndicesUniformly(start, end, indices, seed);
		SampleVerticesFromIndices(start, end, indices, result);
	}

	void SampleVerticesNormalRandomly(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed)
	{
		std::vector<size_t> indices;
		SampleIndicesNormalRandomly(start, end, indices, seed);
		SampleVerticesFromIndices(start, end, indices, result);
	}

	void SampleVerticesWithSoftmaxFeatureDectection(const char* start, const char* end, std::vector<pmp::Point>& result, const std::optional<unsigned int>& seed)
	{
		throw std::runtime_error("SampleVerticesWithSoftmaxFeatureDectection Not implemented\n");
	}

} // anonymous namespace

using namespace IMB;

// --------------------------------------------------------------------------------------------------------

[[nodiscard]] ReconstructionFunction GetReconstructionFunction(const ReconstructionFunctionType& reconstructType)
{
	if (reconstructType == ReconstructionFunctionType::BallPivoting)
		return TriangulatePointsWithBallPivoting;
	else if (reconstructType == ReconstructionFunctionType::Poisson)
		return TriangulatePointsWithPoisson;
	else if (reconstructType == ReconstructionFunctionType::MarchingCubes)
		return TriangulatePointsWithMarchingCubes;
	TriangulatePointsWithLagrangianShrinkWrapping;
}

[[nodiscard]] VertexSelectionFunction GetVertexSelectionFunction(const VertexSelectionType& vertSelType)
{
	if (vertSelType == VertexSelectionType::Sequential)
		return SampleVerticesSequentially;
	else if (vertSelType == VertexSelectionType::UniformRandom)
		return SampleVerticesUniformRandomly;
	else if (vertSelType == VertexSelectionType::NormalRandom)
		return SampleVerticesNormalRandomly;
	SampleVerticesWithSoftmaxFeatureDectection;
}

// --------------------------------------------------------------------------------------------------------


//
// --------------------------------------------------------------------------------------------------------

void IMB::IncrementalMeshBuilder::Init(
	const std::string& fileName, const unsigned int& completionFrequency, 
	const ReconstructionFunctionType& reconstructType, 
	const VertexSelectionType& vertSelType)
{
	m_ReconstructPtCloud = GetReconstructionFunction(reconstructType);
	m_SelectVertices = GetVertexSelectionFunction(vertSelType);
}

void IncrementalMeshBuilder::Triangulate()
{
}

void IncrementalMeshBuilder::SampleVertices()
{
	//m_MeshData.Vertices = m_SelectVertices();
}

//void IncrementalMeshBuilder::SampleVerticesAndNormals()
//{
//}
