#include "PlaneBuilder.h"

namespace Geometry
{
	void PlaneBuilder::BuildBaseData()
	{
		m_BaseResult = std::make_unique<BaseMeshGeometryData>();

		const pmp::Scalar width = m_Settings.Width;
		const pmp::Scalar depth = m_Settings.Depth;
		const auto orig = m_Settings.Origin;

		const size_t nXSegments = m_Settings.nWidthSegments;
		const size_t nYSegments = m_Settings.nDepthSegments;
		const size_t nQuads = nXSegments * nYSegments;
		const size_t nXVerts = nXSegments + 1;
		const size_t nYVerts = nYSegments + 1;
		const size_t nVertices = nXVerts * nYVerts;

		m_BaseResult->Vertices.reserve(nQuads);
		const size_t nPolys = (m_Settings.UseQuads ? 1 : 2) * nQuads;
		m_BaseResult->PolyIndices.reserve(nPolys);
		if (m_Settings.ComputeNormals)
			m_BaseResult->VertexNormals.reserve(nQuads);

		const InsertionFunction polyInsertFn = m_Settings.UseQuads ? EmplaceQuadIndexTuple : EmplaceTriIndexTuples;

		for (unsigned int i = 0; i < nYVerts; i++)
		{
			const pmp::Scalar yParam = orig[1] + static_cast<pmp::Scalar>(i) / static_cast<pmp::Scalar>(nYSegments) * depth;
			for (unsigned int j = 0; j < nXVerts; j++)
			{
				const pmp::Scalar xParam = orig[0] + static_cast<pmp::Scalar>(j) / static_cast<pmp::Scalar>(nXSegments) * width;
				m_BaseResult->Vertices.emplace_back(pmp::vec3{ xParam, yParam, orig[2] });

				if (i > 0 && j > 0)
				{
					polyInsertFn(*m_BaseResult,
						(i - 1) * static_cast<unsigned int>(nXVerts) + j - 1,
						(i - 1) * static_cast<unsigned int>(nXVerts) + j,
						i * static_cast<unsigned int>(nXVerts) + j,
						i * static_cast<unsigned int>(nXVerts) + j - 1
					);
				}
			}
		}

		if (m_Settings.ComputeNormals)
		{
			const auto unitNormal = pmp::vec3{ 0.0, 0.0, 1.0 };
			m_BaseResult->VertexNormals = std::vector(nVertices, unitNormal);
		}
	}
} // namespace Geometry
