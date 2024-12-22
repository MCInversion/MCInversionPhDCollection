#include "TorusBuilder.h"

namespace Geometry
{
	void TorusBuilder::BuildBaseData()
	{
		m_BaseResult = std::make_unique<BaseMeshGeometryData>();

		const pmp::Scalar ringRadius = m_Settings.RingRadius;
		const pmp::Scalar tubeRadius = m_Settings.TubeRadius;

		const size_t nRSegments = m_Settings.RingSegments;
		const size_t nTSegments = m_Settings.TubeSegments;
		const size_t nQuads = nRSegments * nTSegments;

		m_BaseResult->Vertices.reserve(nQuads);
		const size_t nPolys = (m_Settings.UseQuads ? 1 : 2) * nQuads;
		m_BaseResult->PolyIndices.reserve(nPolys);
		if (m_Settings.ComputeNormals)
			m_BaseResult->VertexNormals.reserve(nQuads);

		const InsertionFunction polyInsertFn = m_Settings.UseQuads ? EmplaceQuadIndexTuple : EmplaceTriIndexTuples;

		//pmp::Scalar prevRingParam{ 0.0 };
		//pmp::Scalar prevTubeParam{ 0.0 };
		for (unsigned int i = 0; i < nRSegments; i++)
		{
			const pmp::Scalar ringParam = 2.0 * static_cast<pmp::Scalar>(M_PI * i) / static_cast<pmp::Scalar>(nRSegments);
			for (unsigned int j = 0; j < nTSegments; j++)
			{
				const pmp::Scalar tubeParam = 2.0 * static_cast<pmp::Scalar>(M_PI * j) / static_cast<pmp::Scalar>(nTSegments);

				m_BaseResult->Vertices.emplace_back(pmp::vec3{
					(ringRadius + tubeRadius * cos(tubeParam)) * cos(ringParam),
					(ringRadius + tubeRadius * cos(tubeParam)) * sin(ringParam),
					-tubeRadius * sin(tubeParam)
				});

				if (i > 0 && j > 0)
				{
					polyInsertFn(*m_BaseResult,
						(i - 1) * static_cast<unsigned int>(nTSegments) + j - 1,
						(i - 1) * static_cast<unsigned int>(nTSegments) + j,
						i * static_cast<unsigned int>(nTSegments) + j,
						i * static_cast<unsigned int>(nTSegments) + j - 1
					);
				}
				if (m_Settings.ComputeNormals)
				{
					m_BaseResult->VertexNormals.emplace_back(
						pmp::vec3{ cos(ringParam) * cos(tubeParam), sin(ringParam) * cos(tubeParam), -sin(tubeParam) });
				}
			}
			if (i > 0)
			{
				polyInsertFn(*m_BaseResult,
					(i - 1) * static_cast<unsigned int>(nTSegments) + static_cast<unsigned int>(nTSegments) - 1,
					(i - 1) * static_cast<unsigned int>(nTSegments),
					i * static_cast<unsigned int>(nTSegments),
					i * static_cast<unsigned int>(nTSegments) + static_cast<unsigned int>(nTSegments) - 1
				);
			}
		}
		for (unsigned int j = 1; j < nTSegments; j++)
		{
			//const pmp::Scalar tubeParam = M_PI * (1.0 + 2.0 * static_cast<pmp::Scalar>(j) / static_cast<pmp::Scalar>(nTSegments));
			polyInsertFn(*m_BaseResult,
				(static_cast<unsigned int>(nRSegments) - 1) * static_cast<unsigned int>(nTSegments) + j - 1,
				(static_cast<unsigned int>(nRSegments) - 1) * static_cast<unsigned int>(nTSegments) + j,
				j,
				j - 1
			);
		}
		polyInsertFn(*m_BaseResult,
			(static_cast<unsigned int>(nRSegments) - 1) * static_cast<unsigned int>(nTSegments) + static_cast<unsigned int>(nTSegments) - 1,
			(static_cast<unsigned int>(nRSegments) - 1) * static_cast<unsigned int>(nTSegments),
			0,
			static_cast<unsigned int>(nTSegments) - 1
		);
	}

} // namespace Geometry