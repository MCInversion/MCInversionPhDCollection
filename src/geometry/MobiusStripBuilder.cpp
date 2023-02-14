#include "MobiusStripBuilder.h"

namespace Geometry
{
	/**
	 * \brief A function to compute normal vector from the parametric realization of the Mobius strip.
	 * \param u     ring parameter.
	 * \param v     line parameter.
	 * \param r     ring radius.
	 * \param l     line length.
	 * \return normal vector for the given parameters (u, v, r, l)
	 */
	[[nodiscard]] pmp::vec3 ComputeMobiusNormal(const float& u, const float& v, const float& r, const float& l)
	{
		const float norm = sqrt(l * l * (16 * r * r + 3 * l * l * v * v + 2 * l * v * (8 * r * cos(u / 2) + l * v * cos(u))));
		return pmp::vec3{
			(l * sin(u / 2) * (4 * r * cos(u) - 2 * l * v * sin(u / 2) * sin(u))) / norm,
			(l * (4 * r * sin(u / 2) * sin(u) + l * v * (cos(u) + sin(u) * sin(u)))) / norm,
			(l * (4 * r * cos(u / 2) + l * v * (1 + cos(u)))) / norm
		};
	}

	void MobiusStripBuilder::BuildBaseData()
	{
		m_BaseResult = std::make_unique<BaseMeshGeometryData>();

		const float ringRadius = m_Settings.RingRadius;
		const float sLineLength = m_Settings.SweepLineLength;

		const size_t nRSegments = m_Settings.RingSegments;
		const size_t nLSegments = m_Settings.SweepSegments;
		const size_t nQuads = nRSegments * nLSegments;

		m_BaseResult->Vertices.reserve(nQuads);
		const size_t nPolys = (m_Settings.UseQuads ? 1 : 2) * nQuads;
		m_BaseResult->PolyIndices.reserve(nPolys);
		if (m_Settings.ComputeNormals)
			m_BaseResult->VertexNormals.reserve(nQuads);

		const InsertionFunction polyInsertFn = m_Settings.UseQuads ? EmplaceQuadIndexTuple : EmplaceTriIndexTuples;

		for (unsigned int i = 0; i < nRSegments; i++)
		{
			const float ringParam = 2.0f * static_cast<float>(M_PI * i) / static_cast<float>(nRSegments);
			for (unsigned int j = 0; j < nLSegments; j++)
			{
				const float lineParam = (2.0f * static_cast<float>(j) / static_cast<float>(nLSegments) - 1.0f);
				const float rad = (ringRadius + 0.5f * sLineLength * lineParam * cos(0.5f * ringParam));
				m_BaseResult->Vertices.emplace_back(pmp::vec3{
					rad * cos(ringParam),
					rad * sin(ringParam),
					0.5f * sLineLength * lineParam * sin(0.5f * ringParam)
				});

				if (i > 0 && j > 0)
				{
					polyInsertFn(*m_BaseResult,
						(i - 1) * static_cast<unsigned int>(nLSegments) + j - 1,
						(i - 1) * static_cast<unsigned int>(nLSegments) + j,
						i * static_cast<unsigned int>(nLSegments) + j,
						i * static_cast<unsigned int>(nLSegments) + j - 1
					);
				}
				if (m_Settings.ComputeNormals)
				{
					m_BaseResult->VertexNormals.emplace_back(
						ComputeMobiusNormal(ringParam, lineParam, ringRadius, sLineLength));
				}
			}
		}
		/*for (unsigned int j = 1; j < nLSegments; j++)
		{
			polyInsertFn(*m_BaseResult,
				(static_cast<unsigned int>(nRSegments) - 1) * static_cast<unsigned int>(nLSegments) + j - 1,
				(static_cast<unsigned int>(nRSegments) - 1) * static_cast<unsigned int>(nLSegments) + j,
				j,
				j - 1
			);
		}
		polyInsertFn(*m_BaseResult,
			(static_cast<unsigned int>(nRSegments) - 1) * static_cast<unsigned int>(nLSegments) + static_cast<unsigned int>(nLSegments) - 1,
			(static_cast<unsigned int>(nRSegments) - 1) * static_cast<unsigned int>(nLSegments),
			0,
			static_cast<unsigned int>(nLSegments) - 1
		);*/
	}

} // namespace Geometry