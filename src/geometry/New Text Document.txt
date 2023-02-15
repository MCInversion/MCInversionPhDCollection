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
		const float norm = sqrt(l * l * (8 * r * r + 3 * l * 2 * v * 2 + l * v * (8 * r * cos(u) + l * v * cos(2 * u))));
		constexpr float sqrt2 = static_cast<float>(M_SQRT2);
		return pmp::vec3{
			(sqrt2 * l * sin(u) * (2 * r * cos(u) - l * v * pow(sin(u), 2.0f))) / norm,
			(sqrt2 * l * (l * v * pow(cos(u), 3.0f) + 2 * (r + l * v * cos(u)) * pow(sin(u), 2.0f))) / norm,
			-1.0f * (l * sqrt2 * cos(u) * (2 * r + l * v * cos(u))) / norm
		} * (-1.0f);
	}

	void MobiusStripBuilder::BuildBaseData()
	{
		m_BaseResult = std::make_unique<BaseMeshGeometryData>();

		const float ringRadius = m_Settings.RingRadius;
		const float sLineLength = m_Settings.SweepLineLength;

		const size_t nRSegments = m_Settings.RingSegments;
		const size_t nLSegments = m_Settings.SweepSegments;
		const size_t nLVerts = nLSegments + 1;
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
			for (unsigned int j = 0; j < nLVerts; j++)
			{
				const float lineParam = sLineLength * (static_cast<float>(j) / static_cast<float>(nLVerts) - 0.5f);
				const float rad = (ringRadius + 0.5f * sLineLength * lineParam * cos(ringParam));
				m_BaseResult->Vertices.emplace_back(pmp::vec3{
					rad * cos(ringParam),
					rad * sin(ringParam),
					0.5f * sLineLength * lineParam * sin(ringParam)
				});

				if (i > 0 && j > 0)
				{
					polyInsertFn(*m_BaseResult,
						(i - 1) * static_cast<unsigned int>(nLVerts) + j - 1,
						(i - 1) * static_cast<unsigned int>(nLVerts) + j,
						i * static_cast<unsigned int>(nLVerts) + j,
						i * static_cast<unsigned int>(nLVerts) + j - 1
					);
				}
				if (m_Settings.ComputeNormals)
				{
					m_BaseResult->VertexNormals.emplace_back(
						ComputeMobiusNormal(ringParam, lineParam, ringRadius, sLineLength));
				}
			}
		}
		for (unsigned int j = 1; j < nLVerts; j++)
		{
			polyInsertFn(*m_BaseResult,
				(static_cast<unsigned int>(nRSegments) - 1) * static_cast<unsigned int>(nLVerts) + j - 1,
				(static_cast<unsigned int>(nRSegments) - 1) * static_cast<unsigned int>(nLVerts) + j,
				j,
				j - 1
			);
		}
	}

} // namespace Geometry