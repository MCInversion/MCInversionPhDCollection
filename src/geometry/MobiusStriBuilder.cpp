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
	[[nodiscard]] pmp::vec3 ComputeMobiusNormal(const pmp::Scalar& u, const pmp::Scalar& v, const pmp::Scalar& r, const pmp::Scalar& l)
	{
		const auto norm = static_cast<pmp::Scalar>(sqrt(l * l * (8 * r * r + 3 * l * 2 * v * 2 + l * v * (8 * r * cos(u) + l * v * cos(2 * u)))));
		constexpr auto sqrt2 = static_cast<pmp::Scalar>(M_SQRT2);
		return pmp::vec3{
			static_cast<pmp::Scalar>(sqrt2 * l * sin(u) * (2 * r * cos(u) - l * v * pow(sin(u), 2.0))) / norm,
			static_cast<pmp::Scalar>(sqrt2 * l * (l * v * pow(cos(u), 3.0) + 2 * (r + l * v * cos(u)) * pow(sin(u), 2.0))) / norm,
			(pmp::Scalar)-1.0 * static_cast<pmp::Scalar>(l * sqrt2 * cos(u) * (2 * r + l * v * cos(u))) / norm
		} * (-1.0);
	}

	void MobiusStripBuilder::BuildBaseData()
	{
		m_BaseResult = std::make_unique<BaseMeshGeometryData>();

		const pmp::Scalar ringRadius = m_Settings.RingRadius;
		const pmp::Scalar sLineLength = m_Settings.SweepLineLength;

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
			const pmp::Scalar ringParam = 2.0 * static_cast<pmp::Scalar>(M_PI * i) / static_cast<pmp::Scalar>(nRSegments);
			for (unsigned int j = 0; j < nLVerts; j++)
			{
				const pmp::Scalar lineParam = sLineLength * (static_cast<pmp::Scalar>(j) / static_cast<pmp::Scalar>(nLVerts) - 0.5);
				const pmp::Scalar rad = (ringRadius + 0.5 * sLineLength * lineParam * cos(ringParam));
				m_BaseResult->Vertices.emplace_back(pmp::vec3{
					rad * cos(ringParam),
					rad * sin(ringParam),
					(pmp::Scalar)0.5 * sLineLength * lineParam * sin(ringParam)
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