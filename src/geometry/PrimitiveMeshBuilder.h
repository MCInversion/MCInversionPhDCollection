#pragma once

#include "pmp/SurfaceMesh.h"

#include "GeometryConversionUtils.h"

namespace Geometry
{
	/**
	 * \brief A base builder object for generating parametric meshes.
	 * \class PrimitiveMeshBuilder
	 *
	 * NOTE: consists of two stages:
	 * (1) BaseMeshGeometryData, and
	 * (2) pmp::SurfaceMesh to make use of existing implementation for BaseMeshGeometryData.
	 */
	class PrimitiveMeshBuilder
	{
	public:
		/// \brief Destructor.
		virtual ~PrimitiveMeshBuilder() = default;

		/**
		 * \brief Build BaseMeshGeometryData in *m_BaseResult. Overriden by specific primitive builders.
		 */
		virtual void BuildBaseData() = 0;

		/**
		 * \brief Converts *m_BufferResult to pmp::SurfaceMesh.
		 * \throw std::logic_error if m_BaseResult == nullptr.
		 */
		void BuildPMPSurfaceMesh()
		{
			if (!m_BaseResult)
			{
				throw std::logic_error("PrimitiveMeshBuilder::BuildPMPSurfaceMesh: cannot build pmp::SurfaceMesh without prior PrimitiveMeshBuilder::BuildBaseData()!\n");
			}
			m_Result = std::make_unique<pmp::SurfaceMesh>(ConvertBufferGeomToPMPSurfaceMesh(*m_BaseResult));
		}

		/// \brief base mesh geometry data getter.
		[[nodiscard]] BaseMeshGeometryData&& GetBaseResult() const
		{
			return std::move(*m_BaseResult);
		}

		/// \brief pmp::SurfaceMesh getter.
		[[nodiscard]] pmp::SurfaceMesh&& GetPMPSurfaceMeshResult() const
		{
			return std::move(*m_Result);
		}

	protected:
		std::unique_ptr<BaseMeshGeometryData> m_BaseResult{ nullptr };
		std::unique_ptr<pmp::SurfaceMesh> m_Result{ nullptr };
	};

	/// \brief insertion function for quad vertex indices
	inline void EmplaceQuadIndexTuple(BaseMeshGeometryData& baseData, const unsigned int& i0, const unsigned int& i1, const unsigned int& i2, const unsigned int& i3)
	{
		const std::vector polyIdTuple{ i0, i1, i2, i3 };
		baseData.PolyIndices.push_back(polyIdTuple);
	}

	/// \brief insertion function for triangle pair vertex indices
	inline void EmplaceTriIndexTuples(BaseMeshGeometryData& baseData, const unsigned int& i0, const unsigned int& i1, const unsigned int& i2, const unsigned int& i3)
	{
		const std::vector polyIdTuple0{ i0, i1, i2 };
		const std::vector polyIdTuple1{ i0, i2, i3 };
		baseData.PolyIndices.push_back(polyIdTuple0);
		baseData.PolyIndices.push_back(polyIdTuple1);
	}

	/// \brief an insertion function for vertex indices.
	using InsertionFunction = std::function<void(BaseMeshGeometryData&, const unsigned int&, const unsigned int&, const unsigned int&, const unsigned int&)>;

} // namespace Geometry