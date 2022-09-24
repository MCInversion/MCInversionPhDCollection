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

} // namespace Geometry