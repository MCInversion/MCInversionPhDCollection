#pragma once

#include "PrimitiveMeshBuilder.h"

/**
 * \brief Constants used for estimating vertex counts and edge lengths of ico-sphere meshes.
 */
constexpr unsigned int N_ICO_VERTS_0 = 12; // number of vertices in an icosahedron.
constexpr unsigned int N_ICO_EDGES_0 = 30; // number of edges in an icosahedron.

namespace Geometry
{
	/**
	 * \brief Settings for ico-sphere construction.
	 * \struct IcoSphereSettings
	*/
	struct IcoSphereSettings
	{
		unsigned int SubdivisionLevel{ 0 }; //! subdivision level. Zero corresponds to a basic icosahedron.
		float Radius{ 1.0f }; //! ico-sphere radius
		bool ComputeNormals{ false }; //! whether to compute vertex normals
	};

	/**
	 * \brief A builder object for generating parametric ico-sphere meshes (a 4:1-subdivided icosahedron).
	 * \class IcoSphereBuilder
	 */
	class IcoSphereBuilder : public PrimitiveMeshBuilder
	{
	public:
		/// \brief constructor
		explicit IcoSphereBuilder(const IcoSphereSettings& settings)
			: m_SubdivisionLevel(settings.SubdivisionLevel), m_Radius(settings.Radius), m_ComputeNormals(settings.ComputeNormals)
		{ }

		/// \brief builds BaseMeshGeometryData for an ico-sphere with given settings.
		void BuildBaseData() override;

	private:

		unsigned int m_SubdivisionLevel{ 0 }; //>! subdivision level. Zero corresponds to a basic icosahedron.
		float m_Radius{ 1.0f }; //>! ico-sphere radius
		bool m_ComputeNormals{ false };
	};

} // namespace Geometry