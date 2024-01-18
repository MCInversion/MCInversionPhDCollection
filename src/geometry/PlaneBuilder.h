#pragma once

#include "PrimitiveMeshBuilder.h"

namespace Geometry
{
	/**
	 * \brief Settings for plane construction.
	 * \struct PlaneSettings
	*/
	struct PlaneSettings
	{
		pmp::vec3 Origin{}; //! reference point of this plane
		float Width{ 1.0f }; //! dimension in the x-direction
		float Depth{ 1.0f }; //! dimension in the y-direction
		unsigned int nWidthSegments{ 5 }; //! number of segments in the x-direction
		unsigned int nDepthSegments{ 5 }; //! number of segments in the y-direction
		bool UseQuads{ false }; //! if true the torus surface will be composed of quads instead of triangles.
		bool ComputeNormals{ false }; //! whether to compute vertex normals
	};

	/**
	 * \brief A builder object for generating parametric plane meshes (composed of (triangulated) quads).
	 * \class PlaneBuilder
	 */
	class PlaneBuilder : public PrimitiveMeshBuilder
	{
	public:
		/// \brief constructor
		explicit PlaneBuilder(const PlaneSettings& settings)
			: m_Settings(settings)
		{
			if (m_Settings.nWidthSegments < 1)
			{
				std::cerr << "PlaneBuilder::PlaneBuilder: m_Settings.nWidthSegments < 1, setting m_Settings.nWidthSegments = 1!\n";
				m_Settings.nWidthSegments = 1;
			}
			if (m_Settings.nDepthSegments < 1)
			{
				std::cerr << "PlaneBuilder::PlaneBuilder: m_Settings.nDepthSegments < 4, setting m_Settings.nDepthSegments = 1!\n";
				m_Settings.nDepthSegments = 1;
			}
			if (m_Settings.Width <= 0.0f)
			{
				std::cerr << "PlaneBuilder::PlaneBuilder: m_Settings.Width <= 0.0, setting m_Settings.Width = 1.0!\n";
				m_Settings.Width = 1.0f;
			}
			if (m_Settings.Depth <= 0.0f)
			{
				std::cerr << "PlaneBuilder::PlaneBuilder: m_Settings.Depth <= 0.0, setting m_Settings.Depth = 1.0!\n";
				m_Settings.Depth = 1.0f;
			}
		}

		/// \brief builds BaseMeshGeometryData for an ico-sphere with given settings.
		void BuildBaseData() override;

	private:

		PlaneSettings m_Settings{};
	};

} // namespace Geometry