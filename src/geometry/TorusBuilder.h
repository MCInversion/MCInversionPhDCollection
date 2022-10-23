#pragma once

#include "PrimitiveMeshBuilder.h"

namespace Geometry
{
	/**
	 * \brief Settings for torus construction.
	 * \struct TorusSettings
	*/
	struct TorusSettings
	{
		float RingRadius{ 1.0f }; //! radius of the torus ring.
		float TubeRadius{ 0.3f }; //! radius of the tube of the torus.
		unsigned int RingSegments{ 10 }; //! number of segments per torus ring.
		unsigned int TubeSegments{ 5 }; //! number of segments along torus tube.
		bool UseQuads{ false }; //! if true the torus surface will be composed of quads instead of triangles.
		bool ComputeNormals{ false }; //! whether to compute vertex normals
	};

	/**
	 * \brief A builder object for generating parametric torus meshes.
	 * \class TorusBuilder
	 */
	class TorusBuilder : public PrimitiveMeshBuilder
	{
	public:
		/// \brief constructor
		explicit TorusBuilder(const TorusSettings& settings)
			: m_Settings(settings)
		{
			if (m_Settings.TubeSegments < 4)
			{
				std::cerr << "TorusBuilder::TorusBuilder: m_Settings.TubeSegments < 4, setting m_Settings.TubeSegments = 3!\n";
				m_Settings.TubeSegments = 3;
			}
			if (m_Settings.RingSegments < 4)
			{
				std::cerr << "TorusBuilder::TorusBuilder: m_Settings.RingSegments < 4, setting m_Settings.RingSegments = 3!\n";
				m_Settings.RingSegments = 3;
			}
			if (m_Settings.RingRadius <= 0.0f)
			{
				std::cerr << "TorusBuilder::TorusBuilder: m_Settings.RingRadius <= 0.0, setting m_Settings.RingRadius = 1.0!\n";
				m_Settings.RingRadius = 1.0f;
			}
			if (m_Settings.TubeRadius <= 0.0f)
			{
				std::cerr << "TorusBuilder::TorusBuilder: m_Settings.TubeRadius <= 0.0, setting m_Settings.TubeRadius = 0.3!\n";
				m_Settings.TubeRadius = 0.3f;
			}
		}

		/// \brief builds BaseMeshGeometryData for an torus with given settings.
		void BuildBaseData() override;

	private:

		TorusSettings m_Settings{};
	};

} // namespace Geometry