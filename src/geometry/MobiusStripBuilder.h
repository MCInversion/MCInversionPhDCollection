#pragma once

#include "PrimitiveMeshBuilder.h"

namespace Geometry
{
	/**
	 * \brief Settings for the construction of a Mobius strip using the sweep-line approach.
	 * \struct MobiusStripSettings
	*/
	struct MobiusStripSettings
	{
		float RingRadius{ 1.0f }; //! radius of the Mobius ring.
		float SweepLineLength{ 0.3f }; //! length of the sweep line of the Mobius strip.
		unsigned int RingSegments{ 10 }; //! number of segments per Mobius ring.
		unsigned int SweepSegments{ 5 }; //! number of segments along the swept line.
		bool UseQuads{ false }; //! if true the torus surface will be composed of quads instead of triangles.
		bool ComputeNormals{ false }; //! whether to compute vertex normals
	};

	/**
	 * \brief A builder object for generating parametric Mobius strip meshes.
	 * \class MobiusStripBuilder
	 */
	class MobiusStripBuilder : public PrimitiveMeshBuilder
	{
	public:
		/// \brief constructor
		explicit MobiusStripBuilder(const MobiusStripSettings& settings)
			: m_Settings(settings)
		{
			if (m_Settings.SweepSegments < 1)
			{
				std::cerr << "MobiusStripBuilder::MobiusStripBuilder: m_Settings.SweepSegments < 1, setting m_Settings.SweepSegments = 1!\n";
				m_Settings.SweepSegments = 1;
			}
			if (m_Settings.RingSegments < 4)
			{
				std::cerr << "MobiusStripBuilder::MobiusStripBuilder: m_Settings.RingSegments < 4, setting m_Settings.RingSegments = 3!\n";
				m_Settings.RingSegments = 3;
			}
			if (m_Settings.RingRadius <= 0.0f)
			{
				std::cerr << "MobiusStripBuilder::MobiusStripBuilder: m_Settings.RingRadius <= 0.0, setting m_Settings.RingRadius = 1.0!\n";
				m_Settings.RingRadius = 1.0f;
			}
			if (m_Settings.SweepLineLength <= 0.0f)
			{
				std::cerr << "MobiusStripBuilder::MobiusStripBuilder: m_Settings.SweepLineLength <= 0.0, setting m_Settings.SweepLineLength = 0.3!\n";
				m_Settings.SweepLineLength = 0.3f;
			}
		}

		/// \brief builds BaseMeshGeometryData for a Mobius strip with given settings.
		void BuildBaseData() override;

	private:

		MobiusStripSettings m_Settings{};
	};

} // namespace Geometry