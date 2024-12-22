#pragma once

#include "PrimitiveMeshBuilder.h"

namespace Geometry
{
/**
 * \brief Settings for the construction of the multi torus
 * \struct MultiTorusSettings
 */
struct MultiTorusSettings
{
	unsigned int NumberOfLoops{ 2 }; //! the number of loops which the multi-torus should have
	pmp::BoundingBox BoundingBox{}; //! bounding box where the polygonized scalar field should be generated
	pmp::Scalar RingRadius{ 1.0 }; //! radius of each loop of the multi-torus
	pmp::Scalar TubeRadius{ 0.3 }; //! radius of each tube of each loop of the multi-torus.
	unsigned int RingSegments{ 10 }; //! number of segments per torus ring.
	unsigned int TubeSegments{ 5 }; //! number of segments along torus tube.
	bool UseQuads{ false }; //! if true the torus surface will be composed of quads instead of triangles.
	bool ComputeNormals{ false }; //! whether to compute vertex normals
};

/**
 * \brief A builder object for generating parametric multi-torus meshes.
 * \class MultiTorusBuilder
 */
class MultiTorusBuilder : public PrimitiveMeshBuilder
{
public:
	/// \brief constructor
	explicit MultiTorusBuilder(const MultiTorusSettings& settings)
		: m_Settings(settings)
	{
		if (m_Settings.NumberOfLoops < 1)
		{
			std::cerr << "MultiTorusBuilder::MultiTorusBuilder: m_Settings.NumberOfLoops < 1, setting m_Settings.NumberOfLoops = 1!\n";
			m_Settings.NumberOfLoops = 1;
		}
		if (m_Settings.BoundingBox.is_empty())
		{
			std::cerr << "MultiTorusBuilder::MultiTorusBuilder: m_Settings.BoundingBox.is_empty(), setting m_Settings.BoundingBox to {{}, {1, 1, 1}}!\n";
			m_Settings.BoundingBox = pmp::BoundingBox{ pmp::vec3{}, pmp::vec3{1.0, 1.0, 1.0} };
		}
		if (m_Settings.RingRadius <= 0.0)
		{
			std::cerr << "MultiTorusBuilder::MultiTorusBuilder: m_Settings.RingRadius <= 0.0, setting m_Settings.RingRadius = 1.0!\n";
			m_Settings.RingRadius = 1.0;
		}
		if (m_Settings.TubeRadius <= 0.0)
		{
			std::cerr << "MultiTorusBuilder::MultiTorusBuilder: m_Settings.TubeRadius <= 0.0, setting m_Settings.TubeRadius = 0.3!\n";
			m_Settings.TubeRadius = 0.3;
		}
	}

	/// \brief builds BaseMeshGeometryData for an torus with given settings.
	void BuildBaseData() override;

private:

	MultiTorusSettings m_Settings{};
};

} // namespace Geometry