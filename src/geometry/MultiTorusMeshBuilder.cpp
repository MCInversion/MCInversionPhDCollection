#include "MultiTorusMeshBuilder.h"

namespace Geometry
{
    void MultiTorusBuilder::BuildBaseData()
    {
        m_BaseResult = std::make_unique<BaseMeshGeometryData>();

        const auto& bbox = m_Settings.BoundingBox;

        const auto center = bbox.center();
        const auto& boxMax = bbox.max();
        const auto& boxMin = bbox.min();
        const float maxXYRadius = 0.5 * sqrt((boxMax[0] - boxMin[0]) * (boxMax[0] - boxMin[0]) + (boxMax[1] - boxMin[1]) * (boxMax[1] - boxMin[1]));

        const float ringRadius = m_Settings.RingRadius;
        const float tubeRadius = m_Settings.TubeRadius;

        const size_t nLoops = m_Settings.NumberOfLoops;
        const size_t nRSegments = m_Settings.RingSegments;
        const size_t nTSegments = m_Settings.TubeSegments;

		// TODO
    }

} // namespace Geometry
