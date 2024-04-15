#pragma once
#include "geometry/GeometryConversionUtils.h"

namespace IMB 
{
    /// \brief A functor for rendering a mesh.
    using MeshRenderFunction = std::function<void(const Geometry::BaseMeshGeometryData&)>;

    class MeshUpdateHandler 
    {
    public:
        explicit MeshUpdateHandler(const MeshRenderFunction& renderFunc);
    };
} // namespace IMB