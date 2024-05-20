#pragma once

#include "GeometryConversionUtils.h"

namespace Geometry
{
	/// \brief Generates PolyIndices as triples using Delaunay triangulation from the Fade2D library.
	void TriangulateWithFade2D(BaseMeshGeometryData& data, const std::vector<pmp::Point>& constraintPolyline);

	/// \brief Generates PolyIndices as triples using the ball pivoting algorithm by the VCG library. The ball radius is computed from the 6 nearest neighbor distance.
	void TriangulateWithVCGBPA(BaseMeshGeometryData& data);
	
} // namespace Geometry