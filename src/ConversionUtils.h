#pragma once

#include "geometry/Grid.h"

#include "pmp/SurfaceMesh.h"

/**
 * \brief Exports scalar field to .vti format.
 * \param filename     in output path
 * \param scalarGrid   scalar grid to be exported.
 */
void ExportToVTI(const std::string& filename, const Geometry::ScalarGrid& scalarGrid);