#pragma once

#include "geometry/Grid.h"

#include "pmp/SurfaceMesh.h"

/**
 * \brief Exports scalar field to .vti format.
 * \param filename     in output path.
 * \param scalarGrid   scalar grid to be exported.
 * \throw std::invalid_argument if scalarGrid is invalid or an empty filename is given.
 */
void ExportToVTI(const std::string& filename, const Geometry::ScalarGrid& scalarGrid);

/**
 * \brief Exports vector field to .vtk format as a STRUCTURED_GRID of points with assigned vectors.
 * \param filename      in output path.
 * \param vectorGrid    vector grid to be exported.
 * \throw std::invalid_argument if vectorGrid is invalid or an empty filename is given.
 */
void ExportToVTK(const std::string& filename, const Geometry::VectorGrid& vectorGrid);