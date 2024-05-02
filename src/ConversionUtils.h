#pragma once

#include <Eigen/Sparse>

#include "geometry/GeometryConversionUtils.h"
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


/// \brief identifier for sparse matrix.
using SparseMatrix = Eigen::SparseMatrix<double>;

/// \brief prints out a trilinear system to a file. For debugging purposes.
void DumpMatrixAndRHSToFile(const std::string& filename, const SparseMatrix& A, const Eigen::MatrixXd& b);

/// \brief Nifti import utility (using bet2 functionality).
//[[nodiscard]] Geometry::ScalarGrid ImportNiftiAsScalarGrid(const std::string& fileName);

/// \brief VTI image data import.
[[nodiscard]] Geometry::ScalarGrid ImportVTI(const std::string& fileName);

/// \brief VTK polydata data import.
[[nodiscard]] Geometry::BaseMeshGeometryData ImportVTK(const std::string& fileName);