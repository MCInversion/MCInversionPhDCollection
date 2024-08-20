#pragma once

#include <Eigen/Sparse>

#include "geometry/GeometryConversionUtils.h"
#include "geometry/Grid.h"
#include "pmp/ManifoldCurve2D.h"

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

/// \brief 2D point cloud to VTK polydata data export.
//[[nodiscard]] bool Export2DPointCloudToVTK(const std::vector<pmp::vec2>& points, const std::string& fileName);

/// \brief 2D point cloud to PLY data export.
[[nodiscard]] bool Export2DPointCloudToPLY(const std::vector<pmp::vec2>& points, const std::string& fileName);

/// \brief Export 2D manifold curve to PLY data.
[[nodiscard]] bool ExportManifoldCurve2DToPLY(const pmp::ManifoldCurve2D& curve, const std::string& fileName);

/// \brief Export a 2D scalar field as an RGB PNG image
void ExportScalarGrid2DToPNG(const std::string& filename, const Geometry::ScalarGrid2D& grid,
	const std::function<double(const pmp::vec2&, const Geometry::ScalarGrid2D&)>& interpolate, 
	float nPixelsPerCellX = 1, float nPixelsPerCellY = 1);