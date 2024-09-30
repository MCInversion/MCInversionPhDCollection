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

/// \brief A wrapper for normalized RGB values in a color scheme
struct RGBColor
{
	float r, g, b; // from [0, 1]
};

/// \brief A color map placeholder
using RGBColorScheme = std::map<double, RGBColor>;

/// \brief Default grayscale color scheme
const RGBColorScheme GRAYSCALE_MAP = {
	{ 0.0, { 0.0f, 0.0f, 0.0f } }, // Black
	{ 1.0, { 1.0f, 1.0f, 1.0f } }  // White
};

/// \brief Rainbow to white color scheme
const RGBColorScheme RAINBOW_TO_WHITE_MAP = {
	{ 0.0, { 1.0f, 0.0f, 0.0f } }, // Red
	{ 0.05, { 0.988f, 0.518f, 0.012f } }, // Orange (252, 132, 3)
	{ 0.1, { 0.988f, 0.957f, 0.012f } }, // Yellow (252, 244, 3)
	{ 0.15, { 1.0f, 0.8706f, 0.3490f } }, // Yellowish (255, 222, 89)
	//{ 0.2, {0.5294f, 0.8314f, 0.5059f } }, // Green (135, 212, 129)
	{ 0.225, { 0.498f, 0.859f, 0.961f } }, // Light Blue (127, 219, 245)
	{ 0.3, { 0.298f, 0.459f, 1.0f } }, // Blue
	{ 0.4, { 1.0f, 1.0f, 1.0f } }, // White 1
	{ 1.0, { 1.0f, 1.0f, 1.0f } }  // White 2
};

/// \brief Export a 2D scalar field as an RGB PNG image
void ExportScalarGrid2DToPNG(const std::string& filename, const Geometry::ScalarGrid2D& grid,
	const std::function<double(const pmp::vec2&, const Geometry::ScalarGrid2D&)>& interpolate, 
	float nPixelsPerCellX = 1, float nPixelsPerCellY = 1, const RGBColorScheme& colorMap = GRAYSCALE_MAP);

/// \brief A specialized export of file (*.gdim2d) for containing the proper dimensions of the 2D scalar grid.
void ExportScalarGridDimInfo2D(const std::string& filename, const Geometry::ScalarGrid2D& grid);

/// \brief A specialized export of file (*.gdim2d) for containing the proper dimensions of the 2D vector grid.
void ExportVectorGridDimInfo2D(const std::string& filename, const Geometry::VectorGrid2D& grid);

/// \brief Overload multiplication for RGBColorScheme
RGBColorScheme operator*(const RGBColorScheme& scheme, double factor);

/// \brief Exports an image of the used color scale matched to the given value range.
void ExportColorScaleToPNG(const std::string& filename, 
	const double& minValue, const double& maxValue, 
	const RGBColorScheme& colorMap,
	const unsigned int& imageHeight, const unsigned int& imageWidth);