#pragma once

#include "geometry/GeometryConversionUtils.h"
#include "geometry/Grid.h"

/// \brief 2D point cloud to VTK polydata data export.
//[[nodiscard]] bool Export2DPointCloudToVTK(const std::vector<pmp::vec2>& points, const std::string& fileName);

/// \brief 2D point cloud to PLY data export.
[[nodiscard]] bool Export2DPointCloudToPLY(const std::vector<pmp::vec2>& points, const std::string& fileName);

/// \brief Export 2D manifold curve to PLY data.
[[nodiscard]] bool ExportManifoldCurve2DToPLY(const pmp::ManifoldCurve2D& curve, const std::string& fileName);

/// \brief A PLY exporter for a general curve geometry.
[[nodiscard]] bool ExportBaseCurveGeometryDataToPLY(const Geometry::BaseCurveGeometryData& geomData, const std::string& absFileName);

/// \brief A wrapper for normalized RGB values in a color scheme
struct RGBColor
{
	pmp::Scalar r, g, b; // from [0, 1]
};

/// \brief A color map placeholder
using RGBColorScheme = std::map<double, RGBColor>;

/// \brief Default grayscale color scheme
const RGBColorScheme GRAYSCALE_MAP = {
	{ 0.0, { 0.0, 0.0, 0.0 } }, // Black
	{ 1.0, { 1.0, 1.0, 1.0 } }  // White
};

/// \brief Rainbow to white color scheme
const RGBColorScheme RAINBOW_TO_WHITE_MAP = {
	{ 0.0, { 1.0, 0.0, 0.0 } }, // Red
	{ 0.05, { 0.988, 0.518, 0.012 } }, // Orange (252, 132, 3)
	{ 0.1, { 0.988, 0.957, 0.012 } }, // Yellow (252, 244, 3)
	{ 0.15, { 1.0, 0.8706, 0.3490 } }, // Yellowish (255, 222, 89)
	//{ 0.2, {0.5294, 0.8314, 0.5059 } }, // Green (135, 212, 129)
	{ 0.225, { 0.498, 0.859, 0.961 } }, // Light Blue (127, 219, 245)
	{ 0.3, { 0.298, 0.459, 1.0 } }, // Blue
	{ 0.4, { 1.0, 1.0, 1.0 } }, // White 1
	{ 1.0, { 1.0, 1.0, 1.0 } }  // White 2
};

const RGBColorScheme SIGN_TEMP_MAP = {
	{ 0.0, { 0.498, 0.859, 0.961 } }, // Light Blue (127, 219, 245)
	{ 0.4, {0.0, 0.0, 1.0}}, // Blue
	{ 0.5, { 1.0, 1.0, 1.0 } },  // White
	{ 0.6, { 1.0, 0.0, 0.0 } }, // Red
	{ 1.0, { 0.988, 0.518, 0.012 } } // Orange (252, 132, 3)
};

/// \brief Scales the color map so that zero is in the middle
[[nodiscard]] RGBColorScheme AdjustColorMapForZeroMidpoint(const RGBColorScheme& scheme, double minVal, double maxVal);

/// \brief Export a 2D scalar field as an RGB PNG image
void ExportScalarGrid2DToPNG(const std::string& filename, const Geometry::ScalarGrid2D& grid,
	const std::function<double(const pmp::vec2&, const Geometry::ScalarGrid2D&)>& interpolate,
	pmp::Scalar nPixelsPerCellX = 1, pmp::Scalar nPixelsPerCellY = 1, const RGBColorScheme& colorMap = GRAYSCALE_MAP);

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

/// \brief Exports polylines to ply
void ExportPolyLinesToPLY(const std::vector<std::vector<pmp::Point2>>& polyLines, const std::string& fileName);

/// \brief Imports a png image
///	\param[in] absFileName           absolute filename of the image to be opened.
///	\param[in] normalizationRange    optional range for normalizing the [0, 255] values of the grayscale image.
///	\return optional imageGrid if successful.
[[nodiscard]] std::optional<Geometry::ScalarGrid2D> ImportPNGImageGrayscale(const std::string& absFileName,
	const std::optional<std::pair<pmp::Scalar, pmp::Scalar>>& normalizationRange = std::nullopt);