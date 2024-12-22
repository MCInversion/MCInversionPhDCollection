#include "GeometryIOUtils.h"

#include "geometry/GridUtil.h"
#include "utils/StringUtils.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_utils/stb_image_write.h"

#define STB_TRUETYPE_IMPLEMENTATION
#include "stb_utils/stb_truetype.h"

#include <fstream>
#include <ranges>
#include <iomanip>
#include <sstream>
#include <cmath>

// set up font directory
const std::string fontDirPath = DSYSTEM_FONT_DIR;

//bool Export2DPointCloudToVTK(const std::vector<pmp::vec2>& points, const std::string& fileName)
//{
//	const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(fileName);
//	if (extension != "vtk")
//	{
//		std::cerr << "Export2DPointCloudToVTK: Invalid file extension!" << std::endl;
//		return false;
//	}
//
//	std::ofstream file(fileName);
//	if (!file.is_open())
//	{
//		std::cerr << "Export2DPointCloudToVTK: Failed to open file for writing: " << fileName << std::endl;
//		return false;
//	}
//
//	file << "# vtk DataFile Version 3.0\n";
//	file << "2D Point Cloud\n";
//	file << "ASCII\n";
//	file << "DATASET POLYDATA\n";
//	file << "POINTS " << points.size() << " float\n";
//
//	for (const auto& point : points)
//	{
//		file << point[0] << " " << point[1] << " " << point[2] << "\n";
//	}
//
//	file << "VERTICES " << points.size() << " " << points.size() * 2 << "\n";
//	for (size_t i = 0; i < points.size(); ++i)
//	{
//		file << "1 " << i << "\n";
//	}
//
//	file.close();
//	return true;
//}

bool Export2DPointCloudToPLY(const std::vector<pmp::vec2>& points, const std::string& fileName)
{
	const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(fileName);
	if (extension != "ply")
	{
		std::cerr << "Export2DPointCloudToPLY: Invalid file extension!" << std::endl;
		return false;
	}

	std::ofstream file(fileName);
	if (!file.is_open())
	{
		std::cerr << "Export2DPointCloudToPLY: Failed to open file for writing: " << fileName << std::endl;
		return false;
	}

	file << "ply\n";
	file << "format ascii 1.0\n";
	file << "element vertex " << points.size() << "\n";
	file << "property float x\n";
	file << "property float y\n";
	file << "property float z\n";
	file << "end_header\n";

	for (const auto& point : points)
	{
		file << point[0] << " " << point[1] << " 0\n"; // Z-coordinate set to 0
	}

	file.close();
	return true;
}

bool ExportManifoldCurve2DToPLY(const pmp::ManifoldCurve2D& curve, const std::string& fileName)
{
	return pmp::write_to_ply(curve, fileName);
}

namespace
{
	[[nodiscard]] RGBColor InterpolateColors(const RGBColor& c1, const RGBColor& c2, pmp::Scalar t)
	{
		return {
			c1.r + t * (c2.r - c1.r),
			c1.g + t * (c2.g - c1.g),
			c1.b + t * (c2.b - c1.b)
		};
	}
}


double CalculateSymmetricColorMapScalingFactor(double minVal, double maxVal)
{
	if (minVal >= maxVal) {
		throw std::invalid_argument("CalculateMidpointScalingFactor: minVal must be less than maxVal.\n");
	}

	// Calculate the offset to center the range around zero
	double rangeMidpoint = (minVal + maxVal) / 2.0;
	double centeredMinVal = minVal - rangeMidpoint;
	double centeredMaxVal = maxVal - rangeMidpoint;

	// Calculate the absolute maximum for symmetrical scaling
	double maxAbsVal = std::max(std::abs(centeredMinVal), std::abs(centeredMaxVal));

	// If the maximum value is zero or near zero, return 1.0 to avoid division by zero
	if (maxAbsVal < std::numeric_limits<double>::epsilon()) {
		return 1.0;
	}

	// Calculate scaling factor to map values symmetrically to the [0, 1] range
	double scalingFactor = 1.0 / maxAbsVal;

	return scalingFactor;
}

RGBColorScheme AdjustColorMapForZeroMidpoint(const RGBColorScheme& scheme, double minVal, double maxVal)
{
	if (minVal >= maxVal)
	{
		throw std::invalid_argument("AdjustColorMapForMidpoint: minVal must be less than maxVal.\n");
	}

	// Check if 0.0 is within the range
	if (minVal > 0.0 || maxVal < 0.0)
	{
		// Exit early if 0.0 is not in the range
		std::cerr << "AdjustColorMapForMidpoint: 0.0 is not within the range. No adjustment applied.\n";
		return scheme; // Return the original scheme unmodified
	}

	double rangeMidpoint = 0.0; // Centering around zero
	double maxAbsVal = std::max(std::abs(minVal), std::abs(maxVal));

	if (maxAbsVal < std::numeric_limits<double>::epsilon())
	{
		throw std::invalid_argument("AdjustColorMapForMidpoint: Range is too small to adjust.\n");
	}

	RGBColorScheme adjustedScheme;
	for (const auto& [key, color] : scheme)
	{
		// Adjust keys symmetrically around 0.0 so that 0.0 maps to 0.5
		double adjustedKey;
		if (key < 0.5)
		{
			// Keys below 0.5 (negative side)
			adjustedKey = 0.5 * (key / 0.5) * (std::abs(minVal) / maxAbsVal);
		}
		else
		{
			// Keys above 0.5 (positive side)
			adjustedKey = 0.5 + 0.5 * ((key - 0.5) / 0.5) * (std::abs(maxVal) / maxAbsVal);
		}
		adjustedScheme[adjustedKey] = color;
	}

	return adjustedScheme;
}

void ExportScalarGrid2DToPNG(const std::string& filename, const Geometry::ScalarGrid2D& grid,
	const std::function<double(const pmp::vec2&, const Geometry::ScalarGrid2D&)>& interpolate,
	pmp::Scalar nPixelsPerCellX, pmp::Scalar nPixelsPerCellY, const RGBColorScheme& colorMap)
{
	if (!grid.IsValid())
	{
		throw std::invalid_argument("ExportScalarGrid2DToPNG: grid to be exported is invalid!\n");
	}
	if (filename.empty())
	{
		throw std::invalid_argument("ExportScalarGrid2DToPNG: filename cannot be empty!\n");
	}

	if (colorMap.size() < 2)
	{
		throw std::invalid_argument("ExportScalarGrid2DToPNG: color scheme must contain at least 2 colors!\n");
	}

	const auto& dims = grid.Dimensions();
	const int nx = static_cast<int>(dims.Nx);
	const int ny = static_cast<int>(dims.Ny);
	const auto& orig = grid.Box().min();
	const pmp::Scalar cellSize = grid.CellSize();

	const int imageWidth = static_cast<int>(nx * nPixelsPerCellX);
	const int imageHeight = static_cast<int>(ny * nPixelsPerCellY);

	std::vector<uint8_t> image(imageWidth * imageHeight * 3);

	// Normalize the scalar values to [0, 1]
	const auto minMax = std::minmax_element(grid.Values().begin(), grid.Values().end());
	const double minVal = *minMax.first;
	const double maxVal = *minMax.second;

	for (int i = 0; i < imageWidth; ++i)
	{
		for (int j = 0; j < imageHeight; ++j)
		{
			const auto x = orig[0] + static_cast<pmp::Scalar>(i) / static_cast<pmp::Scalar>(nPixelsPerCellX) * cellSize;
			const auto y = orig[1] + static_cast<pmp::Scalar>(j) / static_cast<pmp::Scalar>(nPixelsPerCellY) * cellSize;

			// Get the corresponding grid value and normalize it
			const double value = interpolate(pmp::vec2(x, y), grid);
			const double normalizedValue = (value - minVal) / (maxVal - minVal);

			// Find the two closest colors in the scheme
			auto upper = colorMap.lower_bound(normalizedValue);
			auto lower = (upper == colorMap.begin()) ? upper : std::prev(upper);

			// Handle edge cases where normalizedValue is less than the first color or greater than the last color
			if (upper == colorMap.end())
			{
				upper = std::prev(upper);
			}
			if (lower == upper && upper == colorMap.begin())
			{
				lower = upper;
			}

			// Interpolate between the two colors
			const double range = upper->first - lower->first;
			const pmp::Scalar t = (range > 0.0) ? static_cast<pmp::Scalar>((normalizedValue - lower->first) / range) : 0.0;
			const auto color = InterpolateColors(lower->second, upper->second, t);

			// Set the pixel value (RGB)
			image[3 * (i + j * imageWidth) + 0] = static_cast<uint8_t>(255.0 * color.r); // R
			image[3 * (i + j * imageWidth) + 1] = static_cast<uint8_t>(255.0 * color.g); // G
			image[3 * (i + j * imageWidth) + 2] = static_cast<uint8_t>(255.0 * color.b); // B
		}
	}

	// Save the image to a file
	if (!stbi_write_png(filename.c_str(), imageWidth, imageHeight, 3, image.data(), imageWidth * 3))
	{
		throw std::runtime_error("ExportToPNG: Failed to save image to file!");
	}
}

void ExportScalarGridDimInfo2D(const std::string& fileName, const Geometry::ScalarGrid2D& grid)
{
	// Validate input grid
	if (!grid.Dimensions().Valid() || grid.CellSize() <= 0.0)
	{
		throw std::invalid_argument("ExportScalarGridDimInfo2D: Invalid grid dimensions or cell size.");
	}

	// Validate filename
	if (fileName.empty())
	{
		throw std::invalid_argument("ExportScalarGridDimInfo2D: Filename cannot be empty.");
	}

	const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(fileName);
	if (extension != "gdim2d")
	{
		throw std::invalid_argument("ExportScalarGridDimInfo2D: Invalid file extension!");
	}

	// Open the output file
	std::ofstream outFile(fileName);
	if (!outFile)
	{
		throw std::runtime_error("ExportScalarGridDimInfo2D: Failed to open file for writing.");
	}

	// Get the grid information
	const auto& bbox = grid.Box();
	const auto& dims = grid.Dimensions();
	const auto cellSize = grid.CellSize();

	// Write the grid information to the file
	outFile << "Grid Dimensions (Nx, Ny): " << dims.Nx << ", " << dims.Ny << "\n";
	outFile << "Bounding Box (min_x, min_y) -> (max_x, max_y): "
		<< "(" << bbox.min()[0] << ", " << bbox.min()[1] << ") -> "
		<< "(" << bbox.max()[0] << ", " << bbox.max()[1] << ")\n";
	outFile << "Cell Size: " << cellSize << "\n";

	// Close the file
	outFile.close();
	if (!outFile)
	{
		throw std::runtime_error("ExportScalarGridDimInfo2D: Error occurred while writing to file.");
	}
}

void ExportVectorGridDimInfo2D(const std::string& filename, const Geometry::VectorGrid2D& grid)
{
	// Validate input grid
	if (!grid.Dimensions().Valid() || grid.CellSize() <= 0.0)
	{
		throw std::invalid_argument("ExportVectorGridDimInfo2D: Invalid grid dimensions or cell size.");
	}

	// Validate filename
	if (filename.empty())
	{
		throw std::invalid_argument("ExportVectorGridDimInfo2D: filename cannot be empty.");
	}

	const auto extension = Utils::ExtractLowercaseFileExtensionFromPath(filename);

	if (extension != "gdim2d")
	{
		throw std::invalid_argument("ExportVectorGridDimInfo2D: Invalid file extension!");
	}

	// Open the output file
	std::ofstream outFile(filename);
	if (!outFile)
	{
		throw std::runtime_error("ExportVectorGridDimInfo2D: Failed to open file for writing.");
	}

	// Get the grid information
	const auto& bbox = grid.Box();
	const auto& dims = grid.Dimensions();
	const auto cellSize = grid.CellSize();

	// Write the grid information to the file
	outFile << "Grid Dimensions (Nx, Ny): " << dims.Nx << ", " << dims.Ny << "\n";
	outFile << "Bounding Box (min_x, min_y) -> (max_x, max_y): "
		<< "(" << bbox.min()[0] << ", " << bbox.min()[1] << ") -> "
		<< "(" << bbox.max()[0] << ", " << bbox.max()[1] << ")\n";
	outFile << "Cell Size: " << cellSize << "\n";

	// Close the file
	outFile.close();
	if (!outFile)
	{
		throw std::runtime_error("ExportVectorGridDimInfo2D: Error occurred while writing to file.");
	}
}

RGBColorScheme operator*(const RGBColorScheme& scheme, double factor)
{
	RGBColorScheme result;

	bool is_first = true;
	for (const auto& [key, color] : scheme) {
		if (is_first) {
			// Keep the first key (0.0) intact
			result[key] = color;
			is_first = false;
		}
		else {
			// Scale the rest of the keys by the factor
			result[key * factor] = color;
		}
	}

	return result;
}

//
// ===================================================================
//

constexpr unsigned int BORDER_SIZE{ 4 };
constexpr unsigned int TICK_WIDTH{ 3 };
constexpr unsigned int N_TICKS{ 4 };
constexpr double PADDING_FRACTION{ 0.12 };  // Padding as a fraction of the smaller dimension
constexpr double TICK_LENGTH_FRACTION{ 0.2 }; // Tick length as a fraction of the smaller dimension

static [[nodiscard]] RGBColor BlendColor(double ratio, const RGBColorScheme& colorMap)
{
	// Handle the case where the value is less than the lowest key
	if (ratio <= colorMap.begin()->first)
	{
		return colorMap.begin()->second;
	}

	// Handle the case where the value is greater than the highest key
	if (ratio >= colorMap.rbegin()->first)
	{
		return colorMap.rbegin()->second;
	}

	// At this point, we know that the value is within the range of the color map keys
	const auto upper = colorMap.upper_bound(ratio);
	const auto lower = std::prev(upper);

	const double range = upper->first - lower->first;
	const double factor = (ratio - lower->first) / range;

	// Linear interpolation between lower and upper colors
	RGBColor interpolated;
	interpolated.r = lower->second.r + factor * (upper->second.r - lower->second.r);
	interpolated.g = lower->second.g + factor * (upper->second.g - lower->second.g);
	interpolated.b = lower->second.b + factor * (upper->second.b - lower->second.b);

	return interpolated;
}

// Define a structure to hold font data and character bitmaps
struct FontData
{
	std::unique_ptr<unsigned char[]> bitmap;  // Dynamically allocated bitmap
	stbtt_bakedchar cdata[96];                // ASCII 32..126 is 95 glyphs

	FontData() : bitmap(std::make_unique<unsigned char[]>(512 * 512)) {} // Allocate bitmap memory on the heap
};

// Loads a font from specified path
static [[nodiscard]] FontData LoadFont(const std::string& fontFilePath, pmp::Scalar fontSize)
{
	FontData fontData;

	// Use std::vector for dynamically allocated buffer on the heap
	std::vector<unsigned char> ttf_buffer(1 << 20);  // 1 MB buffer for font loading

	// Read the font from file
	FILE* fontFile = fopen(fontFilePath.c_str(), "rb");
	if (!fontFile)
	{
		throw std::invalid_argument("LoadFont: Failed to load font file");
	}

	// Load the font into the buffer
	const size_t bytesRead = fread(ttf_buffer.data(), 1, ttf_buffer.size(), fontFile);
	fclose(fontFile);

	if (bytesRead <= 0)
	{
		throw std::runtime_error("LoadFont: Failed to read font data from file");
	}

	// Bake the font bitmap at a given font size
	if (stbtt_BakeFontBitmap(ttf_buffer.data(), 0, fontSize, fontData.bitmap.get(), 512, 512, 32, 96, fontData.cdata) <= 0)
	{
		throw std::runtime_error("LoadFont: Failed to bake font bitmap");
	}

	return fontData;
}

// Function to format the tick value with a specified number of decimal places and rules for scientific notation
static [[nodiscard]] std::string FormatTickValue(double value, unsigned int decimalPlaces, bool forceZero = true)
{
	std::stringstream ss;

	// Case 1: value >= 1000, use scientific notation with decimalPlaces + 1 digits
	if (value >= 1000.0)
	{
		ss << std::scientific << std::setprecision(decimalPlaces + 1) << value;
	}
	// Case 2: value <= 1e-decimalPlaces, use scientific notation with decimalPlaces + 1 digits
	else if (value <= std::pow(10, -static_cast<int>(decimalPlaces)) && !forceZero)
	{
		ss << std::scientific << std::setprecision(decimalPlaces) << value;
	}
	// Case 3: Otherwise, use fixed-point notation with decimalPlaces decimal places
	else
	{
		ss << std::fixed << std::setprecision(decimalPlaces) << value;
	}

	return ss.str();
}

// Function to draw text on the image using the baked font
void DrawTextOnImage(std::vector<unsigned char>& image, unsigned int imgWidth, unsigned int imgHeight, const FontData& fontData, const std::string& text, int posX, int posY)
{
	for (const char ch : text)
	{
		if (ch < 32 || ch > 126) continue;  // Ignore non-printable characters

		const stbtt_bakedchar* b = &fontData.cdata[static_cast<unsigned int>(ch - 32)];
		const int x = posX + static_cast<int>(b->xoff);
		const int y = posY + static_cast<int>(b->yoff);
		const auto w = static_cast<int>(b->x1 - b->x0);
		const auto h = static_cast<int>(b->y1 - b->y0);

		// Blit the character bitmap onto the image buffer
		for (int i = 0; i < h; ++i)
		{
			for (int j = 0; j < w; ++j)
			{
				if (x + j >= 0 && x + j < static_cast<int>(imgWidth) && y + i >= 0 && y + i < static_cast<int>(imgHeight))
				{
					const unsigned char pixel = fontData.bitmap[static_cast<unsigned int>((b->y0 + i) * 512 + (b->x0 + j))];
					if (pixel > 0)
					{
						const int idx = static_cast<int>(((y + i) * imgWidth + (x + j)) * 3);
						image[idx + 0] = 0; // R
						image[idx + 1] = 0; // G
						image[idx + 2] = 0; // B
					}
				}
			}
		}
		posX += static_cast<int>(b->xadvance);
	}
}

void ExportColorScaleToPNG(const std::string& filename, const double& minValue, const double& maxValue, const RGBColorScheme& colorMap, const unsigned int& imageHeight, const unsigned int& imageWidth)
{
	if (minValue > maxValue)
	{
		throw std::invalid_argument("ExportColorScaleToPNG: minValue > maxValue!\n");
	}
	if (imageHeight == 0 || imageWidth == 0)
	{
		throw std::invalid_argument("ExportColorScaleToPNG: imageHeight == 0 || imageWidth == 0!\n");
	}
	if (filename.empty())
	{
		throw std::invalid_argument("ExportColorScaleToPNG: filename cannot be empty!\n");
	}
	if (colorMap.size() < 2)
	{
		throw std::invalid_argument("ExportColorScaleToPNG: color scheme must contain at least 2 colors!\n");
	}

	// Determine the padding and tick length as fractions of the smallest image dimension
	const auto isVertical = imageHeight > imageWidth;
	const auto [smallestDim, largestDim] = std::minmax(imageWidth, imageHeight);
	const auto paddingY = static_cast<unsigned int>(PADDING_FRACTION * smallestDim);
	const auto paddingX = (isVertical ? 3 : 1) * paddingY;
	const int offsetX = isVertical ? -static_cast<int>(paddingX) / 2 : 0;
	const int offsetY = isVertical ? 0 : -static_cast<int>(paddingY) / 2;
	const unsigned int colorBarLargestDim = largestDim - 2 * (BORDER_SIZE + (isVertical ? paddingY : paddingX));
	const unsigned int scalePxStartOffset = (isVertical ? paddingY + offsetY : paddingX + offsetX) + BORDER_SIZE;

	// Allocate image buffer (RGB, so 3 channels)
	std::vector<unsigned char> image(imageWidth * imageHeight * 3, 255); // White background by default

	// Fill the color scale
	for (unsigned int i = scalePxStartOffset; i < scalePxStartOffset + colorBarLargestDim; ++i)
	{
		const double ratio = (i - scalePxStartOffset) / static_cast<double>(colorBarLargestDim);
		const RGBColor color = BlendColor(ratio, colorMap);

		// draw isochromatic lines
		if (isVertical)
		{
			for (unsigned int x = offsetX + paddingX + BORDER_SIZE; x < offsetX + imageWidth - paddingX - BORDER_SIZE; ++x)
			{
				image[(i * imageWidth + x) * 3 + 0] = static_cast<unsigned char>(color.r * 255);
				image[(i * imageWidth + x) * 3 + 1] = static_cast<unsigned char>(color.g * 255);
				image[(i * imageWidth + x) * 3 + 2] = static_cast<unsigned char>(color.b * 255);
			}
		}
		else
		{
			for (unsigned int y = offsetY + paddingY + BORDER_SIZE; y < offsetY + imageHeight - paddingY - BORDER_SIZE; ++y)
			{
				image[(y * imageWidth + i) * 3 + 0] = static_cast<unsigned char>(color.r * 255);
				image[(y * imageWidth + i) * 3 + 1] = static_cast<unsigned char>(color.g * 255);
				image[(y * imageWidth + i) * 3 + 2] = static_cast<unsigned char>(color.b * 255);
			}
		}
	}

	// Draw the black outline
	const auto maxX = offsetX + imageWidth - paddingX - BORDER_SIZE;
	const auto maxY = offsetY + imageHeight - paddingY - BORDER_SIZE;
	for (unsigned int i = 0; i < BORDER_SIZE; ++i)
	{
		// Top and bottom borders (horizontal)
		for (unsigned int x = offsetX + paddingX; x < maxX; ++x)
		{
			image[((i + paddingY) * imageWidth + x) * 3 + 0] = 0;
			image[((i + paddingY) * imageWidth + x) * 3 + 1] = 0;
			image[((i + paddingY) * imageWidth + x) * 3 + 2] = 0;
			image[((imageHeight - 1 - (i + paddingY)) * imageWidth + x) * 3 + 0] = 0;
			image[((imageHeight - 1 - (i + paddingY)) * imageWidth + x) * 3 + 1] = 0;
			image[((imageHeight - 1 - (i + paddingY)) * imageWidth + x) * 3 + 2] = 0;
		}

		// Left and right borders (vertical)
		for (unsigned int y = offsetY + paddingY; y < maxY; ++y)
		{
			image[(y * imageWidth + i + offsetX + paddingX) * 3 + 0] = 0;
			image[(y * imageWidth + i + offsetX + paddingX) * 3 + 1] = 0;
			image[(y * imageWidth + i + offsetX + paddingX) * 3 + 2] = 0;
			image[(y * imageWidth + (imageWidth - 1 - (i - offsetX + paddingX))) * 3 + 0] = 0;
			image[(y * imageWidth + (imageWidth - 1 - (i - offsetX + paddingX))) * 3 + 1] = 0;
			image[(y * imageWidth + (imageWidth - 1 - (i - offsetX + paddingX))) * 3 + 2] = 0;
		}
	}

	// --------------  Draw tick marks ----------------
	const auto tickLength = static_cast<unsigned int>(TICK_LENGTH_FRACTION * smallestDim);

	// Calculate font size dynamically based on padding
	constexpr unsigned int decimalPlaces = 2;
	const auto fontSize = paddingX * 1.2 / static_cast<pmp::Scalar>(decimalPlaces);  // Dynamic font size, smaller than padding
	const FontData fontData = LoadFont(fontDirPath + "\\CENTURY.TTF", fontSize); // Load your TTF font file

	for (unsigned int t = 0; t <= N_TICKS; ++t)
	{
		// Calculate the position of each tick based on the ratio along the color bar
		const double tickRatio = static_cast<double>(t) / N_TICKS;
		const auto tickPos = static_cast<unsigned int>(scalePxStartOffset + tickRatio * colorBarLargestDim);
		const double tickValue = minValue + tickRatio * (maxValue - minValue);

		const std::string formattedTickValue = FormatTickValue(tickValue, decimalPlaces);

		if (isVertical)
		{
			//// Draw vertical ticks on the left edge
			//for (unsigned int y = 0; y < TICK_WIDTH; ++y)
			//{
			//	for (unsigned int x = 0; x < tickLength; ++x)
			//	{
			//		image[((tickPos + y) * imageWidth + (padding + x)) * 3 + 0] = 0;
			//		image[((tickPos + y) * imageWidth + (padding + x)) * 3 + 1] = 0;
			//		image[((tickPos + y) * imageWidth + (padding + x)) * 3 + 2] = 0;
			//	}
			//}

			//// Render the tick value next to the tick mark
			//DrawTextOnImage(image, imageWidth, imageHeight, fontData, formattedTickValue, 5, tickPos - fontSize / 2);

			// Draw vertical ticks on the right edge
			for (unsigned int y = 0; y < TICK_WIDTH; ++y)
			{
				for (unsigned int x = maxX - tickLength - 2; x < maxX + 2; ++x)
				{
					image[((tickPos + y) * imageWidth + x) * 3 + 0] = 0;
					image[((tickPos + y) * imageWidth + x) * 3 + 1] = 0;
					image[((tickPos + y) * imageWidth + x) * 3 + 2] = 0;
				}
			}

			// Render the tick value next to the tick mark
			DrawTextOnImage(image, imageWidth, imageHeight, fontData, formattedTickValue, maxX + 5, tickPos + fontSize / 4);
		}
		else
		{
			// Draw horizontal ticks on the bottom edge
			for (unsigned int x = 0; x < TICK_WIDTH; ++x)
			{
				for (unsigned int y = 0; y < tickLength; ++y)
				{
					image[((imageHeight - paddingY - y - offsetY) * imageWidth + (tickPos + x + offsetX)) * 3 + 0] = 0;
					image[((imageHeight - paddingY - y - offsetY) * imageWidth + (tickPos + x + offsetX)) * 3 + 1] = 0;
					image[((imageHeight - paddingY - y - offsetY) * imageWidth + (tickPos + x + offsetX)) * 3 + 2] = 0;
				}
			}

			// Render the tick value below the tick mark
			DrawTextOnImage(image, imageWidth, imageHeight, fontData, formattedTickValue, tickPos - fontSize / 2, imageHeight - paddingY / 2 + 10);
		}
	}

	// Save the image to a file
	if (!stbi_write_png(filename.c_str(), imageWidth, imageHeight, 3, image.data(), imageWidth * 3))
	{
		throw std::runtime_error("ExportColorScaleToPNG: Failed to save image to file!");
	}
}

void ExportPolyLinesToPLY(const std::vector<std::vector<pmp::Point2>>& polyLines, const std::string& fileName)
{
	// Ensure the file has the correct extension
	const size_t lastSlash = fileName.find_last_of("/\\");
	const size_t dot = fileName.rfind('.');

	if (dot == std::string::npos || (lastSlash != std::string::npos && dot < lastSlash))
	{
		std::cerr << "ExportPolyLinesToPLY [ERROR]: Invalid file extension or path!\n";
		return;
	}

	std::string ext = fileName.substr(dot + 1);
	std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

	if (ext != "ply")
	{
		std::cerr << "ExportPolyLinesToPLY: Invalid file extension. Expected '.ply'!" << std::endl;
		return;
	}

	std::ofstream file(fileName);
	if (!file.is_open())
	{
		std::cerr << "ExportPolyLinesToPLY: Failed to open file for writing: " << fileName << std::endl;
		return;
	}

	// Calculate total number of vertices and edges
	size_t totalVertices = 0;
	size_t totalEdges = 0;
	for (const auto& polyline : polyLines)
	{
		totalVertices += polyline.size();
		if (polyline.size() > 1)
		{
			totalEdges += (polyline.size() - 1);
		}
	}

	// Write the PLY header
	file << "ply" << std::endl;
	file << "format ascii 1.0" << std::endl;
	file << "element vertex " << totalVertices << std::endl;
	file << "property float x" << std::endl;
	file << "property float y" << std::endl;
	file << "property float z" << std::endl;
	file << "element edge " << totalEdges << std::endl;
	file << "property int vertex1" << std::endl;
	file << "property int vertex2" << std::endl;
	file << "end_header" << std::endl;

	// Write vertex data
	size_t vertexIndex = 0;
	std::vector<size_t> vertexIndices; // Store vertex indices for edges
	for (const auto& polyline : polyLines)
	{
		for (const auto& point : polyline)
		{
			file << point[0] << " " << point[1] << " 0.0" << std::endl; // z = 0 for 2D points
			vertexIndices.push_back(vertexIndex++);
		}
	}

	// Write edge data
	size_t indexOffset = 0;
	for (const auto& polyline : polyLines)
	{
		for (size_t i = 1; i < polyline.size(); ++i)
		{
			file << (indexOffset + i - 1) << " " << (indexOffset + i) << std::endl;
		}
		indexOffset += polyline.size();
	}

	// File closes automatically when going out of scope (RAII)
}
