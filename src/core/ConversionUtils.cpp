#include "ConversionUtils.h"

#include <fstream>
#include <ranges>
#include <iomanip>
#include <sstream>
#include <cmath>

#include "geometry/GridUtil.h"
#include "utils/StringUtils.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_utils/stb_image_write.h"

#define STB_TRUETYPE_IMPLEMENTATION
#include "stb_utils/stb_truetype.h"

// set up font directory
const std::string fontDirPath = DSYSTEM_FONT_DIR;

void ExportToVTI(const std::string& filename, const Geometry::ScalarGrid& scalarGrid)
{
    if (!scalarGrid.IsValid())
        throw std::invalid_argument("ExportToVTI: scalarGrid to be exported is invalid!\n");
    if (filename.empty())
        throw std::invalid_argument("ExportToVTI: filename cannot be empty!\n");

    std::fstream vti(filename + ".vti", std::fstream::out);

    const auto& dims = scalarGrid.Dimensions();
    const auto nx = static_cast<unsigned int>(dims.Nx);
    const auto ny = static_cast<unsigned int>(dims.Ny);
    const auto nz = static_cast<unsigned int>(dims.Nz);
    const float dx = scalarGrid.CellSize();

    const pmp::vec3 min = scalarGrid.Box().min();
    //const pmp::vec3 max = min + pmp::vec3(static_cast<float>(nx), static_cast<float>(ny),static_cast<float>(nz)) * dx;
    const pmp::vec3 max = scalarGrid.Box().max();

    vti << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    vti << "	<ImageData WholeExtent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1 << "\" Origin=\"" << min[0] + 0.5f * dx << " " << min[1] + 0.5f * dx << " " << min[2] + 0.5f * dx << "\" Spacing=\"" << dx << " " << dx << " " << dx << "\">\n";
    vti << "		<Piece Extent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1 << "\">\n";
    vti << "			<PointData Scalars=\"Scalars_\">\n";
    vti << "				<DataArray type=\"Float32\" Name=\"Scalars_\" format=\"ascii\" RangeMin=\"" << min << "\" RangeMax=\"" << max << "\">\n";

    for (const auto& val : scalarGrid.Values()) {
        vti << static_cast<float>(val) << "\n";
    }

    vti << "				</DataArray>\n";
    vti << "			</PointData>\n";
    vti << "		<CellData>\n";
    vti << "		</CellData>\n";
    vti << "	</Piece>\n";
    vti << "	</ImageData>\n";
    vti << "</VTKFile>\n";

    vti.close();
}

void ExportToVTK(const std::string& filename, const Geometry::VectorGrid& vectorGrid)
{
    if (!vectorGrid.IsValid())
        throw std::invalid_argument("ExportToVTI: vectorGrid to be exported is invalid!\n");
    if (filename.empty())
        throw std::invalid_argument("ExportToVTI: filename cannot be empty!\n");

    std::fstream vtk(filename + ".vtk", std::fstream::out);

    const auto& [Nx, Ny, Nz] = vectorGrid.Dimensions();
    const auto nx = static_cast<unsigned int>(Nx);
    const auto ny = static_cast<unsigned int>(Ny);
    const auto nz = static_cast<unsigned int>(Nz);
    const float dx = vectorGrid.CellSize();

    const pmp::vec3 orig = vectorGrid.Box().min();

    vtk << "# vtk DataFile Version 3.0\n";
    vtk << "vtk output\n";
    vtk << "ASCII\n";
    vtk << "DATASET STRUCTURED_GRID\n";
    vtk << "DIMENSIONS " << nx << " " << ny << " " << nz <<"\n";

    const size_t nValues = Nx * Ny * Nz;
    vtk << "POINTS " << nValues << " float" << "\n";

    for (unsigned int iz = 0; iz < nz; iz++) {
        for (unsigned int iy = 0; iy < ny; iy++) {
            for (unsigned int ix = 0; ix < nx; ix++) {
                vtk << 
                    orig[0] + static_cast<float>(ix) * dx << " " <<
                    orig[1] + static_cast<float>(iy) * dx << " " <<
                    orig[2] + static_cast<float>(iz) * dx << "\n";
            }
        }
    }

    vtk << "POINT_DATA " << nValues << "\n";
    vtk << "VECTORS Vectors_ double" << "\n";

    const auto& valuesX = vectorGrid.ValuesX();
    const auto& valuesY = vectorGrid.ValuesY();
    const auto& valuesZ = vectorGrid.ValuesZ();

    for (unsigned int i = 0; i < nValues; i++) {
        vtk << valuesX[i] << " " << valuesY[i] << " " << valuesZ[i] << "\n";
    }

    vtk.close();
}

void DumpMatrixAndRHSToFile(const std::string& filename, const SparseMatrix& A, const Eigen::MatrixXd& b)
{
	if (A.cols() != A.rows() || A.rows() != b.rows())
	{
		throw std::invalid_argument("DumpMatrixAndRHSToFile: invalid matrix dimensions!\n");
	}

	std::ofstream fileOStream(filename);
	if (!fileOStream.is_open())
		return;

	std::vector<double> diagElems(A.outerSize());
	std::vector<std::vector<double>> offDiagElems(A.outerSize());
	std::vector<Eigen::Vector3d> rhsElems(b.rows());

	constexpr unsigned int expectedValence = 6;
	unsigned int maxOffDiagRowLength = 0;

	for (unsigned int i = 0; i < A.outerSize(); i++)
	{
		diagElems[i] = A.coeff(i, i);
		rhsElems[i] = { b.coeff(i, 0), b.coeff(i, 1), b.coeff(i, 2) };
		offDiagElems[i].reserve(expectedValence);
		for (SparseMatrix::InnerIterator it(A, i); it; ++it)
		{
			if (it.row() <= it.col())
				continue;

			offDiagElems[i].emplace_back(A.coeff(i, it.col()));
		}

		if (offDiagElems[i].size() > maxOffDiagRowLength)
			maxOffDiagRowLength = offDiagElems[i].size();
	}

	// NOTE: all double values will be written to precision 21 and all columns will need to adjust to that.
	constexpr unsigned int streamDblPrecisionWSpace = 21;
	// Table col names header:

	const std::string diagColName = "A[i][i]"; // 7 letters
	const std::string offDiagColName = "A[i][j != i]"; // 12 letters
	//const std::string 
	const std::string spacesDiag = "       "; // 7 spaces
	const std::string rhsColNames = "|       b[i][0]       ,       b[i][1]       ,       b[i][2]       |"; // 3 * 21 letters + 4 delimiters = 67 letters
	const auto nRowIndexDigits = static_cast<unsigned int>(std::floor(log10(A.outerSize()) + 1));

	for (unsigned int i = 0; i < nRowIndexDigits / 2 - (nRowIndexDigits % 2 == 1 ? 0 : 1); i++) fileOStream << " ";
	fileOStream << "i";
	for (unsigned int i = 0; i < nRowIndexDigits / 2; i++) fileOStream << " ";
	fileOStream << "|" << spacesDiag << diagColName << spacesDiag << "|";

	const unsigned int nLettersBeforeAfterOffDiag = (maxOffDiagRowLength == 0 ? 0 : streamDblPrecisionWSpace * (maxOffDiagRowLength / 2) - (offDiagColName.size() / 2));

	for (unsigned int i = 0; i < nLettersBeforeAfterOffDiag; i++) fileOStream << " ";
	fileOStream << offDiagColName;
	for (unsigned int i = 0; i < nLettersBeforeAfterOffDiag; i++) fileOStream << " ";
	fileOStream << rhsColNames << "\n";

	const auto strRowSizeTotal = nRowIndexDigits + 1 + 2 * spacesDiag.size() + 1 + diagColName.size() + maxOffDiagRowLength * streamDblPrecisionWSpace + rhsColNames.size();
	for (unsigned int i = 0; i < strRowSizeTotal; i++) fileOStream << "=";
	fileOStream << "\n";

	// Values:
	// set stream precision
	fileOStream.precision(streamDblPrecisionWSpace - 2);

	for (unsigned int i = 0; i < A.outerSize(); i++)
	{
		std::cout << "DumpMatrixAndRHSToFile: writing row " << i << "\n";
		// index entry
		const auto iDigits = static_cast<unsigned int>(std::floor(log10(i) + 1));
		fileOStream << i;
		const auto nRemainingDigits = (iDigits == nRowIndexDigits ? 0 : nRowIndexDigits - iDigits);
		for (unsigned int j = 0; j < nRemainingDigits; j++) fileOStream << " ";

		// diag value entry
		fileOStream << "|" << diagElems[i] << "|";

		// off diag entries
		unsigned int symbolCount = 0;
		for (const auto& val : offDiagElems[i])
		{
			fileOStream << val << " ";
			symbolCount += streamDblPrecisionWSpace;
		}
		const auto remainingSymbolCount = streamDblPrecisionWSpace * maxOffDiagRowLength - symbolCount;
		for (unsigned int j = 0; j < remainingSymbolCount; j++) fileOStream << " ";
		fileOStream << "|" << rhsElems[i][0] << " " << rhsElems[i][1] << " " << rhsElems[i][2] << "|\n";
	}


	fileOStream.close();
}

//-----------------------------------------------------------------------------
/*! \brief Verifies whether the given line string begins with a token.
 *  \param[in] line     evaluated line string.
 *  \param[in] token    evaluated token.
*/
//-----------------------------------------------------------------------------
static bool LineBeginsWithToken(const std::string& line, const std::string& token)
{
	if (line.empty())
		return false;

	std::stringstream sStream{ line };
	std::string lineToken;
	sStream >> lineToken;
	return lineToken == token;
}

/**
 * \brief Loads line beginning with a given token.
 * \param line           modifiable line ref.
 * \param fileIStream    file stream.
 * \param token          searched token.
 */
void LoadTokenLine(std::string& line, std::ifstream& fileIStream, const std::string& token)
{
	while (!LineBeginsWithToken(line, token))
	{
		if (!std::getline(fileIStream, line))
		{
			fileIStream.close();
			std::cerr << "LoadTokenLine: Unexpected end of file!\n";
			throw std::runtime_error("LoadTokenLine: Unexpected end of file!\n");
		}
	}
}

/// \brief basic utility for verifying whether a given string is numeric.
bool IsNumber(const std::string& str)
{
	for (const auto& c : str)
	{
		if (std::isdigit(c) == 0 && c != '-' && c != '.') return false;
	}
	return true;
}

/**
 * \brief Parses six extent values from a given extent string.
 * \param valueBuffer    buffer for parsed values.
 * \param extentStr      input extent string.
 */
void ParseExtentValues(std::vector<size_t>& valueBuffer, const std::string& extentStr)
{
	if (extentStr.empty())
		return;

	std::stringstream sStream{ extentStr };
	std::string token;
	constexpr size_t extentTokenCount = 6;
	valueBuffer.reserve(extentTokenCount);
	size_t iter = 0;
	do
	{
		iter++;
		sStream >> token;
		if (!IsNumber(token))
			break;

		valueBuffer.emplace_back(std::stoi(token));
		
	} while (iter < extentTokenCount);
}

/// \brief verifies the validity of extent value buffer.
[[nodiscard]] bool AreExtentValuesValid(const std::vector<size_t>& extentValBuffer)
{
	if (extentValBuffer.size() != 6)
		return false;

	if (extentValBuffer[0] > extentValBuffer[1])
		return false;

	if (extentValBuffer[2] > extentValBuffer[3])
		return false;

	return extentValBuffer[4] <= extentValBuffer[5];
}

/// \brief extracts grid dimensions from extent values.
[[nodiscard]] std::tuple<size_t, size_t, size_t> GetGridDimensions(const std::vector<size_t>& extentValues)
{
	assert(AreExtentValuesValid(extentValues));
	return {
		extentValues[1] - extentValues[0] + 1 /**/,
		extentValues[3] - extentValues[2] + 1 /**/,
		extentValues[5] - extentValues[4] + 1  /**/
	};
}

/**
 * \brief Parses a vector of 3 values.
 * \param valueBuffer     target value buffer.
 * \param str             input string.
 */
void Parse3DPointValues(std::vector<float>& valueBuffer, const std::string& str)
{
	if (str.empty())
		return;

	std::stringstream sStream{ str };
	std::string token;
	constexpr size_t pointTokenCount = 3;
	valueBuffer.reserve(pointTokenCount);
	size_t iter = 0;
	do
	{
		iter++;
		sStream >> token;
		if (!IsNumber(token))
			break;

		valueBuffer.emplace_back(std::stof(token));

	} while (iter < pointTokenCount);
}

/// \brief validates point value buffer size.
[[nodiscard]] bool ArePointValuesValid(const std::vector<float>& valueBuffer)
{
	if (valueBuffer.size() != 3)
		return false;

	return true;
}

/// \brief extracts a cell size from the spacing vector. NOTE: currently only regular grids with cube-spacing are supported.
[[nodiscard]] float GetCellSize(const std::vector<float>& spacingVec)
{
	assert(spacingVec.size() == 3);
	if (std::fabs(spacingVec[0] - spacingVec[1]) < FLT_EPSILON && std::fabs(spacingVec[1] - spacingVec[2]) < FLT_EPSILON)
	{
		return spacingVec[0];
	}

	std::cerr << "ImportVTI::GetCellSize [WARNING]: >>>>>>>>>> unequal spacing not supported! <<<<<<<<<<<<< \n";
	std::cerr << "ImportVTI::GetCellSize [WARNING]: unequal spacing values: " << spacingVec[0] << " " << spacingVec[1] << " " << spacingVec[2] << "!\n";
	std::cerr << "ImportVTI::GetCellSize [WARNING]: only the first spacing value will be chosen for cell size! This may result in unevenly scaled fields!\n";
	return spacingVec[0];
}

Geometry::ScalarGrid ImportVTI(const std::string& fileName)
{
	std::ifstream fileIStream(fileName);
	if (!fileIStream.is_open())
	{
		std::cerr << "ImportVTI: file" + fileName + " could not be opened!\n";
		throw std::invalid_argument("ImportVTI: file" + fileName + " could not be opened!\n");
	}

	const auto extensionOrig = fileName.substr(fileName.find(".") + 1);
	std::string extLower = "";
	std::ranges::transform(extensionOrig, std::back_inserter(extLower), std::tolower);
	if (extLower != "vti")
	{
		std::cerr << "ImportVTI: file" << fileName << " could not be opened!\n";
		throw std::invalid_argument("ImportVTI: file" + fileName + " could not be opened!\n");
	}

	// ===================================================
	// >>>>>>>>>>>>>>>> read header <<<<<<<<<<<<<<<<<<<<<<
    // ===================================================
	std::string line;
	LoadTokenLine(line, fileIStream, "<VTKFile");
	//  type=\"ImageData\"
	const auto imageDataId = line.substr(line.find("<VTKFile") + 9, line.find("type=\"ImageData\"") + 7);
	if (imageDataId != "type=\"ImageData\"")
	{
		std::cerr << "ImportVTI: Header part \"type=\"ImageData\"\" not found!\n";
		throw std::runtime_error("ImportVTI: Header part \"type=\"ImageData\"\" not found!\n");
	}

	LoadTokenLine(line, fileIStream, "<ImageData");

	// ========== Data extent (index dimensions) ===========
	const auto extentBeginId = line.find("WholeExtent=\"") + 13;
	const auto extentEndId = line.find("\" Origin=\"");
	const auto strExtent = line.substr(extentBeginId, extentEndId - extentBeginId);
	std::vector<size_t> extentVals{};
	ParseExtentValues(extentVals, strExtent);
	if (!AreExtentValuesValid(extentVals))
	{
		std::cerr << "ImportVTI: Invalid extent!\n";
		throw std::runtime_error("ImportVTI: Invalid extent!\n");
	}
	const auto [Nx, Ny, Nz] = GetGridDimensions(extentVals);
	// ========== Data origin (point) =====================
	const auto originBeginId = extentEndId + 10;
	const auto originEndId = line.find("\" Spacing=\"");
	const auto strOrigin = line.substr(originBeginId, originEndId - originBeginId);
	std::vector<float> originPtCoords{};
	Parse3DPointValues(originPtCoords, strOrigin);
	if (!ArePointValuesValid(originPtCoords))
	{
		std::cerr << "ImportVTI: Invalid origin!\n";
		throw std::runtime_error("ImportVTI: Invalid origin!\n");
	}

	// ========== Spacing (point) =========================
	const auto spacingBeginId = originEndId + 11;
	const auto spacingEndId = line.find("\"", spacingBeginId + 1);
	const auto strSpacing = line.substr(spacingBeginId, spacingEndId - spacingBeginId);
	std::vector<float> cellSpacing{};
	Parse3DPointValues(cellSpacing, strSpacing);
	if (!ArePointValuesValid(cellSpacing))
	{
		std::cerr << "ImportVTI: Invalid spacing!\n";
		throw std::runtime_error("ImportVTI: Invalid spacing!\n");
	}
	const float cellSize = GetCellSize(cellSpacing);
	// ===================================================

	// >>>>>>> compute box <<<<<<<<<<<<<<
	constexpr float boxEpsilon = 1e-6f; // round-off error for ceil/floor in global grid coords.
	const pmp::vec3 boxMinVec{
		originPtCoords[0] + boxEpsilon,// will be floored
		originPtCoords[1] + boxEpsilon,// will be floored
		originPtCoords[2] + boxEpsilon // will be floored
	};
	const pmp::vec3 boxMaxVec{
		originPtCoords[0] + cellSize * (Nx - 1) - boxEpsilon, // will be ceil-ed
		originPtCoords[1] + cellSize * (Ny - 1) - boxEpsilon, // will be ceil-ed
		originPtCoords[2] + cellSize * (Nz - 1) - boxEpsilon  // will be ceil-ed
	};
	const pmp::BoundingBox gridBox(boxMinVec, boxMaxVec);

	// >>>>>>> initialize grid <<<<<<<<<<<<
	Geometry::ScalarGrid result{ cellSize, gridBox };
	// >>>>>>> load values <<<<<<<<<<<<<<<<
	LoadTokenLine(line, fileIStream, "<DataArray");
	auto& resultValues = result.Values();
	const unsigned int gridExtent = (Nx - 1) * (Ny - 1) * (Nz - 1);
	unsigned int gridPos = 0;
	std::string token;
	while (token != "</DataArray>" && !fileIStream.eof())
	{
		// since cell data is not necessarily divided into individual lines, we proceed by streaming the values
		fileIStream >> token;
		if (!IsNumber(token))
			continue;

		if (gridPos > gridExtent - 1)
		{
			std::cerr << "ImportVTI [WARNING]: gridPos > gridExtent - 1! There are more values than gridExtent!\n";
			while (token != "</DataArray>" && !fileIStream.eof()) { fileIStream >> token; gridPos++; } // iterate towards the actual data extent.
			std::cerr << "ImportVTI [WARNING]: omitting " << gridPos - gridExtent - 1 << " tokens.\n";
			fileIStream.close();
			return result;
		}
		resultValues[gridPos] = std::stod(token);
		++gridPos;
	}
	if (gridExtent > gridPos + 1)
	{
		std::cerr << "ImportVTI [WARNING]: gridExtent > gridPos + 1! Not all values are loaded!\n";
		std::cerr << "ImportVTI [WARNING]: omitting " << gridExtent - gridPos << " values.\n";
	}

	fileIStream.close();

	return result;
}

Geometry::BaseMeshGeometryData ImportVTK(const std::string& fileName)
{
	std::ifstream fileIStream(fileName);
	if (!fileIStream.is_open())
	{
		std::cerr << "ImportVTK: file" + fileName + " could not be opened!\n";
		throw std::invalid_argument("ImportVTK: file" + fileName + " could not be opened!\n");
	}

	const auto extensionOrig = fileName.substr(fileName.find(".") + 1);
	std::string extLower = "";
	std::ranges::transform(extensionOrig, std::back_inserter(extLower), std::tolower);
	if (extLower != "vtk")
	{
		std::cerr << "ImportVTK: file" << fileName << " could not be opened!\n";
		throw std::invalid_argument("ImportVTI: file" + fileName + " could not be opened!\n");
	}

	Geometry::BaseMeshGeometryData meshData;
	std::string line;

	// Parse line by line
	while (std::getline(fileIStream, line)) {
		std::istringstream iss(line);
		std::string token;
		iss >> token;

		// Skip comments and empty lines
		if (token.empty() || token[0] == '#')
			continue;

		// Read vertex positions
		if (token == "POINTS") {
			int nPoints;
			iss >> nPoints;
			meshData.Vertices.reserve(nPoints);
			for (int i = 0; i < nPoints; ++i) {
				pmp::vec3 point;
				fileIStream >> point[0] >> point[1] >> point[2];
				meshData.Vertices.push_back(point);
			}
		}

		// Read polygons
		else if (token == "POLYGONS") {
			int nPolygons, totalIndices;
			iss >> nPolygons >> totalIndices;
			for (int i = 0; i < nPolygons; ++i) {
				int nVerts;
				fileIStream >> nVerts;
				std::vector<unsigned int> indices(nVerts);
				for (int j = 0; j < nVerts; ++j) {
					fileIStream >> indices[j];
				}
				meshData.PolyIndices.push_back(indices);
			}
		}
	}

	return meshData;
}

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
	[[nodiscard]] RGBColor InterpolateColors(const RGBColor& c1, const RGBColor& c2, float t)
	{
		return {
			c1.r + t * (c2.r - c1.r),
			c1.g + t * (c2.g - c1.g),
			c1.b + t * (c2.b - c1.b)
		};
	}	
}


void ExportScalarGrid2DToPNG(const std::string& filename, const Geometry::ScalarGrid2D& grid, 
	const std::function<double(const pmp::vec2&, const Geometry::ScalarGrid2D&)>& interpolate,
	float nPixelsPerCellX, float nPixelsPerCellY, const RGBColorScheme& colorMap)
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
	const float cellSize = grid.CellSize();

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
			const float x = orig[0] + static_cast<float>(i) / static_cast<float>(nPixelsPerCellX) * cellSize;
			const float y = orig[1] + static_cast<float>(j) / static_cast<float>(nPixelsPerCellY) * cellSize;

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
			const float t = (range > 0.0) ? static_cast<float>((normalizedValue - lower->first) / range) : 0.0f;
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
	if (!grid.Dimensions().Valid() || grid.CellSize() <= 0.0f)
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
	if (!grid.Dimensions().Valid() || grid.CellSize() <= 0.0f)
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
static [[nodiscard]] FontData LoadFont(const std::string& fontFilePath, float fontSize)
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
	const auto fontSize = paddingX * 1.2f / static_cast<float>(decimalPlaces);  // Dynamic font size, smaller than padding
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
