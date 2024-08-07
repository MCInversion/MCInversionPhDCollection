#include "ConversionUtils.h"

#include <fstream>

#include "geometry/GridUtil.h"
#include "utils/StringUtils.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image/stb_image_write.h"

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

void ExportScalarGrid2DToPNG(const std::string& filename, const Geometry::ScalarGrid2D& grid, float nPixelsPerCellX, float nPixelsPerCellY)
{
	if (!grid.IsValid())
	{
		throw std::invalid_argument("ExportScalarGrid2DToPNG: grid to be exported is invalid!\n");
	}
	if (filename.empty())
	{
		throw std::invalid_argument("ExportScalarGrid2DToPNG: filename cannot be empty!\n");
	}

	const auto& dims = grid.Dimensions();
	const int nx = static_cast<int>(dims.Nx);
	const int ny = static_cast<int>(dims.Ny);

	const int imageWidth = static_cast<int>(nx * nPixelsPerCellX);
	const int imageHeight = static_cast<int>(ny * nPixelsPerCellY);

	std::vector<uint8_t> image(imageWidth * imageHeight * 3);

	// Normalize the scalar values to [0, 255]
	auto minMax = std::minmax_element(grid.Values().begin(), grid.Values().end());
	double minVal = *minMax.first;
	double maxVal = *minMax.second;

	for (int i = 0; i < imageWidth; ++i)
	{
		for (int j = 0; j < imageHeight; ++j)
		{
			// Map the pixel position to the grid position
			int gridX = static_cast<int>(i / nPixelsPerCellX);
			int gridY = static_cast<int>(j / nPixelsPerCellY);

			// Ensure we do not go out of bounds
			gridX = std::min(gridX, nx - 1);
			gridY = std::min(gridY, ny - 1);

			// Get the corresponding grid value
			double value = grid.Values()[gridX + gridY * nx];
			uint8_t intensity = static_cast<uint8_t>(255.0 * (value - minVal) / (maxVal - minVal));

			// Set the pixel value (RGB)
			image[3 * (i + j * imageWidth) + 0] = intensity; // R
			image[3 * (i + j * imageWidth) + 1] = intensity; // G
			image[3 * (i + j * imageWidth) + 2] = intensity; // B
		}
	}

	// Save the image to a file
	if (!stbi_write_png(filename.c_str(), imageWidth, imageHeight, 3, image.data(), imageWidth * 3))
	{
		throw std::runtime_error("ExportToPNG: Failed to save image to file!");
	}
}
