#include "ConversionUtils.h"

#include "newimage/newimageio.h"

#include <fstream>

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

Geometry::ScalarGrid ImportNiftiAsScalarGrid(const std::string& fileName)
{
	// bet2 implementation of nifti I/O
	NEWIMAGE::volume<float> vol{};
	NEWIMAGE::read_volume(vol, fileName.c_str());

	// TODO: xdim, ydim and zdim are cell dimensions. Now since ScalarGrid does not support non-cube voxels something has to be done about it

	assert(std::fabs(vol.xdim() - vol.ydim()) < FLT_EPSILON && std::fabs(vol.zdim() - vol.ydim()) < FLT_EPSILON);
	assert(vol.xdim() > FLT_EPSILON && vol.ydim() > FLT_EPSILON && vol.zdim() > FLT_EPSILON);
	const float cellSize = vol.xdim(); // assuming cube voxels

	const size_t Nx = vol.xsize();
	const size_t Ny = vol.ysize();
	const size_t Nz = vol.zsize();

	const pmp::vec3 volBoxHalfSize{ cellSize * Nx / 2.0f, cellSize * Ny / 2.0f, cellSize * Nz / 2.0f };

	// because we can't extract NEWIMAGE::volume's origin (it's somehow in integer coordinates)
	// the grid box will be centered at (0,0,0):
	const pmp::BoundingBox volBox(-volBoxHalfSize, volBoxHalfSize);

	Geometry::ScalarGrid result(cellSize, volBox, 0.0);
	auto& resultValues = result.Values();

	for (unsigned int iz = 0; iz < Nz; iz++) {
		for (unsigned int iy = 0; iy < Ny; iy++) {
			for (unsigned int ix = 0; ix < Nx; ix++) {
				const unsigned int gridPos = Nx * Ny * iz + Nx * iy + ix;
				resultValues[gridPos] = vol.value(ix, iy, iz);
			}
		}
	}

	return result;
}
