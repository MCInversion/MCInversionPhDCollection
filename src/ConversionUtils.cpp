#include "ConversionUtils.h"

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
