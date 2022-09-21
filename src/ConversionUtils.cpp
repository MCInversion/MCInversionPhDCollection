#include "ConversionUtils.h"

#include <fstream>

void ExportToVTI(const std::string& filename, const Geometry::ScalarGrid& scalarGrid)
{
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
