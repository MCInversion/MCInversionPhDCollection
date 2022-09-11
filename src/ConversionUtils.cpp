#include "ConversionUtils.h"

#include <fstream>

pmp::SurfaceMesh GetTriangulatedSurfaceMesh(const pmp::SurfaceMesh& mesh)
{
    pmp::SurfaceMesh result(mesh);

    // TODO: Use Poly2Tri

    for (const auto f : result.faces())
    {
        if (result.valence(f) == 3)
            continue;

        const auto vBegin = *result.vertices(f).begin();
        result.split(f, vBegin);
    }

    return result;
}

std::vector<Vec3ui> GetTriangleIndices(const pmp::SurfaceMesh& mesh)
{
    assert(mesh.is_triangle_mesh());
    const size_t nTris = mesh.faces_size();
    std::vector<Vec3ui> result(nTris);

    for (size_t i = 0; i < nTris; i++)
    {
        const pmp::Face f(i);
        std::vector<unsigned int> vIds;
        vIds.reserve(3);
        for (const auto v : mesh.vertices(f))
            vIds.emplace_back(v.idx());
        result[i] = Vec3ui(vIds[0], vIds[1], vIds[2]);
    }

    return result;
}

std::vector<Vec3f> GetVertexPositions(const pmp::SurfaceMesh& mesh)
{
    const size_t nVerts = mesh.vertices_size();
    std::vector<Vec3f> result(nVerts);
    const auto vpointPropData = mesh.get_vertex_property<pmp::Point>("v:point").data();

    for (size_t i = 0; i < nVerts; i++)
    {
        const auto pos = vpointPropData[i];
        result[i] = Vec3f(pos[0], pos[1], pos[2]);
    }

    return result;
}

void ExportToVTI(const std::string& filename, const Vec3f& origin, float dx, int nx, int ny, int nz, Array3f& values)
{
    std::fstream vti(filename + ".vti", std::fstream::out);
    const Vec3f min = origin;
    const Vec3f max = origin + Vec3f(nx, ny, nz) * dx;

    vti << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    vti << "	<ImageData WholeExtent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1 << "\" Origin=\"" << origin[0] + 0.5 * dx << " " << origin[1] + 0.5 * dx << " " << origin[2] + 0.5 * dx << "\" Spacing=\"" << dx << " " << dx << " " << dx << "\">\n";
    vti << "		<Piece Extent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1 << "\">\n";
    vti << "			<PointData Scalars=\"Scalars_\">\n";
    vti << "				<DataArray type=\"Float32\" Name=\"Scalars_\" format=\"ascii\" RangeMin=\"" << min << "\" RangeMax=\"" << max << "\">\n";

    for (const auto& val : values) {
        vti << val << "\n";
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
