// Copyright 2011-2021 the Polygon Mesh Processing Library developers.
// Copyright 2001-2005 by Computer Graphics Group, RWTH Aachen
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "pmp/SurfaceMeshIO.h"

#include <clocale>
#include <cstring>
#include <cctype>

#include <algorithm>
#include <map>
#include <fstream>
#include <limits>

#include <rply.h>

// helper function
template <typename T>
void tfread(FILE* in, const T& t)
{
    size_t n_items = fread((char*)&t, 1, sizeof(t), in);
    PMP_ASSERT(n_items > 0);
}

// helper function
template <typename T>
void tfwrite(FILE* out, const T& t)
{
    size_t n_items = fwrite((char*)&t, 1, sizeof(t), out);
    PMP_ASSERT(n_items > 0);
}

namespace pmp {

void SurfaceMeshIO::read(SurfaceMesh& mesh)
{
    std::setlocale(LC_NUMERIC, "C");

    // clear mesh before reading from file
    mesh.clear();

    // extract file extension
    std::string::size_type dot(filename_.rfind("."));
    if (dot == std::string::npos)
        throw IOException("Could not determine file extension!");
    std::string ext = filename_.substr(dot + 1, filename_.length() - dot - 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);

    mesh.name_ = filename_.substr(filename_.find_last_of("/\\") + 1);

    // extension determines reader
    if (ext == "off")
        read_off(mesh);
    else if (ext == "obj")
        read_obj(mesh);
    else if (ext == "vtk")
        read_vtk(mesh);
    else if (ext == "stl")
        read_stl(mesh);
    else if (ext == "ply")
        read_ply(mesh);
    else if (ext == "pmp")
        read_pmp(mesh);
    else if (ext == "xyz")
        read_xyz(mesh);
    else if (ext == "agi")
        read_agi(mesh);
    else
        throw IOException("Could not find reader for " + filename_);
}

void SurfaceMeshIO::write(const SurfaceMesh& mesh)
{
    // extract file extension
    std::string::size_type dot(filename_.rfind("."));
    if (dot == std::string::npos)
        throw IOException("Could not determine file extension!");
    std::string ext = filename_.substr(dot + 1, filename_.length() - dot - 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);

    // extension determines reader
    if (ext == "off")
        write_off(mesh);
    else if (ext == "obj")
        write_obj(mesh);
    else if (ext == "vtk")
        write_vtk(mesh);
    else if (ext == "stl")
        write_stl(mesh);
    else if (ext == "ply")
        write_ply(mesh);
    else if (ext == "pmp")
        write_pmp(mesh);
    else if (ext == "xyz")
        write_xyz(mesh);
    else
        throw IOException("Could not find writer for " + filename_);
}

void SurfaceMeshIO::read_obj(SurfaceMesh& mesh) const
{
    std::array<char, 200> s;
    float x, y, z;
    std::vector<Vertex> vertices;
    std::vector<TexCoord> all_tex_coords; //individual texture coordinates
    std::vector<int>
        halfedge_tex_idx; //texture coordinates sorted for halfedges
    HalfedgeProperty<TexCoord> tex_coords =
        mesh.halfedge_property<TexCoord>("h:tex");
    bool with_tex_coord = false;

    // open file (in ASCII mode)
    FILE* in = fopen(filename_.c_str(), "r");
    if (!in)
        throw IOException("Failed to open file: " + filename_);

    // clear line once
    memset(s.data(), 0, 200);

    // parse line by line (currently only supports vertex positions & faces
    while (in && !feof(in) && fgets(s.data(), 200, in))
    {
        // comment
        if (s[0] == '#' || isspace(s[0]))
            continue;

        // vertex
        else if (strncmp(s.data(), "v ", 2) == 0)
        {
            if (sscanf(s.data(), "v %f %f %f", &x, &y, &z))
            {
                mesh.add_vertex(Point(x, y, z));
            }
        }

        // normal
        else if (strncmp(s.data(), "vn ", 3) == 0)
        {
            if (sscanf(s.data(), "vn %f %f %f", &x, &y, &z))
            {
                // problematic as it can be either a vertex property when interpolated
                // or a halfedge property for hard edges
            }
        }

        // texture coordinate
        else if (strncmp(s.data(), "vt ", 3) == 0)
        {
            if (sscanf(s.data(), "vt %f %f", &x, &y))
            {
                all_tex_coords.emplace_back(x, y);
            }
        }

        // face
        else if (strncmp(s.data(), "f ", 2) == 0)
        {
            int component(0), nv(0);
            bool end_of_vertex(false);
            char *p0, *p1(s.data() + 1);

            vertices.clear();
            halfedge_tex_idx.clear();

            // skip white-spaces
            while (*p1 == ' ')
                ++p1;

            while (p1)
            {
                p0 = p1;

                // overwrite next separator

                // skip '/', '\n', ' ', '\0', '\r' <-- don't forget Windows
                while (*p1 != '/' && *p1 != '\r' && *p1 != '\n' && *p1 != ' ' &&
                       *p1 != '\0')
                    ++p1;

                // detect end of vertex
                if (*p1 != '/')
                {
                    end_of_vertex = true;
                }

                // replace separator by '\0'
                if (*p1 != '\0')
                {
                    *p1 = '\0';
                    p1++; // point to next token
                }

                // detect end of line and break
                if (*p1 == '\0' || *p1 == '\n')
                {
                    p1 = nullptr;
                }

                // read next vertex component
                if (*p0 != '\0')
                {
                    switch (component)
                    {
                        case 0: // vertex
                        {
                            int idx = atoi(p0);
                            if (idx < 0)
                                idx = mesh.n_vertices() + idx + 1;
                            vertices.emplace_back(idx - 1);
                            break;
                        }
                        case 1: // texture coord
                        {
                            int idx = atoi(p0) - 1;
                            halfedge_tex_idx.push_back(idx);
                            with_tex_coord = true;
                            break;
                        }
                        case 2: // normal
                            break;
                    }
                }

                ++component;

                if (end_of_vertex)
                {
                    component = 0;
                    nv++;
                    end_of_vertex = false;
                }
            }

            Face f = mesh.add_face(vertices);

            // add texture coordinates
            if (with_tex_coord && f.is_valid())
            {
                auto h_fit = mesh.halfedges(f);
                auto h_end = h_fit;
                unsigned v_idx = 0;
                do
                {
                    tex_coords[*h_fit] =
                        all_tex_coords.at(halfedge_tex_idx.at(v_idx));
                    ++v_idx;
                    ++h_fit;
                } while (h_fit != h_end);
            }
        }
        // clear line
        memset(s.data(), 0, 200);
    }

    // if there are no textures, delete texture property!
    if (!with_tex_coord)
    {
        mesh.remove_halfedge_property(tex_coords);
    }

    fclose(in);
}

void SurfaceMeshIO::read_vtk(SurfaceMesh& mesh) const
{
    std::array<char, 1024> s;
    FILE* in = fopen(filename_.c_str(), "r");
    if (!in)
        throw IOException("Failed to open file: " + filename_);

    // Clear line once
    memset(s.data(), 0, s.size());

    // Temporary storage for parsing the file
    std::vector<Vertex> vertices;

    // Parse line by line
    while (in && !feof(in) && fgets(s.data(), s.size(), in))
    {
        // Skip comments and empty lines
        if (s[0] == '#' || isspace(s[0]))
            continue;

        // Read vertex positions
        if (strncmp(s.data(), "POINTS", 6) == 0)
        {
            int nPoints;
            sscanf(s.data() + 6, "%d", &nPoints); // Assumes "POINTS" is followed by the number of points
            vertices.reserve(nPoints);
            for (int i = 0; i < nPoints; ++i)
            {
                float x, y, z;
                if (fscanf(in, "%f %f %f", &x, &y, &z) == 3)
                {
                    auto v = mesh.add_vertex(Point(x, y, z));
                    vertices.push_back(v);
                }
            }
        }

        // Read polygons
        else if (strncmp(s.data(), "POLYGONS", 8) == 0)
        {
            int nPolygons, totalIndices;
            sscanf(s.data() + 8, "%d %d", &nPolygons, &totalIndices); // Assumes "POLYGONS" is followed by the number of polygons and total indices

            for (int i = 0; i < nPolygons; ++i)
            {
                int nVerts;
                fscanf(in, "%d", &nVerts); // Read the number of vertices for this polygon

                std::vector<Vertex> faceVertices;
                faceVertices.reserve(nVerts);
                for (int j = 0; j < nVerts; ++j)
                {
                    int idx;
                    fscanf(in, "%d", &idx); // Read vertex index
                    faceVertices.push_back(vertices.at(idx));
                }
                mesh.add_face(faceVertices);
            }
        }

        // Clear line for next read
        memset(s.data(), 0, s.size());
    }

    fclose(in);
}

void SurfaceMeshIO::write_obj(const SurfaceMesh& mesh) const
{
    FILE* out = fopen(filename_.c_str(), "w");
    if (!out)
        throw IOException("Failed to open file: " + filename_);

    // comment
    fprintf(out, "# OBJ export from PMP\n");

    // write vertices
    auto points = mesh.get_vertex_property<Point>("v:point");
    for (auto v : mesh.vertices())
    {
        const Point& p = points[v];
        fprintf(out, "v %.10f %.10f %.10f\n", p[0], p[1], p[2]);
    }

    // write normals
    auto normals = mesh.get_vertex_property<Normal>("v:normal");
    if (normals)
    {
        for (auto v : mesh.vertices())
        {
            const Normal& n = normals[v];
            fprintf(out, "vn %.10f %.10f %.10f\n", n[0], n[1], n[2]);
        }
    }

    // write texture coordinates
    auto tex_coords = mesh.get_halfedge_property<TexCoord>("h:tex");
    if (tex_coords)
    {
        for (auto h : mesh.halfedges())
        {
            const TexCoord& pt = tex_coords[h];
            fprintf(out, "vt %.10f %.10f\n", pt[0], pt[1]);
        }
    }

    // write faces
    for (auto f : mesh.faces())
    {
        fprintf(out, "f");

        auto h = mesh.halfedges(f);
        for (auto v : mesh.vertices(f))
        {
            auto idx = v.idx() + 1;
            if (tex_coords)
            {
                // write vertex index, texCoord index and normal index
                fprintf(out, " %d/%d/%d", idx, (*h).idx() + 1, idx);
                ++h;
            }
            else
            {
                // write vertex index and normal index
                fprintf(out, " %d//%d", idx, idx);
            }
        }
        fprintf(out, "\n");
    }

    fclose(out);
}

//-----------------------------------------------------------------------------
/*! \brief Generates a polydata points header based on point count.
 *  \param[in] mesh           exported mesh.
 *  \return polydata points header string
*/
//-----------------------------------------------------------------------------
[[nodiscard]] std::string GenerateVTKPolydataPointsHeaderFromData(const SurfaceMesh& mesh)
{
    std::string result;
    result += "\nPOINTS " + std::to_string(mesh.n_vertices()) + " double\n";
    return result;
}

//-----------------------------------------------------------------------------
/*! \brief Generates a polydata polygons header based on vertex indices counts.
 *  \param[in] mesh           exported mesh
 *  \return polydata polygons header string
*/
//-----------------------------------------------------------------------------
[[nodiscard]] std::string GenerateVTKPolydataPolygonsHeaderFromData(const SurfaceMesh& mesh)
{
    std::string result;
    size_t vertIdsCount = 0;
    for (const auto f : mesh.faces()) 
    {
        vertIdsCount += std::distance(mesh.vertices(f).begin(), mesh.vertices(f).end()) + 1;
    }

    // each row is (vertIdsCount + 1) entries long
    result += "\nPOLYGONS " + std::to_string(mesh.n_faces()) + " " + std::to_string(vertIdsCount) + "\n";
    return result;
}

constexpr Scalar MAX_ALLOWED_SCALAR_VAL{ 1e+10 };

/// \brief A simple verification whether the vertex property contains values larger than MAX_ALLOWED_SCALAR_VAL
[[nodiscard]] bool HasInvalidValues(const VertexProperty<Scalar>& scalarVProp)
{
    return std::ranges::any_of(scalarVProp.vector(),
        [](const auto val) { return val > MAX_ALLOWED_SCALAR_VAL; });
}

//!> \brief Header string for VTK polydata file.
const auto VTK_Polydata_Header_Str =
	std::string("# vtk DataFile Version 4.2\n") +
	"vtk output\n" +
	"ASCII\n" +
	"DATASET POLYDATA\n";

void SurfaceMeshIO::write_vtk(const SurfaceMesh& mesh) const
{
    FILE* out = fopen(filename_.c_str(), "w");
    if (!out)
        throw IOException("Failed to open file: " + filename_);

    // header
    fprintf(out, "%s", VTK_Polydata_Header_Str.c_str());

    // points header
    fprintf(out, "%s", GenerateVTKPolydataPointsHeaderFromData(mesh).c_str());

    // write vertices
    auto points = mesh.get_vertex_property<Point>("v:point");
    for (const auto v : mesh.vertices())
    {
        const Point& p = points[v];
        fprintf(out, "%.10f %.10f %.10f\n", p[0], p[1], p[2]);
    }

    // polygons header
    fprintf(out, "%s", GenerateVTKPolydataPolygonsHeaderFromData(mesh).c_str());

    // write faces
    for (const auto f : mesh.faces())
    {
        const auto nFaceVerts = static_cast<IndexType>(std::distance(mesh.vertices(f).begin(), mesh.vertices(f).end()));
        fprintf(out, "%d", nFaceVerts);
        for (auto v : mesh.vertices(f))
        {
            const auto idx = v.idx();
            fprintf(out, " %d", idx);
        }
        fprintf(out, "\n");
    }

    if (mesh.vertex_properties().size() > 3)
    {
	    // write additional vertex properties as scalar lookup table

        fprintf(out, "\nPOINT_DATA %d\n", static_cast<int>(mesh.n_vertices()));

        for (const auto& propName : mesh.vertex_properties())
        {
            if (propName == "v:point")
                continue;
            if (propName == "v:connectivity")
                continue;
            if (propName == "v:normal")
                continue;
            if (propName == "v:deleted")
                continue;

            const auto vPropScalar = mesh.get_vertex_property<Scalar>(propName);
            if (!vPropScalar)
                continue;

            if (HasInvalidValues(vPropScalar))
                continue;

            auto propNameToExport = propName.substr(propName.find(":") + 1);
            propNameToExport[0] = std::toupper(propNameToExport[0]); // make first letter upper case

            fprintf(out, "\nSCALARS %s double 1\n", propNameToExport.c_str());
            fprintf(out, "LOOKUP_TABLE default\n");
 
            for (const auto& propVal : vPropScalar.vector())
            {
                fprintf(out, "%f\n", propVal);
            }
        }
    }

    fclose(out);
}

void read_off_ascii(SurfaceMesh& mesh, FILE* in, const bool has_normals,
                    const bool has_texcoords, const bool has_colors)
{
    std::array<char, 1000> line;
    int nc;
    unsigned int i, j, items, idx;
    unsigned int nv, nf, ne;
    float x, y, z, r, g, b;
    Vertex v;

    // properties
    VertexProperty<Normal> normals;
    VertexProperty<TexCoord> texcoords;
    VertexProperty<Color> colors;
    if (has_normals)
        normals = mesh.vertex_property<Normal>("v:normal");
    if (has_texcoords)
        texcoords = mesh.vertex_property<TexCoord>("v:tex");
    if (has_colors)
        colors = mesh.vertex_property<Color>("v:color");

    // #Vertice, #Faces, #Edges
    items = fscanf(in, "%d %d %d\n", (int*)&nv, (int*)&nf, (int*)&ne);
    PMP_ASSERT(items);

    mesh.reserve(nv, std::max(3 * nv, ne), nf);

    // read vertices: pos [normal] [color] [texcoord]
    for (i = 0; i < nv && !feof(in); ++i)
    {
        // read line
        auto lp = fgets(line.data(), 1000, in);
        lp = line.data();

        // position
        items = sscanf(lp, "%f %f %f%n", &x, &y, &z, &nc);
        assert(items == 3);
        v = mesh.add_vertex(Point(x, y, z));
        lp += nc;

        // normal
        if (has_normals)
        {
            if (sscanf(lp, "%f %f %f%n", &x, &y, &z, &nc) == 3)
            {
                normals[v] = Normal(x, y, z);
            }
            lp += nc;
        }

        // color
        if (has_colors)
        {
            if (sscanf(lp, "%f %f %f%n", &r, &g, &b, &nc) == 3)
            {
                if (r > 1.0 || g > 1.0 || b > 1.0)
                {
                    r /= 255.0;
                    g /= 255.0;
                    b /= 255.0;
                }
                colors[v] = Color(r, g, b);
            }
            lp += nc;
        }

        // tex coord
        if (has_texcoords)
        {
            items = sscanf(lp, "%f %f%n", &x, &y, &nc);
            assert(items == 2);
            texcoords[v][0] = x;
            texcoords[v][1] = y;
            lp += nc;
        }
    }

    // read faces: #N v[1] v[2] ... v[n-1]
    std::vector<Vertex> vertices;
    for (i = 0; i < nf; ++i)
    {
        // read line
        auto lp = fgets(line.data(), 1000, in);
        lp = line.data();

        // #vertices
        items = sscanf(lp, "%d%n", (int*)&nv, &nc);
        assert(items == 1);
        vertices.resize(nv);
        lp += nc;

        // indices
        for (j = 0; j < nv; ++j)
        {
            items = sscanf(lp, "%d%n", (int*)&idx, &nc);
            assert(items == 1);
            vertices[j] = Vertex(idx);
            lp += nc;
        }
        mesh.add_face(vertices);
    }
}

void read_off_binary(SurfaceMesh& mesh, FILE* in, const bool has_normals,
                     const bool has_texcoords, const bool has_colors)
{
    IndexType i, j, idx(0);
    IndexType nv(0), nf(0), ne(0);
    Point p, n;
    vec2 t;
    Vertex v;

    // binary cannot (yet) read colors
    if (has_colors)
        throw IOException("Colors not supported for binary OFF file.");

    // properties
    VertexProperty<Normal> normals;
    VertexProperty<TexCoord> texcoords;
    if (has_normals)
        normals = mesh.vertex_property<Normal>("v:normal");
    if (has_texcoords)
        texcoords = mesh.vertex_property<TexCoord>("v:tex");

    // #Vertice, #Faces, #Edges
    tfread(in, nv);
    tfread(in, nf);
    tfread(in, ne);
    mesh.reserve(nv, std::max(3 * nv, ne), nf);

    // read vertices: pos [normal] [color] [texcoord]
    for (i = 0; i < nv && !feof(in); ++i)
    {
        // position
        tfread(in, p);
        v = mesh.add_vertex((Point)p);

        // normal
        if (has_normals)
        {
            tfread(in, n);
            normals[v] = (Normal)n;
        }

        // tex coord
        if (has_texcoords)
        {
            tfread(in, t);
            texcoords[v][0] = t[0];
            texcoords[v][1] = t[1];
        }
    }

    // read faces: #N v[1] v[2] ... v[n-1]
    std::vector<Vertex> vertices;
    for (i = 0; i < nf; ++i)
    {
        tfread(in, nv);
        vertices.resize(nv);
        for (j = 0; j < nv; ++j)
        {
            tfread(in, idx);
            vertices[j] = Vertex(idx);
        }
        mesh.add_face(vertices);
    }
}

void SurfaceMeshIO::write_off_binary(const SurfaceMesh& mesh) const
{
    FILE* out = fopen(filename_.c_str(), "w");
    if (!out)
        throw IOException("Failed to open file: " + filename_);

    fprintf(out, "OFF BINARY\n");
    fclose(out);
    auto nv = (IndexType)mesh.n_vertices();
    auto nf = (IndexType)mesh.n_faces();
    auto ne = IndexType{};

    out = fopen(filename_.c_str(), "ab");
    tfwrite(out, nv);
    tfwrite(out, nf);
    tfwrite(out, ne);
    auto points = mesh.get_vertex_property<Point>("v:point");
    for (auto v : mesh.vertices())
    {
        const Point& p = points[v];
        tfwrite(out, p);
    }

    for (auto f : mesh.faces())
    {
        IndexType valence = mesh.valence(f);
        tfwrite(out, valence);
        for (auto fv : mesh.vertices(f))
            tfwrite(out, (IndexType)fv.idx());
    }
    fclose(out);
}

void SurfaceMeshIO::read_off(SurfaceMesh& mesh) const
{
    std::array<char, 200> line;
    bool has_texcoords = false;
    bool has_normals = false;
    bool has_colors = false;
    bool has_hcoords = false;
    bool has_dim = false;
    bool is_binary = false;

    // open file (in ASCII mode)
    FILE* in = fopen(filename_.c_str(), "r");
    if (!in)
        throw IOException("Failed to open file: " + filename_);

    // read header: [ST][C][N][4][n]OFF BINARY
    auto c = fgets(line.data(), 200, in);
    assert(c != nullptr);
    c = line.data();
    if (c[0] == 'S' && c[1] == 'T')
    {
        has_texcoords = true;
        c += 2;
    }
    if (c[0] == 'C')
    {
        has_colors = true;
        ++c;
    }
    if (c[0] == 'N')
    {
        has_normals = true;
        ++c;
    }
    if (c[0] == '4')
    {
        has_hcoords = true;
        ++c;
    }
    if (c[0] == 'n')
    {
        has_dim = true;
        ++c;
    }
    if (strncmp(c, "OFF", 3) != 0)
    {
        fclose(in);
        throw IOException("Failed to parse OFF header");
    }
    if (strncmp(c + 4, "BINARY", 6) == 0)
        is_binary = true;

    if (has_hcoords)
    {
        fclose(in);
        throw IOException("Error: Homogeneous coordinates not supported.");
    }
    if (has_dim)
    {
        fclose(in);
        throw IOException("Error: vertex dimension != 3 not supported");
    }

    // if binary: reopen file in binary mode
    if (is_binary)
    {
        fclose(in);
        in = fopen(filename_.c_str(), "rb");
        c = fgets(line.data(), 200, in);
        assert(c != nullptr);
    }

    // read as ASCII or binary
    if (is_binary)
        read_off_binary(mesh, in, has_normals, has_texcoords, has_colors);
    else
        read_off_ascii(mesh, in, has_normals, has_texcoords, has_colors);

    fclose(in);
}

void SurfaceMeshIO::write_off(const SurfaceMesh& mesh) const
{
    if (flags_.use_binary)
    {
        write_off_binary(mesh);
        return;
    }

    FILE* out = fopen(filename_.c_str(), "w");
    if (!out)
        throw IOException("Failed to open file: " + filename_);

    bool has_normals = false;
    bool has_texcoords = false;
    bool has_colors = false;

    auto normals = mesh.get_vertex_property<Normal>("v:normal");
    auto texcoords = mesh.get_vertex_property<TexCoord>("v:tex");
    auto colors = mesh.get_vertex_property<Color>("v:color");

    if (normals && flags_.use_vertex_normals)
        has_normals = true;
    if (texcoords && flags_.use_vertex_texcoords)
        has_texcoords = true;
    if (colors && flags_.use_vertex_colors)
        has_colors = true;

    // header
    if (has_texcoords)
        fprintf(out, "ST");
    if (has_colors)
        fprintf(out, "C");
    if (has_normals)
        fprintf(out, "N");
    fprintf(out, "OFF\n%zu %zu 0\n", mesh.n_vertices(), mesh.n_faces());

    // vertices, and optionally normals and texture coordinates
    VertexProperty<Point> points = mesh.get_vertex_property<Point>("v:point");
    for (auto v : mesh.vertices())
    {
        const Point& p = points[v];
        fprintf(out, "%.10f %.10f %.10f", p[0], p[1], p[2]);

        if (has_normals)
        {
            const Normal& n = normals[v];
            fprintf(out, " %.10f %.10f %.10f", n[0], n[1], n[2]);
        }

        if (has_colors)
        {
            const Color& c = colors[v];
            fprintf(out, " %.10f %.10f %.10f", c[0], c[1], c[2]);
        }

        if (has_texcoords)
        {
            const TexCoord& t = texcoords[v];
            fprintf(out, " %.10f %.10f", t[0], t[1]);
        }

        fprintf(out, "\n");
    }

    // faces
    for (auto f : mesh.faces())
    {
        auto nv = mesh.valence(f);
        fprintf(out, "%zu", nv);
        auto fv = mesh.vertices(f);
        auto fvend = fv;
        do
        {
            fprintf(out, " %d", (*fv).idx());
        } while (++fv != fvend);
        fprintf(out, "\n");
    }

    fclose(out);
}

void SurfaceMeshIO::read_pmp(SurfaceMesh& mesh) const
{
    // open file (in binary mode)
    FILE* in = fopen(filename_.c_str(), "rb");
    if (!in)
        throw IOException("Failed to open file: " + filename_);

    // how many elements?
    unsigned int nv(0), ne(0), nh(0), nf(0);
    tfread(in, nv);
    tfread(in, ne);
    tfread(in, nf);
    nh = 2 * ne;

    // texture coordinates?
    bool has_htex(false);
    tfread(in, has_htex);

    // resize containers
    mesh.vprops_.resize(nv);
    mesh.hprops_.resize(nh);
    mesh.eprops_.resize(ne);
    mesh.fprops_.resize(nf);

    // get properties
    auto vconn =
        mesh.vertex_property<SurfaceMesh::VertexConnectivity>("v:connectivity");
    auto hconn = mesh.halfedge_property<SurfaceMesh::HalfedgeConnectivity>(
        "h:connectivity");
    auto fconn =
        mesh.face_property<SurfaceMesh::FaceConnectivity>("f:connectivity");
    auto point = mesh.vertex_property<Point>("v:point");

    // read properties from file
    size_t nvc = fread((char*)vconn.data(),
                       sizeof(SurfaceMesh::VertexConnectivity), nv, in);
    size_t nhc = fread((char*)hconn.data(),
                       sizeof(SurfaceMesh::HalfedgeConnectivity), nh, in);
    size_t nfc = fread((char*)fconn.data(),
                       sizeof(SurfaceMesh::FaceConnectivity), nf, in);
    size_t np = fread((char*)point.data(), sizeof(Point), nv, in);
    PMP_ASSERT(nvc == nv);
    PMP_ASSERT(nhc == nh);
    PMP_ASSERT(nfc == nf);
    PMP_ASSERT(np == nv);

    // read texture coordiantes
    if (has_htex)
    {
        auto htex = mesh.halfedge_property<TexCoord>("h:tex");
        size_t nhtc = fread((char*)htex.data(), sizeof(TexCoord), nh, in);
        PMP_ASSERT(nhtc == nh);
    }

    fclose(in);
}

void SurfaceMeshIO::read_xyz(SurfaceMesh& mesh) const
{
    // open file (in ASCII mode)
    FILE* in = fopen(filename_.c_str(), "r");
    if (!in)
        throw IOException("Failed to open file: " + filename_);

    // add normal property
    // \todo this adds property even if no normals present. change it.
    auto vnormal = mesh.vertex_property<Normal>("v:normal");

    std::array<char, 200> line;
    float x, y, z;
    float nx, ny, nz;
    int n;
    Vertex v;

    // read data
    while (in && !feof(in) && fgets(line.data(), 200, in))
    {
        n = sscanf(line.data(), "%f %f %f %f %f %f", &x, &y, &z, &nx, &ny, &nz);
        if (n >= 3)
        {
            v = mesh.add_vertex(Point(x, y, z));
            if (n >= 6)
            {
                vnormal[v] = Normal(nx, ny, nz);
            }
        }
    }

    fclose(in);
}

// \todo remove duplication with read_xyz
void SurfaceMeshIO::read_agi(SurfaceMesh& mesh) const
{
    // open file (in ASCII mode)
    FILE* in = fopen(filename_.c_str(), "r");
    if (!in)
        throw IOException("Failed to open file: " + filename_);

    // add normal property
    auto normal = mesh.vertex_property<Normal>("v:normal");
    auto color = mesh.vertex_property<Color>("v:color");

    std::array<char, 200> line;
    float x, y, z;
    float nx, ny, nz;
    float r, g, b;
    int n;
    Vertex v;

    // read data
    while (in && !feof(in) && fgets(line.data(), 200, in))
    {
        n = sscanf(line.data(), "%f %f %f %f %f %f %f %f %f", &x, &y, &z, &r,
                   &g, &b, &nx, &ny, &nz);
        if (n == 9)
        {
            v = mesh.add_vertex(Point(x, y, z));
            normal[v] = Normal(nx, ny, nz);
            color[v] = Color(r / 255.0, g / 255.0, b / 255.0);
        }
    }

    fclose(in);
}

void SurfaceMeshIO::write_pmp(const SurfaceMesh& mesh) const
{
    // open file (in binary mode)
    FILE* out = fopen(filename_.c_str(), "wb");
    if (!out)
        throw IOException("Failed to open file: " + filename_);

    // get properties
    auto vconn = mesh.get_vertex_property<SurfaceMesh::VertexConnectivity>(
        "v:connectivity");
    auto hconn = mesh.get_halfedge_property<SurfaceMesh::HalfedgeConnectivity>(
        "h:connectivity");
    auto fconn =
        mesh.get_face_property<SurfaceMesh::FaceConnectivity>("f:connectivity");
    auto point = mesh.get_vertex_property<Point>("v:point");
    auto htex = mesh.get_halfedge_property<TexCoord>("h:tex");

    // how many elements?
    unsigned int nv, ne, nh, nf;
    nv = mesh.n_vertices();
    ne = mesh.n_edges();
    nh = mesh.n_halfedges();
    nf = mesh.n_faces();

    // write header
    tfwrite(out, nv);
    tfwrite(out, ne);
    tfwrite(out, nf);
    tfwrite(out, (bool)htex);

    // write properties to file
    fwrite((char*)vconn.data(), sizeof(SurfaceMesh::VertexConnectivity), nv,
           out);
    fwrite((char*)hconn.data(), sizeof(SurfaceMesh::HalfedgeConnectivity), nh,
           out);
    fwrite((char*)fconn.data(), sizeof(SurfaceMesh::FaceConnectivity), nf, out);
    fwrite((char*)point.data(), sizeof(Point), nv, out);

    // texture coordinates
    if (htex)
        fwrite((char*)htex.data(), sizeof(TexCoord), nh, out);

    fclose(out);
}

// helper to assemble vertex data
static int vertexCallback(p_ply_argument argument)
{
    long idx;
    void* pdata;
    ply_get_argument_user_data(argument, &pdata, &idx);

    auto* mesh = (pmp::SurfaceMesh*)pdata;
    auto point = mesh->get_object_property<pmp::Point>("g:point");
    point[0][idx] = ply_get_argument_value(argument);

    if (idx == 2)
        mesh->add_vertex(point[0]);

    return 1;
}

// helper to assemble face data
static int faceCallback(p_ply_argument argument)
{
    long length, value_index;
    void* pdata;
    long idata;
    ply_get_argument_user_data(argument, &pdata, &idata);
    ply_get_argument_property(argument, nullptr, &length, &value_index);

    auto* mesh = (pmp::SurfaceMesh*)pdata;
    auto vertices =
        mesh->get_object_property<std::vector<pmp::Vertex>>("g:vertices");

    if (value_index == 0)
        vertices[0].clear();

    auto idx = (IndexType)ply_get_argument_value(argument);
    vertices[0].push_back(pmp::Vertex(idx));

    if (value_index == length - 1)
        mesh->add_face(vertices[0]);

    return 1;
}

void SurfaceMeshIO::read_ply(SurfaceMesh& mesh) const
{
    // add object properties to hold temporary data
    auto point = mesh.add_object_property<Point>("g:point");
    auto vertices = mesh.add_object_property<std::vector<Vertex>>("g:vertices");

    // open file, read header
    p_ply ply = ply_open(filename_.c_str(), nullptr, 0, nullptr);

    if (!ply)
        throw IOException("Failed to open file: " + filename_);

    if (!ply_read_header(ply))
        throw IOException("Failed to read PLY header!");

    // setup callbacks for basic properties
    ply_set_read_cb(ply, "vertex", "x", vertexCallback, &mesh, 0);
    ply_set_read_cb(ply, "vertex", "y", vertexCallback, &mesh, 1);
    ply_set_read_cb(ply, "vertex", "z", vertexCallback, &mesh, 2);

    ply_set_read_cb(ply, "face", "vertex_indices", faceCallback, &mesh, 0);

    // read the data
    if (!ply_read(ply))
        throw IOException("Failed to read PLY data!");

    ply_close(ply);

    // clean-up properties
    mesh.remove_object_property(point);
    mesh.remove_object_property(vertices);
}

void SurfaceMeshIO::write_ply(const SurfaceMesh& mesh) const
{
    e_ply_storage_mode mode = flags_.use_binary ? PLY_LITTLE_ENDIAN : PLY_ASCII;
    p_ply ply = ply_create(filename_.c_str(), mode, nullptr, 0, nullptr);

    ply_add_comment(ply, "File written with pmp-library");
    ply_add_element(ply, "vertex", mesh.n_vertices());
    ply_add_scalar_property(ply, "x", PLY_FLOAT);
    ply_add_scalar_property(ply, "y", PLY_FLOAT);
    ply_add_scalar_property(ply, "z", PLY_FLOAT);
    ply_add_element(ply, "face", mesh.n_faces());
    ply_add_property(ply, "vertex_indices", PLY_LIST, PLY_UCHAR, PLY_INT);
    ply_write_header(ply);

    // write vertices
    auto points = mesh.get_vertex_property<Point>("v:point");
    for (auto v : mesh.vertices())
    {
        ply_write(ply, points[v][0]);
        ply_write(ply, points[v][1]);
        ply_write(ply, points[v][2]);
    }

    // write faces
    for (auto f : mesh.faces())
    {
        ply_write(ply, mesh.valence(f));
        for (auto fv : mesh.vertices(f))
            ply_write(ply, fv.idx());
    }

    ply_close(ply);
}

// helper class for STL reader
class CmpVec
{
public:
    CmpVec(Scalar eps = std::numeric_limits<Scalar>::min()) : eps_(eps) {}

    bool operator()(const vec3& v0, const vec3& v1) const
    {
        if (fabs(v0[0] - v1[0]) <= eps_)
        {
            if (fabs(v0[1] - v1[1]) <= eps_)
            {
                return (v0[2] < v1[2] - eps_);
            }
            else
                return (v0[1] < v1[1] - eps_);
        }
        else
            return (v0[0] < v1[0] - eps_);
    }

private:
    Scalar eps_;
};

void SurfaceMeshIO::read_stl(SurfaceMesh& mesh) const
{
    std::array<char, 100> line;
    unsigned int i, nT(0);
    vec3 p;
    Vertex v;
    std::vector<Vertex> vertices(3);
    size_t n_items(0);

    CmpVec comp(std::numeric_limits<Scalar>::min());
    std::map<vec3, Vertex, CmpVec> vMap(comp);
    std::map<vec3, Vertex, CmpVec>::iterator vMapIt;

    // open file (in ASCII mode)
    FILE* in = fopen(filename_.c_str(), "r");
    if (!in)
        throw IOException("Failed to open file: " + filename_);

    // ASCII or binary STL?
    auto c = fgets(line.data(), 6, in);
    PMP_ASSERT(c != nullptr);
    const bool binary = ((strncmp(line.data(), "SOLID", 5) != 0) &&
                         (strncmp(line.data(), "solid", 5) != 0));

    // parse binary STL
    if (binary)
    {
        // re-open file in binary mode
        fclose(in);
        in = fopen(filename_.c_str(), "rb");
        if (!in)
            throw IOException("Failed to open file: " + filename_);

        // skip dummy header
        n_items = fread(line.data(), 1, 80, in);
        PMP_ASSERT(n_items > 0);

        // read number of triangles
        tfread(in, nT);

        // read triangles
        while (nT)
        {
            // skip triangle normal
            n_items = fread(line.data(), 1, 12, in);
            PMP_ASSERT(n_items > 0);
            // triangle's vertices
            for (i = 0; i < 3; ++i)
            {
                tfread(in, p);

                // has vector been referenced before?
                if ((vMapIt = vMap.find(p)) == vMap.end())
                {
                    // No : add vertex and remember idx/vector mapping
                    v = mesh.add_vertex((Point)p);
                    vertices[i] = v;
                    vMap[p] = v;
                }
                else
                {
                    // Yes : get index from map
                    vertices[i] = vMapIt->second;
                }
            }

            // Add face only if it is not degenerated
            if ((vertices[0] != vertices[1]) && (vertices[0] != vertices[2]) &&
                (vertices[1] != vertices[2]))
                mesh.add_face(vertices);

            n_items = fread(line.data(), 1, 2, in);
            PMP_ASSERT(n_items > 0);
            --nT;
        }
    }

    // parse ASCII STL
    else
    {
        // parse line by line
        while (in && !feof(in) && fgets(line.data(), 100, in))
        {
            // skip white-space
            for (c = line.data(); isspace(*c) && *c != '\0'; ++c)
            {
            };

            // face begins
            if ((strncmp(c, "outer", 5) == 0) || (strncmp(c, "OUTER", 5) == 0))
            {
                // read three vertices
                for (i = 0; i < 3; ++i)
                {
                    // read line
                    c = fgets(line.data(), 100, in);
                    PMP_ASSERT(c != nullptr);

                    // skip white-space
                    for (c = line.data(); isspace(*c) && *c != '\0'; ++c)
                    {
                    };

                    // read x, y, z
#if PMP_SCALAR_TYPE_64
                    sscanf(c + 6, "%lf %lf %lf", &p[0], &p[1], &p[2]);
#else
                    sscanf(c + 6, "%f %f %f", &p[0], &p[1], &p[2]);
#endif

                    // has vector been referenced before?
                    if ((vMapIt = vMap.find(p)) == vMap.end())
                    {
                        // No : add vertex and remember idx/vector mapping
                        v = mesh.add_vertex((Point)p);
                        vertices[i] = v;
                        vMap[p] = v;
                    }
                    else
                    {
                        // Yes : get index from map
                        vertices[i] = vMapIt->second;
                    }
                }

                // Add face only if it is not degenerated
                if ((vertices[0] != vertices[1]) &&
                    (vertices[0] != vertices[2]) &&
                    (vertices[1] != vertices[2]))
                    mesh.add_face(vertices);
            }
        }
    }

    fclose(in);
}

void SurfaceMeshIO::write_stl(const SurfaceMesh& mesh) const
{
    if (!mesh.is_triangle_mesh())
    {
        auto what = "SurfaceMeshIO::write_stl: Not a triangle mesh.";
        throw InvalidInputException(what);
    }

    auto fnormals = mesh.get_face_property<Normal>("f:normal");
    if (!fnormals)
    {
        auto what = "SurfaceMeshIO::write_stl: No face normals present.";
        throw InvalidInputException(what);
    }

    std::ofstream ofs(filename_.c_str());
    auto points = mesh.get_vertex_property<Point>("v:point");

    ofs << "solid stl" << std::endl;
    Normal n;
    Point p;

    for (auto f : mesh.faces())
    {
        n = fnormals[f];
        ofs << "  facet normal ";
        ofs << n[0] << " " << n[1] << " " << n[2] << std::endl;
        ofs << "    outer loop" << std::endl;
        for (auto v : mesh.vertices(f))
        {
            p = points[v];
            ofs << "      vertex ";
            ofs << p[0] << " " << p[1] << " " << p[2] << std::endl;
        }
        ofs << "    endloop" << std::endl;
        ofs << "  endfacet" << std::endl;
    }
    ofs << "endsolid" << std::endl;
    ofs.close();
}

void SurfaceMeshIO::write_xyz(const SurfaceMesh& mesh) const
{
    std::ofstream ofs(filename_);
    if (!ofs)
        throw IOException("Failed to open file: " + filename_);

    auto vnormal = mesh.get_vertex_property<Normal>("v:normal");
    for (auto v : mesh.vertices())
    {
        ofs << mesh.position(v);
        ofs << " ";
        if (vnormal)
        {
            ofs << vnormal[v];
        }
        ofs << std::endl;
    }

    ofs.close();
}

} // namespace pmp
