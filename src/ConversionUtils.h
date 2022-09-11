#pragma once

#include "sdfgen/vec.h"
#include "sdfgen/array3.h"

#include "pmp/SurfaceMesh.h"

/**
 * \brief Converts polygons from pmp::SurfaceMesh to triangles.
 * \param mesh     Input surface mesh.
 * \return a triangulated pmp::SurfaceMesh.
 */
[[nodiscard]] pmp::SurfaceMesh GetTriangulatedSurfaceMesh(const pmp::SurfaceMesh& mesh);

/**
 * \brief Extracts triangulated polygons from pmp::SurfaceMesh as a Vec3ui array of vertex index triples.
 * \param mesh     Input surface mesh.
 * \return a std::vector<Vec3ui> array of vertex index triples.
 */
[[nodiscard]] std::vector<Vec3ui> GetTriangleIndices(const pmp::SurfaceMesh& mesh);

/**
 * \brief Extracts vertex positions from pmp::SurfaceMesh.
 * \param mesh     Input surface mesh.
 * \return a std::vector<Vec3f> array of vertex positions.
 */
[[nodiscard]] std::vector<Vec3f> GetVertexPositions(const pmp::SurfaceMesh& mesh);

/**
 * \brief Exports scalar field to .vti format.
 * \param filename     in output path
 * \param origin       grid origin point
 * \param dx           cell size
 * \param nx           number of cells in x
 * \param ny           number of cells in y
 * \param nz           number of cells in z
 * \param values       array of values
 */
void ExportToVTI(const std::string& filename, const Vec3f& origin, float dx, int nx, int ny, int nz, Array3f& values);