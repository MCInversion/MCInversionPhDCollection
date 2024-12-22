// Copyright 2011-2021 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#pragma once

#include "pmp/Config.h"
#include "pmp/MatVec.h"

//! \def PMP_ASSERT(x)
//! Custom assert macro that allows to silence unused variable warnings with no
//! overhead. Generates no code in release mode since if the argument to
//! sizeof() is an expression it is not evaluated. In debug mode we just fall
//! back to the default assert().
#ifdef NDEBUG
#define PMP_ASSERT(x)    \
    do                   \
    {                    \
        (void)sizeof(x); \
    } while (0)
#else
#define PMP_ASSERT(x) assert(x)
#endif

//! The pmp-library namespace
namespace pmp {

//! \addtogroup core
//! @{

//! Overriding scalar type aliases
#if PMP_SCALAR_TYPE_64
    using vec2 = dvec2;  // Use double-based 2D vector
    using vec3 = dvec3;  // Use double-based 3D vector
#else
    using vec2 = Vector<float, 2>;  // Use float-based 2D vector
    using vec3 = Vector<float, 3>;  // Use float-based 3D vector
#endif

//! Point type
using Point = vec3;
using Point2 = vec2;

//! Normal type
using Normal = vec3;
using Normal2 = vec2;

//! Color type
//! \details RGB values in the range of [0,1]
using Color = vec3;

//! Texture coordinate type
using TexCoord = vec2;

//! Common IO flags for reading and writing
struct IOFlags
{
    bool use_binary = false;             //!< read / write binary format
    bool use_vertex_normals = false;     //!< read / write vertex normals
    bool use_vertex_colors = false;      //!< read / write vertex colors
    bool use_vertex_texcoords = false;   //!< read / write vertex texcoords
    bool use_face_normals = false;       //!< read / write face normals
    bool use_face_colors = false;        //!< read / write face colors
    bool use_halfedge_texcoords = false; //!< read / write halfedge texcoords
};

//! @}

//! \defgroup core core
//! \brief Core data structure and utilities.

//! \defgroup algorithms algorithms
//! \brief Geometry processing algorithms.

//! \defgroup visualization visualization
//! \brief Visualization tools using OpenGL.

} // namespace pmp
