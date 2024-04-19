/*
======================================================================
Copyright 2018, The ilastik development team

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided
with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived
from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.
======================================================================
 */

#pragma once

#include <map>
#include <cstddef>

namespace MarchingCubes
{
	typedef float Point[3];

	///	==========================================================================
	/// \brief Point identifier
	///	\struct IdPoint
	///	==========================================================================
	struct IdPoint
	{
		size_t id;
		float x, y, z;
	};

	typedef std::map<size_t, IdPoint> PointIdMapping;

	///	==========================================================================
	/// \brief Triangle point indices.
	///	\struct Triangle
	///	==========================================================================
	struct Triangle
	{
		size_t pointId[3];
	};

	///	==========================================================================
	/// \brief the primary struct used to pass around the components of a mesh
	///	\struct MC_Mesh
	///	==========================================================================
	struct MC_Mesh
	{
		size_t vertexCount; //! the number of vertices/normals
		Point* vertices; //! the vertex positions as an array of points
		Point* normals; //! the normal direction of each vertex as an array of points
		size_t faceCount; //! the number of faces
		size_t* faces; //! the faces given by 3 vertex indices (length = faceCount * 3)

		MC_Mesh(size_t, Point*, Point*, size_t, size_t*);
		MC_Mesh() = default;
	};

	/**
	 * \brief The marching cubes algorithm as described here: http://paulbourke.net/geometry/polygonise/
	 * \param volume      contains the data (size = xDim * yDim * zDim).
	 * \param xDim        the x dimension of the grid.
	 * \param yDim        the y dimension of the grid.
	 * \param zDim        the z dimension of the grid.
	 * \param isoLevel    the minimum isoLevel, all values >= isoLevel will contribute to the mesh.
	 * \return the mesh is returned, the caller takes ownership over the pointers.
	 */
	template<typename T>
	MC_Mesh GetMarchingCubesMesh(const T* volume, size_t xDim, size_t yDim, size_t zDim, T isoLevel);
	
} // namespace MarchingCubes


