#pragma once

#include <vector>
#include "BoundaryElement.h"

namespace MAT
{
	class BoundaryGenerator 
	{
	public:
		// returns the boundaryElements vector
		std::vector<BoundaryElement> getBoundaryElements(unsigned char shape);

		// Overload to return boundary elements for custom shapes
		std::vector<BoundaryElement> getBoundaryElementsFromCurve(const std::vector<Vector2d>& customShapePoints);
	private:
		// populates the boundaryElements vector with edges and concave vertices
		void generateShape(const Vector2d *shape, int vertexCount);

		// member variable
		std::vector<BoundaryElement> boundaryElements;
	};
} // namespace MAT