#include "GeometryUtil.h"

#include "pmp/MatVec.h"

namespace Geometry
{
	// ====== Helper macros for vectors (to increase speed) ============

	#define CROSS(dest, v1, v2)						     \
	          dest[0] = v1[1] * v2[2] - v1[2] * v2[1];   \
	          dest[1] = v1[2] * v2[0] - v1[0] * v2[2];	 \
	          dest[2] = v1[0] * v2[1] - v1[1] * v2[0]

	#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

	#define SUB(dest, v1, v2)						\
	          dest[0] = v1[0] - v2[0];				\
	          dest[1] = v1[1] - v2[1];				\
	          dest[2] = v1[2] - v2[2]

	#define FINDMINMAX(x0, x1, x2, min, max)		\
			  min = max = x0;						\
			  if (x1 < min) min = x1;				\
			  if (x1 > max) max = x1;				\
			  if (x2 < min) min = x2;				\
			  if (x2 > max) max = x2

	double GetDistanceToTriangleSq(const std::vector<pmp::vec3>& vertices, const pmp::vec3& point)
	{
		assert(vertices.size() == 3); // only vertex triples allowed

		pmp::vec3 diff = point - vertices[0];
		pmp::vec3 edge0 = vertices[1] - vertices[0];
		pmp::vec3 edge1 = vertices[2] - vertices[0];
		const double a00 = DOT(edge0, edge0);
		const double a01 = DOT(edge0, edge1);
		const double a11 = DOT(edge1, edge1);
		const double b0 = -DOT(diff, edge0);
		const double b1 = -DOT(diff, edge1);
		constexpr double zero = 0.0;
		constexpr double one = 1.0;
		const double det = a00 * a11 - a01 * a01;
		double t0 = a01 * b1 - a11 * b0;
		double t1 = a01 * b0 - a00 * b1;

		if (t0 + t1 <= det) {
			if (t0 < zero) {
				if (t1 < zero) { // region 4			
					if (b0 < zero) {
						t1 = zero;
						if (-b0 >= a00) { // V1					
							t0 = one;
						}
						else { // E01					
							t0 = -b0 / a00;
						}
					}
					else {
						t0 = zero;
						if (b1 >= zero) { // V0					
							t1 = zero;
						}
						else if (-b1 >= a11) { // V2					
							t1 = one;
						}
						else { // E20					
							t1 = -b1 / a11;
						}
					}
				}
				else { // region 3			
					t0 = zero;
					if (b1 >= zero) { // V0				
						t1 = zero;
					}
					else if (-b1 >= a11) { // V2				
						t1 = one;
					}
					else { // E20				
						t1 = -b1 / a11;
					}
				}
			}
			else if (t1 < zero) { // region 5		
				t1 = zero;
				if (b0 >= zero) { // V0			
					t0 = zero;
				}
				else if (-b0 >= a00) { // V1			
					t0 = one;
				}
				else { // E01			
					t0 = -b0 / a00;
				}
			}
			else { // region 0, interior		
				const double invDet = one / det;
				t0 *= invDet;
				t1 *= invDet;
			}
		}
		else {
			double tmp0, tmp1, numer, denom;

			if (t0 < zero) { // region 2		
				tmp0 = a01 + b0;
				tmp1 = a11 + b1;
				if (tmp1 > tmp0) {
					numer = tmp1 - tmp0;
					denom = a00 - 2.0 * a01 + a11;
					if (numer >= denom) { // V1				
						t0 = one;
						t1 = zero;
					}
					else { // E12				
						t0 = numer / denom;
						t1 = one - t0;
					}
				}
				else {
					t0 = zero;
					if (tmp1 <= zero) { // V2				
						t1 = one;
					}
					else if (b1 >= zero) { // V0				
						t1 = zero;
					}
					else { // E20				
						t1 = -b1 / a11;
					}
				}
			}
			else if (t1 < zero) { // region 6		
				tmp0 = a01 + b1;
				tmp1 = a00 + b0;
				if (tmp1 > tmp0) {
					numer = tmp1 - tmp0;
					denom = a00 - 2.0 * a01 + a11;
					if (numer >= denom) { // V2				
						t1 = one;
						t0 = zero;
					}
					else { // E12				
						t1 = numer / denom;
						t0 = one - t1;
					}
				}
				else {
					t1 = zero;
					if (tmp1 <= zero) { // V1				
						t0 = one;
					}
					else if (b0 >= zero) { // V0				
						t0 = zero;
					}
					else { // E01				
						t0 = -b0 / a00;
					}
				}
			}
			else { // region 1		
				numer = a11 + b1 - a01 - b0;
				if (numer <= zero) { // V2			
					t0 = zero;
					t1 = one;
				}
				else {
					denom = a00 - 2.0 * a01 + a11;
					if (numer >= denom) { // V1				
						t0 = one;
						t1 = zero;
					}
					else { // 12				
						t0 = numer / denom;
						t1 = one - t0;
					}
				}
			}
		}

		pmp::vec3 closest = vertices[0] + t0 * edge0 + t1 * edge1;
		SUB(diff, point, closest);

		return DOT(diff, diff);
	}

	/**
	 * \brief Utility to compute intersection between plane (defined by normal and ref point) and a box (defined only by its max point).
	 * \param normal     plane normal.
	 * \param refPt      plane reference point.
	 * \param boxMax     max point of a box.
	 * \return true if the plane intersects the box.
	 */
	[[nodiscard]] bool PlaneIntersectsBox(const pmp::vec3& normal, const pmp::vec3& refPt, const pmp::vec3& boxMax) {
		pmp::vec3 vmin, vmax;

		for (int q = 0; q <= 2; q++) 
		{
			if (normal[q] > 0.0f) 
			{
				vmin[q] = -boxMax[q] - refPt[q];
				vmax[q] = boxMax[q] - refPt[q];
				continue;
			}

			vmin[q] = boxMax[q] - refPt[q];
			vmax[q] = -boxMax[q] - refPt[q];
		}

		if (DOT(normal, vmin) > 0.0f) return false;
		if (DOT(normal, vmax) >= 0.0f) return true;
		return false;
	}

	// =================== Helper macros for axis tests ======================
	// X-tests:

	#define AXISTEST_X01(a, b, fa, fb)									       \
			p0 = (a) * v0[1] - (b) * v0[2];			       				       \
			p2 = (a) * v2[1] - (b) * v2[2];			       					   \
	        if (p0 < p2) {min = p0; max = p2;} else {min = p2; max = p0;}      \
			rad = (fa) * boxHalfSize[1] + (fb) * boxHalfSize[2];		       \
			if (min > rad || max < -rad) return false

	#define AXISTEST_X2(a, b, fa, fb)									       \
			p0 = (a) * v0[1] - (b) * v0[2];									   \
			p1 = (a) * v1[1] - (b) * v1[2];		       						   \
	        if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;}      \
			rad = (fa) * boxHalfSize[1] + (fb) * boxHalfSize[2];			   \
			if (min > rad || max < -rad) return false

	// Y-tests:

	#define AXISTEST_Y02(a, b, fa, fb)									       \
			p0 = -(a) * v0[0] + (b) * v0[2];				   				   \
			p2 = -(a) * v2[0] + (b) * v2[2];				       			   \
			if (p0 < p2) {min = p0; max = p2;} else {min = p2; max = p0;}      \
			rad = (fa) * boxHalfSize[0] + (fb) * boxHalfSize[2];		       \
			if (min > rad || max < -rad) return false

	#define AXISTEST_Y1(a, b, fa, fb)							    	       \
			p0 = -(a) * v0[0] + (b) * v0[2];		      				       \
			p1 = -(a) * v1[0] + (b) * v1[2];					               \
			if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;}      \
			rad = (fa) * boxHalfSize[0] + (fb) * boxHalfSize[2];			   \
			if (min > rad || max < -rad) return false

	// Z-tests:

	#define AXISTEST_Z12(a, b, fa, fb)									       \
			p1 = (a) * v1[0] - (b) * v1[1];			 			               \
			p2 = (a) * v2[0] - (b) * v2[1];			       	                   \
			if (p2 < p1) {min = p2; max = p1;} else {min = p1; max = p2;}      \
			rad = (fa) * boxHalfSize[0] + (fb) * boxHalfSize[1];			   \
			if (min > rad || max < -rad) return false

	#define AXISTEST_Z0(a, b, fa, fb)									       \
			p0 = (a) * v0[0] - (b) * v0[1];									   \
			p1 = (a) * v1[0] - (b) * v1[1];									   \
			if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;}      \
			rad = (fa) * boxHalfSize[0] + (fb) * boxHalfSize[1];			   \
			if (min > rad || max < -rad) return false

	// =====================================================

	bool TriangleIntersectsBox(const std::vector<pmp::vec3>& vertices, const pmp::vec3& boxCenter, const pmp::vec3& boxHalfSize)
	{
		assert(vertices.size() == 3); // only vertex triples allowed

		pmp::vec3 v0, v1, v2;
		float min, max, p0, p1, p2, rad, fex, fey, fez;
		pmp::vec3 normal, e0, e1, e2;

		SUB(v0, vertices[0], boxCenter);
		SUB(v1, vertices[1], boxCenter);
		SUB(v2, vertices[2], boxCenter);

		// tri edges:
		SUB(e0, v1, v0);
		SUB(e1, v2, v1);
		SUB(e2, v0, v2);

		// 9 axis tests:
		fex = fabsf(e0[0]);
		fey = fabsf(e0[1]);
		fez = fabsf(e0[2]);

		AXISTEST_X01(e0[2], e0[1], fez, fey);
		AXISTEST_Y02(e0[2], e0[0], fez, fex);
		AXISTEST_Z12(e0[1], e0[0], fey, fex);

		fex = fabsf(e1[0]);
		fey = fabsf(e1[1]);
		fez = fabsf(e1[2]);

		AXISTEST_X01(e1[2], e1[1], fez, fey);
		AXISTEST_Y02(e1[2], e1[0], fez, fex);
		AXISTEST_Z0(e1[1], e1[0], fey, fex);

		fex = fabsf(e2[0]);
		fey = fabsf(e2[1]);
		fez = fabsf(e2[2]);

		AXISTEST_X2(e2[2], e2[1], fez, fey);
		AXISTEST_Y1(e2[2], e2[0], fez, fex);
		AXISTEST_Z12(e2[1], e2[0], fey, fex);

		// test for AABB overlap in x, y, and z:
		// test in x:
		FINDMINMAX(v0[0], v1[0], v2[0], min, max);
		if (min > boxHalfSize[0] || max < -boxHalfSize[0]) return false;

		// test in y:
		FINDMINMAX(v0[1], v1[1], v2[1], min, max);
		if (min > boxHalfSize[1] || max < -boxHalfSize[1]) return false;

		// test in z:
		FINDMINMAX(v0[2], v1[2], v2[2], min, max);
		if (min > boxHalfSize[2] || max < -boxHalfSize[2]) return false;

		// test if the box intersects the triangle plane  dot(normal, x) + d = 0

		CROSS(normal, e0, e1);

		if (!PlaneIntersectsBox(normal, v0, boxHalfSize)) return false;

		return true;
	}

	// TODO: Implement hitCount!
	bool RayIntersectsTriangle(
		const pmp::vec3& rayStart, const pmp::vec3& rayDir,
		const std::vector<pmp::vec3>& triVertices, const float& minParam, const float& maxParam)
	{
		constexpr float eps = 1e-8f;
		pmp::vec3 h, s, q;
		float a, f, u, v;
		pmp::vec3 edge1 = triVertices[1] - triVertices[0];
		pmp::vec3 edge2 = triVertices[2] - triVertices[0];
		CROSS(h, rayDir, edge2);
		a = DOT(edge1, h);
		if (a > -eps && a < eps)
		{
			return false;    // This ray is parallel to this triangle.
		}

		f = 1.0f / a;
		s = rayStart - triVertices[0];
		u = f * DOT(s, h);
		if (u < 0.0 || u > 1.0)
		{
			return false;
		}
		q = cross(s, edge1);
		v = f * DOT(rayDir, q);
		if (v < 0.0 || u + v > 1.0) 
		{
			return false;
		}

		// At this stage we can compute t to find out where the intersection point is on the line.
		float t = f * DOT(edge2, q);
		//bool inOrNoRange = t < minParam;
		//bool outOrNoRange = t > maxParam;
		if (t > eps && t < (1.0f / eps)) 
		{
			return true;
		}

		// This means that there is a line intersection but not a ray intersection.
		return false;
	}
} // namespace Geometry