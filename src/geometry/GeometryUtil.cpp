#include "GeometryUtil.h"
#include "CollisionKdTree.h"
#include "IcoSphereBuilder.h"

#include "pmp/BoundingBox.h"
#include "pmp/SurfaceMesh.h"
#include "pmp/algorithms/DifferentialGeometry.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <optional>
#include <random>
#include <vector>

namespace Geometry
{
	// ====== Helper macros for vectors (to increase speed) ============

	#define CROSS(dest, v1, v2)						     \
	          dest[0] = v1[1] * v2[2] - v1[2] * v2[1];   \
	          dest[1] = v1[2] * v2[0] - v1[0] * v2[2];	 \
	          dest[2] = v1[0] * v2[1] - v1[1] * v2[0];

	#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

	#define SUB(dest, v1, v2)						\
	          dest[0] = v1[0] - v2[0];				\
	          dest[1] = v1[1] - v2[1];				\
	          dest[2] = v1[2] - v2[2];

	#define FINDMINMAX(x0, x1, x2, min, max)		\
			  min = max = x0;						\
			  if (x1 < min) min = x1;				\
			  if (x1 > max) max = x1;				\
			  if (x2 < min) min = x2;				\
			  if (x2 > max) max = x2;

	#define SCALAR(dest,alpha,v) dest[0] = alpha * v[0]; \
	                             dest[1] = alpha * v[1]; \
	                             dest[2] = alpha * v[2];

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
	[[nodiscard]] bool PlaneIntersectsBox(const pmp::vec3& normal, const pmp::vec3& refPt, const pmp::vec3& boxMax) 
	{
		pmp::vec3 vmin, vmax;

		for (int q = 0; q <= 2; q++) 
		{
			if (normal[q] > 0.0) 
			{
				vmin[q] = -boxMax[q] - refPt[q];
				vmax[q] = boxMax[q] - refPt[q];
				continue;
			}

			vmin[q] = boxMax[q] - refPt[q];
			vmax[q] = -boxMax[q] - refPt[q];
		}

		if (DOT(normal, vmin) > 0.0) return false;
		if (DOT(normal, vmax) >= 0.0) return true;
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
		pmp::Scalar min, max, p0, p1, p2, rad, fex, fey, fez;
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

	// ====== Helper tri-tri macros for intersection test functions ============
	// Source: Contours by benardp, https://github.com/benardp/contours, freestyle/view_map/triangle_triangle_intersection.c

#define ORIENT_2D(a, b, c)  ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]))


#define INTERSECTION_TEST_VERTEX(P1, Q1, R1, P2, Q2, R2) { \
  if (ORIENT_2D(R2,P2,Q1) >= 0.0)                         \
    if (ORIENT_2D(R2,Q2,Q1) <= 0.0)                       \
      if (ORIENT_2D(P1,P2,Q1) > 0.0) {                    \
	if (ORIENT_2D(P1,Q2,Q1) <= 0.0) return 1;             \
	else return 0;} else {                                 \
	if (ORIENT_2D(P1,P2,R1) >= 0.0)                       \
	  if (ORIENT_2D(Q1,R1,P2) >= 0.0) return 1;           \
	  else return 0;                                       \
	else return 0;}                                        \
    else                                                   \
      if (ORIENT_2D(P1,Q2,Q1) <= 0.0)                     \
	if (ORIENT_2D(R2,Q2,R1) <= 0.0)                       \
	  if (ORIENT_2D(Q1,R1,Q2) >= 0.0) return 1;           \
	  else return 0;                                       \
	else return 0;                                         \
      else return 0;                                       \
  else                                                     \
    if (ORIENT_2D(R2,P2,R1) >= 0.0)                       \
      if (ORIENT_2D(Q1,R1,R2) >= 0.0)                     \
	if (ORIENT_2D(P1,P2,R1) >= 0.0) return 1;             \
	else return 0;                                         \
      else                                                 \
	if (ORIENT_2D(Q1,R1,Q2) >= 0.0) {                     \
	  if (ORIENT_2D(R2,R1,Q2) >= 0.0) return 1;           \
	  else return 0; }                                     \
	else return 0;                                         \
    else  return 0;                                        \
 }
	

#define INTERSECTION_TEST_EDGE(P1, Q1, R1, P2, Q2, R2) { \
  if (ORIENT_2D(R2,P2,Q1) >= 0.0) {\
    if (ORIENT_2D(P1,P2,Q1) >= 0.0) { \
        if (ORIENT_2D(P1,Q1,R2) >= 0.0) return 1; \
        else return 0;} else { \
      if (ORIENT_2D(Q1,R1,P2) >= 0.0){ \
	if (ORIENT_2D(R1,P1,P2) >= 0.0) return 1; else return 0;} \
      else return 0; } \
  } else {\
    if (ORIENT_2D(R2,P2,R1) >= 0.0) {\
      if (ORIENT_2D(P1,P2,R1) >= 0.0) {\
	if (ORIENT_2D(P1,R1,R2) >= 0.0) return 1;  \
	else {\
	  if (ORIENT_2D(Q1,R1,R2) >= 0.0) return 1; else return 0;}}\
      else  return 0; }\
    else return 0; }}

	// ====== Helper functions for tri-tri intersection test functions ============
	// Source: Contours by benardp, https://github.com/benardp/contours, freestyle/view_map/triangle_triangle_intersection.c

	int ccw_tri_tri_intersection_2d(
		const pmp::Scalar p1[2], const pmp::Scalar q1[2], const pmp::Scalar r1[2],
		const pmp::Scalar p2[2], const pmp::Scalar q2[2], const pmp::Scalar r2[2])
	{
		if (ORIENT_2D(p2, q2, p1) >= 0.0)
		{
			if (ORIENT_2D(q2, r2, p1) >= 0.0)
			{
				if (ORIENT_2D(r2, p2, p1) >= 0.0) return 1;
				INTERSECTION_TEST_EDGE(p1, q1, r1, p2, q2, r2)
			}
			if (ORIENT_2D(r2, p2, p1) >= 0.0)
				INTERSECTION_TEST_EDGE(p1, q1, r1, r2, p2, q2)
				INTERSECTION_TEST_VERTEX(p1, q1, r1, p2, q2, r2)
		}
		if (ORIENT_2D(q2, r2, p1) >= 0.0)
		{
			if (ORIENT_2D(r2, p2, p1) >= 0.0)
				INTERSECTION_TEST_EDGE(p1, q1, r1, q2, r2, p2)
				INTERSECTION_TEST_VERTEX(p1, q1, r1, q2, r2, p2)
		}
		else INTERSECTION_TEST_VERTEX(p1, q1, r1, r2, p2, q2)
	}

	int tri_tri_overlap_test_2d(
		const pmp::Scalar p1[2], const pmp::Scalar q1[2], const pmp::Scalar r1[2], 
		const pmp::Scalar p2[2], const pmp::Scalar q2[2], const pmp::Scalar r2[2])
	{
		if (ORIENT_2D(p1, q1, r1) < 0.0)
		{
			if (ORIENT_2D(p2, q2, r2) < 0.0)
			{
				return ccw_tri_tri_intersection_2d(p1, r1, q1, p2, r2, q2);
			}
			return ccw_tri_tri_intersection_2d(p1, r1, q1, p2, q2, r2);
		}
		if (ORIENT_2D(p2, q2, r2) < 0.0)
		{
			return ccw_tri_tri_intersection_2d(p1, q1, r1, p2, r2, q2);
		}

		return ccw_tri_tri_intersection_2d(p1, q1, r1, p2, q2, r2);
	}

	int coplanar_tri_tri3d(
		const pmp::Scalar p1[3], const pmp::Scalar q1[3], const pmp::Scalar r1[3],
		const pmp::Scalar p2[3], const pmp::Scalar q2[3], const pmp::Scalar r2[3],
		const pmp::Scalar normal_1[3])
	{

		pmp::Scalar P1[2], Q1[2], R1[2];
		pmp::Scalar P2[2], Q2[2], R2[2];

		pmp::Scalar n_x, n_y, n_z;

		n_x = ((normal_1[0] < 0) ? -normal_1[0] : normal_1[0]);
		n_y = ((normal_1[1] < 0) ? -normal_1[1] : normal_1[1]);
		n_z = ((normal_1[2] < 0) ? -normal_1[2] : normal_1[2]);

		/* Projection of the triangles in 3D onto 2D such that the area of
		   the projection is maximized. */
		if ((n_x > n_z) && (n_x >= n_y)) {
			// Project onto plane YZ

			P1[0] = q1[2]; P1[1] = q1[1];
			Q1[0] = p1[2]; Q1[1] = p1[1];
			R1[0] = r1[2]; R1[1] = r1[1];

			P2[0] = q2[2]; P2[1] = q2[1];
			Q2[0] = p2[2]; Q2[1] = p2[1];
			R2[0] = r2[2]; R2[1] = r2[1];

		}
		else if ((n_y > n_z) && (n_y >= n_x)) {
			// Project onto plane XZ

			P1[0] = q1[0]; P1[1] = q1[2];
			Q1[0] = p1[0]; Q1[1] = p1[2];
			R1[0] = r1[0]; R1[1] = r1[2];

			P2[0] = q2[0]; P2[1] = q2[2];
			Q2[0] = p2[0]; Q2[1] = p2[2];
			R2[0] = r2[0]; R2[1] = r2[2];

		}
		else {
			// Project onto plane XY

			P1[0] = p1[0]; P1[1] = p1[1];
			Q1[0] = q1[0]; Q1[1] = q1[1];
			R1[0] = r1[0]; R1[1] = r1[1];

			P2[0] = p2[0]; P2[1] = p2[1];
			Q2[0] = q2[0]; Q2[1] = q2[1];
			R2[0] = r2[0]; R2[1] = r2[1];
		}

		return tri_tri_overlap_test_2d(P1, Q1, R1, P2, Q2, R2);
	}

	/* Permutation in a canonical form of T2's vertices */

#define CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2) {\
	  SUB(v1,p2,q1);\
	  SUB(v2,p1,q1);\
	  CROSS(N1,v1,v2);\
	  SUB(v1,q2,q1);\
	  if (DOT(v1,N1) > 0.0) return 0;\
	  SUB(v1,p2,p1);\
	  SUB(v2,r1,p1);\
	  CROSS(N1,v1,v2);\
	  SUB(v1,r2,p1);\
	  if (DOT(v1,N1) > 0.0) return 0;\
	  else return 1; }

#define TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
	  if (dp2 > 0.0) { \
	     if (dq2 > 0.0) CHECK_MIN_MAX(p1,r1,q1,r2,p2,q2) \
	     else if (dr2 > 0.0) CHECK_MIN_MAX(p1,r1,q1,q2,r2,p2)\
	     else CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2) }\
	  else if (dp2 < 0.0) { \
	    if (dq2 < 0.0) CHECK_MIN_MAX(p1,q1,r1,r2,p2,q2)\
	    else if (dr2 < 0.0) CHECK_MIN_MAX(p1,q1,r1,q2,r2,p2)\
	    else CHECK_MIN_MAX(p1,r1,q1,p2,q2,r2)\
	  } else { \
	    if (dq2 < 0.0) { \
	      if (dr2 >= 0.0)  CHECK_MIN_MAX(p1,r1,q1,q2,r2,p2)\
	      else CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2)\
	    } \
	    else if (dq2 > 0.0) { \
	      if (dr2 > 0.0) CHECK_MIN_MAX(p1,r1,q1,p2,q2,r2)\
	      else  CHECK_MIN_MAX(p1,q1,r1,q2,r2,p2)\
	    } \
	    else  { \
	      if (dr2 > 0.0) CHECK_MIN_MAX(p1,q1,r1,r2,p2,q2)\
	      else if (dr2 < 0.0) CHECK_MIN_MAX(p1,r1,q1,r2,p2,q2)\
	      else return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1);\
	     }}}

	static [[nodiscard]] int tri_tri_overlap_test_3d(
		const pmp::Scalar p1[3], const pmp::Scalar q1[3], const pmp::Scalar r1[3],
		const pmp::Scalar p2[3], const pmp::Scalar q2[3], const pmp::Scalar r2[3])
	{
		pmp::Scalar v1[3], v2[3];
		pmp::Scalar N1[3], N2[3];

		/* Compute distance signs  of p1, q1 and r1 to the plane of
		   triangle(p2,q2,r2) */
		SUB(v1, p2, r2)
		SUB(v2, q2, r2)
		CROSS(N2, v1, v2)

		SUB(v1, p1, r2)
		pmp::Scalar dp1 = DOT(v1, N2);
		SUB(v1, q1, r2)
		pmp::Scalar dq1 = DOT(v1, N2);
		SUB(v1, r1, r2)
		pmp::Scalar dr1 = DOT(v1, N2);

		if (((dp1 * dq1) > 0.0) && ((dp1 * dr1) > 0.0))  return 0;
		/* Compute distance signs  of p2, q2 and r2 to the plane of
		   triangle(p1,q1,r1) */
		SUB(v1, q1, p1)
		SUB(v2, r1, p1)
		CROSS(N1, v1, v2)

		SUB(v1, p2, r1)
		pmp::Scalar dp2 = DOT(v1, N1);
		SUB(v1, q2, r1)
		pmp::Scalar dq2 = DOT(v1, N1);
		SUB(v1, r2, r1)
		pmp::Scalar dr2 = DOT(v1, N1);
		if (((dp2 * dq2) > 0.0) && ((dp2 * dr2) > 0.0)) return 0;

		/* Permutation in a canonical form of T1's vertices */
		if (dp1 > 0.0) 
		{
			if (dq1 > 0.0) TRI_TRI_3D(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2)
			if (dr1 > 0.0) TRI_TRI_3D(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2)
			TRI_TRI_3D(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2)
		}
		if (dp1 < 0.0)
		{
			if (dq1 < 0.0) TRI_TRI_3D(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2)
			if (dr1 < 0.0) TRI_TRI_3D(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2)
			TRI_TRI_3D(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2)
		}
		if (dq1 < 0.0) 
		{
			if (dr1 >= 0.0) TRI_TRI_3D(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2)
			TRI_TRI_3D(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2)
		}
		if (dq1 > 0.0)
		{
			if (dr1 > 0.0) TRI_TRI_3D(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2)
			TRI_TRI_3D(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2)
		}
		if (dr1 > 0.0) TRI_TRI_3D(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2)
		if (dr1 < 0.0) TRI_TRI_3D(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2)
		return coplanar_tri_tri3d(p1, q1, r1, p2, q2, r2, N1);
	}

	bool TriangleIntersectsTriangle(const std::vector<pmp::vec3>& vertices0, const std::vector<pmp::vec3>& vertices1)
	{
		return tri_tri_overlap_test_3d(
			vertices0[0].data(), vertices0[1].data(),vertices0[2].data(),
			vertices1[0].data(), vertices1[1].data(), vertices1[2].data()) > 0;
	}


	/*
	* =========================================================================
	*  Three-dimensional Triangle-Triangle Intersection
	* =========================================================================
	*/

	// actual intersection function for debugging instead of the macro.
	int construct_intersection(
		const pmp::Scalar p1[3], const pmp::Scalar q1[3], const pmp::Scalar r1[3],
		const pmp::Scalar p2[3], const pmp::Scalar q2[3], const pmp::Scalar r2[3],
		const pmp::Scalar N1Orig[3], const pmp::Scalar N2Orig[3],
		pmp::Scalar startPt[3], pmp::Scalar endPt[3])
	{
		pmp::Scalar v1[3], v2[3], v[3];
		pmp::Scalar N1[3] = { N1Orig[0], N1Orig[1], N1Orig[2] }, N2[3] = { N2Orig[0], N2Orig[1], N2Orig[2] }, N[3];
		pmp::Scalar alpha;

		SUB(v1, q1, p1)
		SUB(v2, r2, p1)
		CROSS(N, v1, v2)
		SUB(v, p2, p1)
		if (DOT(v, N) > 0.0) 
		{
			SUB(v1, r1, p1)
			CROSS(N, v1, v2)
			if (DOT(v, N) <= 0.0)
			{
				SUB(v2, q2, p1)
				CROSS(N,v1,v2)
				if (DOT(v, N) > 0.0)
				{
					SUB(v1, p1, p2)
					SUB(v2, p1, r1)
					alpha = DOT(v1, N2) / DOT(v2, N2);
					SCALAR(v1, alpha, v2)
					SUB(startPt, p1, v1)
					SUB(v1, p2, p1)
					SUB(v2, p2, r2)
					alpha = DOT(v1, N1) / DOT(v2, N1);
					SCALAR(v1, alpha, v2)
					SUB(endPt, p2, v1)
					return 1;
				}

				SUB(v1, p2, p1)
				SUB(v2, p2, q2)
				alpha = DOT(v1, N1) / DOT(v2, N1);
				SCALAR(v1, alpha, v2)
				SUB(startPt, p2, v1)
				SUB(v1, p2, p1)
				SUB(v2, p2, r2)
				alpha = DOT(v1, N1) / DOT(v2, N1);
				SCALAR(v1, alpha, v2)
				SUB(endPt, p2, v1)
				return 1;
			}
			
			return 0;
		}

		SUB(v2, q2, p1)
		CROSS(N, v1, v2)
		if (DOT(v, N) < 0.0)
		{
			return 0;
		}
		SUB(v1, r1, p1)
		CROSS(N, v1, v2)
		if (DOT(v, N) >= 0.0)
		{
			SUB(v1, p1, p2)
			SUB(v2, p1, r1)
			alpha = DOT(v1, N2) / DOT(v2, N2);
			SCALAR(v1, alpha, v2)
			SUB(startPt, p1, v1)
			SUB(v1, p1, p2)
			SUB(v2, p1, q1)
			alpha = DOT(v1, N2) / DOT(v2, N2);
			SCALAR(v1, alpha, v2)
			SUB(endPt, p1, v1)
			return 1;
		}

		SUB(v1, p2, p1)
		SUB(v2, p2, q2)
		alpha = DOT(v1, N1) / DOT(v2, N1);
		SCALAR(v1, alpha, v2)
		SUB(startPt, p2, v1)
		SUB(v1, p1, p2)
		SUB(v2, p1, q1)
		alpha = DOT(v1, N1) / DOT(v2, N1);
		SCALAR(v1, alpha, v2)
		SUB(endPt, p1, v1)
		return 1;
	}

	/*
	   This macro is called when the triangles surely intersect
	   It constructs the segment of intersection of the two triangles
	   if they are not coplanar.
	*/
	#define CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) { \
	  SUB(v1,q1,p1) \
	  SUB(v2,r2,p1) \
	  CROSS(N,v1,v2) \
	  SUB(v,p2,p1) \
	  if (DOT(v,N) > 0.0) {\
	    SUB(v1,r1,p1) \
	    CROSS(N,v1,v2) \
	    if (DOT(v,N) <= 0.0) { \
	      SUB(v2,q2,p1) \
	      CROSS(N,v1,v2) \
	      if (DOT(v,N) > 0.0) { \
		SUB(v1,p1,p2) \
		SUB(v2,p1,r1) \
		alpha = DOT(v1,N2) / DOT(v2,N2); \
		SCALAR(v1,alpha,v2) \
		SUB(source,p1,v1) \
		SUB(v1,p2,p1) \
		SUB(v2,p2,r2) \
		alpha = DOT(v1,N1) / DOT(v2,N1); \
		SCALAR(v1,alpha,v2) \
		SUB(target,p2,v1) \
		return 1; \
	      } else { \
		SUB(v1,p2,p1) \
		SUB(v2,p2,q2) \
		alpha = DOT(v1,N1) / DOT(v2,N1); \
		SCALAR(v1,alpha,v2) \
		SUB(source,p2,v1) \
		SUB(v1,p2,p1) \
		SUB(v2,p2,r2) \
		alpha = DOT(v1,N1) / DOT(v2,N1); \
		SCALAR(v1,alpha,v2) \
		SUB(target,p2,v1) \
		return 1; \
	      } \
	    } else { \
	      return 0; \
	    } \
	  } else { \
	    SUB(v2,q2,p1) \
	    CROSS(N,v1,v2) \
	    if (DOT(v,N) < 0.0) { \
	      return 0; \
	    } else { \
	      SUB(v1,r1,p1) \
	      CROSS(N,v1,v2) \
	      if (DOT(v,N) >= 0.0) { \
		SUB(v1,p1,p2) \
		SUB(v2,p1,r1) \
		alpha = DOT(v1,N2) / DOT(v2,N2); \
		SCALAR(v1,alpha,v2) \
		SUB(source,p1,v1) \
		SUB(v1,p1,p2) \
		SUB(v2,p1,q1) \
		alpha = DOT(v1,N2) / DOT(v2,N2); \
		SCALAR(v1,alpha,v2) \
		SUB(target,p1,v1) \
		return 1; \
	      } else { \
		SUB(v1,p2,p1) \
		SUB(v2,p2,q2) \
		alpha = DOT(v1,N1) / DOT(v2,N1); \
		SCALAR(v1,alpha,v2) \
		SUB(source,p2,v1) \
		SUB(v1,p1,p2) \
		SUB(v2,p1,q1) \
		alpha = DOT(v1,N2) / DOT(v2,N2); \
		SCALAR(v1,alpha,v2) \
		SUB(target,p1,v1) \
		return 1; \
	      }}}} 

	  // #define TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
	  //if (dp2 > 0.0) { \
	  //   if (dq2 > 0.0) return construct_intersection(p1,r1,q1,r2,p2, q2, N1, N2, source, target); \
	  //   else if (dr2 > 0.0) return construct_intersection(p1,r1,q1,q2,r2,p2, N1, N2, source, target);\
	  //   else return construct_intersection(p1,q1,r1,p2,q2,r2, N1, N2, source, target); }\
	  //else if (dp2 < 0.0) { \
	  //  if (dq2 < 0.0) return construct_intersection(p1,q1,r1,r2,p2,q2, N1, N2, source, target);\
	  //  else if (dr2 < 0.0) return construct_intersection(p1,q1,r1,q2,r2,p2, N1, N2, source, target);\
	  //  else return construct_intersection(p1,r1,q1,p2,q2,r2, N1, N2, source, target);\
	  //} else { \
	  //  if (dq2 < 0.0) { \
	  //    if (dr2 >= 0.0)  return construct_intersection(p1,r1,q1,q2,r2,p2, N1, N2, source, target);\
	  //    else return construct_intersection(p1,q1,r1,p2,q2,r2, N1, N2, source, target);\
	  //  } \
	  //  else if (dq2 > 0.0) { \
	  //    if (dr2 > 0.0) return construct_intersection(p1,r1,q1,p2,q2,r2, N1, N2, source, target);\
	  //    else  return construct_intersection(p1,q1,r1,q2,r2,p2, N1, N2, source, target);\
	  //  } \
	  //  else  { \
	  //    if (dr2 > 0.0) return construct_intersection(p1,q1,r1,r2,p2,q2, N1, N2, source, target);\
	  //    else if (dr2 < 0.0) return construct_intersection(p1,r1,q1,r2,p2,q2, N1, N2, source, target);\
	  //    else { \
      //    *coplanar = 1; \
	  //    return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1);\
	  //   } \
	  //}} }
	#define TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
	  if (dp2 > 0.0) { \
	     if (dq2 > 0.0) CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2) \
	     else if (dr2 > 0.0) CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
	     else CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) }\
	  else if (dp2 < 0.0) { \
	    if (dq2 < 0.0) CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
	    else if (dr2 < 0.0) CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
	    else CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
	  } else { \
	    if (dq2 < 0.0) { \
	      if (dr2 >= 0.0)  CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
	      else CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2)\
	    } \
	    else if (dq2 > 0.0) { \
	      if (dr2 > 0.0) CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
	      else  CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
	    } \
	    else  { \
	      if (dr2 > 0.0) CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
	      else if (dr2 < 0.0) CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2)\
	      else { \
       		*coplanar = 1; \
		return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1);\
	     } \
	  }} }


	int tri_tri_intersection_test_3d(
		const pmp::Scalar p1[3], const pmp::Scalar q1[3], const pmp::Scalar r1[3],
		const pmp::Scalar p2[3], const pmp::Scalar q2[3], const pmp::Scalar r2[3],
		int* coplanar, pmp::Scalar source[3], pmp::Scalar target[3])
	{
		pmp::Scalar v1[3], v2[3];
		pmp::Scalar N1[3], N2[3];
		pmp::Scalar v[3], N[3];
		pmp::Scalar alpha;

		// Compute distance signs  of p1, q1 and r1 
		// to the plane of triangle(p2,q2,r2)
		SUB(v1, p2, r2)
		SUB(v2, q2, r2)
		CROSS(N2, v1, v2)

		SUB(v1, p1, r2)
		pmp::Scalar dp1 = DOT(v1, N2);
		SUB(v1, q1, r2)
		pmp::Scalar dq1 = DOT(v1, N2);
		SUB(v1, r1, r2)
		pmp::Scalar dr1 = DOT(v1, N2);

		if (((dp1 * dq1) > 0.0) && ((dp1 * dr1) > 0.0))  return 0;

		// Compute distance signs  of p2, q2 and r2 
		// to the plane of triangle(p1,q1,r1)
		SUB(v1, q1, p1)
		SUB(v2, r1, p1)
		CROSS(N1, v1, v2)

		SUB(v1, p2, r1)
		pmp::Scalar dp2 = DOT(v1, N1);
		SUB(v1, q2, r1)
		pmp::Scalar dq2 = DOT(v1, N1);
		SUB(v1, r2, r1)
		pmp::Scalar dr2 = DOT(v1, N1);

		if (((dp2 * dq2) > 0.0) && ((dp2 * dr2) > 0.0)) return 0;
		// Permutation in a canonical form of T1's vertices

		//  printf("d1 = [%f %f %f], d2 = [%f %f %f]\n", dp1, dq1, dr1, dp2, dq2, dr2);
		/*
		// added by Aaron
		if (ZERO_TEST(dp1) || ZERO_TEST(dq1) ||ZERO_TEST(dr1) ||ZERO_TEST(dp2) ||ZERO_TEST(dq2) ||ZERO_TEST(dr2))
		  {
			coplanar = 1;
			return 0;
		  }
		*/
		if (dp1 > 0.0) 
		{
			if (dq1 > 0.0) TRI_TRI_INTER_3D(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2)
			if (dr1 > 0.0) TRI_TRI_INTER_3D(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2)
			TRI_TRI_INTER_3D(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2)
		}
		if (dp1 < 0.0)
		{
			if (dq1 < 0.0) TRI_TRI_INTER_3D(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2)
			if (dr1 < 0.0) TRI_TRI_INTER_3D(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2)
			TRI_TRI_INTER_3D(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2)
		}
		if (dq1 < 0.0) 
		{
			if (dr1 >= 0.0) TRI_TRI_INTER_3D(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2)
			TRI_TRI_INTER_3D(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2)
		}
		if (dq1 > 0.0) 
		{
			if (dr1 > 0.0) TRI_TRI_INTER_3D(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2)
			TRI_TRI_INTER_3D(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2)
		}
		if (dr1 > 0.0) TRI_TRI_INTER_3D(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2)
		if (dr1 < 0.0) TRI_TRI_INTER_3D(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2)
		// triangles are co-planar
		* coplanar = 1;
		return coplanar_tri_tri3d(p1, q1, r1, p2, q2, r2, N1);
	}

	std::optional<std::pair<pmp::vec3, pmp::vec3>> ComputeTriangleTriangleIntersectionLine(const std::vector<pmp::vec3>& vertices0, const std::vector<pmp::vec3>& vertices1)
	{
		int coplanar{ 0 };
		pmp::vec3 startPt;
		pmp::vec3 endPt;
		if (tri_tri_intersection_test_3d(
			vertices0[0].data(), vertices0[1].data(), vertices0[2].data(),
			vertices1[0].data(), vertices1[1].data(), vertices1[2].data(), 
			&coplanar, startPt.data(), endPt.data()) == 0)
		{
			return {};
		}

		return std::pair{ startPt, endPt };
	}

	/// \brief intersection tolerance for Moller-Trumbore algorithm.
	constexpr pmp::Scalar MT_INTERSECTION_EPSILON = 1e-5;

	bool RayIntersectsTriangle(Ray& ray, const std::vector<pmp::vec3>& triVertices)
	{
		pmp::vec3 edge1 = triVertices[1] - triVertices[0];
		pmp::vec3 edge2 = triVertices[2] - triVertices[0];
		pmp::vec3 cross1;
		CROSS(cross1, ray.Direction, edge2);
		const pmp::Scalar det = DOT(edge1, cross1);
		if (std::fabs(det) < MT_INTERSECTION_EPSILON)
		{
			return false; // Ray is parallel to triangle
		}

		const pmp::Scalar invDet = 1.0 / det;
		const auto startToTri0 = ray.StartPt - triVertices[0];
		const pmp::Scalar u = invDet * DOT(startToTri0, cross1);
		if (u < -MT_INTERSECTION_EPSILON || u > 1.0 + MT_INTERSECTION_EPSILON)
		{
			return false;
		}

		pmp::vec3 cross2;
		CROSS(cross2, startToTri0, edge1);
		const pmp::Scalar v = invDet * DOT(ray.Direction, cross2);
		if (v < -MT_INTERSECTION_EPSILON || u + v > 1.0 + MT_INTERSECTION_EPSILON)
		{
			return false;
		}

		const pmp::Scalar t = invDet * DOT(edge2, cross2);
		if (t > MT_INTERSECTION_EPSILON && t < 1.0 / MT_INTERSECTION_EPSILON &&
			t >= ray.ParamMin && t <= ray.ParamMax)
		{
			ray.HitParam = t;
			return true;
		}
		return false; // No valid intersection
	}

	//
	// =========================================================================
	//

	/**
	 * \brief Computes the preference for ray abs direction during ray-triangle intersection.
	 * \param absDir    evaluated direction vector.
	 * \return index of the preferred axis.
	 */
	[[nodiscard]] unsigned int MaxDim(const pmp::vec3& absDir)
	{
		if ((absDir[0] > absDir[1]) && (absDir[0] > absDir[2]))
		{
			return 0; // X-axis;
		}

		if (absDir[1] > absDir[2])
		{
			return 1; //  Y-axis;
		}

		return 2; // Z-axis;
	}

	Ray::Ray(const pmp::vec3& startPt, const pmp::vec3& dir)
		: StartPt(startPt), Direction(dir)
	{
		// dir vector must be normalized.
		const pmp::Scalar dirLenSq = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2];
		if (dirLenSq >= 1.0 + FLT_EPSILON || dirLenSq <= 1.0 - FLT_EPSILON)
		{
			throw std::logic_error("Ray::Ray: ||dir|| != 1 ! Ray direction vector must be normalized!\n");
		}

		InvDirection[0] = 1.0 / dir[0];
		InvDirection[1] = 1.0 / dir[1];
		InvDirection[2] = 1.0 / dir[2];

		// calculate id of max dimension of ray direction.
		const pmp::vec3 absDir{ std::fabs(dir[0]), std::fabs(dir[1]), std::fabs(dir[2]) };
		kz = MaxDim(absDir);
		//kx = kz + 1; if (kx == 3) kx = 0;
		//ky = kx + 1; if (ky == 3) ky = 0;
		kx = (kz + 1) % 3;
		ky = (kx + 1) % 3;
		// swap kx and ky dims to preserve winding direction of triangles.
		if (dir[kz] < 0.0) std::swap(kx, ky);
		// calculate shear constants
		Sx = dir[kx] / dir[kz];
		Sy = dir[ky] / dir[kz];
		Sz = 1.0 / dir[kz];
	}

	/// \brief if true, we use backface-culling (i.e. skipping triangles whose normals point away from the ray).
#define BACKFACE_CULLING false
	
	bool RayIntersectsTriangleWatertight(Ray& ray, const std::vector<pmp::vec3>& triVertices)
	{
		// actual alg
		const pmp::vec3 A = triVertices[0] - ray.StartPt;
		const pmp::vec3 B = triVertices[1] - ray.StartPt;
		const pmp::vec3 C = triVertices[2] - ray.StartPt;
		const pmp::Scalar Ax = A[ray.kx] - ray.Sx * A[ray.kz];
		const pmp::Scalar Ay = A[ray.ky] - ray.Sy * A[ray.kz];
		const pmp::Scalar Bx = B[ray.kx] - ray.Sx * B[ray.kz];
		const pmp::Scalar By = B[ray.ky] - ray.Sy * B[ray.kz];
		const pmp::Scalar Cx = C[ray.kx] - ray.Sx * C[ray.kz];
		const pmp::Scalar Cy = C[ray.ky] - ray.Sy * C[ray.kz];
		pmp::Scalar U = Cx * By - Cy * Bx;
		pmp::Scalar V = Ax * Cy - Ay * Cx;
		pmp::Scalar W = Bx * Ay - By * Ax;
		if (U == 0.0 || V == 0.0 || W == 0.0) 
		{
			const double CxBy = static_cast<double>(Cx) * static_cast<double>(By);
			const double CyBx = static_cast<double>(Cy) * static_cast<double>(Bx);
			U = static_cast<pmp::Scalar>(CxBy - CyBx);
			const double AxCy = static_cast<double>(Ax) * static_cast<double>(Cy);
			const double AyCx = static_cast<double>(Ay) * static_cast<double>(Cx);
			V = static_cast<pmp::Scalar>(AxCy - AyCx);
			const double BxAy = static_cast<double>(Bx) * static_cast<double>(Ay);
			const double ByAx = static_cast<double>(By) * static_cast<double>(Ax);
			W = static_cast<pmp::Scalar>(BxAy - ByAx);
		}
#if BACKFACE_CULLING
		if (U < 0.0 || V < 0.0 || W < 0.0) return false;
#else
		if ((U < 0.0 || V < 0.0 || W < 0.0) &&
		   (U > 0.0 || V > 0.0 || W > 0.0)) return false;
#endif
		const pmp::Scalar det = U + V + W;
		if (det == 0.0) return false;
		const pmp::Scalar Az = ray.Sz * A[ray.kz];
		const pmp::Scalar Bz = ray.Sz * B[ray.kz];
		const pmp::Scalar Cz = ray.Sz * C[ray.kz];
		const pmp::Scalar T = U * Az + V * Bz + W * Cz;
		if (T < 0.0 || T > ray.HitParam * det)
			return false;
		
		const pmp::Scalar rcpDet = 1.0 / det; // reciprocal det
		/*const pmp::Scalar hitBCoordU = U * rcpDet;
		const pmp::Scalar hitBCoordV = V * rcpDet;
		const pmp::Scalar hitBCoordW = W * rcpDet;*/
		ray.HitParam = T * rcpDet;
		
		return true;
	}

	// =========================================================================

	constexpr unsigned int idVec[3] = { 0, 1, 2 };

	// accelerated roundoff for watertightness [Woop, Benthin, Wald, 2013, p. 70]
	pmp::Scalar p = 1.0 + 2e-23;
	pmp::Scalar m = 1.0 - 2e-23;
	[[nodiscard]] pmp::Scalar up(const pmp::Scalar a) { return a > 0.0 ? a * p : a * m; }
	[[nodiscard]] pmp::Scalar dn(const pmp::Scalar a) { return a > 0.0 ? a * m : a * p; }
	[[nodiscard]] pmp::Scalar Up(const pmp::Scalar a) { return a * p; }
	[[nodiscard]] pmp::Scalar Dn(const pmp::Scalar a) { return a * m; }
	constexpr pmp::Scalar eps = 5.0 * 2e-24;

	bool RayIntersectsABox(const Ray& ray, const pmp::BoundingBox& box)
	{
		int nearX = static_cast<int>(idVec[ray.kx]), farX = static_cast<int>(idVec[ray.kx]);
		int nearY = static_cast<int>(idVec[ray.ky]), farY = static_cast<int>(idVec[ray.ky]);
		int nearZ = static_cast<int>(idVec[ray.kz]), farZ = static_cast<int>(idVec[ray.kz]);
		if (ray.Direction[ray.kx] < 0.0) std::swap(nearX, farX);
		if (ray.Direction[ray.ky] < 0.0) std::swap(nearY, farY);
		if (ray.Direction[ray.kz] < 0.0) std::swap(nearZ, farZ);

		const pmp::vec3 absDMin{
			std::fabs(ray.StartPt[0] - box.min()[0]),
			std::fabs(ray.StartPt[1] - box.min()[1]),
			std::fabs(ray.StartPt[2] - box.min()[2])
		};
		const pmp::vec3 absDMax{
			std::fabs(ray.StartPt[0] - box.max()[0]),
			std::fabs(ray.StartPt[1] - box.max()[1]),
			std::fabs(ray.StartPt[2] - box.max()[2])
		};
		pmp::vec3 lower{Dn(absDMin[0]),	Dn(absDMin[1]),	Dn(absDMin[2])};
		pmp::vec3 upper{Up(absDMax[0]),	Up(absDMax[1]),	Up(absDMax[2])};
		const pmp::Scalar max_z = std::max(lower[ray.kz], upper[ray.kz]);
		const pmp::Scalar err_near_x = Up(lower[ray.kx] + max_z);
		const pmp::Scalar err_near_y = Up(lower[ray.ky] + max_z);
		pmp::Scalar start_near_x = up(ray.StartPt[ray.kx] + Up(eps * err_near_x));
		pmp::Scalar start_near_y = up(ray.StartPt[ray.ky] + Up(eps * err_near_y));
		const pmp::Scalar start_near_z = ray.StartPt[ray.kz];
		const pmp::Scalar err_far_x = Up(upper[ray.kx] + max_z);
		const pmp::Scalar err_far_y = Up(upper[ray.ky] + max_z);
		pmp::Scalar start_far_x = dn(ray.StartPt[ray.kx] - Up(eps * err_far_x));
		pmp::Scalar start_far_y = dn(ray.StartPt[ray.ky] - Up(eps * err_far_y));
		const pmp::Scalar start_far_z = ray.StartPt[ray.kz];
		if (ray.Direction[ray.kx] < 0.0) std::swap(start_near_x, start_far_x);
		if (ray.Direction[ray.ky] < 0.0) std::swap(start_near_y, start_far_y);
		const pmp::Scalar rdir_near_x = Dn(Dn(ray.InvDirection[ray.kx]));
		const pmp::Scalar rdir_near_y = Dn(Dn(ray.InvDirection[ray.ky]));
		const pmp::Scalar rdir_near_z = Dn(Dn(ray.InvDirection[ray.kz]));
		const pmp::Scalar rdir_far_x = Up(Up(ray.InvDirection[ray.kx]));
		const pmp::Scalar rdir_far_y = Up(Up(ray.InvDirection[ray.ky]));
		const pmp::Scalar rdir_far_z = Up(Up(ray.InvDirection[ray.kz]));
		pmp::Scalar tNearX = (box.min()[nearX] - start_near_x) * rdir_near_x;
		pmp::Scalar tNearY = (box.min()[nearY] - start_near_y) * rdir_near_y;
		pmp::Scalar tNearZ = (box.min()[nearZ] - start_near_z) * rdir_near_z;
		pmp::Scalar tFarX = (box.max()[farX] - start_far_x) * rdir_far_x;
		pmp::Scalar tFarY = (box.max()[farY] - start_far_y) * rdir_far_y;
		pmp::Scalar tFarZ = (box.max()[farZ] - start_far_z) * rdir_far_z;
		const pmp::Scalar tNear = std::max<pmp::Scalar>({ tNearX, tNearY, tNearZ, ray.ParamMin });
		const pmp::Scalar tFar = std::min<pmp::Scalar>({ tFarX, tFarY, tFarZ, ray.ParamMax });
		return tNear <= tFar;
	}

	// =====================================================================================

	bool CircleIntersectsCircle2D(const pmp::Point2& center1, const pmp::Scalar& radius1, const pmp::Point2& center2, const pmp::Scalar& radius2)
	{
		// Calculate the distance between the centers of the two circles
		const pmp::Scalar distance = norm(center1 - center2);

		// Check if the circles intersect
		return distance < (radius1 + radius2);
	}

	bool SphereIntersectsSphere3D(const pmp::Point& center1, const pmp::Scalar& radius1, const pmp::Point& center2, const pmp::Scalar& radius2)
	{
		// Calculate the distance between the centers of the two spheres
		const pmp::Scalar distance = norm(center1 - center2);

		// Check if the spheres intersect
		return distance < (radius1 + radius2);
	}

	// =====================================================================================


	double GetDistanceToLine2DSq(const std::vector<pmp::Point2>& vertices, const pmp::vec2& point)
	{
		if (vertices.size() != 2)
		{
			throw std::invalid_argument("GetDistanceToLine2DSq: The vertices vector must contain exactly two points!\n");
		}

		// Extract vertices
		const pmp::Point2& p1 = vertices[0];
		const pmp::Point2& p2 = vertices[1];

		// Vector from p1 to p2
		pmp::vec2 v(p2[0] - p1[0], p2[1] - p1[1]);

		// Vector from p1 to the point
		pmp::vec2 w(point[0] - p1[0], point[1] - p1[1]);

		// Calculate the projection of w onto v to find the closest point on the line segment
		double c1 = w[0] * v[0] + w[1] * v[1];
		double c2 = v[0] * v[0] + v[1] * v[1];

		// If the line segment is degenerate (i.e., p1 == p2), return the distance to p1
		if (std::abs(c2) < 1e-6) 
		{
			return (w[0] * w[0] + w[1] * w[1]);
		}

		double t = c1 / c2;

		// Clamp t to the range [0, 1] to ensure the closest point is on the segment
		t = std::max(0.0, std::min(1.0, t));

		// Calculate the closest point on the line segment
		pmp::vec2 closestPoint(p1[0] + t * v[0], p1[1] + t * v[1]);

		// Vector from the closest point to the point
		pmp::vec2 distVec(point[0] - closestPoint[0], point[1] - closestPoint[1]);

		// Return the squared distance
		return distVec[0] * distVec[0] + distVec[1] * distVec[1];
	}

	bool Line2DIntersectsBox(const std::vector<pmp::Point2>& vertices, const pmp::Point2& boxCenter, const pmp::vec2& boxHalfSize)
	{
		if (vertices.size() != 2)
		{
			throw std::invalid_argument("Line2DIntersectsBox: The vertices vector must contain exactly two points!\n");
		}

		const pmp::Point2& p1 = vertices[0];
		const pmp::Point2& p2 = vertices[1];

		// direction vector
		pmp::vec2 d(p2[0] - p1[0], p2[1] - p1[1]);

		// Compute the half extents of the box
		const pmp::vec2& e = boxHalfSize;

		// Compute the center of the line segment
		pmp::vec2 lineCenter((p1[0] + p2[0]) * 0.5, (p1[1] + p2[1]) * 0.5);

		// Translate the box and the line segment to the origin
		pmp::vec2 T(lineCenter[0] - boxCenter[0], lineCenter[1] - boxCenter[1]);

		// Compute the absolute values of the direction vector components
		pmp::vec2 absD(std::fabs(d[0]), std::fabs(d[1]));

		// Perform the SAT test
		// Test the box's x and y axes
		if (std::fabs(T[0]) > e[0] + 0.5 * absD[0]) return false;
		if (std::fabs(T[1]) > e[1] + 0.5 * absD[1]) return false;

		// Test the line segment's direction vector
		pmp::Scalar r = e[1] * std::fabs(d[0]) + e[0] * std::fabs(d[1]);
		if (std::fabs(T[0] * d[1] - T[1] * d[0]) > r) return false;

		// No separating axis found, the line segment intersects the box
		return true;
	}

	/// \brief a helper func for the 2x2 determinant
	inline double Det(double a, double b, double c, double d)
	{
		return a * d - b * c;
	}

	bool Line2DIntersectsLine2D(const std::vector<pmp::Point2>& vertices0, const std::vector<pmp::Point2>& vertices1)
	{
		if (vertices0.size() < 2 || vertices1.size() < 2)
			throw std::invalid_argument("Line2DIntersectsLine2D requires exactly two points per line.");

		// Extract coordinates for the first line
		pmp::Scalar x1 = vertices0[0][0], y1 = vertices0[0][1];
		pmp::Scalar x2 = vertices0[1][0], y2 = vertices0[1][1];

		// Extract coordinates for the second line
		pmp::Scalar x3 = vertices1[0][0], y3 = vertices1[0][1];
		pmp::Scalar x4 = vertices1[1][0], y4 = vertices1[1][1];

		// Compute determinants
		pmp::Scalar detL1 = Det(x1, y1, x2, y2);
		pmp::Scalar detL2 = Det(x3, y3, x4, y4);
		pmp::Scalar x1mx2 = x1 - x2;
		pmp::Scalar x3mx4 = x3 - x4;
		pmp::Scalar y1my2 = y1 - y2;
		pmp::Scalar y3my4 = y3 - y4;

		pmp::Scalar xnom = Det(detL1, x1mx2, detL2, x3mx4);
		pmp::Scalar ynom = Det(detL1, y1my2, detL2, y3my4);
		pmp::Scalar denom = Det(x1mx2, y1my2, x3mx4, y3my4);

		if (std::abs(denom) < FLT_EPSILON) // Lines are parallel or collinear
		{
			return false;
		}

		// Calculate intersection point
		pmp::Scalar ix = xnom / denom;
		pmp::Scalar iy = ynom / denom;

		if (!std::isfinite(ix) || !std::isfinite(iy)) // Check for numerical stability
		{
			return false;
		}

		// Ensure the intersection point is within both line segments
		bool onSegment1 = (std::min(x1, x2) <= ix && ix <= std::max(x1, x2)) && (std::min(y1, y2) <= iy && iy <= std::max(y1, y2));
		bool onSegment2 = (std::min(x3, x4) <= ix && ix <= std::max(x3, x4)) && (std::min(y3, y4) <= iy && iy <= std::max(y3, y4));

		return onSegment1 && onSegment2;
	}

	bool IsPointLeftOfLine2D(const pmp::Point2& point, const std::pair<pmp::Point2, pmp::Point2>& line)
	{
		return ((line.second[0] - line.first[0]) * (point[1] - line.first[1])) >
			((line.second[1] - line.first[1]) * (point[0] - line.first[0]));
	}

	static [[nodiscard]] pmp::Scalar Sign2D(const pmp::Point2& p1, const pmp::Point2& p2, const pmp::Point2& p3)
	{
		return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1]);
	}

	bool IsPointInTriangle2D(const pmp::Point2& point, const std::vector<pmp::Point2>& triPts)
	{
		if (triPts.size() < 3)
		{
			std::cerr << "Geometry::IsPointInTriangle2D: riPts.size() < 3!\n";
			return false;
		}

		const bool b1 = Sign2D(point, triPts[0], triPts[1]);
		const bool b2 = Sign2D(point, triPts[1], triPts[2]);
		const bool b3 = Sign2D(point, triPts[2], triPts[0]);

		return ((b1 == b2) && (b2 == b3));
	}

	namespace
	{

		/// \brief Calculates the parametric distance at which two 2D rays intersect.
		[[nodiscard]] pmp::Scalar CalculateIntersectionParametricDistance(const Ray2D& ray1, const Ray2D& ray2)
		{
			// Represent rays as parametric equations: P1 = StartPt1 + t1 * Direction1, P2 = StartPt2 + t2 * Direction2
			pmp::vec2 dirCross = pmp::vec2(-ray2.Direction[1], ray2.Direction[0]); // Perpendicular direction to ray2
			pmp::Scalar det = dirCross[0] * ray1.Direction[0] + dirCross[1] * ray1.Direction[1];

			// If det is zero, rays are parallel and do not intersect
			if (std::abs(det) < INTERSECTION_EPSILON)
			{
				return FLT_MAX;
			}

			pmp::vec2 startDiff = ray2.StartPt - ray1.StartPt;
			return (startDiff[0] * dirCross[0] + startDiff[1] * dirCross[1]) / det;
		}

	} // anonymous namespace

	Ray2D& Ray2D::operator+=(const Ray2D& other)
	{
		// Calculate intersection point between the two rays
		const auto t1 = CalculateIntersectionParametricDistance(*this, other);

		// Check if the intersection is within valid parametric range
		if (t1 > ParamMin && t1 < ParamMax)
		{
			// Update HitParam and ParamMax if the intersection is closer
			HitParam = std::min<pmp::Scalar>(HitParam, t1);
			ParamMax = HitParam;
		}

		return *this;
	}

	pmp::vec2 Ray2D::GetVector() const
	{
		// Calculate the starting offset using ParamMin
		pmp::vec2 offsetStart = StartPt + ParamMin * Direction;

		// Calculate the vector representation using ParamMax and offsetStart
		return offsetStart + (ParamMax - ParamMin) * Direction;
	}

	pmp::Point2 Ray2D::GetMin() const
	{
		// Calculate the point at the start offset using ParamMin
		return StartPt + ParamMin * Direction;
	}

	pmp::Point2 Ray2D::GetMax() const
	{
		// Calculate the endpoint of the ray using ParamMax
		return StartPt + ParamMax * Direction;
	}

	bool RayBoxIntersection2D(const pmp::Point2& startPt, const pmp::vec2& direction, const pmp::BoundingBox2& box, pmp::Scalar& tMinOut, pmp::Scalar& tMaxOut)
	{
		// Initialize parameters for the intersection distances
		pmp::Scalar tMin = 0.0;
		pmp::Scalar tMax = std::numeric_limits<pmp::Scalar>::infinity();

		// Get the min and max points of the bounding box
		pmp::Point2 boxMin = box.min();
		pmp::Point2 boxMax = box.max();

		// Check for intersection with each axis (x and y)
		for (int i = 0; i < 2; ++i) // 0 for x-axis, 1 for y-axis
		{
			if (std::abs(direction[i]) > INTERSECTION_EPSILON) // Avoid division by zero
			{
				// Calculate intersection distances with the slab boundaries
				pmp::Scalar t1 = (boxMin[i] - startPt[i]) / direction[i];
				pmp::Scalar t2 = (boxMax[i] - startPt[i]) / direction[i];

				// Ensure t1 is the min and t2 is the max distance
				if (t1 > t2) std::swap(t1, t2);

				// Update tMin and tMax to include the current axis
				tMin = std::max(tMin, t1);
				tMax = std::min(tMax, t2);

				// Check for no intersection
				if (tMin > tMax) return false;
			}
			else
			{
				// The ray is parallel to the current axis; check if it's outside the box
				if (startPt[i] < boxMin[i] || startPt[i] > boxMax[i]) return false;
			}
		}

		// If we reach here, the ray intersects the bounding box
		tMinOut = tMin;
		tMaxOut = tMax;
		return true;
	}

	pmp::SurfaceMesh ConstructIcoSphere(const Sphere3D& sphere, unsigned int subdiv)
	{
		Geometry::IcoSphereBuilder icoBuilder({ subdiv, sphere.Radius });
		icoBuilder.BuildBaseData();
		icoBuilder.BuildPMPSurfaceMesh();
		if (sphere.Center == pmp::Point(0, 0, 0))
			return icoBuilder.GetPMPSurfaceMeshResult();

		auto mesh = icoBuilder.GetPMPSurfaceMeshResult();
		const auto translationMatrix = translation_matrix(sphere.Center);
		mesh *= translationMatrix;
		return mesh;
	}

	pmp::Scalar CalculateNeighborhoodRingArea(const pmp::Point& center, std::vector<pmp::Point>& neighbors)
	{
		if (neighbors.size() < 3)
		{
			return 0.0;
		}

		// compute vectors from center to each neighbor
		std::vector<pmp::vec3> vecs;
		vecs.reserve(neighbors.size());
		for (auto const& nb : neighbors)
		{
			pmp::vec3 v = nb - center;
			vecs.push_back(v);
		}

		// find a non-degenerate normal by crossing the first vector with another
		pmp::vec3 normal{ 0.0, 0.0, 0.0 };
		for (size_t i = 1; i < vecs.size(); ++i)
		{
			normal = pmp::cross(vecs[0], vecs[i]);
			if (pmp::norm(normal) > std::numeric_limits<pmp::Scalar>::epsilon())
			{
				normal /= pmp::norm(normal);
				break;
			}
			if (i + 1 == vecs.size())
			{
				return 0.0;
			}
		}

		// build tangent basis t1, t2
		double proj = pmp::dot(vecs[0], normal);
		pmp::vec3 t1 = vecs[0] - proj * normal;
		double t1len = pmp::norm(t1);
		if (t1len < std::numeric_limits<pmp::Scalar>::epsilon())
		{
			return 0.0;
		}
		t1 /= t1len;
		pmp::vec3 t2 = pmp::cross(normal, t1);

		// compute angles for sorting
		struct AngledVec { pmp::Scalar theta; pmp::vec3 v; };
		std::vector<AngledVec> angled;
		angled.reserve(vecs.size());
		for (auto const& v : vecs)
		{
			double x = pmp::dot(v, t1);
			double y = pmp::dot(v, t2);
			double theta = std::atan2(y, x);
			angled.push_back({ theta, v });
		}
		std::ranges::sort(angled, [](auto const& a, auto const& b) { return a.theta < b.theta; });

		// sum triangle areas between consecutive neighbors
		pmp::Scalar totalArea = 0.0;
		for (size_t i = 0; i < angled.size(); ++i)
		{
			const pmp::vec3& vi = angled[i].v;
			const pmp::vec3& vj = angled[(i + 1) % angled.size()].v;
			// triangle_area(center, center+vi, center+vj) == 0.5 * norm(cross(vi, vj))
			pmp::Scalar area = pmp::triangle_area(
				pmp::Point{ 0.0, 0.0, 0.0 },
				pmp::Point{ vi[0], vi[1], vi[2] },
				pmp::Point{ vj[0], vj[1], vj[2] });
			totalArea += area;
		}

		return totalArea;
	}

} // namespace Geometry