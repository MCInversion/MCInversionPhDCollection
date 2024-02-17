#include "GeometryUtil.h"

#include "sdf/CollisionKdTree.h"

#include "pmp/MatVec.h"
#include "pmp/BoundingBox.h"

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
	[[nodiscard]] bool PlaneIntersectsBox(const pmp::vec3& normal, const pmp::vec3& refPt, const pmp::vec3& boxMax) 
	{
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

	// ====== Helper tri-tri macros for intersection test functions ============
	// Source: Contours by benardp, https://github.com/benardp/contours, freestyle/view_map/triangle_triangle_intersection.c

#define ORIENT_2D(a, b, c)  ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]))


#define INTERSECTION_TEST_VERTEX(P1, Q1, R1, P2, Q2, R2) { \
  if (ORIENT_2D(R2,P2,Q1) >= 0.0f)                         \
    if (ORIENT_2D(R2,Q2,Q1) <= 0.0f)                       \
      if (ORIENT_2D(P1,P2,Q1) > 0.0f) {                    \
	if (ORIENT_2D(P1,Q2,Q1) <= 0.0f) return 1;             \
	else return 0;} else {                                 \
	if (ORIENT_2D(P1,P2,R1) >= 0.0f)                       \
	  if (ORIENT_2D(Q1,R1,P2) >= 0.0f) return 1;           \
	  else return 0;                                       \
	else return 0;}                                        \
    else                                                   \
      if (ORIENT_2D(P1,Q2,Q1) <= 0.0f)                     \
	if (ORIENT_2D(R2,Q2,R1) <= 0.0f)                       \
	  if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) return 1;           \
	  else return 0;                                       \
	else return 0;                                         \
      else return 0;                                       \
  else                                                     \
    if (ORIENT_2D(R2,P2,R1) >= 0.0f)                       \
      if (ORIENT_2D(Q1,R1,R2) >= 0.0f)                     \
	if (ORIENT_2D(P1,P2,R1) >= 0.0f) return 1;             \
	else return 0;                                         \
      else                                                 \
	if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) {                     \
	  if (ORIENT_2D(R2,R1,Q2) >= 0.0f) return 1;           \
	  else return 0; }                                     \
	else return 0;                                         \
    else  return 0;                                        \
 }
	

#define INTERSECTION_TEST_EDGE(P1, Q1, R1, P2, Q2, R2) { \
  if (ORIENT_2D(R2,P2,Q1) >= 0.0f) {\
    if (ORIENT_2D(P1,P2,Q1) >= 0.0f) { \
        if (ORIENT_2D(P1,Q1,R2) >= 0.0f) return 1; \
        else return 0;} else { \
      if (ORIENT_2D(Q1,R1,P2) >= 0.0f){ \
	if (ORIENT_2D(R1,P1,P2) >= 0.0f) return 1; else return 0;} \
      else return 0; } \
  } else {\
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) {\
      if (ORIENT_2D(P1,P2,R1) >= 0.0f) {\
	if (ORIENT_2D(P1,R1,R2) >= 0.0f) return 1;  \
	else {\
	  if (ORIENT_2D(Q1,R1,R2) >= 0.0f) return 1; else return 0;}}\
      else  return 0; }\
    else return 0; }}

	// ====== Helper functions for tri-tri intersection test functions ============
	// Source: Contours by benardp, https://github.com/benardp/contours, freestyle/view_map/triangle_triangle_intersection.c

	int ccw_tri_tri_intersection_2d(
		const pmp::Scalar p1[2], const pmp::Scalar q1[2], const pmp::Scalar r1[2],
		const pmp::Scalar p2[2], const pmp::Scalar q2[2], const pmp::Scalar r2[2])
	{
		if (ORIENT_2D(p2, q2, p1) >= 0.0f)
		{
			if (ORIENT_2D(q2, r2, p1) >= 0.0f)
			{
				if (ORIENT_2D(r2, p2, p1) >= 0.0f) return 1;
				INTERSECTION_TEST_EDGE(p1, q1, r1, p2, q2, r2)
			}
			if (ORIENT_2D(r2, p2, p1) >= 0.0f)
				INTERSECTION_TEST_EDGE(p1, q1, r1, r2, p2, q2)
				INTERSECTION_TEST_VERTEX(p1, q1, r1, p2, q2, r2)
		}
		if (ORIENT_2D(q2, r2, p1) >= 0.0f)
		{
			if (ORIENT_2D(r2, p2, p1) >= 0.0f)
				INTERSECTION_TEST_EDGE(p1, q1, r1, q2, r2, p2)
				INTERSECTION_TEST_VERTEX(p1, q1, r1, q2, r2, p2)
		}
		else INTERSECTION_TEST_VERTEX(p1, q1, r1, r2, p2, q2)
	}

	int tri_tri_overlap_test_2d(
		const pmp::Scalar p1[2], const pmp::Scalar q1[2], const pmp::Scalar r1[2], 
		const pmp::Scalar p2[2], const pmp::Scalar q2[2], const pmp::Scalar r2[2])
	{
		if (ORIENT_2D(p1, q1, r1) < 0.0f)
		{
			if (ORIENT_2D(p2, q2, r2) < 0.0f)
			{
				return ccw_tri_tri_intersection_2d(p1, r1, q1, p2, r2, q2);
			}
			return ccw_tri_tri_intersection_2d(p1, r1, q1, p2, q2, r2);
		}
		if (ORIENT_2D(p2, q2, r2) < 0.0f)
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
	  if (DOT(v1,N1) > 0.0f) return 0;\
	  SUB(v1,p2,p1);\
	  SUB(v2,r1,p1);\
	  CROSS(N1,v1,v2);\
	  SUB(v1,r2,p1);\
	  if (DOT(v1,N1) > 0.0f) return 0;\
	  else return 1; }

#define TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
	  if (dp2 > 0.0f) { \
	     if (dq2 > 0.0f) CHECK_MIN_MAX(p1,r1,q1,r2,p2,q2) \
	     else if (dr2 > 0.0f) CHECK_MIN_MAX(p1,r1,q1,q2,r2,p2)\
	     else CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2) }\
	  else if (dp2 < 0.0f) { \
	    if (dq2 < 0.0f) CHECK_MIN_MAX(p1,q1,r1,r2,p2,q2)\
	    else if (dr2 < 0.0f) CHECK_MIN_MAX(p1,q1,r1,q2,r2,p2)\
	    else CHECK_MIN_MAX(p1,r1,q1,p2,q2,r2)\
	  } else { \
	    if (dq2 < 0.0f) { \
	      if (dr2 >= 0.0f)  CHECK_MIN_MAX(p1,r1,q1,q2,r2,p2)\
	      else CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2)\
	    } \
	    else if (dq2 > 0.0f) { \
	      if (dr2 > 0.0f) CHECK_MIN_MAX(p1,r1,q1,p2,q2,r2)\
	      else  CHECK_MIN_MAX(p1,q1,r1,q2,r2,p2)\
	    } \
	    else  { \
	      if (dr2 > 0.0f) CHECK_MIN_MAX(p1,q1,r1,r2,p2,q2)\
	      else if (dr2 < 0.0f) CHECK_MIN_MAX(p1,r1,q1,r2,p2,q2)\
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
		SUB(v1, p2, r2);
		SUB(v2, q2, r2);
		CROSS(N2, v1, v2);

		SUB(v1, p1, r2);
		pmp::Scalar dp1 = DOT(v1, N2);
		SUB(v1, q1, r2);
		pmp::Scalar dq1 = DOT(v1, N2);
		SUB(v1, r1, r2);
		pmp::Scalar dr1 = DOT(v1, N2);

		if (((dp1 * dq1) > 0.0f) && ((dp1 * dr1) > 0.0f))  return 0;
		/* Compute distance signs  of p2, q2 and r2 to the plane of
		   triangle(p1,q1,r1) */
		SUB(v1, q1, p1);
		SUB(v2, r1, p1);
		CROSS(N1, v1, v2);

		SUB(v1, p2, r1);
		pmp::Scalar dp2 = DOT(v1, N1);
		SUB(v1, q2, r1);
		pmp::Scalar dq2 = DOT(v1, N1);
		SUB(v1, r2, r1);
		pmp::Scalar dr2 = DOT(v1, N1);
		if (((dp2 * dq2) > 0.0f) && ((dp2 * dr2) > 0.0f)) return 0;

		/* Permutation in a canonical form of T1's vertices */
		if (dp1 > 0.0f) 
		{
			if (dq1 > 0.0f) TRI_TRI_3D(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2)
			if (dr1 > 0.0f) TRI_TRI_3D(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2)
			TRI_TRI_3D(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2)
		}
		if (dp1 < 0.0f)
		{
			if (dq1 < 0.0f) TRI_TRI_3D(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2)
			if (dr1 < 0.0f) TRI_TRI_3D(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2)
			TRI_TRI_3D(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2)
		}
		if (dq1 < 0.0f) 
		{
			if (dr1 >= 0.0f) TRI_TRI_3D(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2)
			TRI_TRI_3D(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2)
		}
		if (dq1 > 0.0f)
		{
			if (dr1 > 0.0f) TRI_TRI_3D(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2)
			TRI_TRI_3D(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2)
		}
		if (dr1 > 0.0f) TRI_TRI_3D(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2)
		if (dr1 < 0.0f) TRI_TRI_3D(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2)
		return coplanar_tri_tri3d(p1, q1, r1, p2, q2, r2, N1);
	}

	bool TriangleIntersectsTriangle(const std::vector<pmp::vec3>& vertices0, const std::vector<pmp::vec3>& vertices1)
	{
		return tri_tri_overlap_test_3d(
			vertices0[0].data(), vertices0[1].data(),vertices0[2].data(),
			vertices1[0].data(), vertices1[1].data(), vertices1[2].data()) > 0;
	}

	/// \brief intersection tolerance for Moller-Trumbore algorithm.
	constexpr float MT_INTERSECTION_EPSILON = 1e-6f;

	bool RayIntersectsTriangle(Ray& ray, const std::vector<pmp::vec3>& triVertices)
	{
		pmp::vec3 edge1 = triVertices[1] - triVertices[0];
		pmp::vec3 edge2 = triVertices[2] - triVertices[0];
		pmp::vec3 cross1;
		CROSS(cross1, ray.Direction, edge2);
		const float det = DOT(edge1, cross1);
		if (det > -MT_INTERSECTION_EPSILON && det < MT_INTERSECTION_EPSILON)
		{
			return false; // This ray is parallel to this triangle.
		}

		const float invDet = 1.0f / det;
		const auto startToTri0 = ray.StartPt - triVertices[0];
		const float u = invDet * DOT(startToTri0, cross1);
		if (u < 0.0f || u > 1.0f)
		{
			return false;
		}
		pmp::vec3 cross2;
		CROSS(cross2, startToTri0, edge1);
		const float v = invDet * DOT(ray.Direction, cross2);
		if (v < 0.0f || u + v > 1.0f) 
		{
			return false;
		}

		// At this stage we can compute t to find out where the intersection point is on the line.
		const float t = invDet * DOT(edge2, cross2);
		if (t > MT_INTERSECTION_EPSILON && t < 1.0f / MT_INTERSECTION_EPSILON &&
			t >= ray.ParamMin && t <= ray.ParamMax)
		{
			ray.HitParam = t;
			return true;
		}
		// there is a line intersection but not a ray intersection.
		return false;
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
		const float dirLenSq = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2];
		if (dirLenSq >= 1.0f + FLT_EPSILON || dirLenSq <= 1.0f - FLT_EPSILON)
		{
			throw std::logic_error("Ray::Ray: ||dir|| != 1 ! Ray direction vector must be normalized!\n");
		}

		InvDirection[0] = 1.0f / dir[0];
		InvDirection[1] = 1.0f / dir[1];
		InvDirection[2] = 1.0f / dir[2];

		// calculate id of max dimension of ray direction.
		const pmp::vec3 absDir{ std::fabs(dir[0]), std::fabs(dir[1]), std::fabs(dir[2]) };
		kz = MaxDim(absDir);
		//kx = kz + 1; if (kx == 3) kx = 0;
		//ky = kx + 1; if (ky == 3) ky = 0;
		kx = (kz + 1) % 3;
		ky = (kx + 1) % 3;
		// swap kx and ky dims to preserve winding direction of triangles.
		if (dir[kz] < 0.0f) std::swap(kx, ky);
		// calculate shear constants
		Sx = dir[kx] / dir[kz];
		Sy = dir[ky] / dir[kz];
		Sz = 1.0f / dir[kz];
	}

	/// \brief if true, we use backface-culling (i.e. skipping triangles whose normals point away from the ray).
#define BACKFACE_CULLING false
	
	bool RayIntersectsTriangleWatertight(Ray& ray, const std::vector<pmp::vec3>& triVertices)
	{
		// actual alg
		const pmp::vec3 A = triVertices[0] - ray.StartPt;
		const pmp::vec3 B = triVertices[1] - ray.StartPt;
		const pmp::vec3 C = triVertices[2] - ray.StartPt;
		const float Ax = A[ray.kx] - ray.Sx * A[ray.kz];
		const float Ay = A[ray.ky] - ray.Sy * A[ray.kz];
		const float Bx = B[ray.kx] - ray.Sx * B[ray.kz];
		const float By = B[ray.ky] - ray.Sy * B[ray.kz];
		const float Cx = C[ray.kx] - ray.Sx * C[ray.kz];
		const float Cy = C[ray.ky] - ray.Sy * C[ray.kz];
		float U = Cx * By - Cy * Bx;
		float V = Ax * Cy - Ay * Cx;
		float W = Bx * Ay - By * Ax;
		if (U == 0.0f || V == 0.0f || W == 0.0f) 
		{
			const double CxBy = static_cast<double>(Cx) * static_cast<double>(By);
			const double CyBx = static_cast<double>(Cy) * static_cast<double>(Bx);
			U = static_cast<float>(CxBy - CyBx);
			const double AxCy = static_cast<double>(Ax) * static_cast<double>(Cy);
			const double AyCx = static_cast<double>(Ay) * static_cast<double>(Cx);
			V = static_cast<float>(AxCy - AyCx);
			const double BxAy = static_cast<double>(Bx) * static_cast<double>(Ay);
			const double ByAx = static_cast<double>(By) * static_cast<double>(Ax);
			W = static_cast<float>(BxAy - ByAx);
		}
#if BACKFACE_CULLING
		if (U < 0.0f || V < 0.0f || W < 0.0f) return false;
#else
		if ((U < 0.0f || V < 0.0f || W < 0.0f) &&
		   (U > 0.0f || V > 0.0f || W > 0.0f)) return false;
#endif
		const float det = U + V + W;
		if (det == 0.0f) return false;
		const float Az = ray.Sz * A[ray.kz];
		const float Bz = ray.Sz * B[ray.kz];
		const float Cz = ray.Sz * C[ray.kz];
		const float T = U * Az + V * Bz + W * Cz;
		if (T < 0.0f || T > ray.HitParam * det)
			return false;
		
		const float rcpDet = 1.0f / det; // reciprocal det
		/*const float hitBCoordU = U * rcpDet;
		const float hitBCoordV = V * rcpDet;
		const float hitBCoordW = W * rcpDet;*/
		ray.HitParam = T * rcpDet;
		
		return true;
	}

	// =========================================================================

	constexpr unsigned int idVec[3] = { 0, 1, 2 };

	// accelerated roundoff for watertightness [Woop, Benthin, Wald, 2013, p. 70]
	float p = 1.0f + 2e-23f;
	float m = 1.0f - 2e-23f;
	[[nodiscard]] float up(const float a) { return a > 0.0f ? a * p : a * m; }
	[[nodiscard]] float dn(const float a) { return a > 0.0f ? a * m : a * p; }
	[[nodiscard]] float Up(const float a) { return a * p; }
	[[nodiscard]] float Dn(const float a) { return a * m; }
	constexpr float eps = 5.0f * 2e-24f;

	bool RayIntersectsABox(const Ray& ray, const pmp::BoundingBox& box)
	{
		int nearX = static_cast<int>(idVec[ray.kx]), farX = static_cast<int>(idVec[ray.kx]);
		int nearY = static_cast<int>(idVec[ray.ky]), farY = static_cast<int>(idVec[ray.ky]);
		int nearZ = static_cast<int>(idVec[ray.kz]), farZ = static_cast<int>(idVec[ray.kz]);
		if (ray.Direction[ray.kx] < 0.0f) std::swap(nearX, farX);
		if (ray.Direction[ray.ky] < 0.0f) std::swap(nearY, farY);
		if (ray.Direction[ray.kz] < 0.0f) std::swap(nearZ, farZ);

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
		const float max_z = std::max(lower[ray.kz], upper[ray.kz]);
		const float err_near_x = Up(lower[ray.kx] + max_z);
		const float err_near_y = Up(lower[ray.ky] + max_z);
		float start_near_x = up(ray.StartPt[ray.kx] + Up(eps * err_near_x));
		float start_near_y = up(ray.StartPt[ray.ky] + Up(eps * err_near_y));
		const float start_near_z = ray.StartPt[ray.kz];
		const float err_far_x = Up(upper[ray.kx] + max_z);
		const float err_far_y = Up(upper[ray.ky] + max_z);
		float start_far_x = dn(ray.StartPt[ray.kx] - Up(eps * err_far_x));
		float start_far_y = dn(ray.StartPt[ray.ky] - Up(eps * err_far_y));
		const float start_far_z = ray.StartPt[ray.kz];
		if (ray.Direction[ray.kx] < 0.0f) std::swap(start_near_x, start_far_x);
		if (ray.Direction[ray.ky] < 0.0f) std::swap(start_near_y, start_far_y);
		const float rdir_near_x = Dn(Dn(ray.InvDirection[ray.kx]));
		const float rdir_near_y = Dn(Dn(ray.InvDirection[ray.ky]));
		const float rdir_near_z = Dn(Dn(ray.InvDirection[ray.kz]));
		const float rdir_far_x = Up(Up(ray.InvDirection[ray.kx]));
		const float rdir_far_y = Up(Up(ray.InvDirection[ray.ky]));
		const float rdir_far_z = Up(Up(ray.InvDirection[ray.kz]));
		float tNearX = (box.min()[nearX] - start_near_x) * rdir_near_x;
		float tNearY = (box.min()[nearY] - start_near_y) * rdir_near_y;
		float tNearZ = (box.min()[nearZ] - start_near_z) * rdir_near_z;
		float tFarX = (box.max()[farX] - start_far_x) * rdir_far_x;
		float tFarY = (box.max()[farY] - start_far_y) * rdir_far_y;
		float tFarZ = (box.max()[farZ] - start_far_z) * rdir_far_z;
		const float tNear = std::max({ tNearX, tNearY, tNearZ, ray.ParamMin });
		const float tFar = std::min({ tFarX, tFarY, tFarZ, ray.ParamMax });
		return tNear <= tFar;
	}

} // namespace Geometry