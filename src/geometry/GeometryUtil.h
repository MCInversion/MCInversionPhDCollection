#pragma once

#include <optional>

#include "pmp/MatVec.h"
#include "pmp/Types.h"

// forward declarations
namespace pmp
{
	class BoundingBox;
}

namespace SDF
{
	struct Ray;
}

namespace Geometry
{
	/**
	 * \brief Computes the squared distance from point to a triangle.
	 * \param vertices     list of (three) vertices of a triangle.
	 * \param point        point from which the distance is to be computed.
	 * \return squared distance from point to triangle.
	 */
	[[nodiscard]] double GetDistanceToTriangleSq(const std::vector<pmp::vec3>& vertices, const pmp::vec3& point);

	/**
	 * \brief An intersection test between a triangle and a box.
	 * \param vertices     list of (three) vertices of a triangle.
	 * \param boxCenter    center of the box to be queried.
	 * \param boxHalfSize  half-size of the box to be queried.
	 * \return true if the box intersects the triangle.
	 */
	[[nodiscard]] bool TriangleIntersectsBox(const std::vector<pmp::vec3>& vertices, const pmp::vec3& boxCenter, const pmp::vec3& boxHalfSize);

	/**
	 * \brief Computes the squared distance from point to a 2D line segment.
	 * \param vertices     list of (two) vertices of a 2D line segment.
	 * \param point        point from which the distance is to be computed.
	 * \return squared distance from point to the 2D line segment.
	 */
	[[nodiscard]] double GetDistanceToLine2DSq(const std::vector<pmp::Point2>& vertices, const pmp::vec2& point);

	/**
	 * \brief An intersection test between a 2D line segment and a box.
	 * \param vertices     list of (three) vertices of a  2D line segment.
	 * \param boxCenter    center of the box to be queried.
	 * \param boxHalfSize  half-size of the box to be queried.
	 * \return true if the box intersects the 2D line segment.
	 */
	[[nodiscard]] bool Line2DIntersectsBox(const std::vector<pmp::Point2>& vertices, const pmp::Point2& boxCenter, const pmp::vec2& boxHalfSize);

	/**
	 * \brief An intersection test between a two triangles
	 * \param vertices0    first triangle vertices list.
	 * \param vertices1    second triangle vertices list.
	 * \return true if the triangles intersect.
	 */
	[[nodiscard]] bool TriangleIntersectsTriangle(const std::vector<pmp::vec3>& vertices0, const std::vector<pmp::vec3>& vertices1);

	/**
	 * \brief An intersection test between a two 2D line segments
	 * \param vertices0    first line segment endpoints list.
	 * \param vertices1    second line segment endpoints list.
	 * \return true if the line segments intersect.
	 */
	[[nodiscard]] bool Line2DIntersectsLine2D(const std::vector<pmp::Point2>& vertices0, const std::vector<pmp::Point2>& vertices1);

	/**
	 * \brief A utility that returns an intersector line between two triangles.
	 * \param vertices0    first triangle vertices list.
	 * \param vertices1    second triangle vertices list.
	 * \return optional pair { start pt, end pt } of the intersection line. If std::nullopt, triangles do not intersect.
	 */
	[[nodiscard]] std::optional<std::pair<pmp::vec3, pmp::vec3>> ComputeTriangleTriangleIntersectionLine(const std::vector<pmp::vec3>& vertices0, const std::vector<pmp::vec3>& vertices1);

	// ======================================================================

	/// \brief a wrapper for the parameters of a ray intersecting KD-tree boxes.
	struct Ray
	{
		pmp::vec3 StartPt{};
		pmp::vec3 Direction{ 1.0f, 0.0f, 0.0f };
		pmp::vec3 InvDirection{ 1.0f, FLT_MAX, FLT_MAX };
		float ParamMin{ 0.0f };
		float ParamMax{ FLT_MAX };
		float HitParam{ FLT_MAX };

		// dir vector dimension indices:
		unsigned int kx{ 0 };
		unsigned int ky{ 0 };
		unsigned int kz{ 0 };

		// shear constants:
		float Sx{ 0.0f };
		float Sy{ 0.0f };
		float Sz{ 0.0f };

		/**
		 * \brief Constructor. Construct from startPt and direction vector as in [Woop, Benthin, Wald, 2013].
		 * \param startPt    start point of the ray.
		 * \param dir        (normalized) direction vector of the ray.
		 * \throw std::logic_error if dir is not normalized.
		 */
		Ray(const pmp::vec3& startPt, const pmp::vec3& dir);
	};

	/**
	 * \brief An intersection test between a triangle and a ray. [Moller, Trumbore, 1997].
	 * \param ray           intersecting ray (with modifiable hit param).
	 * \param triVertices   list of (three) vertices of a triangle.
	 * \return true if the ray intersects the triangle.
	 */
	[[nodiscard]] bool RayIntersectsTriangle(Ray& ray, const std::vector<pmp::vec3>& triVertices);

	/**
	 * \brief A watertight intersection test between a triangle and a ray. [Woop, Benthin, Wald, 2013].
	 * \param ray           intersecting ray (with modifiable hit param).
	 * \param triVertices   list of (three) vertices of a triangle.
	 * \return true if the ray intersects the triangle.
	 */
	[[nodiscard]] bool RayIntersectsTriangleWatertight(Ray& ray, const std::vector<pmp::vec3>& triVertices);

	/**
	 * \brief Ray-box intersection test.
	 * \param ray           intersecting ray.
	 * \param box           intersected box.
	 * \return true if the ray intersects the box.
	 */
	[[nodiscard]] bool RayIntersectsABox(const Ray& ray, const pmp::BoundingBox& box);

	/**
	 * \brief Verifies whether two planar circles intersect.
	 * \param center1       center of the first circle.
	 * \param radius1       radius of the first circle.
	 * \param center2       center of the second circle.
	 * \param radius2       radius of the second circle.
	 * \return true if circles intersect.
	 */
	[[nodiscard]] bool CircleIntersectsCircle2D(const pmp::Point2& center1, const pmp::Scalar& radius1, const pmp::Point2& center2, const pmp::Scalar& radius2);

	/**
	 * \brief Verifies whether two spheres intersect.
	 * \param center1       center of the first sphere.
	 * \param radius1       radius of the first sphere.
	 * \param center2       center of the second sphere.
	 * \param radius2       radius of the second sphere.
	 * \return true if spheres intersect.
	 */
	[[nodiscard]] bool SphereIntersectsSphere3D(const pmp::Point& center1, const pmp::Scalar& radius1, const pmp::Point& center2, const pmp::Scalar& radius2);

} // namespace Geometry