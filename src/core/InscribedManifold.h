#pragma once

#include "pmp/Types.h"
#include "pmp/MatVec.h"
#include "pmp/ManifoldCurve2D.h"
#include "pmp/SurfaceMesh.h"
#include "pmp/algorithms/CurveFactory.h"

#include "geometry/Grid.h"
#include "geometry/IcoSphereBuilder.h"

/// \brief A wrapper for input data for the InscribedCircleBuilder
struct InscribedCircleInputData
{
	std::vector<pmp::Point2> Points{}; //>! the evaluated point cloud.
	std::shared_ptr<Geometry::ScalarGrid2D> DistanceField{ nullptr }; //>! a pointer to a pre-computed distance field to Points. 
};

/// \brief Parameters of a 2D circle.
struct Circle2D
{
	pmp::Point2 Center{};
	pmp::Scalar Radius{ -1.0f };
};

/// \brief A base utility to calculate the centers and radii of circles inscribed to a point cloud.
class InscribedCircleCalculator
{
public:
	virtual ~InscribedCircleCalculator() = default;

	/**
	 * \brief Estimates inscribed circles to a point cloud.
	 * \param data        input point cloud data.
	 * \return a vector of circles.
	 */
	virtual [[nodiscard]] std::vector<Circle2D> Calculate(const InscribedCircleInputData& data) = 0;
};

/// \brief Calculates the centers and radii of circles inscribed to a point cloud using the naive approach:
/// Use the center of the bounding box of the input point cloud, and the distance to the closest point as radius.
class NaiveInscribedCircleCalculator : public InscribedCircleCalculator
{
public:
	/// \copydoc InscribedCircleCalculator::Calculate
	[[nodiscard]] std::vector<Circle2D> Calculate(const InscribedCircleInputData& data) override;
};

/// \brief Calculates the centers and radii of circles inscribed to a point cloud using the distance-field-based approach:
/// Find the approximate locations of local maxima of the point cloud distance field which will serve as centers, and the distance to the closest point as radii.
class DistanceFieldInscribedCircleCalculator : public InscribedCircleCalculator
{
public:
	/// \copydoc InscribedCircleCalculator::Calculate
	[[nodiscard]] std::vector<Circle2D> Calculate(const InscribedCircleInputData& data) override;
};

/// \brief Calculates the centers and radii of circles inscribed to a point cloud using the distance-field-based approach:
/// Find the approximate locations of local maxima of the point cloud distance field which will serve as centers, and the distance to the closest point as radii.
/// In this calculator, we use a quadtree-based approach instead of analyzing all grid points
class HierarchicalDistanceFieldInscribedCircleCalculator : public InscribedCircleCalculator
{
public:
	/// \copydoc InscribedCircleCalculator::Calculate
	[[nodiscard]] std::vector<Circle2D> Calculate(const InscribedCircleInputData& data) override;
};

/// \brief Calculates the centers and radii of circles inscribed to a point cloud using the distance-field-based approach:
/// Find the approximate locations of local maxima of the point cloud distance field which will serve as centers, and the distance to the closest point as radii.
/// In this calculator, we use a particle-based approach: simulating the movement of points across the grid towards the local maxima.
class ParticleSwarmDistanceFieldInscribedCircleCalculator : public InscribedCircleCalculator
{
public:
	/// \copydoc InscribedCircleCalculator::Calculate
	[[nodiscard]] std::vector<Circle2D> Calculate(const InscribedCircleInputData& data) override;
};

/**
 * \brief Tessellates the resulting curve from the computed radius and center.
 * \param circle        A parametric circle to be reconstructed.
 * \param nSegments     the number of segments to be used in the reconstruction.
 * \return The tessellation of the resulting circle (must be called after CalculateInnerCircleRadiusAndCenter).
 */
inline [[nodiscard]] pmp::ManifoldCurve2D ConstructCircle(const Circle2D& circle, size_t nSegments = 32)
{
	return pmp::CurveFactory::circle(circle.Center, circle.Radius, nSegments);
}

// =====================================================================================
//                                  3D Utils
// -------------------------------------------------------------------------------------

/// \brief A wrapper for input data for the construction of an inscribed sphere
struct InscribedSphereInputData
{
	std::vector<pmp::Point> Points{}; //>! the evaluated point cloud.
	std::shared_ptr<Geometry::ScalarGrid> DistanceField{ nullptr }; //>! a pointer to a pre-computed distance field to Points. 
};

/// \brief Parameters of a 3D sphere.
struct Sphere3D
{
	pmp::Point Center{};
	pmp::Scalar Radius{ -1.0f };
};

/// \brief A base utility to calculate the centers and radii of spheres inscribed to a point cloud.
class InscribedSphereCalculator
{
public:
	virtual ~InscribedSphereCalculator() = default;

	/**
	 * \brief Estimates inscribed spheres to a point cloud.
	 * \param data        input point cloud data.
	 * \return a vector of spheres.
	 */
	virtual [[nodiscard]] std::vector<Sphere3D> Calculate(const InscribedSphereInputData& data) = 0;
};

/// \brief Calculates the centers and radii of spheres inscribed to a point cloud using the naive approach:
/// Use the center of the bounding box of the input point cloud, and the distance to the closest point as radius.
class NaiveInscribedSphereCalculator : public InscribedSphereCalculator
{
public:
	/// \copydoc InscribedSphereCalculator::Calculate
	[[nodiscard]] std::vector<Sphere3D> Calculate(const InscribedSphereInputData& data) override;
};

/// \brief Calculates the centers and radii of spheres inscribed to a point cloud using the distance-field-based approach:
/// Find the approximate locations of local maxima of the point cloud distance field which will serve as centers, and the distance to the closest point as radii.
class DistanceFieldInscribedSphereCalculator : public InscribedSphereCalculator
{
public:
	/// \copydoc InscribedSphereCalculator::Calculate
	[[nodiscard]] std::vector<Sphere3D> Calculate(const InscribedSphereInputData& data) override;
};

/// \brief Calculates the centers and radii of spheres inscribed to a point cloud using the distance-field-based approach:
/// Find the approximate locations of local maxima of the point cloud distance field which will serve as centers, and the distance to the closest point as radii.
/// In this calculator, we use a octree-based approach instead of analyzing all grid points
class HierarchicalDistanceFieldInscribedSphereCalculator : public InscribedSphereCalculator
{
public:
	/// \copydoc InscribedSphereCalculator::Calculate
	[[nodiscard]] std::vector<Sphere3D> Calculate(const InscribedSphereInputData& data) override;
};

/// \brief Calculates the centers and radii of spheres inscribed to a point cloud using the distance-field-based approach:
/// Find the approximate locations of local maxima of the point cloud distance field which will serve as centers, and the distance to the closest point as radii.
/// In this calculator, we use a particle-based approach: simulating the movement of points across the grid towards the local maxima.
class ParticleSwarmDistanceFieldInscribedSphereCalculator : public InscribedSphereCalculator
{
public:
	/// \copydoc InscribedSphereCalculator::Calculate
	[[nodiscard]] std::vector<Sphere3D> Calculate(const InscribedSphereInputData& data) override;
};

/**
 * \brief Tessellates the resulting curve from the computed radius and center.
 * \param sphere        A parametric sphere to be reconstructed.
 * \param subdiv        the level of subdivision used for construction
 * \return The tessellation of the resulting Sphere (must be called after CalculateInnerSphereRadiusAndCenter).
 */
inline [[nodiscard]] pmp::SurfaceMesh ConstructSphere(const Sphere3D& sphere, unsigned int subdiv = 2)
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
