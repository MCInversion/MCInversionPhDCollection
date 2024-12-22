#pragma once

#include "pmp/Types.h"
#include "pmp/SurfaceMesh.h"
#include "pmp/ManifoldCurve2D.h"

#include "GeometryConversionUtils.h"

#include <optional>


namespace Geometry
{
	// forward declarations
	class ScalarGrid;
	class CollisionKdTree;

	/// \brief Computes minimum internal angle per triangle averaged for each vertex over adjacent triangles
	///        & stores the values as vertex scalar data.
	///        This metric takes values from [0, M_PI / 3]. Preferred values are from [M_PI / 6, M_PI / 3].
	[[nodiscard]] bool ComputeTriangleMinAngleVertexValues(pmp::SurfaceMesh& mesh);

	/// \brief Computes maximum internal angle per triangle averaged for each vertex over adjacent triangles
	///        & stores the values as vertex scalar data.
	///        This metric takes values from [M_PI / 3, M_PI]. Preferred values are from [M_PI / 3, M_PI / 2].
	[[nodiscard]] bool ComputeTriangleMaxAngleVertexValues(pmp::SurfaceMesh& mesh);

	/// \brief Computes the condition number of each triangle's Jacobian averaged for each vertex over adjacent triangles
	///        & stores the values as vertex scalar data.
	///        This metric takes values from [1, infinity). Preferred values are from [1, 1.3].
	[[nodiscard]] bool ComputeTriangleJacobianConditionNumberVertexValues(pmp::SurfaceMesh& mesh);

	/// \brief Computes the condition number of each triangle's equilateral Jacobian averaged for each vertex over adjacent triangles
	///        & stores the values as vertex scalar data.
	///        This metric takes values from [1, infinity). Preferred values are from [1, 1.3].
	[[nodiscard]] bool ComputeEquilateralTriangleJacobianConditionNumbers(pmp::SurfaceMesh& mesh);

	/// \brief Computes the conditioning of the stiffness matrix of a spring model of triangles according to [ch. 3, Schewchuk, 2002]
	///        & stores the values as vertex scalar data.
	[[nodiscard]] bool ComputeStiffnessMatrixConditioningVertexValues(pmp::SurfaceMesh& mesh);

	// ============================================================================================================
	//                                      Exposed quality functions
	// ------------------------------------------------------------------------------------------------------------

	/// \brief Computes the condition number of a equilateral Jacobian of this triangle.
	///        Note: The Jacobian corresponds to the planar transformation from triangle (-0.5, 0), (0.5, 0), (0, 1) to the evaluated triangle.
	[[nodiscard]] pmp::Scalar GetConditionNumberOfEquilateralTriangleJacobian(const pmp::SurfaceMesh& mesh, const pmp::Face& face);

	/// \brief Preferred min and max values for the equilateral Jacobian condition number.
	constexpr pmp::Scalar JACOBIAN_COND_MIN = 1.0;
	constexpr pmp::Scalar JACOBIAN_COND_MAX = 1.5;

	// ============================================================================================================

	/// \brief verifies whether a metric with a given name is registered.
	[[nodiscard]] bool IsMetricRegistered(const std::string& name);

	/// \brief triangle metric computation function.
	using TriMetricFunction = std::function<bool(pmp::SurfaceMesh&)>;

	/// \brief provides a metric function for a given metricName.
	[[nodiscard]] TriMetricFunction IdentifyMetricFunction(const std::string& metricName);

	/// \brief Computes dihedral angle for mesh vertices, averaged from mesh edges.
	void ComputeEdgeDihedralAngles(pmp::SurfaceMesh& mesh);

	/// \brief Computes vertex principal curvatures and mean curvature + properties related to feature detection.
	void ComputeVertexCurvaturesAndRelatedProperties(pmp::SurfaceMesh& mesh, const pmp::Scalar& principalCurvatureFactor = 2.0);

	/// \brief Computes z-level elevation values for each vertex.
	void ComputeZLevelElevations(pmp::SurfaceMesh& mesh);

	/// \brief Computes edge length range and average.
	[[nodiscard]] std::tuple<pmp::Scalar, pmp::Scalar, pmp::Scalar> ComputeEdgeLengthMinAverageAndMax(const pmp::SurfaceMesh& mesh);

	/// \brief Counts all unreferenced vertices in BaseMeshGeometryData.
	[[nodiscard]] size_t CountUnreferencedVertices(const BaseMeshGeometryData& data);

	/// \brief Computes saliency as vertex property according to [Lee, et al., 2005]
	///	\param[in] mesh               input mesh.
	///	\param[in] forcedVariance     if positive, this value will be used as basis for the sigmas in saliency evaluation.
	///	\param[in] normalizeValues    if true values will be normalized.
	[[nodiscard]] bool EvaluatePMPSurfaceMeshSaliency(pmp::SurfaceMesh& mesh, const double& forcedVariance = -1.0, const bool& normalizeValues = false);

	/// \brief Prints the evaluated histogram data.
	void PrintHistogramResultData(const std::pair<std::pair<pmp::Scalar, pmp::Scalar>, std::vector<unsigned int>>& histData, std::ostream& os);

	/// \brief A test function for subdivision mesh counts estimation.
	/// See: "Cavarga, Mesh Primitive Counting Formula for Subdivision Surfaces, SCG 2023".
	[[nodiscard]] std::pair<std::vector<size_t>, std::vector<size_t>> GetEdgeVertCountsTheoreticalEstimate(const pmp::SurfaceMesh& mesh, const size_t& maxSubdivLevel, const bool& evalOutput = false);

	/// \brief Counts self-intersecting faces of the input mesh.
	/// \param[in] mesh               the evaluated SurfaceMesh.
	///	\param[in] setFaceProperty    if true, the evaluated mesh will be given a face property identifying the self-intersecting faces.
	///	\return the total number of faces which intersect another face.
	[[nodiscard]] size_t CountPMPSurfaceMeshSelfIntersectingFaces(pmp::SurfaceMesh& mesh, const bool& setFaceProperty = false);

	/// \brief A fast verification for the presence of self-intersecting faces of the input mesh.
	[[nodiscard]] bool PMPSurfaceMeshHasSelfIntersections(const pmp::SurfaceMesh& mesh);

	/// \brief Converts a face property to an averaged scalar vertex property.
	void ConvertPMPSurfaceMeshBoolFacePropertyToScalarVertexProperty(pmp::SurfaceMesh& mesh, const std::string& propName);

	/// \brief Extracts self-intersection multimap of each face to all intersecting faces.: f -> {f1, f2, ... , fN_intersecting}
	[[nodiscard]] std::unordered_multimap<unsigned int, unsigned int> ExtractPMPSurfaceMeshFaceIntersectionMultimap(const pmp::SurfaceMesh& mesh);

	/// \brief Computes self-intersection polylines from intersection lines for each triangle-triangle intersection.
	[[nodiscard]] std::vector<std::vector<pmp::vec3>> ComputeSurfaceMeshSelfIntersectionPolylines(const pmp::SurfaceMesh& mesh);

	/// \brief A fast verification for the presence of self-intersecting faces of the input manifold curve.
	[[nodiscard]] bool PMPManifoldCurve2DHasSelfIntersections(const pmp::ManifoldCurve2D& curve);

	/// \brief A ray-casting verification whether a given 2D point is inside a ManifoldCurve2D.
	[[nodiscard]] bool IsPointInsidePMPManifoldCurve(const pmp::Point2& point, const pmp::ManifoldCurve2D& curve);

	/// \brief A ray-casting verification whether a given 3D point is inside a SurfaceMesh.
	[[nodiscard]] bool IsPointInsidePMPSurfaceMesh(const pmp::Point& point, const pmp::SurfaceMesh& mesh);

	/// \brief A ray-casting verification whether a given 3D point is inside a SurfaceMesh. The mesh is already converted to a kd-tree.
	[[nodiscard]] bool IsPointInsidePMPSurfaceMesh(const pmp::Point& point, const std::shared_ptr<CollisionKdTree>& meshKdTree);

	/// \brief Uses an approximation via a tangent quadric (ellipsoid/hyperboloid) for determining the deviations of triangle faces adjacent to vertex v.
	[[nodiscard]] pmp::Scalar CalculateQuadricApproximationErrorAtVertex(const pmp::SurfaceMesh& mesh, pmp::Vertex v);

	/// \brief Uses circular approximation for determining the deviations of edges adjacent to vertex v.
	[[nodiscard]] pmp::Scalar CalculateCircularApproximationErrorAtVertex(const pmp::ManifoldCurve2D& curve, pmp::Vertex v);

	/// \brief Calculates the signed area of input curve.
	/// \throw std::invalid_argument if curve has self-intersections, is open, or empty.
	[[nodiscard]] pmp::Scalar CalculateSignedAreaOfASimpleClosedCurve(const pmp::ManifoldCurve2D& curve);

	/// \brief face quality metric computation function.
	using FaceQualityFunction = std::function<pmp::Scalar(const pmp::SurfaceMesh&, const pmp::Face&)>;

	/// \brief A simple wrapper for the range of values for FaceQualityFunction.
	struct FaceQualityRange
	{
		pmp::Scalar Min{ 0.0 };
		pmp::Scalar Max{ 1.0 };

		bool operator() (const pmp::Scalar& val) const
		{
			return val > Min && val < Max;
		}
	};

	/// \brief A debug utility to print curve values in topological order (based on connectivity)
	void PrintCurveValuesInTopologicalOrder(const pmp::ManifoldCurve2D& curve, const std::vector<pmp::Scalar>& values, std::ostream& os);

	/// \brief A utility for calculating the minimum dimension of planar curve (used primarily for determining distance field cell size).
	[[nodiscard]] pmp::Scalar GetCurveBoundsMinDimension(const pmp::ManifoldCurve2D& curve);

	/// \brief A utility for calculating the maximum dimension of planar curve.
	[[nodiscard]] pmp::Scalar GetCurveBoundsMaxDimension(const pmp::ManifoldCurve2D& curve);

	/// \brief A utility for calculating the minimum and maximum dimensions of planar curve.
	[[nodiscard]] std::pair<pmp::Scalar, pmp::Scalar> GetCurveBoundsMinMaxDimensions(const pmp::ManifoldCurve2D& curve);

} // namespace Geometry