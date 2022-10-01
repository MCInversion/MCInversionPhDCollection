#include "MeshAnalysis.h"

#include <set>

#include "pmp/algorithms/DifferentialGeometry.h"

namespace Geometry
{

	/// \brief Computes minimum angle of a triangle
	[[nodiscard]] float GetMinAngleOfTriangle(const pmp::SurfaceMesh& mesh, const pmp::Face& face)
	{
		const auto nVerts = std::distance(mesh.vertices(face).begin(), mesh.vertices(face).end());
		if (nVerts != 3)
			return -1.0f; // Error, face with three vertices is needed!

		float result = FLT_MAX;
		for (const auto h : mesh.halfedges(face))
		{	
			const auto hPrev = mesh.prev_halfedge(h);
			const auto v = mesh.from_vertex(h); // angle at this vertex

			const auto v0 = mesh.from_vertex(hPrev);
			const auto v1 = mesh.to_vertex(h);

			const auto e0 = mesh.position(v0) - mesh.position(v);
			const auto e1 = mesh.position(v1) - mesh.position(v);
			const float l0 = pmp::norm(e0);
			const float l1 = pmp::norm(e1);
			if (l0 < FLT_EPSILON || l1 < FLT_EPSILON)
				return 0.0f;

			const float dot01 = pmp::dot(e0, e1);
			const float vertAngle = acos(dot01 / (l0 * l1));

			if (vertAngle < result)
				result = vertAngle;
		}

		return result;
	}

	/// \brief Computes maximum angle of a triangle
	[[nodiscard]] float GetMaxAngleOfTriangle(const pmp::SurfaceMesh& mesh, const pmp::Face& face)
	{
		const auto nVerts = std::distance(mesh.vertices(face).begin(), mesh.vertices(face).end());
		if (nVerts != 3)
			return -1.0f; // Error, face with three vertices is needed!

		float result = -FLT_MAX;
		for (const auto h : mesh.halfedges(face))
		{
			const auto hPrev = mesh.prev_halfedge(h);
			const auto v = mesh.from_vertex(h); // angle at this vertex

			const auto v0 = mesh.from_vertex(hPrev);
			const auto v1 = mesh.to_vertex(h);

			const auto e0 = mesh.position(v0) - mesh.position(v);
			const auto e1 = mesh.position(v1) - mesh.position(v);
			const float l0 = pmp::norm(e0);
			const float l1 = pmp::norm(e1);
			if (l0 < FLT_EPSILON || l1 < FLT_EPSILON)
				return 0.0f;

			const float dot01 = pmp::dot(e0, e1);
			const float vertAngle = acos(dot01 / (l0 * l1));

			if (vertAngle > result)
				result = vertAngle;
		}

		return result;
	}

	/// \brief Computes the condition number of a Jacobian of this triangle.
	///        Note: The Jacobian corresponds to the planar transformation from unit triangle (0, 0), (1, 0), (0, 1) to the evaluated triangle.
	[[nodiscard]] float GetConditionNumberOfTriangleJacobian(const pmp::SurfaceMesh& mesh, const pmp::Face& face)
	{
		const auto nVerts = std::distance(mesh.vertices(face).begin(), mesh.vertices(face).end());
		if (nVerts != 3)
			return -1.0f; // Error, face with three vertices is needed!

		const auto h = mesh.halfedge(face);
		const auto hPrev = mesh.prev_halfedge(h);

		const auto v = mesh.from_vertex(h);
		const auto v0 = mesh.from_vertex(hPrev);
		const auto v1 = mesh.to_vertex(h);

		const auto e0 = mesh.position(v0) - mesh.position(v);
		const auto e1 = mesh.position(v1) - mesh.position(v);
		const auto e2 = pmp::cross(e0, e1); // tri plane normal

		const auto xVector = pmp::normalize(e0);
		const auto zVector = pmp::normalize(e2);
		const auto yVector = pmp::perp(xVector, zVector);		

		const float e0X = pmp::norm(e0);
		constexpr float e0Y = 0.0f;

		const float e1X = pmp::dot(e1, xVector);
		const float e1Y = pmp::dot(e1, yVector);

		const float detJ = e0X * e1Y - e1X * e0Y;
		if (std::fabs(detJ) < 1e-5f)
			return FLT_MAX; // singular Jacobian has an infinite condition number

		const pmp::mat3 J{
			e0X, e1X, 0.0f,
			e0Y, e1Y, 0.0f,
			0.0f, 0.0f, 1.0f
		};
		const auto JInv = pmp::inverse(J);

		/*const float jNorm1 = norm(J);
		const float jInvNorm1 = norm(JInv);*/

		// 1-norm of Jacobian and its inverse
		float jNorm1 = -FLT_MAX;
		float jInvNorm1 = -FLT_MAX;
		for (unsigned int i = 0; i < 2; i++)
		{
			float jColSum = 0.0f;
			float jInvColSum = 0.0f;
			for (unsigned int j = 0; j < 2; j++)
			{
				const unsigned int elementId = 3 * i + j;
				jColSum += std::fabs(J[elementId]);
				jInvColSum += std::fabs(JInv[elementId]);
			}
			if (jColSum > jNorm1)
				jNorm1 = jColSum;
			if (jInvColSum > jInvNorm1)
				jInvNorm1 = jInvColSum;
		}

		return jNorm1 * jInvNorm1;
	}

	/// \brief Computes Schewchuk's "stiffness matrix conditioning" number of a given triangle face. For more information see chapter 3 of [Schewchuk, 2002].
	[[nodiscard]] float ComputeStiffnessMatrixConditioningForTriangle(const pmp::SurfaceMesh& mesh, const pmp::Face& face)
	{
		const auto nVerts = std::distance(mesh.vertices(face).begin(), mesh.vertices(face).end());
		if (nVerts != 3)
			return -1.0f; // Error, face with three vertices is needed!

		const auto area_notNorm = pmp::triangle_area(mesh, face);
		if (area_notNorm < 1e-6f)
			return FLT_MAX;

		// TODO: normalization still produces complex roots, what to do about it?

		auto hc = mesh.halfedges(face);
		const auto e0 = mesh.edge(*hc);
		const auto e1 = mesh.edge(*(++hc));
		const auto e2 = mesh.edge(*(++hc));

		const auto l0Sq_notNorm = mesh.edge_length_sq(e0);
		const auto l1Sq_notNorm = mesh.edge_length_sq(e1);
		const auto l2Sq_notNorm = mesh.edge_length_sq(e2);

		const auto maxLenSq = std::max({ l0Sq_notNorm , l1Sq_notNorm, l2Sq_notNorm });

		const auto l0Sq = l0Sq_notNorm / maxLenSq;
		const auto l1Sq = l1Sq_notNorm / maxLenSq;
		const auto l2Sq = l2Sq_notNorm / maxLenSq;
		const auto area = area_notNorm / maxLenSq;

		const auto lSqSumSq = (l0Sq + l1Sq + l2Sq) * (l0Sq + l1Sq + l2Sq);

		return (l0Sq + l1Sq + l2Sq + sqrt(lSqSumSq - 48.0f * area)) / (8.0f * area);
	}

	using FaceMetricFunction = std::function<float(const pmp::SurfaceMesh&, const pmp::Face&)>;

	/**
	 * \brief Computes interpolated triangle metric for each mesh vertex & stores it as vertex property.
	 * \param mesh         input mesh to be analyzed.
	 * \param metricFunc   metric function to be used.
	 * \param metricName   name of the computed vertex property.
	 * \return true if computation was successful.
	 */
	static bool ComputeInterpolatedTrianglesMetric(pmp::SurfaceMesh& mesh, const FaceMetricFunction& metricFunc, const std::string& metricName)
	{
		if (!mesh.is_triangle_mesh())
			return false; // unable to process non-triangle meshes

		auto vProp = mesh.vertex_property<pmp::Scalar>("v:" + metricName, 0.0f);

		for (const auto v : mesh.vertices())
		{
			float mean = 0.0;
			for (const auto f : mesh.faces(v))
			{
				const auto val = metricFunc(mesh, f);
				if (val < 0.0f)
					return false; // Error, invalid triangle
				mean += val;
			}
			const auto nAdjacentFaces = std::distance(mesh.faces(v).begin(), mesh.faces(v).end());
			mean /= static_cast<float>(nAdjacentFaces);
			vProp[v] = mean;
		}
		return true;
	}

	bool ComputeTriangleMinAngleVertexValues(pmp::SurfaceMesh& mesh)
	{
		return ComputeInterpolatedTrianglesMetric(mesh, GetMinAngleOfTriangle, "minAngle");
	}

	bool ComputeTriangleMaxAngleVertexValues(pmp::SurfaceMesh& mesh)
	{
		return ComputeInterpolatedTrianglesMetric(mesh, GetMaxAngleOfTriangle, "maxAngle");
	}

	bool ComputeTriangleJacobianConditionNumberVertexValues(pmp::SurfaceMesh& mesh)
	{
		return ComputeInterpolatedTrianglesMetric(mesh, GetConditionNumberOfTriangleJacobian, "jacobianConditionNumber");
	}

	bool ComputeStiffnessMatrixConditioningVertexValues(pmp::SurfaceMesh& mesh)
	{
		return ComputeInterpolatedTrianglesMetric(mesh, ComputeStiffnessMatrixConditioningForTriangle, "stiffnessMatrixConditioning");
	}

	std::set<std::string> REGISTERED_METRIC_FUNCTION_NAMES{
		"minAngle", "maxAngle", "jacobianConditionNumber", "stiffnessMatrixConditioning"
	};

	bool IsMetricRegistered(const std::string& name)
	{
		if (REGISTERED_METRIC_FUNCTION_NAMES.contains(name))
			return true;

		return false;
	}

	TriMetricFunction IdentifyMetricFunction(const std::string& metricName)
	{
		if (metricName == "minAngle")
			return ComputeTriangleMinAngleVertexValues;

		if (metricName == "maxAngle")
			return ComputeTriangleMaxAngleVertexValues;

		if (metricName == "jacobianConditionNumber")
			return ComputeTriangleJacobianConditionNumberVertexValues;

		if (metricName == "stiffnessMatrixConditioning")
			return ComputeStiffnessMatrixConditioningVertexValues;

		return {};
	}

} // namespace Geometry