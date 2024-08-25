#include "MeshAnalysis.h"

#include "CollisionKdTree.h"
#include "GeometryUtil.h"
#include "GridUtil.h"
#include "MeshSelfIntersection.h"

#include "pmp/algorithms/Curvature.h"
#include "pmp/algorithms/DifferentialGeometry.h"
#include "pmp/algorithms/Normals.h"
#include "pmp/algorithms/Features.h"
#include "pmp/algorithms/TriangleKdTree.h"

#include <set>
#include <unordered_set>
#include <ranges>
#include <iomanip> // for std::setprecision

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

		// Jacobian is weighed by uniform scaling of the triangle. We normalize the edge vectors to assume the triangle is a unit triangle.
		const float l0 = pmp::norm(e0);
		const float l1 = pmp::norm(e1);
		const float lNorm = pmp::norm(e2);
		if (l0 < FLT_EPSILON || l1 < FLT_EPSILON || lNorm < FLT_EPSILON)
			return FLT_MAX; // singular Jacobian has an infinite condition number

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
				const unsigned int elementId = 3 * j + i;
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

	/// \brief Computes the condition number of a equilateral Jacobian of this triangle.
	///        Note: The Jacobian corresponds to the planar transformation from triangle (-0.5, 0), (0.5, 0), (0, 1) to the evaluated triangle.
	[[nodiscard]] float GetConditionNumberOfEquilateralTriangleJacobian(const pmp::SurfaceMesh& mesh, const pmp::Face& face)
	{
		const auto nVerts = std::distance(mesh.vertices(face).begin(), mesh.vertices(face).end());
		if (nVerts != 3)
			return -1.0f; // Error, face with three vertices is needed!

		const auto h = mesh.halfedge(face);
		const auto hPrev = mesh.prev_halfedge(h);

		const auto v = mesh.from_vertex(h);
		const auto v0 = mesh.from_vertex(hPrev);
		const auto v1 = mesh.to_vertex(h);
		const auto e0Midpoint = 0.5f * (mesh.position(v) + mesh.position(v0));

		const auto e0 = mesh.position(v0) - e0Midpoint;
		const auto e1 = mesh.position(v1) - e0Midpoint;
		const auto e2 = pmp::cross(e0, e1); // tri plane normal

		// Jacobian is weighed by uniform scaling of the triangle. We normalize the edge vectors to assume the triangle is a unit triangle.
		const float l0 = pmp::norm(e0);
		const float l1 = pmp::norm(e1);
		const float lNorm = pmp::norm(e2);
		if (l0 < FLT_EPSILON || l1 < FLT_EPSILON || lNorm < FLT_EPSILON)
			return FLT_MAX; // singular Jacobian has an infinite condition number

		const auto xVector = pmp::normalize(e0);
		const auto zVector = pmp::normalize(e2);
		const auto yVector = pmp::perp(xVector, zVector);

		const float e0X = 2.0f * pmp::norm(e0);
		constexpr float e0Y = 0.0f;

		const float e1X = pmp::dot(e1, xVector);
		const float e1Y = pmp::dot(e1, yVector);

		const float detJ = (e0X * e1Y - e1X * e0Y);
		if (std::fabs(detJ) < FLT_EPSILON)
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
				const unsigned int elementId = 3 * j + i;
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
		const auto valReal = (l0Sq + l1Sq + l2Sq) / (8.0f * area);
		const auto valImag = sqrt(48.0f * area - lSqSumSq) / (8.0f * area);
		return sqrt(valReal * valReal + valImag * valImag);
	}

	using FaceMetricFunction = std::function<float(const pmp::SurfaceMesh&, const pmp::Face&)>;
	//using VertexMetricFunction = std::function<float(const pmp::SurfaceMesh&, const pmp::Vertex&)>;

	constexpr float METRIC_MAX_VAL = 1e+12;

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

				if (std::fabsf(val) > METRIC_MAX_VAL)
					continue; // skipping triangle
				mean += val;
			}
			if (mean < FLT_EPSILON)
				return false;
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

	bool ComputeEquilateralTriangleJacobianConditionNumbers(pmp::SurfaceMesh& mesh)
	{
		return ComputeInterpolatedTrianglesMetric(mesh, GetConditionNumberOfEquilateralTriangleJacobian, "equilateralJacobianCondition");
	}

	bool ComputeStiffnessMatrixConditioningVertexValues(pmp::SurfaceMesh& mesh)
	{
		return ComputeInterpolatedTrianglesMetric(mesh, ComputeStiffnessMatrixConditioningForTriangle, "stiffnessMatrixConditioning");
	}

	std::set<std::string> REGISTERED_METRIC_FUNCTION_NAMES{
		"minAngle", "maxAngle", "jacobianConditionNumber", "equilateralJacobianCondition", "stiffnessMatrixConditioning"
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

		if (metricName == "equilateralJacobianCondition")
			return ComputeEquilateralTriangleJacobianConditionNumbers;

		if (metricName == "stiffnessMatrixConditioning")
			return ComputeStiffnessMatrixConditioningVertexValues;

		return {};
	}

	void ComputeEdgeDihedralAngles(pmp::SurfaceMesh& mesh)
	{
		auto eProp = mesh.edge_property<pmp::Scalar>("e:dihedralAngle", 2.0f * M_PI);

		for (const auto e : mesh.edges())
		{
			if (mesh.is_boundary(e))
				continue;

			const auto f0 = mesh.face(mesh.halfedge(e, 0));
			const auto f1 = mesh.face(mesh.halfedge(e, 1));

			const pmp::Normal n0 = pmp::Normals::compute_face_normal(mesh, f0);
			const pmp::Normal n1 = pmp::Normals::compute_face_normal(mesh, f1);

			const auto angleBetweenNormals = angle(n0, n1);
			eProp[e] = angleBetweenNormals + M_PI_2;
		}

		auto vProp = mesh.vertex_property<pmp::Scalar>("v:dihedralAngle", 2.0f * M_PI);
		for (const auto v : mesh.vertices())
		{
			if (mesh.is_boundary(v))
				continue;

			pmp::Scalar dihedralAngleSum = 0.0f;
			for (const auto w : mesh.vertices(v))
			{
				const auto hFrom = mesh.find_halfedge(v, w);
				const auto e = mesh.edge(hFrom);
				dihedralAngleSum += eProp[e];
			}
			const auto valence = mesh.valence(v);
			if (valence == 0)
				continue;

			vProp[v] = dihedralAngleSum / static_cast<pmp::Scalar>(valence);
		}
	}

	void ComputeVertexCurvaturesAndRelatedProperties(pmp::SurfaceMesh& mesh, const float& principalCurvatureFactor)
	{
		pmp::Curvature curvAlg{ mesh };
		curvAlg.analyze_tensor(1);

		auto vMinCurvature = mesh.vertex_property<pmp::Scalar>("v:minCurvature");
		auto vMaxCurvature = mesh.vertex_property<pmp::Scalar>("v:maxCurvature");
		auto vMeanCurvature = mesh.vertex_property<pmp::Scalar>("v:meanCurvature");
		auto vGaussianCurvature = mesh.vertex_property<pmp::Scalar>("v:GaussianCurvature");
		auto vIsCDSVal = mesh.vertex_property<pmp::Scalar>("v:isCDS", -1.0f);

		for (const auto v : mesh.vertices())
		{
			vMinCurvature[v] = curvAlg.min_curvature(v);
			vMaxCurvature[v] = curvAlg.max_curvature(v);
			vMeanCurvature[v] = vMinCurvature[v] + vMaxCurvature[v];
			vGaussianCurvature[v] = vMinCurvature[v] * vMaxCurvature[v];
			vIsCDSVal[v] = pmp::IsConvexDominantSaddle(vMinCurvature[v], vMaxCurvature[v], principalCurvatureFactor) ? 1.0f : -1.0f;
		}
	}

	void ComputeZLevelElevations(pmp::SurfaceMesh& mesh)
	{
		auto vProp = mesh.vertex_property<pmp::Scalar>("v:zLevelElevation", 0.0f);
		for (const auto v : mesh.vertices())
		{
			const auto& vPos = mesh.position(v);
			vProp[v] = vPos[2];
		}
	}

	std::tuple<pmp::Scalar, pmp::Scalar, pmp::Scalar> ComputeEdgeLengthMinAverageAndMax(const pmp::SurfaceMesh& mesh)
	{
		pmp::Scalar lengthSqMin = FLT_MAX;
		pmp::Scalar lengthSqMax = -FLT_MAX;
		pmp::Scalar lengthSqMean= 0.0f;
		for (const auto e : mesh.edges())
		{
			const pmp::Scalar edgeLengthSq = mesh.edge_length_sq(e);
			if (edgeLengthSq > lengthSqMax)
				lengthSqMax = edgeLengthSq;
			if (edgeLengthSq < lengthSqMin)
				lengthSqMin = edgeLengthSq;
			lengthSqMean += sqrt(edgeLengthSq);
		}
		lengthSqMean /= static_cast<pmp::Scalar>(mesh.n_edges());
		return { sqrt(lengthSqMin), lengthSqMean, sqrt(lengthSqMax) };
	}

	size_t CountUnreferencedVertices(const BaseMeshGeometryData& data)
	{
		if (data.Vertices.empty())
			return 0;
		std::vector visited(data.Vertices.size(), false);
		for (const auto& polyIds : data.PolyIndices)
		{
			for (const auto& vId : polyIds)
			{
				visited[vId] = true;
			}
		}
		const size_t nVertsTotal = data.Vertices.size();
		return nVertsTotal - std::ranges::count_if(visited, [](const auto& vis) { return vis; });
	}

	//
	// ==========================================================================
	// ...................... Saliency utils ...................................
	//

	static void ComputeVertexCurvatures(pmp::SurfaceMesh& mesh)
	{
		pmp::Curvature curvAlg{ mesh };
		curvAlg.analyze_tensor(1); // Use tensor-based analysis to get mean curvature

		auto vMinCurvature = mesh.vertex_property<pmp::Scalar>("v:minCurvature");
		auto vMaxCurvature = mesh.vertex_property<pmp::Scalar>("v:maxCurvature");
		auto vMeanCurvature = mesh.vertex_property<pmp::Scalar>("v:meanCurvature");

		for (const auto v : mesh.vertices())
		{
			vMinCurvature[v] = curvAlg.min_curvature(v);
			vMaxCurvature[v] = curvAlg.max_curvature(v);
			vMeanCurvature[v] = (vMinCurvature[v] + vMaxCurvature[v]) * 0.5f;
		}
	}

	static [[nodiscard]] double GaussianWeight(const pmp::SurfaceMesh& mesh, const pmp::Vertex& v, const double& sigma)
	{
		if (!mesh.has_vertex_property("v:meanCurvature"))
		{
			throw std::invalid_argument("Geometry::GaussianWeight: input mesh has no vertex property called \"v:meanCurvature\"\n");
		}

		auto vCurvature = mesh.get_vertex_property<pmp::Scalar>("v:meanCurvature");
		auto positions = mesh.get_vertex_property<pmp::Point>("v:point");
		double weightSum = 0.0;
		double curvatureSum = 0.0;
		const pmp::Point p = positions[v];

		for (const auto u : mesh.vertices(v))
		{
			pmp::Point q = positions[u];
			const pmp::Scalar distance = pmp::norm(p - q);
			const pmp::Scalar weight = exp(-pow(distance, 2) / (2 * pow(sigma, 2)));
			curvatureSum += vCurvature[u] * weight;
			weightSum += weight;
		}

		return curvatureSum / weightSum;
	}

	static void ComputeSaliency(pmp::SurfaceMesh& mesh, const std::vector<double>& sigmas)
	{
		auto saliency = mesh.vertex_property<pmp::Scalar>("v:saliency", 0.0f);

		for (auto v : mesh.vertices()) 
		{
			pmp::Scalar saliencyValue = 0.0f;
			for (double sigma : sigmas) 
			{
				const double fine = GaussianWeight(mesh, v, sigma);
				const double coarse = GaussianWeight(mesh, v, 2 * sigma);
				saliencyValue += std::abs(fine - coarse);
			}
			saliency[v] = saliencyValue;
		}
	}

	static void NormalizeSaliency(pmp::SurfaceMesh& mesh)
	{
		if (!mesh.has_vertex_property("v:saliency"))
		{
			throw std::invalid_argument("Geometry::NormalizeSaliency: input mesh has no vertex property called \"v:saliency\"\n");
		}

		auto saliency = mesh.get_vertex_property<pmp::Scalar>("v:saliency");
		double maxSaliency = *std::ranges::max_element(saliency.vector());

		for (auto v : mesh.vertices())
		{
			saliency[v] = pow(saliency[v] / maxSaliency, 2); // Square the ratio to enhance high values
		}
	}

	bool EvaluatePMPSurfaceMeshSaliency(pmp::SurfaceMesh& mesh, const double& forcedVariance, const bool& normalizeValues)
	{
		if (mesh.is_empty())
		{
			std::cerr << "Geometry::EvaluatePMPSurfaceMeshSaliency: invalid input mesh!\n";
			return false;
		}

		ComputeVertexCurvatures(mesh);

		// Dynamic evaluation of sigma based on mesh properties

		double averageEdgeLength;
		if (forcedVariance > 0.0)
		{
			averageEdgeLength = forcedVariance;
		}
		else
		{
			averageEdgeLength = 0.0;
			for (const auto e : mesh.edges()) 
			{
				averageEdgeLength += static_cast<double>(mesh.edge_length(e));
			}
			averageEdgeLength /= static_cast<double>(mesh.n_edges());
		}
		const std::vector sigmas = {
			averageEdgeLength * 0.5, averageEdgeLength, averageEdgeLength * 1.5,
			averageEdgeLength * 2, averageEdgeLength * 2.5
		};

		ComputeSaliency(mesh, sigmas);
		if (normalizeValues)
		{
			NormalizeSaliency(mesh);
		}
		return true;
	}

	void PrintHistogramResultData(const std::pair<std::pair<float, float>, std::vector<unsigned int>>& histData, std::ostream& os)
	{
		if (std::abs(histData.first.first - histData.first.second) < FLT_EPSILON || histData.second.empty())
		{
			std::cerr << "PrintHistogramResultData: INVALID HISTOGRAM DATA!\n";
			return;
		}
		const auto& [range, bins] = histData;
		const auto& [minDistVal, maxDistVal] = range;
		const float binSize = (maxDistVal - minDistVal) / static_cast<float>(bins.size());

		// Print the first bin with -inf lower bound
		os << "( -inf.. " << std::fixed << std::setprecision(7) << (minDistVal + binSize) << ") : " << bins.front() << '\n';

		// Print the intermediate bins
		for (size_t i = 1; i < bins.size() - 1; ++i) 
		{
			os << "[ "
				<< std::fixed << std::setprecision(7) << (minDistVal + binSize * i) << ".. "
				<< std::fixed << std::setprecision(7) << (minDistVal + binSize * (i + 1)) << ") : "
				<< bins[i] << '\n';
		}

		// Print the last bin with +inf upper bound
		if (bins.size() > 1) 
		{
			os << "[ "
				<< std::fixed << std::setprecision(7) << (maxDistVal - binSize) << ".. +inf) : "
				<< bins.back() << '\n';
		}
	}

	std::pair<std::vector<size_t>, std::vector<size_t>> GetEdgeVertCountsTheoreticalEstimate(const pmp::SurfaceMesh& mesh, const size_t& maxSubdivLevel, const bool& evalOutput)
	{
		const auto nBdEdges0 = mesh.n_boundary_edges();
		const size_t sMax = maxSubdivLevel;

		if (nBdEdges0 == 0)
		{
			// mesh is watertight
			const auto nEdges0 = mesh.n_edges();
			const auto nVerts0 = mesh.n_vertices();

			if (evalOutput)
			{
				std::cout << "............................................................\n";
				std::cout << "GetEdgeVertCountsTheoreticalEstimate:\n";
				std::cout << "nEdges0 = " << nEdges0 << "\n";
				std::cout << "nVerts0 = " << nVerts0 << "\n";
				std::cout << "............................................................\n";
			}

			const auto edgeCountEstimate = [&nEdges0](const size_t& s) { return (static_cast<size_t>(pow(4, s)) * nEdges0); };
			const auto vertCountEstimate = [&nEdges0, &nVerts0](const size_t& s) { return ((nEdges0 * static_cast<size_t>(pow(4, s) - 1) + 3 * nVerts0) / 3); };

			std::vector<size_t> edgeCounts(sMax);
			std::vector<size_t> vertCounts(sMax);

			for (size_t s = 0; s < sMax; s++)
			{
				edgeCounts[s] = edgeCountEstimate(s);
				vertCounts[s] = vertCountEstimate(s);
			}

			return { edgeCounts, vertCounts };
		}

		// mesh is not watertight
		const auto nEdges0 = mesh.n_edges();
		const auto nVerts0 = mesh.n_vertices();

		const auto nIntEdges0 = nEdges0 - nBdEdges0;

		if (evalOutput)
		{
			std::cout << "............................................................\n";
			std::cout << "GetEdgeVertCountsTheoreticalEstimate:\n";
			std::cout << "nIntEdges0 = " << nIntEdges0 << ", nBdEdges = " << nBdEdges0 << "\n";
			std::cout << "nVerts0 = " << nVerts0 << "\n";
			std::cout << "............................................................\n";
		}

		const auto intEdgeCountEstimate = [&nIntEdges0, &nBdEdges0](const size_t& s) -> size_t
		{
			return (pow(2, s - 1)) * ((pow(2, s) - 1) * nBdEdges0 + nIntEdges0 * pow(2, s + 1));
		};
		const auto bdEdgeCountEstimate = [&nBdEdges0](const size_t& s) -> size_t
		{
			return (pow(2, s)) * nBdEdges0;
		};
		const auto vertCountEstimate = [&nIntEdges0, &nBdEdges0, &nVerts0](const size_t& s) -> size_t
		{
			return nBdEdges0 * (pow(4, s) - 4 + 3 * pow(2, s)) / 6 + nIntEdges0 * (pow(4, s) - 1) / 3 + nVerts0;
		};

		std::vector<size_t> edgeCounts(sMax);
		std::vector<size_t> vertCounts(sMax);

		for (size_t s = 0; s < sMax; s++)
		{
			edgeCounts[s] = intEdgeCountEstimate(s) + bdEdgeCountEstimate(s);
			vertCounts[s] = vertCountEstimate(s);
		}

		return { edgeCounts, vertCounts };
	}

	size_t CountPMPSurfaceMeshSelfIntersectingFaces(pmp::SurfaceMesh& mesh, const bool& setFaceProperty)
	{
		if (!mesh.is_triangle_mesh())
		{
			throw std::invalid_argument("CountPMPSurfaceMeshSelfIntersectingFaces: non-triangle SurfaceMesh not supported for this function!\n");
		}

		pmp::FaceProperty<bool> fIsSelfIntersecting;
		if (setFaceProperty)
		{
			fIsSelfIntersecting = mesh.add_face_property<bool>("f:isSelfIntersecting", false);
		}

		const PMPSurfaceMeshAdapter meshAdapter(std::make_shared<pmp::SurfaceMesh>(mesh));
		const auto ptrMeshCollisionKdTree = std::make_unique<CollisionKdTree>(meshAdapter, CenterSplitFunction);
		size_t nSelfIntFaceCountResult = 0;
		for (const auto f : mesh.faces())
		{
			pmp::BoundingBox fBBox;
			std::vector<pmp::Point> vertices0;
			vertices0.reserve(3);
			std::unordered_set<unsigned int> neighboringFaceIds;
			for (const auto v : mesh.vertices(f))
			{
				for (const auto nf : mesh.faces(v))
					neighboringFaceIds.insert(nf.idx());
				const auto& vPos = mesh.position(v);
				vertices0.push_back(vPos);
				fBBox += vPos;
			}

			// Query the kd-tree for candidates
			std::vector<unsigned int> candidateIds;
			ptrMeshCollisionKdTree->GetTrianglesInABox(fBBox, candidateIds);

			for (const auto ci : candidateIds)
			{
				const auto cf = pmp::Face(ci);
				if (cf == f || neighboringFaceIds.contains(ci))
				{
					continue; // Skip self and neighboring faces
				}

				std::vector<pmp::Point> vertices1;
				vertices1.reserve(3);
				for (const auto cv : mesh.vertices(cf))
				{
					vertices1.push_back(mesh.position(cv));
				}

				if (TriangleIntersectsTriangle(vertices0, vertices1))
				{
					++nSelfIntFaceCountResult; // Found an intersection
					if (setFaceProperty) 
					{
						fIsSelfIntersecting[f] = true;
						fIsSelfIntersecting[cf] = true;
					}
					break; // Only count once per face
				}
			}
		}

		return nSelfIntFaceCountResult;
	}

	bool PMPSurfaceMeshHasSelfIntersections(const pmp::SurfaceMesh& mesh)
	{
		if (!mesh.is_triangle_mesh())
		{
			throw std::invalid_argument("PMPSurfaceMeshHasSelfIntersections: non-triangle SurfaceMesh not supported for this function!\n");
		}

		const PMPSurfaceMeshAdapter meshAdapter(std::make_shared<pmp::SurfaceMesh>(mesh));
		const auto ptrMeshCollisionKdTree = std::make_unique<CollisionKdTree>(meshAdapter, CenterSplitFunction);

		for (const auto f : mesh.faces()) 
		{
			pmp::BoundingBox fBBox;
			std::vector<pmp::Point> vertices0;
			vertices0.reserve(3);
			std::unordered_set<unsigned int> neighboringFaceIds;
			for (const auto v : mesh.vertices(f))
			{
				for (const auto nf : mesh.faces(v))
					neighboringFaceIds.insert(nf.idx());
				const auto& vPos = mesh.position(v);
				vertices0.push_back(vPos);
				fBBox += vPos;
			}

			// Query the kd-tree for candidates
			std::vector<unsigned int> candidateIds;
			ptrMeshCollisionKdTree->GetTrianglesInABox(fBBox, candidateIds);

			for (const auto ci : candidateIds)
			{
				const auto cf = pmp::Face(ci);
				if (cf == f || neighboringFaceIds.contains(ci))
				{
					continue; // Skip self and neighboring faces
				}

				std::vector<pmp::Point> vertices1;
				vertices1.reserve(3);
				for (const auto cv : mesh.vertices(cf))
				{
					vertices1.push_back(mesh.position(cv));
				}

				if (TriangleIntersectsTriangle(vertices0, vertices1)) 
				{
					return true; // Found an intersection, return immediately
				}
			}
		}

		return false; // No intersections found
	}

	void ConvertPMPSurfaceMeshBoolFacePropertyToScalarVertexProperty(pmp::SurfaceMesh& mesh, const std::string& propName)
	{
		if (!mesh.has_face_property(propName))
		{
			std::cerr << "ConvertPMPSurfaceMeshBoolFacePropertyToScalarVertexProperty: face property \"" << propName << "\" not found! Aborting...\n";
			return;
		}

		auto fProp = mesh.get_face_property<bool>(propName);
		auto vProp = mesh.vertex_property<pmp::Scalar>("v:" + propName, 0.0f);

		for (const auto v : mesh.vertices())
		{
			pmp::Scalar mean = 0.0f;
			for (const auto f : mesh.faces(v))
			{
				const auto val = (fProp[f] ? 1.0f : -1.0f);
				mean += val;
			}
			const auto nAdjacentFaces = std::distance(mesh.faces(v).begin(), mesh.faces(v).end());
			mean /= static_cast<float>(nAdjacentFaces);
			vProp[v] = mean;
		}
	}

	std::unordered_multimap<unsigned int, unsigned int> ExtractPMPSurfaceMeshFaceIntersectionMultimap(const pmp::SurfaceMesh& mesh)
	{
		if (!mesh.is_triangle_mesh())
		{
			throw std::invalid_argument("PMPSurfaceMeshHasSelfIntersections: non-triangle SurfaceMesh not supported for this function!\n");
		}

		std::unordered_multimap<unsigned int, unsigned int> faceToIntersectingFaces;
		const PMPSurfaceMeshAdapter meshAdapter(std::make_shared<pmp::SurfaceMesh>(mesh));
		const auto ptrMeshCollisionKdTree = std::make_unique<CollisionKdTree>(meshAdapter, CenterSplitFunction);

		for (const auto f : mesh.faces()) 
		{
			pmp::BoundingBox fBBox;
			std::vector<pmp::Point> vertices0;
			vertices0.reserve(3);
			std::unordered_set<unsigned int> neighboringFaceIds;

			for (const auto v : mesh.vertices(f)) 
			{
				for (const auto nf : mesh.faces(v)) 
				{
					neighboringFaceIds.insert(nf.idx());
				}
				const auto& vPos = mesh.position(v);
				vertices0.push_back(vPos);
				fBBox += vPos;
			}

			std::vector<unsigned int> candidateIds;
			ptrMeshCollisionKdTree->GetTrianglesInABox(fBBox, candidateIds);

			for (const auto& ci : candidateIds)
			{
				if (ci == f.idx() || neighboringFaceIds.contains(ci))
				{
					continue; // Skip self and neighboring faces
				}

				std::vector<pmp::Point> vertices1;
				vertices1.reserve(3);
				const auto cf = pmp::Face(ci);
				for (const auto cv : mesh.vertices(cf))
				{
					vertices1.push_back(mesh.position(cv));
				}

				if (TriangleIntersectsTriangle(vertices0, vertices1))
				{
					faceToIntersectingFaces.emplace(f.idx(), ci);
				}
			}
		}

		return faceToIntersectingFaces;
	}

	std::vector<std::vector<pmp::vec3>> ComputeSurfaceMeshSelfIntersectionPolylines(const pmp::SurfaceMesh& mesh)
	{
		MeshSelfIntersectionBucketCollector intersectionDataCollector(mesh);
		const auto faceDataBuckets = intersectionDataCollector.Retrieve();
		std::vector<std::vector<pmp::vec3>> resultPolylines;
		resultPolylines.reserve(faceDataBuckets.size());

		for (const auto& bucket : faceDataBuckets)
		{
			if (bucket.Empty())
				continue;

			std::vector<pmp::vec3> currentPolyline;
			currentPolyline.reserve(bucket.Size());
			const auto bucketOrientation = bucket.GetOrientation();
			pmp::Point barycenter;
			auto facePtData = ExtractFaceData(bucket, barycenter);
			const auto& refPt = facePtData[0].second;
			const auto refVec = refPt - barycenter;

			// sort face points in CCW orientation in the plane defined by barycenter and bucketOrientation.
			// DISCLAIMER: this sorting does NOT account for polylines with negative winding segments, this will require a subroutine to account for that.
			std::ranges::sort(facePtData, [&](const auto& a, const auto& b) {
				const auto vecA = a.second - barycenter;
				const auto vecB = b.second - barycenter;
				pmp::Scalar angleA = atan2(dot(cross(refVec, vecA), bucketOrientation), dot(refVec, vecA));
				pmp::Scalar angleB = atan2(dot(cross(refVec, vecB), bucketOrientation), dot(refVec, vecB));
				angleA = angleA < 0.0f ? (2.0f * static_cast<pmp::Scalar>(M_PI) + angleA) : angleA;
				angleB = angleB < 0.0f ? (2.0f * static_cast<pmp::Scalar>(M_PI) + angleB) : angleB;
				if (std::abs(angleA - angleB) < FLT_MIN)
					return norm(vecA) < norm(vecB);
				return angleA < angleB;
			});

			for (unsigned int i = 0; i < facePtData.size(); ++i)
			{
				const auto& p = facePtData[i].second;
				if (i > 0 && norm(p - facePtData[i - 1].second) < 1e-6f)
				{
					continue; // duplicate pts
				}
				currentPolyline.push_back(p);
			}
			currentPolyline.push_back(refPt); // close the polyline
			resultPolylines.push_back(currentPolyline);
		}

		return resultPolylines;
	}

	bool PMPManifoldCurve2DHasSelfIntersections(const pmp::ManifoldCurve2D& curve)
	{
		if (curve.is_empty())
		{
			throw std::invalid_argument("PMPManifoldCurve2DHasSelfIntersections: cannot process empty curves.\n");
		}

		const ManifoldCurve2DAdapter curveAdapter(std::make_shared<pmp::ManifoldCurve2D>(curve));
		const auto ptrCurveCollisionKdTree = std::make_unique<Collision2DTree>(curveAdapter, CenterSplitFunction2D);

		for (const auto e : curve.edges())
		{
			const auto [v0, v1] = curve.vertices(e);
			pmp::BoundingBox2 eBBox({curve.position(v0), curve.position(v1)});
			const auto ePrev = curve.edge_to(v0);
			const auto eNext = curve.edge_from(v1);

			// Query the kd-tree for candidates
			std::vector<unsigned int> candidateIds;
			ptrCurveCollisionKdTree->GetEdgesInABox(eBBox, candidateIds);

			for (const auto ci : candidateIds)
			{
				const auto ce = pmp::Edge(ci);
				if (ce == e || ce == ePrev || ce == eNext)
				{
					continue; // Skip self and neighboring edges
				}

				const auto [v0Candidate, v1Candidate] = curve.vertices(ce);
				if (Line2DIntersectsLine2D({ curve.position(v0), curve.position(v1) }, { curve.position(v0Candidate), curve.position(v1Candidate) }))
				{
					return true; // Found an intersection, return immediately
				}
			}
		}

		return false; // No intersections found
	}

	bool IsPointInsidePMPManifoldCurve(const pmp::Point2& point, const pmp::ManifoldCurve2D& curve)
	{
		const auto& positions = curve.positions();
		const size_t n = positions.size();
		bool inside = false;

		// Ray-casting algorithm
		for (size_t i = 0, j = n - 1; i < n; j = i++)
		{
			const auto& pi = positions[i];
			const auto& pj = positions[j];

			if (((pi[1] > point[1]) != (pj[1] > point[1])) &&
				(point[0] < (pj[0] - pi[0]) * (point[1] - pi[1]) / (pj[1] - pi[1]) + pi[0]))
			{
				inside = !inside;
			}
		}

		return inside;
	}

	bool IsPointInsidePMPSurfaceMesh(const pmp::Point& point, const pmp::SurfaceMesh& mesh)
	{
		if (mesh.is_empty())
			return false;

		if (!mesh.is_triangle_mesh())
		{
			throw std::invalid_argument("IsPointInsidePMPSurfaceMesh: non-triangle SurfaceMesh not supported for this function!\n");
		}

		const PMPSurfaceMeshAdapter meshAdapter(std::make_shared<pmp::SurfaceMesh>(mesh));
		const auto ptrMeshCollisionKdTree = std::make_unique<CollisionKdTree>(meshAdapter, CenterSplitFunction);

		// Cast a ray in an arbitrary direction (e.g., +X direction)
		Ray ray(point, pmp::vec3(1.0f, 0.0f, 0.0f));
		return ptrMeshCollisionKdTree->IsRayStartPointInsideTriangleMesh(ray);
	}

	bool IsPointInsidePMPSurfaceMesh(const pmp::Point& point, const std::shared_ptr<CollisionKdTree>& meshKdTree)
	{
		if (!meshKdTree)
			return false;

		// Cast a ray in an arbitrary direction (e.g., +X direction)
		Ray ray(point, pmp::vec3(1.0f, 0.0f, 0.0f));
		return meshKdTree->IsRayStartPointInsideTriangleMesh(ray);
	}

	namespace
	{
		struct TangentQuadricParams
		{
			pmp::Scalar a{ 1.0f }, b{ 1.0f }, c{ 1.0f };
			pmp::Point center{};
		};

		[[nodiscard]] pmp::Scalar DistanceToTangentEllipsoid(const pmp::Point& point, const TangentQuadricParams& params)
		{
			pmp::vec3 P(point[0] - params.center[0], point[1] - params.center[1], point[2] - params.center[2]);

			const float k0 = std::sqrt((P[0] / params.a) * (P[0] / params.a) + (P[1] / params.b) * (P[1] / params.b) + (P[2] / params.c) * (P[2] / params.c));
			const float k1 = std::sqrt((P[0] / (params.a * params.a)) * (P[0] / (params.a * params.a)) +
				(P[1] / (params.b * params.b)) * (P[1] / (params.b * params.b)) +
				(P[2] / (params.c * params.c)) * (P[2] / (params.c * params.c)));
			return std::max(k0 * (k0 - 1.f) / k1, 0.f);
		}

		[[nodiscard]] pmp::Scalar DistanceToTangentHyperboloid(const pmp::Point& point, const TangentQuadricParams& params)
		{
			pmp::vec3 P(point[0] - params.center[0], point[1] - params.center[1], point[2] - params.center[2]);
			// Calculate the one-sheet hyperboloid function value at the point
			const float hyperboloidValue = (P[0] / params.a) * (P[0] / params.a) +
				(P[1] / params.b) * (P[1] / params.b) -
				(P[2] / params.c) * (P[2] / params.c) - 1.0f;
			return std::sqrt(std::abs(hyperboloidValue)); // need to make the distance positive
		}

		struct TangentCircleParams
		{
			pmp::Scalar r{ 1.0f };
			pmp::Point2 center{};
		};

		[[nodiscard]] pmp::Scalar DistanceToTangentCircle(const pmp::Point2& point, const TangentCircleParams& cParams)
		{
			const pmp::Scalar distanceToCenter = std::sqrt((point[0] - cParams.center[0]) * (point[0] - cParams.center[0]) +
				(point[1] - cParams.center[1]) * (point[1] - cParams.center[1]));
			return std::abs(distanceToCenter - cParams.r);
		}

	} // anonymous namespace

	constexpr pmp::Scalar CURVATURE_EPSILON{ 1e-6f };

	pmp::Scalar CalculateQuadricApproximationErrorAtVertex(const pmp::SurfaceMesh& mesh, pmp::Vertex v)
	{
		const auto curvature = vertex_curvature(mesh, v);

		if (std::fabs(curvature.gauss) < CURVATURE_EPSILON) // curvature.gauss = curvature.min * curvature.max			
		{
			// no quadric of reasonable size can be fitted to the mesh vertex surroundings because the surface is flat.
			return 0.0f;
		}

		TangentQuadricParams qParams;

		if (curvature.gauss > 0)
		{
			// Estimate tangent ellipsoid parameters (a, b, c) from the curvature
			qParams.a = 1.0f / std::sqrt(std::max(curvature.max, CURVATURE_EPSILON));
			qParams.b = 1.0f / std::sqrt(std::max(curvature.mean, CURVATURE_EPSILON));
			qParams.c = 1.0f / std::sqrt(std::max(curvature.min, CURVATURE_EPSILON));
		}
		else
		{
			// Estimate tangent single-sheet hyperboloid parameters (a, b, c) from the curvature
			qParams.a = 1.0f / std::sqrt(std::max(std::fabs(curvature.max), CURVATURE_EPSILON));
			qParams.b = 1.0f / std::sqrt(std::max(std::fabs(curvature.mean), CURVATURE_EPSILON));
			qParams.c = 1.0f / std::sqrt(std::max(std::fabs(curvature.min), CURVATURE_EPSILON));
		}

		const auto distanceFromQuadric	= (curvature.gauss > 0) ? DistanceToTangentEllipsoid : DistanceToTangentHyperboloid;

		// Use the vertex position as the ellipsoid center approximation
		const auto normal = pmp::Normals::compute_vertex_normal(mesh, v);
		const auto standardNormal = pmp::vec3{ 0, 0, 1 };
		const auto reorientToStandardNormal = rotation_matrix(normal, standardNormal);
		const auto translateToCurrentPosition = translation_matrix(-mesh.position(v));
		const auto affineTransformToCurrentPos = translateToCurrentPosition * reorientToStandardNormal;
		qParams.center = standardNormal * curvature.mean;

		pmp::Scalar maxDeviation = 0.0;
		for (const auto f : mesh.faces(v))
		{
			pmp::Point faceCentroid = affine_transform(affineTransformToCurrentPos, centroid(mesh, f));
			const auto deviation = distanceFromQuadric(faceCentroid, qParams);
			maxDeviation = std::max(maxDeviation, deviation);
		}

		return maxDeviation;
	}

	pmp::Scalar CalculateCircularApproximationErrorAtVertex(const pmp::ManifoldCurve2D& curve, pmp::Vertex v)
	{
		const auto curvature = vertex_curvature(curve, v);

		if (curvature < CURVATURE_EPSILON)
		{
			// no circle of reasonable size can be fitted to the curve vertex surroundings because the curve is flat.
			return 0.0;
		}

		TangentCircleParams cParams;
		cParams.r = 1.0f / std::sqrt(std::max(curvature, CURVATURE_EPSILON));

		const auto normal = pmp::Normals2::compute_vertex_normal(curve, v);

		return pmp::Scalar();
	}

} // namespace Geometry