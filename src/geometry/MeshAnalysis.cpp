#include "MeshAnalysis.h"

#include <set>
#include <stack>
#include <unordered_set>

#include "CollisionKdTree.h"
#include "GeometryUtil.h"

#include "pmp/algorithms/Curvature.h"
#include "pmp/algorithms/DifferentialGeometry.h"
#include "pmp/algorithms/Normals.h"
#include "pmp/algorithms/Features.h"
#include "pmp/algorithms/TriangleKdTree.h"

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

	void ComputeVertexCurvatures(pmp::SurfaceMesh& mesh, const float& principalCurvatureFactor)
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

		const auto ptrMeshCollisionKdTree = std::make_unique<CollisionKdTree>(mesh, CenterSplitFunction);
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

		const auto ptrMeshCollisionKdTree = std::make_unique<CollisionKdTree>(mesh, CenterSplitFunction);

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
		const auto ptrMeshCollisionKdTree = std::make_unique<CollisionKdTree>(mesh, CenterSplitFunction);

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

	constexpr float POLYLINE_END_DISTANCE_TOLERANCE = 1e-6f;

	// TODO: find a faster way to do this by integrating the intersection computation after successful face querying
	std::vector<std::vector<pmp::vec3>> ComputeSurfaceMeshSelfIntersectionPolylines(const pmp::SurfaceMesh& mesh)
	{
		auto faceIntersections = ExtractPMPSurfaceMeshFaceIntersectionMultimap(mesh);
		std::vector<std::vector<pmp::vec3>> resultPolylines;

		auto fIt = faceIntersections.begin();
		std::stack<pmp::Halfedge> heStack;
		while (!faceIntersections.empty())
		{
			//if (faceIntersections.size() < 90)
			//	std::cout << "Here\n";

			if (fIt == faceIntersections.end())
				fIt = faceIntersections.begin();
			std::vector<pmp::vec3> currentPolyline;

			auto intersectingFaceRange = faceIntersections.equal_range(fIt->first);
			const auto baseFace = pmp::Face(fIt->first);
			std::vector<pmp::vec3> vertices0;
			vertices0.reserve(3);
			for (const auto v : mesh.vertices(baseFace))
			{
				vertices0.push_back(mesh.position(v));
			}

			// process all faces intersecting baseFace
			for (auto rangeIt = intersectingFaceRange.first; rangeIt != intersectingFaceRange.second; ++rangeIt)
			{
				const auto intersectingFace = pmp::Face(rangeIt->second);
				std::vector<pmp::vec3> vertices1;
				vertices1.reserve(3);
				for (const auto v : mesh.vertices(intersectingFace))
				{
					vertices1.push_back(mesh.position(v));
				}

				auto intersectionLineOpt = ComputeTriangleTriangleIntersectionLine(vertices0, vertices1);
				if (!intersectionLineOpt.has_value())
					continue;

				const auto& intersectionLine = intersectionLineOpt.value();
				if (norm(intersectionLine.second - intersectionLine.first) < POLYLINE_END_DISTANCE_TOLERANCE)
					continue; // degenerate intersection line

				if (currentPolyline.empty())
					currentPolyline.push_back(intersectionLine.first); // first point in the polyline

				currentPolyline.push_back(intersectionLine.second);
			}

			// fill the half-edge stack with current baseFace's half-edges
			for (const auto he : mesh.halfedges(baseFace))
				heStack.push(he);

			// search for neighbors logged as self-intersecting faces
			bool neighborFound = false;
			do
			{
				const auto he = heStack.top();
				heStack.pop();
				const auto oppHe = mesh.opposite_halfedge(he);
				if (!oppHe.is_valid())
					continue;

				const auto fNeighbor = mesh.face(oppHe);
				const auto fItNew = faceIntersections.find(fNeighbor.idx());
				if (fItNew == faceIntersections.cend())
					continue;

				neighborFound = true;
				faceIntersections.erase(fIt->first); // erasing all entries for this face
				fIt = fItNew;
				break;
			}
			while (!heStack.empty());

			if (!neighborFound)
			{
				// this is the last segment of currentPolyline
				resultPolylines.push_back(currentPolyline);
				currentPolyline = std::vector<pmp::vec3>();
				faceIntersections.erase(fIt->first); // erasing all entries for this face
				fIt = intersectingFaceRange.second;
				heStack = {};
			}
		}

		return resultPolylines;
	}

} // namespace Geometry