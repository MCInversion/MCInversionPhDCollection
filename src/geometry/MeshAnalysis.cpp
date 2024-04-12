#include "MeshAnalysis.h"

#include "CollisionKdTree.h"
#include "GeometryUtil.h"
#include "GridUtil.h"
#include "MeshSelfIntersection.h"

#include "sdf/SDF.h"

#include "pmp/algorithms/Curvature.h"
#include "pmp/algorithms/DifferentialGeometry.h"
#include "pmp/algorithms/Normals.h"
#include "pmp/algorithms/Features.h"
#include "pmp/algorithms/TriangleKdTree.h"

#include <set>
#include <unordered_set>
#include <ranges>
#include <limits>
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
		lengthSqMean /= mesh.n_edges();
		return { sqrt(lengthSqMin), lengthSqMean, sqrt(lengthSqMax) };
	}

	constexpr unsigned int DEFAULT_N_VOXELS_PER_MIN_DIMENSION = 40;

	std::optional<std::pair<std::pair<float, float>, std::vector<unsigned int>>> ComputeMeshDistanceToPointCloudPerVertexHistogram(const pmp::SurfaceMesh& mesh, const std::vector<pmp::vec3>& ptCloud, const unsigned int& nBins)
	{
		// Check for empty mesh or point cloud
		if (mesh.n_vertices() == 0 || ptCloud.empty())
		{
			return {};
		}

		// Compute distance field for the point cloud
		const pmp::BoundingBox ptCloudBBox(ptCloud);
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float cellSize = minSize / DEFAULT_N_VOXELS_PER_MIN_DIMENSION;
		constexpr float volExpansionFactor = 1.0f;
		const SDF::PointCloudDistanceFieldSettings dfSettings{
			cellSize,
			volExpansionFactor,
			Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::BlurPostprocessingType::None
		};

		auto df = SDF::PointCloudDistanceFieldGenerator::Generate(ptCloud, dfSettings);

		// Prepare for histogram computation
		std::vector<unsigned int> histogram(nBins, 0);
		float maxDistVal = std::numeric_limits<float>::lowest();
		float minDistVal = std::numeric_limits<float>::max();

		// Iterate through vertices to compute distances and populate the histogram
		for (const auto& v : mesh.vertices())
		{
			const auto vPos = mesh.position(v);
			const double vDistanceToTarget = TrilinearInterpolateScalarValue(vPos, df);
			maxDistVal = std::max(maxDistVal, static_cast<float>(vDistanceToTarget));
			minDistVal = std::min(minDistVal, static_cast<float>(vDistanceToTarget));
		}

		// Calculate histogram bin size and populate histogram
		const float binSize = (maxDistVal - minDistVal) / nBins;
		for (const auto& v : mesh.vertices())
		{
			const auto vPos = mesh.position(v);
			const double vDistanceToTarget = TrilinearInterpolateScalarValue(vPos, df);
			int binIndex = static_cast<int>((vDistanceToTarget - minDistVal) / binSize);
			binIndex = std::min(binIndex, static_cast<int>(nBins) - 1); // Ensure binIndex is within range
			histogram[binIndex]++;
		}

		return { {{minDistVal, maxDistVal}, histogram} };
	}

	std::optional<double> ComputeMeshToPointCloudHausdorffDistance(const pmp::SurfaceMesh& mesh, const std::vector<pmp::Point>& ptCloud, const unsigned int& nVoxelsPerMinDimension)
	{
		if (mesh.n_vertices() == 0 || ptCloud.empty())
		{
			return {};
		}

		// Compute distance field for the point cloud
		const pmp::BoundingBox ptCloudBBox(ptCloud);
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const float ptCloudMinSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float ptCloudCellSize = ptCloudMinSize / static_cast<float>(nVoxelsPerMinDimension);
		const SDF::PointCloudDistanceFieldSettings ptCloudDfSettings{
			ptCloudCellSize,
			1.0f, // volExpansionFactor
			DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::BlurPostprocessingType::None
		};
		const auto ptCloudDf = SDF::PointCloudDistanceFieldGenerator::Generate(ptCloud, ptCloudDfSettings);

		// Compute distance field for the mesh
		const auto meshBBox = mesh.bounds();
		const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
		const float meshMinSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		const float meshCellSize = meshMinSize / static_cast<float>(nVoxelsPerMinDimension);
		const SDF::DistanceFieldSettings meshDfSettings{
			meshCellSize,
			1.0f, // volExpansionFactor
			DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::KDTreeSplitType::Center,
			SDF::SignComputation::None, // Unsigned distance field
			SDF::BlurPostprocessingType::None,
			SDF::PreprocessingType::Octree
		};
		const auto meshDf = SDF::DistanceFieldGenerator::Generate(mesh, meshDfSettings);

		double maxDistMeshToPointCloud = std::numeric_limits<double>::lowest();
		double maxDistPointCloudToMesh = std::numeric_limits<double>::lowest();

		// Mesh to Point Cloud: Compute max distance using the distance field
		for (const auto v : mesh.vertices())
		{
			const auto vPos = mesh.position(v);
			const double vDistanceToPointCloud = TrilinearInterpolateScalarValue(vPos, ptCloudDf);
			maxDistMeshToPointCloud = std::max(maxDistMeshToPointCloud, vDistanceToPointCloud);
		}

		// Point Cloud to Mesh: Compute max distance using the distance field
		for (const auto& p : ptCloud)
		{
			const double pDistanceToMesh = TrilinearInterpolateScalarValue(p, meshDf);
			maxDistPointCloudToMesh = std::max(maxDistPointCloudToMesh, pDistanceToMesh);
		}

		// Compute Hausdorff Distance as the maximum of these two distances
		return std::max(maxDistMeshToPointCloud, maxDistPointCloudToMesh);
	}

	std::optional<double> ComputeMeshToPointCloudHausdorffDistance(const pmp::SurfaceMesh& mesh, const std::vector<pmp::Point>& ptCloud, const ScalarGrid& ptCloudDf, const unsigned int& nVoxelsPerMinDimension)
	{
		if (mesh.n_vertices() == 0)
		{
			return {};
		}

		// Compute distance field for the mesh
		const auto meshBBox = mesh.bounds();
		const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
		const float meshMinSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
		const float meshCellSize = meshMinSize / static_cast<float>(nVoxelsPerMinDimension);
		const SDF::DistanceFieldSettings meshDfSettings{
			meshCellSize,
			1.0f, // volExpansionFactor
			DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::KDTreeSplitType::Center,
			SDF::SignComputation::None, // Unsigned distance field
			SDF::BlurPostprocessingType::None,
			SDF::PreprocessingType::Octree
		};
		const auto meshDf = SDF::DistanceFieldGenerator::Generate(mesh, meshDfSettings);

		double maxDistMeshToPointCloud = std::numeric_limits<double>::lowest();
		double maxDistPointCloudToMesh = std::numeric_limits<double>::lowest();

		// Mesh to Point Cloud: Compute max distance using the distance field
		for (const auto v : mesh.vertices())
		{
			const auto vPos = mesh.position(v);
			const double vDistanceToPointCloud = TrilinearInterpolateScalarValue(vPos, ptCloudDf);
			maxDistMeshToPointCloud = std::max(maxDistMeshToPointCloud, vDistanceToPointCloud);
		}

		// Point Cloud to Mesh: Compute max distance using the distance field
		for (const auto& p : ptCloud)
		{
			const double pDistanceToMesh = TrilinearInterpolateScalarValue(p, meshDf);
			maxDistPointCloudToMesh = std::max(maxDistPointCloudToMesh, pDistanceToMesh);
		}

		// Compute Hausdorff Distance as the maximum of these two distances
		return std::max(maxDistMeshToPointCloud, maxDistPointCloudToMesh);
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

} // namespace Geometry