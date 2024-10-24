// ------------------------- NewEvolverInnerOuterLSW ------------------------------
// Summer/Fall of 2024. MDPI Mathematics Journal.
// ..................................................................................
/*
 * This work focuses on extending the framework used while publishing [Cavarga, 2024].
 * Its core idea is to solve the stopping problem, i.e.: what is the right time to stop
 * Lagrangian shrink-wrapping (LSW)? Evolving meshes flowing into gaps in the target
 * point cloud prevent the use of a straight-forward heuristic of: if it doesn't move,
 * it's done! So we counterbalance LSW by introducing an inner manifold.
 *
 * DISCLAIMER: This is work-in-progress!
 *
 * [Cavarga, 2024]
 * Cavarga, M. (2024).
 * Automated Surface Extraction: Adaptive Remeshing Meets Lagrangian Shrink-Wrapping.
 * Journal of WSCG 2024. 32, 21-30.
 */
 // --------------------------------------------------------------------------------

#include "pmp/algorithms/CurveFactory.h"
#include "pmp/algorithms/Normals.h"

#include "geometry/GridUtil.h"
#include "geometry/PlaneBuilder.h"

#include "sdf/SDF.h"

#include "utils/TimingUtils.h"

#include "core/ConversionUtils.h"
#include "core/EvolverUtilsCommon.h"
#include "core/InscribedManifold.h"
#include "core/ManifoldEvolver.h"
#include "core/SurfaceEvolver.h"

#include "IOEnvironment.h"
#include "../Experiments.h"

#include <random>

void DistanceFieldHashTest()
{
	//
	// ====== Full tetra ========
	//
	const std::vector tetraVertices = { pmp::Point(0, 0, 0), pmp::Point(1, 0, 0), pmp::Point(0, 1, 0), pmp::Point(0, 0, 1) };
	Geometry::BaseMeshGeometryData meshData;
	meshData.Vertices = tetraVertices;
	meshData.PolyIndices = { {0, 1, 3}, {1, 0, 2}, {0, 3, 2}, {1, 2, 3} };
	const Geometry::BaseMeshAdapter meshAdapter(std::make_shared<Geometry::BaseMeshGeometryData>(meshData));

	const auto meshBBox = meshAdapter.GetBounds();
	const auto meshBBoxSize = meshBBox.max() - meshBBox.min();
	const float minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
	const float cellSize = minSize / 10.0f;

	const SDF::DistanceFieldSettings sdfSettings{
		cellSize,
		1.0f,
		DBL_MAX,
		SDF::KDTreeSplitType::Center,
		SDF::SignComputation::VoxelFloodFill,
		SDF::BlurPostprocessingType::None,
		SDF::PreprocessingType::Octree
	};
	const auto sdf = SDF::DistanceFieldGenerator::Generate(meshAdapter, sdfSettings);

	const auto sdfHash = Geometry::HashScalarGrid(sdf);
	std::cout << "Hash of fullTetra_SDF: " << sdfHash << "\n";
	ExportToVTI(dataOutPath + "fullTetra_SDF", sdf);

	//
	// ============== Tetra without bottom face ===========
	//

	Geometry::BaseMeshGeometryData meshData2;
	meshData2.Vertices = tetraVertices;
	meshData2.PolyIndices = { {0, 1, 3}, {0, 3, 2}, {1, 2, 3} };
	const Geometry::BaseMeshAdapter mesh2Adapter(std::make_shared<Geometry::BaseMeshGeometryData>(meshData2));
	const SDF::DistanceFieldSettings sdfSettings2{
		cellSize,
		1.0f,
		DBL_MAX,
		SDF::KDTreeSplitType::Center,
		SDF::SignComputation::None,
		SDF::BlurPostprocessingType::None,
		SDF::PreprocessingType::Octree
	};
	const auto sdf2 = SDF::DistanceFieldGenerator::Generate(mesh2Adapter, sdfSettings2);

	const auto sdfHash2 = Geometry::HashScalarGrid(sdf2);
	std::cout << "Hash of tetraWithoutBottomFace_SDF: " << sdfHash2 << "\n";
	ExportToVTI(dataOutPath + "tetraWithoutBottomFace_SDF", sdf2);

	//
	// ============ Tetra points =========================
	//

	const SDF::PointCloudDistanceFieldSettings sdfSettings3{
			cellSize,
			1.0f,
			DBL_MAX,
			SDF::BlurPostprocessingType::None
	};

	const auto sdf3 = SDF::PointCloudDistanceFieldGenerator::Generate(tetraVertices, sdfSettings3);

	const auto sdfHash3 = Geometry::HashScalarGrid(sdf3);
	std::cout << "Hash of tetraPoints_SDF: " << sdfHash3 << "\n";
	ExportToVTI(dataOutPath + "tetraPoints_SDF", sdf3);
}


void DistanceField2DHashTest()
{
	constexpr double truncationFactor = DBL_MAX; //  0.75

	//
	// ================= Simple closed square =======================
	//

	{
		const std::vector curveVertices = { pmp::Point2(0, 0), pmp::Point2(1, 0), pmp::Point2(1, 1), pmp::Point2(0, 1) };
		Geometry::BaseCurveGeometryData curveData;
		curveData.Vertices = curveVertices;
		curveData.EdgeIndices = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
		const Geometry::BaseCurveAdapter curveAdapter(std::make_shared<Geometry::BaseCurveGeometryData>(curveData));

		const auto curveBBox = curveAdapter.GetBounds();
		const auto curveBBoxSize = curveBBox.max() - curveBBox.min();
		const float minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
		const float cellSize = minSize / 10.0f;

		const SDF::DistanceField2DSettings sdfSettings{
			cellSize,
			1.0f,
			truncationFactor,
			SDF::KDTreeSplitType::Center,
			SDF::SignComputation2D::PixelFloodFill,
			SDF::PreprocessingType2D::Quadtree
		};
		const auto sdf = SDF::PlanarDistanceFieldGenerator::Generate(curveAdapter, sdfSettings);

		const auto sdfHash = Geometry::HashScalarGrid(sdf);
		std::cout << "Hash of square_SDF: " << sdfHash << "\n";
		ExportScalarGrid2DToPNG(dataOutPath + "square_SDF.png", sdf, Geometry::BilinearInterpolateScalarValue, 10, 10);
	}

	//
	// ================= Simple open square =======================
	//

	{
		const std::vector curveVertices = { pmp::Point2(0, 0), pmp::Point2(1, 0), pmp::Point2(1, 1), pmp::Point2(0, 1) };
		Geometry::BaseCurveGeometryData curveData;
		curveData.Vertices = curveVertices;
		curveData.EdgeIndices = { {0, 1}, {1, 2}, {2, 3} };
		const Geometry::BaseCurveAdapter curveAdapter(std::make_shared<Geometry::BaseCurveGeometryData>(curveData));

		const auto curveBBox = curveAdapter.GetBounds();
		const auto curveBBoxSize = curveBBox.max() - curveBBox.min();
		const float minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
		const float cellSize = minSize / 10.0f;

		const SDF::DistanceField2DSettings sdfSettings{
			cellSize,
			1.0f,
			truncationFactor,
			SDF::KDTreeSplitType::Center,
			SDF::SignComputation2D::None,
			SDF::PreprocessingType2D::Quadtree
		};
		const auto sdf = SDF::PlanarDistanceFieldGenerator::Generate(curveAdapter, sdfSettings);

		const auto sdfHash = Geometry::HashScalarGrid(sdf);
		std::cout << "Hash of squareOpen_SDF: " << sdfHash << "\n";
		ExportScalarGrid2DToPNG(dataOutPath + "squareOpen_SDF.png", sdf, Geometry::BilinearInterpolateScalarValue, 10, 10);
	}

	//
	// ============= Square vertices pt cloud ==============
	//

	{
		const std::vector points = { pmp::Point2(0, 0), pmp::Point2(1, 0), pmp::Point2(1, 1), pmp::Point2(0, 1) };
		const auto pointBBox = pmp::BoundingBox2(points);
		const auto pointBBoxSize = pointBBox.max() - pointBBox.min();
		const float minSize = std::min(pointBBoxSize[0], pointBBoxSize[1]);
		const float cellSize = minSize / 10.0f;

		const SDF::PointCloudDistanceField2DSettings sdfSettings{
			cellSize,
			1.0f,
			truncationFactor
		};

		const auto sdf = SDF::PlanarPointCloudDistanceFieldGenerator::Generate(points, sdfSettings);

		const auto sdfHash = Geometry::HashScalarGrid(sdf);
		std::cout << "Hash of squarePts_SDF: " << sdfHash << "\n";
		ExportScalarGrid2DToPNG(dataOutPath + "squarePts_SDF.png", sdf, Geometry::BilinearInterpolateScalarValue, 10, 10);
	}

	//
	// ============= Ellipse ==============
	//

	{
		auto ellipse = pmp::CurveFactory::circle(pmp::Point2(0, 0), 1.0f, 32);
		ellipse *= scaling_matrix_2d(pmp::vec2(2.0f, 1.0f));

		const auto pointBBox = pmp::BoundingBox2(ellipse.positions());
		const auto pointBBoxSize = pointBBox.max() - pointBBox.min();
		const float minSize = std::min(pointBBoxSize[0], pointBBoxSize[1]);
		const float cellSize = minSize / 20.0f;

		const SDF::PointCloudDistanceField2DSettings sdfSettings{
			cellSize,
			1.0f,
			truncationFactor
		};

		const auto sdf = SDF::PlanarPointCloudDistanceFieldGenerator::Generate(ellipse.positions(), sdfSettings);

		const auto sdfHash = Geometry::HashScalarGrid(sdf);
		std::cout << "Hash of ellipse_SDF: " << sdfHash << "\n";
		ExportScalarGrid2DToPNG(dataOutPath + "ellipse_SDF.png", sdf, Geometry::BilinearInterpolateScalarValue, 10, 10);
	}
}


void ManifoldCurve2DTests()
{
	const auto closedArc = pmp::CurveFactory::circle(pmp::Point2(0, 0), 1.0, 32);

	if (!ExportManifoldCurve2DToPLY(closedArc, dataOutPath + "closedArc.ply"))
	{
		std::cerr << "ExportManifoldCurve2DToPLY: internal error!\n";
	}

	const auto openArc = pmp::CurveFactory::circle(pmp::Point2(0, 0), 1.0, 16, 0, M_PI);

	if (!ExportManifoldCurve2DToPLY(openArc, dataOutPath + "openArc.ply"))
	{
		std::cerr << "ExportManifoldCurve2DToPLY: internal error!\n";
	}
}


void PointCloud2DSliceTests()
{
	// extract 2D point clouds
	const std::vector<std::string> meshForPtCloudNames{
		"armadillo",
		"bunny",
		"maxPlanck",
		"nefertiti"
	};

	const std::map<std::string, pmp::Point> slicingPlaneRefPts{
		{"armadillo", pmp::Point{-0.10348621158928151f, 21.427067319905646f, 9.79369240592005f}},
		{"bunny", pmp::Point{-0.01684039831161499f, 0.11015420407056808f, 0.0012007840834242693f} },
		{"maxPlanck", pmp::Point{30.59686279296875f, -18.105804443359375f, 82.29149055480957f} },
		{"nefertiti", pmp::Point{0.0f, 0.0f, 0.0f} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"armadillo", pmp::vec3{-0.03070969905335075f, 0.12876712096541565f, 0.9911992448253433f}},
		{"bunny", pmp::vec3{0.0f, 0.0f, 1.0f} },
		{"maxPlanck", pmp::vec3{1.0f, 0.0f, 0.0f} },
		{"nefertiti", pmp::vec3{1.0f, 0.0f, 0.0f} }
	};

	for (const auto& meshName : meshForPtCloudNames)
	{
		constexpr size_t samplingLevel = 2;
		constexpr size_t nSamplings = 10;
		constexpr size_t minVerts = 9; // Minimum number of vertices to sample

		constexpr unsigned int seed = 5000; // seed for the pt cloud sampling RNG

		std::cout << "==================================================================\n";
		std::cout << "Mesh To Pt Cloud 2D: " << meshName << ".obj -> " << meshName << "_Pts_2D.ply\n";
		std::cout << "------------------------------------------------------------------\n";
		const auto baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", false);
		if (!baseDataOpt.has_value())
		{
			std::cerr << "baseDataOpt == nullopt!\n";
			break;
		}
		std::cout << meshName << ".obj" << " imported as BaseMeshGeometryData.\n";

		const auto& baseData = baseDataOpt.value();
		const size_t maxVerts = baseData.Vertices.size(); // Maximum number of vertices available
		size_t nVerts = minVerts + (maxVerts - minVerts) * samplingLevel / (nSamplings - 1);
		nVerts = std::max(minVerts, std::min(nVerts, maxVerts));

		std::cout << "Sampling " << nVerts << "/" << maxVerts << " vertices...\n";

		const auto pts3D = SamplePointsFromMeshData(baseData, nVerts, seed);
		if (pts3D.empty())
		{
			std::cerr << "SamplePointsFromMeshData sampled no points for mesh " << meshName << "!\n";
			continue;
		}

		const pmp::BoundingBox ptCloudBBox(pts3D);
		const auto center = ptCloudBBox.center();
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		//const float maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float distTolerance = 0.01f * minSize;

		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0f, 0.0f, 0.0f });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(pts3D, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}

		if (!Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
		{
			std::cerr << "Export2DPointCloudToPLY: internal error during export!\n";
			continue;
		}
	}
}


void CurveReorientTests()
{
	auto circleCurve = pmp::CurveFactory::circle(pmp::Point2{}, 1.0f, 25);
	circleCurve.negate_orientation();
	if (!pmp::write_to_ply(circleCurve, dataOutPath + "circleCurve_Flipped.ply"))
		std::cerr << "Error writing file!\n";
}


void MeshReorientTests()
{
	const std::vector<std::string> meshNames{
		"ico",
		"armadillo",
		"blub",
		"bunny",
		"maxPlanck",
		"nefertiti",
		"ogre",
		"spot"
	};

	for (const auto& meshName : meshNames)
	{
		std::cout << "----------------------------------------------\n";
		std::cout << " flipping " << meshName << ".obj ...\n";
		pmp::SurfaceMesh mesh;
		mesh.read(dataDirPath + meshName + ".obj");
		mesh.negate_orientation();

		//// Function to test forward and backward circulation and check topological consistency of halfedges
		//auto test_halfedge_circulation = [&](pmp::Face f) {
		//	pmp::Halfedge h0 = mesh.halfedge(f); // Start with the first halfedge of the face
		//	pmp::Halfedge h = h0;

		//	bool testPassed = true;

		//	// Check forward and backward circulation and test key methods
		//	do
		//	{
		//		pmp::Halfedge next_h = mesh.next_halfedge(h);
		//		pmp::Halfedge prev_h = mesh.prev_halfedge(h);
		//		pmp::Vertex from_v = mesh.from_vertex(h);
		//		pmp::Vertex to_v = mesh.to_vertex(h);

		//		// 1. Test if next and previous halfedge form a consistent loop
		//		if (mesh.next_halfedge(prev_h) != h)
		//		{
		//			//std::cout << "ERR [0]: next_halfedge(prev_halfedge(h)) != h for halfedge " << h.idx() << "\n";
		//			testPassed = false;
		//		}
		//		if (mesh.prev_halfedge(next_h) != h)
		//		{
		//			//std::cout << "ERR [1]: prev_halfedge(next_halfedge(h)) != h for halfedge " << h.idx() << "\n";
		//			testPassed = false;
		//		}

		//		// 2. Check vertex connectivity between next/prev
		//		if (mesh.from_vertex(next_h) != to_v)
		//		{
		//			//std::cout << "ERR [2]: to_vertex of current halfedge does not match from_vertex of next_halfedge at halfedge " << h.idx() << "\n";
		//			testPassed = false;
		//		}
		//		if (mesh.to_vertex(prev_h) != from_v)
		//		{
		//			//std::cout << "ERR [3]: from_vertex of current halfedge does not match to_vertex of prev_halfedge at halfedge " << h.idx() << "\n";
		//			testPassed = false;
		//		}

		//		// 4. Output debug info for vertices and positions
		//		std::cout << "Halfedge " << h.idx() << " connects vertex " << from_v.idx() << " to vertex " << to_v.idx() << "\n";

		//		h = next_h;
		//	} while (h != h0); // Continue until we've looped back to the start

		//	return testPassed;
		//};

		//// Iterate through faces and check halfedge circulation and topological consistency
		//for (auto f : mesh.faces())
		//{
		//	if (!test_halfedge_circulation(f))
		//	{
		//		std::cout << "Circulation test failed for face " << f.idx() << "\n";
		//	}
		//}

		pmp::Normals::compute_vertex_normals(mesh);
		auto vNormalsProp = mesh.get_vertex_property<pmp::vec3>("v:normal");

		int seed = 9999;
		std::mt19937 rng(seed); // Random number generator with a seed
		std::uniform_int_distribution<size_t> dist(0, mesh.n_vertices() - 1); // Uniform distribution

		// Randomly sample a small subset of vertex normals
		size_t num_samples = std::min(mesh.n_vertices(), static_cast<size_t>(10)); // Adjust number of samples
		std::cout << "Randomly sampled vertex normals:\n";
		for (size_t i = 0; i < num_samples; ++i)
		{
			size_t random_index = dist(rng); // Get a random index
			pmp::Vertex v(random_index);     // Create a vertex handle using the random index

			// Output the normal of the sampled vertex
			auto normal = vNormalsProp[v];
			std::cout << "Vertex " << random_index << ": Normal = ("
				<< normal[0] << ", " << normal[1] << ", " << normal[2] << ")\n";
		}

		mesh.write(dataOutPath + meshName + "_Flipped.obj");
	}
}


void PropertyPerformanceTests()
{
	const std::vector<std::string> meshNames{
		"BentChair",
		"blub",
		"bunny",
		"maxPlanck",
		//"nefertiti", // takes ~16.5 sec per 80-step evolution!
		"ogre",
		"spot"
	};

	size_t nMainLoopIters = 80;

	/// ====== The old evolver mock ======

	class OrigSurfaceEvolverMock
	{
	public:

		OrigSurfaceEvolverMock(const std::string& meshName, size_t nSteps)
			: m_NSteps(nSteps)
		{
			m_Mesh = std::make_shared<pmp::SurfaceMesh>();
			m_Mesh->read(dataDirPath + meshName + ".obj");
		}

		void Evolve()
		{
			auto vProp1 = m_Mesh->add_vertex_property<pmp::Scalar>("v:prop1");
			auto vProp2 = m_Mesh->add_vertex_property<pmp::Scalar>("v:prop2");

			pmp::VertexProperty<pmp::Point> vNormalsProp{};

			auto aFunctionThatNeedsNormals = [&]()
			{
				pmp::vec3 vNormalSum{};
				for (const auto v : m_Mesh->vertices())
				{
					vNormalSum += vNormalsProp[v];
				}
				vNormalSum;
			};

			for (const auto v : m_Mesh->vertices())
			{
				vProp1[v] = static_cast<pmp::Scalar>(v.idx());
				vProp2[v] = static_cast<pmp::Scalar>(m_Mesh->n_vertices() - v.idx());
			}

			for (size_t i = 0; i < m_NSteps; ++i)
			{
				pmp::Normals::compute_vertex_normals(*m_Mesh);
				vNormalsProp = m_Mesh->vertex_property<pmp::Point>("v:normal");

				aFunctionThatNeedsNormals();

				for (const auto v : m_Mesh->vertices())
				{
					vProp1[v] = static_cast<pmp::Scalar>(v.idx());
					vProp2[v] = static_cast<pmp::Scalar>(m_Mesh->n_vertices() - v.idx());
				}
			}
		}

	private:

		size_t m_NSteps{ 10 };
		std::shared_ptr<pmp::SurfaceMesh> m_Mesh{ nullptr };
	};

	/// ====== These objects simulate the new evolver strategy approach ======

	class NewEvolverStrategyMock
	{
	public:
		NewEvolverStrategyMock(const std::string& meshName)
		{
			m_Mesh = std::make_shared<pmp::SurfaceMesh>();
			m_Mesh->read(dataDirPath + meshName + ".obj");
		}

		void Preprocess()
		{
			auto vProp1 = m_Mesh->add_vertex_property<pmp::Scalar>("v:prop1");
			auto vProp2 = m_Mesh->add_vertex_property<pmp::Scalar>("v:prop2");
		}

		void PerformEvolutionStep(size_t step)
		{
			auto vProp1 = m_Mesh->vertex_property<pmp::Scalar>("v:prop1");
			auto vProp2 = m_Mesh->vertex_property<pmp::Scalar>("v:prop2");

			pmp::Normals::compute_vertex_normals(*m_Mesh);
			AFunctionThatNeedsNormals();

			for (const auto v : m_Mesh->vertices())
			{
				vProp1[v] = static_cast<pmp::Scalar>(v.idx());
				vProp2[v] = static_cast<pmp::Scalar>(m_Mesh->n_vertices() - v.idx());
			}
		}

	private:

		void AFunctionThatNeedsNormals()
		{
			auto vNormalsProp = m_Mesh->vertex_property<pmp::Point>("v:normal");
			pmp::vec3 vNormalSum{};
			for (const auto v : m_Mesh->vertices())
			{
				vNormalSum += vNormalsProp[v];
			}
			vNormalSum;
		}

		std::shared_ptr<pmp::SurfaceMesh> m_Mesh{ nullptr };
	};

	class NewEvolverMock
	{
	public:

		NewEvolverMock(size_t nSteps, std::shared_ptr<NewEvolverStrategyMock> strategy)
			: m_NSteps(nSteps), m_Strategy(std::move(strategy))
		{
		}

		void Evolve() const
		{
			for (size_t i = 0; i < m_NSteps; ++i)
			{
				m_Strategy->PerformEvolutionStep(i);
			}
		}
	private:

		size_t m_NSteps{ 10 };
		std::shared_ptr<NewEvolverStrategyMock> m_Strategy{ nullptr };
	};

	constexpr size_t nRuns = 20;

	for (const auto& name : meshNames)
	{
		std::cout << "-------------------------------------------------------------------------------------------\n";
		std::cout << "OrigSurfaceEvolverMock vs {NewEvolverStrategyMock, NewEvolverMock} for \"" << name << "\":\n";
		AVERAGE_TIMING(OrigEvolverPropertiesApproach, nRuns, {
			OrigSurfaceEvolverMock se(name, nMainLoopIters);
			se.Evolve();
			}, true);

		AVERAGE_TIMING(ProposedEncapsulatedApproach, nRuns, {
			auto strategy = std::make_shared<NewEvolverStrategyMock>(name);
			NewEvolverMock se(nMainLoopIters, strategy);
			se.Evolve();
			}, true);
	}
}

void DeletePointsContainedInCircles(std::vector<pmp::Point2>& pointsToFilter, const std::vector<Circle2D>& cutCircles)
{
	auto isPointInCircle = [](const pmp::Point2& point, const Circle2D& circle) {
		return sqrnorm(point - circle.Center) <= (circle.Radius * circle.Radius);
	};

	auto isInAnyCircle = [&cutCircles, &isPointInCircle](const pmp::Point2& point) {
		for (const auto& circle : cutCircles) {
			if (isPointInCircle(point, circle)) return true;
		}
		return false;
	};

	pointsToFilter.erase(std::ranges::remove_if(pointsToFilter, isInAnyCircle).begin(), pointsToFilter.end());
}

void DeletePointsContainedInSpheres(std::vector<pmp::Point>& pointsToFilter, const std::vector<Sphere3D>& cutSpheres)
{
	auto isPointInSphere = [](const pmp::Point& point, const Sphere3D& sphere) {
		return sqrnorm(point - sphere.Center) <= (sphere.Radius * sphere.Radius);
	};

	auto isInAnySphere = [&cutSpheres, &isPointInSphere](const pmp::Point& point) {
		for (const auto& sphere : cutSpheres) {
			if (isPointInSphere(point, sphere)) return true;
		}
		return false;
	};

	pointsToFilter.erase(std::ranges::remove_if(pointsToFilter, isInAnySphere).begin(), pointsToFilter.end());
}

void DeletePointsWithNormalsContainedInSpheres(std::vector<std::pair<pmp::Point, pmp::vec3>>& pointsToFilter, const std::vector<Sphere3D>& cutSpheres)
{
	auto isPointInSphere = [](const pmp::Point& point, const Sphere3D& sphere) {
		return sqrnorm(point - sphere.Center) <= (sphere.Radius * sphere.Radius);
	};

	auto isInAnySphere = [&cutSpheres, &isPointInSphere](const std::pair<pmp::Point, pmp::vec3>& point) {
		for (const auto& sphere : cutSpheres) {
			if (isPointInSphere(point.first, sphere)) return true;
		}
		return false;
	};

	pointsToFilter.erase(std::ranges::remove_if(pointsToFilter, isInAnySphere).begin(), pointsToFilter.end());
}


void OldVsNewLSWTests()
{
	const std::vector<std::string> meshForPtCloudNames{
		//"armadillo",
		//"blub",
		"bunny",
		//"maxPlanck",
		//"nefertiti",
		//"ogre",
		//"spot"
	};
	const std::map<std::string, double> timeStepSizesForPtClouds{
		{"armadillo", 0.05 },
		{ "blub", 0.05 },
		{ "bunny", 0.05 },
		{ "maxPlanck", 0.05 },
		{ "nefertiti", 0.05 },
		{ "ogre", 0.05 },
		{ "spot", 0.05 }
	};
	const std::map<std::string, double> isoLevelOffsetFactors{
		{"armadillo", 0.5 },
		{ "blub", 0.5 },
		{ "bunny", 0.5 },
		{ "maxPlanck", 0.5 },
		{ "nefertiti", 0.5 },
		{ "ogre", 0.5 },
		{ "spot", 0.5 }
	};

	const std::map<std::string, pmp::Point> slicingPlaneRefPts{
		{"armadillo", pmp::Point{-0.10348621158928151f, 21.427067319905646f, 9.79369240592005f}},
		{"bunny", pmp::Point{-0.01684039831161499f, 0.11015420407056808f, 0.0012007840834242693f} },
		{"maxPlanck", pmp::Point{30.59686279296875f, -18.105804443359375f, 82.29149055480957f} },
		{"nefertiti", pmp::Point{0.0f, 0.0f, 0.0f} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"armadillo", pmp::vec3{-0.03070969905335075f, 0.12876712096541565f, 0.9911992448253433f}},
		{"bunny", pmp::vec3{0.0f, 0.0f, 1.0f} },
		{"maxPlanck", pmp::vec3{1.0f, 0.0f, 0.0f} },
		{"nefertiti", pmp::vec3{1.0f, 0.0f, 0.0f} }
	};

	const std::map<std::string, Circle2D> outerCircles{
		{"armadillo", Circle2D{pmp::Point2{0.372234f, 16.6515f}, 121.558f} },
		{"bunny", Circle2D{pmp::Point2{-0.0155906f, 0.102261f}, 0.142831f} },
		{"maxPlanck", Circle2D{pmp::Point2{-17.82f, 82.5006f}, 292.263f} },
		{"nefertiti", Circle2D{pmp::Point2{0.178497f, -0.0410004f}, 441.436f} }
	};

	const std::map<std::string, Sphere3D> outerSpheres{
		{"armadillo", Sphere3D{pmp::Point{0.0122509f, 21.4183f, -0.000249863f}, 136.963f} },
		{"bunny", Sphere3D{pmp::Point{-0.0168297f, 0.110217f, -0.0015718f}, 0.141622f} },
		{"maxPlanck", Sphere3D{pmp::Point{30.658f, -17.9765f, 82.2885f}, 271.982f} },
		{"nefertiti", Sphere3D{pmp::Point{0.0144997f, -0.00499725f, -0.0215073f}, 392.184f} }
	};

	constexpr unsigned int nVoxelsPerMinDimension = 50;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;

	constexpr size_t samplingLevel = 3;
	constexpr size_t nSamplings = 10;
	constexpr size_t minVerts = 9; // Minimum number of vertices to sample

	constexpr unsigned int seed = 5000; // seed for the pt cloud sampling RNG

	//SetRemeshingAdjustmentTimeIndices({}); // no remeshing adjustment
	SetRemeshingAdjustmentTimeIndices({ /*3, 10,*/ 30 /*, 50 , 100, 120, 140, 145*/ });

	constexpr unsigned int NTimeSteps = 400;

	constexpr bool executeOldEvolver = false;
	constexpr bool executeNewSurfaceEvolver = false;
	constexpr bool executeNewCurveEvolver = true;
	constexpr bool executeNewSurfaceCustomEvolver = false;
	constexpr bool executeNewCurveCustomEvolver = false;

	for (const auto& meshName : meshForPtCloudNames)
	{
		// =======================================================================
		//   - - - - - - - - - - - - -   Data   Prep   - - - - - - - - - - - -
		// -----------------------------------------------------------------------

		std::cout << "==================================================================\n";
		std::cout << "Mesh To Pt Cloud: " << meshName << ".obj -> " << meshName << "Pts_" << samplingLevel << ".ply\n";
		std::cout << "------------------------------------------------------------------\n";
		const auto baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", false);
		if (!baseDataOpt.has_value())
		{
			std::cerr << "baseDataOpt == nullopt!\n";
			break;
		}
		std::cout << "meshName.obj" << " imported as BaseMeshGeometryData.\n";
		const auto& baseData = baseDataOpt.value();
		const size_t maxVerts = baseData.Vertices.size(); // Maximum number of vertices available
		size_t nVerts = minVerts + (maxVerts - minVerts) * samplingLevel / (nSamplings - 1);
		nVerts = std::max(minVerts, std::min(nVerts, maxVerts));

		std::cout << "Sampling " << nVerts << "/" << maxVerts << " vertices...\n";

		std::string filename = dataOutPath + meshName + "Pts_" + std::to_string(samplingLevel) + ".ply";
		if (!ExportSampledVerticesToPLY(baseData, nVerts, filename, seed))
		{
			std::cerr << "ExportSampledVerticesToPLY failed!\n";
			break;
		}

		const auto ptCloudName = meshName + "Pts_" + std::to_string(samplingLevel);
		const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataOutPath + ptCloudName + ".ply", true);
		if (!ptCloudOpt.has_value())
		{
			std::cerr << "ptCloudOpt == nullopt!\n";
			break;
		}

		const auto& ptCloud = ptCloudOpt.value();
		const pmp::BoundingBox ptCloudBBox(ptCloud);
		const auto center = ptCloudBBox.center();
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float cellSize = minSize / nVoxelsPerMinDimension;
		constexpr float volExpansionFactor = 1.0f;
		const SDF::PointCloudDistanceFieldSettings dfSettings{
			cellSize,
			volExpansionFactor,
			Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::BlurPostprocessingType::None
		};

		const double isoLvlOffsetFactor = (timeStepSizesForPtClouds.contains(ptCloudName) ? isoLevelOffsetFactors.at(ptCloudName) : defaultOffsetFactor);
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		const double tau = (timeStepSizesForPtClouds.contains(ptCloudName) ? timeStepSizesForPtClouds.at(ptCloudName) : defaultTimeStep); // time step

		// ==========================================================================
		// - - - - - - - - - - -    old SurfaceEvolver   - - - - - - - - - - - -
		// ==========================================================================

		if (executeOldEvolver)
		{
			std::cout << "==================================================================\n";
			std::cout << "Pt Cloud to DF: " << ptCloudName << ".ply -> " << ptCloudName << "_DF_" << nVoxelsPerMinDimension << "voxPerMinDim.vti\n";
			std::cout << "------------------------------------------------------------------\n";

			const auto startSDF = std::chrono::high_resolution_clock::now();
			auto sdf = SDF::PointCloudDistanceFieldGenerator::Generate(ptCloud, dfSettings);
			const auto endSDF = std::chrono::high_resolution_clock::now();
			SDF::ReportOutput(sdf, std::cout);
			const std::chrono::duration<double> timeDiff = endSDF - startSDF;
			std::cout << "DF Time: " << timeDiff.count() << " s\n";
			std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";

			std::cout << "Setting up SurfaceEvolutionSettings.\n";

			MeshTopologySettings topoParams;
			topoParams.MinEdgeMultiplier = 0.14f;
			topoParams.UseBackProjection = false;
			topoParams.PrincipalCurvatureFactor = 3.2f;
			topoParams.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);
			topoParams.EdgeLengthDecayFactor = 0.7f;
			topoParams.ExcludeEdgesWithoutBothFeaturePts = true; // this one is internally set to true for ManifoldEvolver by default
			topoParams.FeatureType = FeatureDetectionType::MeanCurvature;

			AdvectionDiffusionParameters adParams{
				1.0, 1.0,
				1.0, 1.0
			};

			SurfaceEvolutionSettings seSettings{
				ptCloudName,
				NTimeSteps,
				tau,
				fieldIsoLevel,
				3, // IcoSphereSubdivisionLevel
				adParams,
				topoParams,
				minSize, maxSize,
				ptCloudBBox.center(),
				true, false,
				dataOutPath,
				MeshLaplacian::Voronoi,
				{"minAngle", "maxAngle", "jacobianConditionNumber", "equilateralJacobianCondition",/* "stiffnessMatrixConditioning" */},
				0.05f,
				true,
				false
			};
			ReportInput(seSettings, std::cout);

			std::cout << "Setting up SurfaceEvolver.\n";

			SurfaceEvolver oldEvolver(sdf, volExpansionFactor, seSettings);

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				oldEvolver.Evolve();
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > SurfaceEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}

			std::cout << "done.\n";
		}

		// ==========================================================================
		// - - - - - - - - -  New Manifold Evolver (Surface)   - - - - - - - - - - - 
		// ==========================================================================

		if (executeNewSurfaceEvolver)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = false;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "newEvol_Pts" + std::to_string(samplingLevel);
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			std::cout << "Setting up ManifoldSurfaceEvolutionStrategy.\n";

			// Set up the evolution strategy
			auto surfaceStrategy = std::make_shared<ManifoldSurfaceEvolutionStrategy>(
				strategySettings, MeshLaplacian::Voronoi,
				std::make_shared<std::vector<pmp::Point>>(ptCloud));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver newEvolver(globalSettings, std::move(surfaceStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				newEvolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}

		// ==========================================================================
		// - - - - - - - - -  New Manifold Evolver (Curve)  - - - - - - - - - - - - 
		// ==========================================================================

		const float distTolerance = 0.01f * minSize;
		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0f, 0.0f, 0.0f });
		auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}
		const std::vector cutCircles{ 
			Circle2D{pmp::Point2{-0.02f, 0.04f}, 0.015f },
			Circle2D{pmp::Point2{0.03f, 0.04f}, 0.015f } };
		DeletePointsContainedInCircles(pts2D, cutCircles);

		if (!Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
		{
			std::cerr << "Export2DPointCloudToPLY: internal error during export!\n";
			continue;
		}

		// -----------------------
		if (executeNewCurveEvolver)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = false;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 0.5 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 4;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "newEvol_Pts2D" + std::to_string(samplingLevel);
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			auto curveStrategy = std::make_shared<ManifoldCurveEvolutionStrategy>(
				strategySettings, std::make_shared<std::vector<pmp::Point2>>(pts2D));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}

		if (executeNewSurfaceCustomEvolver)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "newEvol_Pts" + std::to_string(samplingLevel);
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerSpheres.contains(meshName))
			{
				std::cerr << "!outerSpheres.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			const auto outerSphere = outerSpheres.at(meshName);
			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;

			// construct outer ico-sphere
			Geometry::IcoSphereBuilder icoBuilder({ strategySettings.LevelOfDetail, outerSphere.Radius });
			icoBuilder.BuildBaseData();
			icoBuilder.BuildPMPSurfaceMesh();
			auto outerSurface = icoBuilder.GetPMPSurfaceMeshResult();
			const pmp::mat4 transfMatrixGeomMove{
				1.0f, 0.0f, 0.0f, outerSphere.Center[0],
				0.0f, 1.0f, 0.0f, outerSphere.Center[1],
				0.0f, 0.0f, 1.0f, outerSphere.Center[2],
				0.0f, 0.0f, 0.0f, 1.0f
			};
			outerSurface *= transfMatrixGeomMove;

			std::vector<pmp::SurfaceMesh> innerSurfaces;
			auto surfaceStrategy = std::make_shared<CustomManifoldSurfaceEvolutionStrategy>(
				strategySettings, MeshLaplacian::Voronoi,
				outerSurface, innerSurfaces,
				std::make_shared<std::vector<pmp::Point>>(ptCloud));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(surfaceStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}

		if (executeNewCurveCustomEvolver)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "newEvol_Pts2D" + std::to_string(samplingLevel);
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerCircles.contains(meshName))
			{
				std::cerr << "!outerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			const auto outerCircle = outerCircles.at(meshName);
			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);
			std::vector<pmp::ManifoldCurve2D> innerCurves;
			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, outerCurve, innerCurves,
				std::make_shared<std::vector<pmp::Point2>>(pts2D));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	}
}


void PairedLSWTests()
{
	const std::vector<std::string> meshForPtCloudNames{
		//"armadillo",
		//"blub",
		//"bunny",
		//"maxPlanck",
		"nefertiti",
		//"ogre",
		//"spot"
	};
	const std::map<std::string, double> timeStepSizesForPtClouds{
		{"armadillo", 0.05 },
		{ "blub", 0.05 },
		{ "bunny", 0.05 },
		{ "maxPlanck", 0.05 },
		{ "nefertiti", 0.05 },
		{ "ogre", 0.05 },
		{ "spot", 0.05 }
	};
	const std::map<std::string, double> isoLevelOffsetFactors{
		{"armadillo", 0.5 },
		{ "blub", 0.5 },
		{ "bunny", 0.5 },
		{ "maxPlanck", 0.5 },
		{ "nefertiti", 0.5 },
		{ "ogre", 0.5 },
		{ "spot", 0.5 }
	};

	const std::map<std::string, pmp::Point> slicingPlaneRefPts{
		{"armadillo", pmp::Point{-0.10348621158928151f, 21.427067319905646f, 9.79369240592005f}},
		{"bunny", pmp::Point{-0.01684039831161499f, 0.11015420407056808f, 0.0012007840834242693f} },
		{"maxPlanck", pmp::Point{30.59686279296875f, -18.105804443359375f, 82.29149055480957f} },
		{"nefertiti", pmp::Point{0.0f, 0.0f, 0.0f} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"armadillo", pmp::vec3{-0.03070969905335075f, 0.12876712096541565f, 0.9911992448253433f}},
		{"bunny", pmp::vec3{0.0f, 0.0f, 1.0f} },
		{"maxPlanck", pmp::vec3{1.0f, 0.0f, 0.0f} },
		{"nefertiti", pmp::vec3{1.0f, 0.0f, 0.0f} }
	};

	const std::map<std::string, Circle2D> outerCircles{
		{"armadillo", Circle2D{pmp::Point2{0.372234f, 16.6515f}, 121.558f} },
		{"bunny", Circle2D{pmp::Point2{-0.0155906f, 0.102261f}, 0.142831f} },
		{"maxPlanck", Circle2D{pmp::Point2{-17.82f, 82.5006f}, 292.263f} },
		{"nefertiti", Circle2D{pmp::Point2{0.178497f, -0.0410004f}, 441.436f} }
	};
	const std::map<std::string, Circle2D> innerCircles{
		{"armadillo", Circle2D{pmp::Point2{-3.0f, 52.0f}, 20.0f}},
		{"bunny", Circle2D{pmp::Point2{-0.025f, 0.08f}, 0.025f}},
		{"maxPlanck", Circle2D{pmp::Point2{8.0f, 85.0f}, 50.0f}},
		{"nefertiti", Circle2D{pmp::Point2{-20.0f, 100.0f}, 55.0f}}
	};

	const std::map<std::string, Sphere3D> outerSpheres{
		{"armadillo", Sphere3D{pmp::Point{0.0122509f, 21.4183f, -0.000249863f}, 136.963f} },
		{"bunny", Sphere3D{pmp::Point{-0.0168297f, 0.110217f, -0.0015718f}, 0.141622f} },
		{"maxPlanck", Sphere3D{pmp::Point{30.658f, -17.9765f, 82.2885f}, 271.982f} },
		{"nefertiti", Sphere3D{pmp::Point{0.0144997f, -0.00499725f, -0.0215073f}, 392.184f} }
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;

	constexpr size_t samplingLevel = 3;
	constexpr size_t nSamplings = 10;
	constexpr size_t minVerts = 9; // Minimum number of vertices to sample

	constexpr unsigned int seed = 5000; // seed for the pt cloud sampling RNG

	//SetRemeshingAdjustmentTimeIndices({}); // no remeshing adjustment
	SetRemeshingAdjustmentTimeIndices({ /*3, 10,*/ 30 /*, 50 , 100, 120, 140, 145*/ });

	constexpr unsigned int NTimeSteps = 180;

	constexpr bool executeCustomCurveEvolver = true;
	constexpr bool executeCustomSurfaceEvolver = false;

	for (const auto& meshName : meshForPtCloudNames)
	{
		// =======================================================================
		//   - - - - - - - - - - - - -   Data   Prep   - - - - - - - - - - - -
		// -----------------------------------------------------------------------

		std::cout << "==================================================================\n";
		std::cout << "Mesh To Pt Cloud: " << meshName << ".obj -> " << meshName << "Pts_" << samplingLevel << ".ply\n";
		std::cout << "------------------------------------------------------------------\n";
		const auto baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", false);
		if (!baseDataOpt.has_value())
		{
			std::cerr << "baseDataOpt == nullopt!\n";
			break;
		}
		std::cout << "meshName.obj" << " imported as BaseMeshGeometryData.\n";
		const auto& baseData = baseDataOpt.value();
		const size_t maxVerts = baseData.Vertices.size(); // Maximum number of vertices available
		size_t nVerts = minVerts + (maxVerts - minVerts) * samplingLevel / (nSamplings - 1);
		nVerts = std::max(minVerts, std::min(nVerts, maxVerts));

		std::cout << "Sampling " << nVerts << "/" << maxVerts << " vertices...\n";

		std::string filename = dataOutPath + meshName + "Pts_" + std::to_string(samplingLevel) + ".ply";
		if (!ExportSampledVerticesToPLY(baseData, nVerts, filename, seed))
		{
			std::cerr << "ExportSampledVerticesToPLY failed!\n";
			break;
		}

		const auto ptCloudName = meshName + "Pts_" + std::to_string(samplingLevel);
		const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataOutPath + ptCloudName + ".ply", true);
		if (!ptCloudOpt.has_value())
		{
			std::cerr << "ptCloudOpt == nullopt!\n";
			break;
		}

		const auto& ptCloud = ptCloudOpt.value();
		const pmp::BoundingBox ptCloudBBox(ptCloud);
		const auto center = ptCloudBBox.center();
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		// const float maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float cellSize = minSize / nVoxelsPerMinDimension;

		const double isoLvlOffsetFactor = (isoLevelOffsetFactors.contains(ptCloudName) ? isoLevelOffsetFactors.at(ptCloudName) : defaultOffsetFactor);
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		const double tau = (timeStepSizesForPtClouds.contains(ptCloudName) ? timeStepSizesForPtClouds.at(ptCloudName) : defaultTimeStep); // time step

		// ==========================================================================
		// - - - - - - - - -  New Manifold Evolver (Curve)  - - - - - - - - - - - - 
		// ==========================================================================

		const float distTolerance = 0.01f * minSize;
		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0f, 0.0f, 0.0f });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}

		if (!Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
		{
			std::cerr << "Export2DPointCloudToPLY: internal error during export!\n";
			continue;
		}

		// -----------------------
		if (executeCustomSurfaceEvolver)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "newEvol_Pts" + std::to_string(samplingLevel);
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerSpheres.contains(meshName))
			{
				std::cerr << "!outerSpheres.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			const auto outerSphere = outerSpheres.at(meshName);
			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;

			// construct outer ico-sphere
			Geometry::IcoSphereBuilder icoBuilder({ strategySettings.LevelOfDetail, outerSphere.Radius });
			icoBuilder.BuildBaseData();
			icoBuilder.BuildPMPSurfaceMesh();
			auto outerSurface = icoBuilder.GetPMPSurfaceMeshResult();
			const pmp::mat4 transfMatrixGeomMove{
				1.0f, 0.0f, 0.0f, outerSphere.Center[0],
				0.0f, 1.0f, 0.0f, outerSphere.Center[1],
				0.0f, 0.0f, 1.0f, outerSphere.Center[2],
				0.0f, 0.0f, 0.0f, 1.0f
			};
			outerSurface *= transfMatrixGeomMove;

			std::vector<pmp::SurfaceMesh> innerSurfaces;
			auto surfaceStrategy = std::make_shared<CustomManifoldSurfaceEvolutionStrategy>(
				strategySettings, MeshLaplacian::Voronoi,
				outerSurface, innerSurfaces,
				std::make_shared<std::vector<pmp::Point>>(ptCloud));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(surfaceStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}

		if (executeCustomCurveEvolver)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				return 0.0 * TRIVIAL_EPSILON(distance);
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 1.5 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 4;
			strategySettings.TangentialVelocityWeight = 0.05;

			//strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.22f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "newEvol_Pts2D" + std::to_string(samplingLevel) + "_Repulsionless";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerCircles.contains(meshName))
			{
				std::cerr << "!outerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}
			if (!innerCircles.contains(meshName))
			{
				std::cerr << "!innerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			const auto outerCircle = outerCircles.at(meshName);
			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);

			const auto innerCircle = innerCircles.at(meshName);
			const auto nInnerSegments = static_cast<unsigned int>(static_cast<pmp::Scalar>(nSegments) * innerCircle.Radius / outerCircle.Radius * 2);
			auto innerCurve = pmp::CurveFactory::circle(innerCircle.Center, innerCircle.Radius, nInnerSegments);
			innerCurve.negate_orientation();
			std::vector innerCurves{ innerCurve };

			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, outerCurve, innerCurves,
				std::make_shared<std::vector<pmp::Point2>>(pts2D));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	}
}


void PairedLSWRepulsionTests()
{
	const std::vector<std::string> meshForPtCloudNames{
		//"armadillo",
		//"blub",
		"bunny",
		//"maxPlanck",
		//"nefertiti",
		//"ogre",
		//"spot"
	};
	const std::map<std::string, double> timeStepSizesForPtClouds{
		{"armadillo", 0.05 },
		{ "blub", 0.05 },
		{ "bunny", 0.05 },
		{ "maxPlanck", 0.05 },
		{ "nefertiti", 0.05 },
		{ "ogre", 0.05 },
		{ "spot", 0.05 }
	};
	const std::map<std::string, double> isoLevelOffsetFactors{
		{"armadillo", 0.5 },
		{ "blub", 0.5 },
		{ "bunny", 0.5 },
		{ "maxPlanck", 0.5 },
		{ "nefertiti", 0.5 },
		{ "ogre", 0.5 },
		{ "spot", 0.5 }
	};

	const std::map<std::string, pmp::Point> slicingPlaneRefPts{
		{"armadillo", pmp::Point{-0.10348621158928151f, 21.427067319905646f, 9.79369240592005f}},
		{"bunny", pmp::Point{-0.01684039831161499f, 0.11015420407056808f, 0.0012007840834242693f} },
		{"maxPlanck", pmp::Point{30.59686279296875f, -18.105804443359375f, 82.29149055480957f} },
		{"nefertiti", pmp::Point{0.0f, 0.0f, 0.0f} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"armadillo", pmp::vec3{-0.03070969905335075f, 0.12876712096541565f, 0.9911992448253433f}},
		{"bunny", pmp::vec3{0.0f, 0.0f, 1.0f} },
		{"maxPlanck", pmp::vec3{1.0f, 0.0f, 0.0f} },
		{"nefertiti", pmp::vec3{1.0f, 0.0f, 0.0f} }
	};

	const std::map<std::string, Circle2D> outerCircles{
		{"armadillo", Circle2D{pmp::Point2{0.372234f, 16.6515f}, 121.558f} },
		{"bunny", Circle2D{pmp::Point2{-0.0155906f, 0.102261f}, 0.142831f} },
		{"maxPlanck", Circle2D{pmp::Point2{-17.82f, 82.5006f}, 292.263f} },
		{"nefertiti", Circle2D{pmp::Point2{0.178497f, -0.0410004f}, 441.436f} }
	};
	const std::map<std::string, Circle2D> innerCircles{
		{"armadillo", Circle2D{pmp::Point2{-3.0f, 52.0f}, 20.0f}},
		{"bunny", Circle2D{pmp::Point2{-0.025f, 0.08f}, 0.025f}},
		{"maxPlanck", Circle2D{pmp::Point2{8.0f, 85.0f}, 50.0f}},
		{"nefertiti", Circle2D{pmp::Point2{-20.0f, 100.0f}, 55.0f}}
	};

	const std::map<std::string, Sphere3D> outerSpheres{
		{"armadillo", Sphere3D{pmp::Point{0.0122509f, 21.4183f, -0.000249863f}, 136.963f} },
		{"bunny", Sphere3D{pmp::Point{-0.0168297f, 0.110217f, -0.0015718f}, 0.141622f} },
		{"maxPlanck", Sphere3D{pmp::Point{30.658f, -17.9765f, 82.2885f}, 271.982f} },
		{"nefertiti", Sphere3D{pmp::Point{0.0144997f, -0.00499725f, -0.0215073f}, 392.184f} }
	};

	const std::map<std::string, double> criticalDistances{
		{"armadillo", 5.0},
		{"bunny", 0.005},
		{"maxPlanck", 15.0},
		{"nefertiti", 20.0}
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;

	constexpr size_t samplingLevel = 3;
	constexpr size_t nSamplings = 10;
	constexpr size_t minVerts = 9; // Minimum number of vertices to sample

	constexpr unsigned int seed = 5000; // seed for the pt cloud sampling RNG

	//SetRemeshingAdjustmentTimeIndices({}); // no remeshing adjustment
	SetRemeshingAdjustmentTimeIndices({ /*3, 10,*/ 30 /*, 50 , 100, 120, 140, 145*/ });

	constexpr unsigned int NTimeSteps = 180;

	constexpr bool executeCustomCurveEvolver = true;
	constexpr bool executeCustomSurfaceEvolver = false;

	for (const auto& meshName : meshForPtCloudNames)
	{
		// =======================================================================
		//   - - - - - - - - - - - - -   Data   Prep   - - - - - - - - - - - -
		// -----------------------------------------------------------------------

		std::cout << "==================================================================\n";
		std::cout << "Mesh To Pt Cloud: " << meshName << ".obj -> " << meshName << "Pts_" << samplingLevel << ".ply\n";
		std::cout << "------------------------------------------------------------------\n";
		const auto baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", false);
		if (!baseDataOpt.has_value())
		{
			std::cerr << "baseDataOpt == nullopt!\n";
			break;
		}
		std::cout << "meshName.obj" << " imported as BaseMeshGeometryData.\n";
		const auto& baseData = baseDataOpt.value();
		const size_t maxVerts = baseData.Vertices.size(); // Maximum number of vertices available
		size_t nVerts = minVerts + (maxVerts - minVerts) * samplingLevel / (nSamplings - 1);
		nVerts = std::max(minVerts, std::min(nVerts, maxVerts));

		std::cout << "Sampling " << nVerts << "/" << maxVerts << " vertices...\n";

		std::string filename = dataOutPath + meshName + "Pts_" + std::to_string(samplingLevel) + ".ply";
		if (!ExportSampledVerticesToPLY(baseData, nVerts, filename, seed))
		{
			std::cerr << "ExportSampledVerticesToPLY failed!\n";
			break;
		}

		const auto ptCloudName = meshName + "Pts_" + std::to_string(samplingLevel);
		const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataOutPath + ptCloudName + ".ply", true);
		if (!ptCloudOpt.has_value())
		{
			std::cerr << "ptCloudOpt == nullopt!\n";
			break;
		}

		const auto& ptCloud = ptCloudOpt.value();
		const pmp::BoundingBox ptCloudBBox(ptCloud);
		const auto center = ptCloudBBox.center();
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float cellSize = minSize / nVoxelsPerMinDimension;

		const double isoLvlOffsetFactor = (timeStepSizesForPtClouds.contains(ptCloudName) ? isoLevelOffsetFactors.at(ptCloudName) : defaultOffsetFactor);
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		const double tau = (timeStepSizesForPtClouds.contains(ptCloudName) ? timeStepSizesForPtClouds.at(ptCloudName) : defaultTimeStep); // time step

		// ==========================================================================
		// - - - - - - - - -  New Manifold Evolver (Curve)  - - - - - - - - - - - - 
		// ==========================================================================

		const float distTolerance = 0.01f * minSize;
		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0f, 0.0f, 0.0f });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}

		if (!Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
		{
			std::cerr << "Export2DPointCloudToPLY: internal error during export!\n";
			continue;
		}

		if (!criticalDistances.contains(meshName))
		{
			std::cerr << "!criticalDistances.contains(\"" << meshName << "\") ... skipping.\n";
			continue;
		}
		const double criticalDistance = criticalDistances.at(meshName);

		// -----------------------
		if (executeCustomSurfaceEvolver)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldRepulsion = [&criticalDistance](double distance)
			{
				return 1.0 * (1.0 / (criticalDistance + 0.2 * criticalDistance) - 1.0 / (distance + 0.2 * criticalDistance));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "newEvol_Pts" + std::to_string(samplingLevel);
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerSpheres.contains(meshName))
			{
				std::cerr << "!outerSpheres.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			const auto outerSphere = outerSpheres.at(meshName);
			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;

			// construct outer ico-sphere
			Geometry::IcoSphereBuilder icoBuilder({ strategySettings.LevelOfDetail, outerSphere.Radius });
			icoBuilder.BuildBaseData();
			icoBuilder.BuildPMPSurfaceMesh();
			auto outerSurface = icoBuilder.GetPMPSurfaceMeshResult();
			const pmp::mat4 transfMatrixGeomMove{
				1.0f, 0.0f, 0.0f, outerSphere.Center[0],
				0.0f, 1.0f, 0.0f, outerSphere.Center[1],
				0.0f, 0.0f, 1.0f, outerSphere.Center[2],
				0.0f, 0.0f, 0.0f, 1.0f
			};
			outerSurface *= transfMatrixGeomMove;

			std::vector<pmp::SurfaceMesh> innerSurfaces;
			auto surfaceStrategy = std::make_shared<CustomManifoldSurfaceEvolutionStrategy>(
				strategySettings, MeshLaplacian::Voronoi,
				outerSurface, innerSurfaces,
				std::make_shared<std::vector<pmp::Point>>(ptCloud));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(surfaceStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}

		if (executeCustomCurveEvolver)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.OuterManifoldRepulsion = [&criticalDistance](double distance)
			{
				return 0.001 * (1.0 / (criticalDistance + 0.5 * criticalDistance) - 1.0 / (distance + 0.5 * criticalDistance));
			};
			strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				return 0.5 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 0.5 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.InnerManifoldRepulsion = [&criticalDistance](double distance)
			{
				return 0.001 * (1.0 / (criticalDistance + 0.5 * criticalDistance) - 1.0 / (distance + 0.5 * criticalDistance));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "newEvol_Pts2D" + std::to_string(samplingLevel) + "_WithRepulsion";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerCircles.contains(meshName))
			{
				std::cerr << "!outerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}
			if (!innerCircles.contains(meshName))
			{
				std::cerr << "!innerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			const auto outerCircle = outerCircles.at(meshName);
			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);

			const auto innerCircle = innerCircles.at(meshName);
			const auto nInnerSegments = static_cast<unsigned int>(static_cast<pmp::Scalar>(nSegments) * innerCircle.Radius / outerCircle.Radius * 2);
			auto innerCurve = pmp::CurveFactory::circle(innerCircle.Center, innerCircle.Radius, nInnerSegments);
			innerCurve.negate_orientation();
			std::vector<pmp::ManifoldCurve2D> innerCurves{ innerCurve };

			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, outerCurve, innerCurves,
				std::make_shared<std::vector<pmp::Point2>>(pts2D));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	}
}


void OutwardEvolvingInnerCircleTest()
{
	const Circle2D innerTestCircle{ pmp::Point2{-3.0f, 52.0f}, 20.0f };

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;
	constexpr unsigned int NTimeSteps = 180;

	const pmp::BoundingBox2 bbox{
			pmp::Point2{
				innerTestCircle.Center[0] - innerTestCircle.Radius,
				innerTestCircle.Center[1] - innerTestCircle.Radius
			},
			pmp::Point2{
				innerTestCircle.Center[0] + innerTestCircle.Radius,
				innerTestCircle.Center[1] + innerTestCircle.Radius
			}
	};
	const auto bboxSize = bbox.max() - bbox.min();
	const float minSize = std::min(bboxSize[0], bboxSize[1]);
	const float maxSize = std::max(bboxSize[0], bboxSize[1]);
	const float cellSize = minSize / nVoxelsPerMinDimension;

	const double isoLvlOffsetFactor = defaultOffsetFactor;
	const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

	{
		std::cout << "Setting up ManifoldEvolutionSettings.\n";

		ManifoldEvolutionSettings strategySettings;
		strategySettings.UseInnerManifolds = true;
		strategySettings.InnerManifoldEpsilon = [](double distance)
		{
			return -0.05 * TRIVIAL_EPSILON(distance);
		};
		strategySettings.TimeStep = defaultTimeStep;
		strategySettings.LevelOfDetail = 3;
		strategySettings.TangentialVelocityWeight = 0.05;

		strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.22f;
		strategySettings.RemeshingSettings.UseBackProjection = false;

		strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
		strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

		strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
		strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

		std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

		GlobalManifoldEvolutionSettings globalSettings;
		globalSettings.NSteps = NTimeSteps;
		globalSettings.DoRemeshing = true;
		globalSettings.DetectFeatures = false;
		globalSettings.ExportPerTimeStep = true;
		globalSettings.ExportTargetDistanceFieldAsImage = true;
		globalSettings.ProcedureName = "singleInnerCircleTestPhaseTwo";
		globalSettings.OutputPath = dataOutPath;
		globalSettings.ExportResult = false;

		globalSettings.RemeshingResizeFactor = 0.7f;
		globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

		const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
		auto innerCurve = pmp::CurveFactory::circle(innerTestCircle.Center, innerTestCircle.Radius, nSegments);
		innerCurve.negate_orientation();
		std::vector<pmp::ManifoldCurve2D> innerCurves{ innerCurve };

		auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
			strategySettings, std::nullopt, innerCurves, nullptr);

		std::cout << "Setting up ManifoldEvolver.\n";

		ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

		std::cout << "ManifoldEvolver::Evolve ... ";

		try
		{
			evolver.Evolve();
		}
		catch (std::invalid_argument& ex)
		{
			std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
		}
		catch (std::runtime_error& ex)
		{
			std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
		}
		catch (...)
		{
			std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
		}
	}
}


void AdditionalCurveShapesTests()
{
	//
	// ==============  pmp::CurveFactory::rectangle =============
	//
	// Test case 1: Full rectangle without chamfering
	pmp::ManifoldCurve2D rectCurve1 = pmp::CurveFactory::rectangle(pmp::Point2{ 0, 0 }, 10, 5, 10, false);
	if (!pmp::write_to_ply(rectCurve1, dataOutPath + "rectangleCurve1_noChamfer.ply"))
		std::cerr << "Error writing rectangleCurve1_noChamfer.ply!\n";

	// Test case 2: Full rectangle with chamfering
	pmp::ManifoldCurve2D rectCurve2 = pmp::CurveFactory::rectangle(pmp::Point2{ 0, 0 }, 10, 5, 10, true);
	if (!pmp::write_to_ply(rectCurve2, dataOutPath + "rectangleCurve2_withChamfer.ply"))
		std::cerr << "Error writing rectangleCurve2_withChamfer.ply!\n";

	// Test case 3: Small rectangle without chamfering
	pmp::ManifoldCurve2D rectCurve3 = pmp::CurveFactory::rectangle(pmp::Point2{ 0, 0 }, 4, 2, 10, false);
	if (!pmp::write_to_ply(rectCurve3, dataOutPath + "rectangleCurve3_noChamfer_small.ply"))
		std::cerr << "Error writing rectangleCurve3_noChamfer_small.ply!\n";

	// Test case 4: Small rectangle with chamfering
	pmp::ManifoldCurve2D rectCurve4 = pmp::CurveFactory::rectangle(pmp::Point2{ 0, 0 }, 4, 2, 10, true);
	if (!pmp::write_to_ply(rectCurve4, dataOutPath + "rectangleCurve4_withChamfer_small.ply"))
		std::cerr << "Error writing rectangleCurve4_withChamfer_small.ply!\n";

	//
	// ==============  pmp::CurveFactory::sampled_polygon =============
	//
	// Test case 1: equilateral triangle without chamfering
	pmp::ManifoldCurve2D polyCurve1 = pmp::CurveFactory::sampled_polygon(pmp::BASE_EQUILATERAL_TRIANGLE_VERTS, 20, false);
	if (!pmp::write_to_ply(polyCurve1, dataOutPath + "triangle_noChamfer.ply"))
		std::cerr << "Error writing triangle_noChamfer.ply!\n";

	// Test case 2: equilateral triangle with chamfering
	pmp::ManifoldCurve2D polyCurve2 = pmp::CurveFactory::sampled_polygon(pmp::BASE_EQUILATERAL_TRIANGLE_VERTS, 20, true);
	if (!pmp::write_to_ply(polyCurve2, dataOutPath + "triangle_withChamfer.ply"))
		std::cerr << "Error writing triangle_withChamfer.ply!\n";

	const std::vector squareVertices = {
		pmp::Point2{0.0f, 0.0f},
		pmp::Point2{1.0f, 0.0f},
		pmp::Point2{1.0f, 1.0f},
		pmp::Point2{0.0f, 1.0f}
	};
	// Test case 3: unit square without chamfering
	pmp::ManifoldCurve2D polyCurve3 = pmp::CurveFactory::sampled_polygon(squareVertices, 20, false);
	if (!pmp::write_to_ply(polyCurve3, dataOutPath + "square_noChamfer.ply"))
		std::cerr << "Error writing square_noChamfer.ply!\n";

	// Test case 4: unit square with chamfering
	pmp::ManifoldCurve2D polyCurve4 = pmp::CurveFactory::sampled_polygon(squareVertices, 20, true);
	if (!pmp::write_to_ply(polyCurve4, dataOutPath + "square_withChamfer.ply"))
		std::cerr << "Error writing square_withChamfer.ply!\n";

	const std::vector nonClosedPolyline = {
		pmp::Point2{0.0f, 0.0f},
		pmp::Point2{1.0f, 0.0f},
		pmp::Point2{1.5f, 0.5f},
		pmp::Point2{2.0f, 0.0f}
	};
	// Test case 5: non-closed polyline without chamfering
	pmp::ManifoldCurve2D polyCurve5 = pmp::CurveFactory::sampled_polygon(nonClosedPolyline, 20, false, false);
	if (!pmp::write_to_ply(polyCurve5, dataOutPath + "nonClosedPolyline_noChamfer.ply"))
		std::cerr << "Error writing nonClosedPolyline_noChamfer.ply!\n";
	// Test case 6: non-closed polyline with chamfering
	pmp::ManifoldCurve2D polyCurve6 = pmp::CurveFactory::sampled_polygon(nonClosedPolyline, 20, true, false);
	if (!pmp::write_to_ply(polyCurve6, dataOutPath + "nonClosedPolyline_withChamfer.ply"))
		std::cerr << "Error writing nonClosedPolyline_withChamfer.ply!\n";

	const std::vector hexagonVertices = {
		pmp::Point2{0.0f, 1.0f},
		pmp::Point2{0.866f, 0.5f},
		pmp::Point2{0.866f, -0.5f},
		pmp::Point2{0.0f, -1.0f},
		pmp::Point2{-0.866f, -0.5f},
		pmp::Point2{-0.866f, 0.5f}
	};
	// Test case 7: hexagon without chamfering
	pmp::ManifoldCurve2D polyCurve7 = pmp::CurveFactory::sampled_polygon(hexagonVertices, 30, false, true);
	if (!pmp::write_to_ply(polyCurve7, dataOutPath + "hexagon_noChamfer.ply"))
		std::cerr << "Error writing hexagon_noChamfer.ply!\n";
	// Test case 8: hexagon with chamfering
	pmp::ManifoldCurve2D polyCurve8 = pmp::CurveFactory::sampled_polygon(hexagonVertices, 30, true, true);
	if (!pmp::write_to_ply(polyCurve8, dataOutPath + "hexagon_withChamfer.ply"))
		std::cerr << "Error writing hexagon_withChamfer.ply!\n";
}


void AdvectionDrivenInnerCircleTests()
{
	const std::map<std::string, pmp::ManifoldCurve2D> targetCurves{
		{ "circle", pmp::CurveFactory::circle(pmp::Point2{-3.0f, 52.0f}, 35.0f, 25, 0.0, 2.0 * M_PI)},
		{ "incompleteCircle", pmp::CurveFactory::circle(pmp::Point2{-3.0f, 52.0f}, 35.0f, 25, M_PI_2, 2.0 * M_PI)},
		{ "sineDeformedCircle", pmp::CurveFactory::sine_deformed_circle(pmp::Point2{-3.0f, 52.0f}, 35.0f, 25, 7.0f, 4.0f, 0.0, 2.0 * M_PI)},
		{ "sineDeformedIncompleteCircle", pmp::CurveFactory::sine_deformed_circle(pmp::Point2{-3.0f, 52.0f}, 35.0f, 25, 7.0f, 4.0f, M_PI_2, 2.0 * M_PI)},
		{ "chamferedRectangle", pmp::CurveFactory::rectangle(pmp::Point2{-3.0f, 52.0f}, 60.0f, 70.0f, 15, true)},
		{ "incompleteChamferedRectangle", pmp::CurveFactory::sampled_polygon({
			pmp::Point2{-30.0f, -35.0f} + pmp::Point2{-3.0f, 52.0f},
			pmp::Point2{30.0f, -35.0f} + pmp::Point2{-3.0f, 52.0f},
			pmp::Point2{30.0f, 35.0f} + pmp::Point2{-3.0f, 52.0f},
			pmp::Point2{-30.0f, 35.0f} + pmp::Point2{-3.0f, 52.0f}}, 30, true, false)},
		{ "chamferedTriangle", pmp::CurveFactory::sampled_polygon({
			pmp::Point2{-0.5f, -sqrtf(3.0f) / 6.0f} *120.0f + pmp::Point2{-3.0f, 52.0f},
			pmp::Point2{0.5f, -sqrtf(3.0f) / 6.0f} *120.0f + pmp::Point2{-3.0f, 52.0f},
			pmp::Point2{0.0f, sqrtf(3.0f) / 3.0f} *120.0f + pmp::Point2{-3.0f, 52.0f}}, 30, true)},
		{ "incompleteChamferedTriangle", pmp::CurveFactory::sampled_polygon({
			pmp::Point2{-0.5f, -sqrtf(3.0f) / 6.0f} *120.0f + pmp::Point2{-3.0f, 52.0f},
			pmp::Point2{0.5f, -sqrtf(3.0f) / 6.0f} *120.0f + pmp::Point2{-3.0f, 52.0f},
			pmp::Point2{0.0f, sqrtf(3.0f) / 3.0f} *120.0f + pmp::Point2{-3.0f, 52.0f}}, 30, true, false)}
	};
	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;
	constexpr unsigned int NTimeSteps = 180;
	const auto innerTestCircle = Circle2D{ pmp::Point2{-3.0f, 52.0f}, 20.0f };

	for (const auto& [ptCloudName, curve] : targetCurves)
	{
		std::cout << "========================================================\n";
		std::cout << "         inner circle LSW for: " << ptCloudName << " ... \n";
		std::cout << " ------------------------------------------------------ \n";
		const auto& pts2D = curve.positions();
		const pmp::BoundingBox2 bbox{ pts2D };

		const auto bboxSize = bbox.max() - bbox.min();
		const float minSize = std::min(bboxSize[0], bboxSize[1]);
		//const float maxSize = std::max(bboxSize[0], bboxSize[1]);
		const float cellSize = minSize / nVoxelsPerMinDimension;

		const double isoLvlOffsetFactor = defaultOffsetFactor;
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				return 0.001 * TRIVIAL_EPSILON(distance);
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 1.5 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = defaultTimeStep;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.22f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = "innerCircle_" + ptCloudName + "_PtsTest";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto innerCurve = pmp::CurveFactory::circle(innerTestCircle.Center, innerTestCircle.Radius, nSegments);
			innerCurve.negate_orientation();
			std::vector<pmp::ManifoldCurve2D> innerCurves{ innerCurve };

			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, std::nullopt, innerCurves,
				std::make_shared<std::vector<pmp::Point2>>(pts2D));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	}
}


void ConcentricCirclesTests()
{
	// Define the inner and outer circle pairs directly
	const std::vector<std::pair<Circle2D, Circle2D>> circlePairs{
		{Circle2D{pmp::Point2{-3.0f, 52.0f}, 100.0f}, Circle2D{pmp::Point2{-3.0f, 52.0f}, 121.558f}},
		{Circle2D{pmp::Point2{-25.0f, 8.0f}, 0.055f}, Circle2D{pmp::Point2{-25.0f, 8.0f}, 0.142831f}},
		{Circle2D{pmp::Point2{8.0f, 85.0f}, 50.0f}, Circle2D{pmp::Point2{8.0f, 85.0f}, 292.263f}},
		{Circle2D{pmp::Point2{-20.0f, 90.0f}, 55.0f}, Circle2D{pmp::Point2{-20.0f, 90.0f}, 441.436f}}
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 0.25; // 1.5;
	constexpr unsigned int NTimeSteps = 180;

	for (unsigned int curveId = 0; const auto & circlePair : circlePairs)
	{
		std::cout << " ====================================================== \n";
		std::cout << "Circles " << std::to_string(curveId) << ": \n";
		std::cout << " ------------------------------------------------------ \n";

		// Retrieve the inner and outer circles from the pair
		const auto& innerCircle = circlePair.first;
		const auto& outerCircle = circlePair.second;

		const pmp::BoundingBox2 bbox{
			pmp::Point2{
				std::min(innerCircle.Center[0] - innerCircle.Radius, outerCircle.Center[0] - outerCircle.Radius),
				std::min(innerCircle.Center[1] - innerCircle.Radius, outerCircle.Center[1] - outerCircle.Radius)
			},
			pmp::Point2{
				std::max(innerCircle.Center[0] + innerCircle.Radius, outerCircle.Center[0] + outerCircle.Radius),
				std::max(innerCircle.Center[1] + innerCircle.Radius, outerCircle.Center[1] + outerCircle.Radius)
			}
		};
		const auto bboxSize = bbox.max() - bbox.min();
		const float minSize = std::min(bboxSize[0], bboxSize[1]);
		const float cellSize = minSize / nVoxelsPerMinDimension;

		//if (executeCustomCurveEvolver)
		{
			const double isoLvlOffsetFactor = defaultOffsetFactor;
			const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.AdvectionInteractWithOtherManifolds = true;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				return 0.001 * TRIVIAL_EPSILON(distance);
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 2.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = defaultTimeStep;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			strategySettings.ExportVariableScalarFieldsDimInfo = true;
			strategySettings.ExportVariableVectorFieldsDimInfo = true;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = "concentricCircles" + std::to_string(curveId) + "_Repulsionless";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);

			const auto nInnerSegments = static_cast<unsigned int>(static_cast<pmp::Scalar>(nSegments) * innerCircle.Radius / outerCircle.Radius * 2);
			auto innerCurve = pmp::CurveFactory::circle(innerCircle.Center, innerCircle.Radius, nInnerSegments);
			innerCurve.negate_orientation();
			std::vector innerCurves{ innerCurve };

			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, outerCurve, innerCurves, nullptr);

			std::cout << "Setting up ManifoldEvolver.\n";


			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}

		curveId++;
	}
}


void NonConcentricCirclesTest()
{
	// Define the inner and outer circle pairs directly
	const std::vector<std::pair<Circle2D, Circle2D>> circlePairs{
		{Circle2D{pmp::Point2{-3.0f, 52.0f}, 20.0f}, Circle2D{pmp::Point2{0.372234f, 16.6515f}, 121.558f}},
		{Circle2D{pmp::Point2{-0.025f, 0.08f}, 0.025f}, Circle2D{pmp::Point2{-0.0155906f, 0.102261f}, 0.142831f}},
		{Circle2D{pmp::Point2{8.0f, 85.0f}, 50.0f}, Circle2D{pmp::Point2{-17.82f, 82.5006f}, 292.263f}},
		{Circle2D{pmp::Point2{-20.0f, 90.0f}, 55.0f}, Circle2D{pmp::Point2{0.178497f, -0.0410004f}, 441.436f}}
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 0.25; //1.5;
	constexpr unsigned int NTimeSteps = 180;

	for (unsigned int curveId = 0; const auto & circlePair : circlePairs)
	{
		std::cout << " ====================================================== \n";
		std::cout << "Circles " << std::to_string(curveId) << ": \n";
		std::cout << " ------------------------------------------------------ \n";

		// Retrieve the inner and outer circles from the pair
		const auto& innerCircle = circlePair.first;
		const auto& outerCircle = circlePair.second;

		const pmp::BoundingBox2 bbox{
			pmp::Point2{
				std::min(innerCircle.Center[0] - innerCircle.Radius, outerCircle.Center[0] - outerCircle.Radius),
				std::min(innerCircle.Center[1] - innerCircle.Radius, outerCircle.Center[1] - outerCircle.Radius)
			},
			pmp::Point2{
				std::max(innerCircle.Center[0] + innerCircle.Radius, outerCircle.Center[0] + outerCircle.Radius),
				std::max(innerCircle.Center[1] + innerCircle.Radius, outerCircle.Center[1] + outerCircle.Radius)
			}
		};
		const auto bboxSize = bbox.max() - bbox.min();
		const float minSize = std::min(bboxSize[0], bboxSize[1]);
		//const float maxSize = std::max(bboxSize[0], bboxSize[1]);
		const float cellSize = minSize / nVoxelsPerMinDimension;

		const double isoLvlOffsetFactor = defaultOffsetFactor;
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		//if (executeCustomCurveEvolver)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.AdvectionInteractWithOtherManifolds = true;
			// use PreComputeAdvectionDiffusionParams?
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				return 0.001 * TRIVIAL_EPSILON(distance);
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 2.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = defaultTimeStep;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "fieldIsoLevel: " << fieldIsoLevel << "\n";

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = "curve" + std::to_string(curveId) + "_Repulsionless";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);

			const auto nInnerSegments = static_cast<unsigned int>(static_cast<pmp::Scalar>(nSegments) * innerCircle.Radius / outerCircle.Radius * 2);
			auto innerCurve = pmp::CurveFactory::circle(innerCircle.Center, innerCircle.Radius, nInnerSegments);
			innerCurve.negate_orientation();
			std::vector innerCurves{ innerCurve };

			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, outerCurve, innerCurves, nullptr);

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}

		curveId++;
	}
}

void EquilibriumPairedManifoldTests()
{
	const auto center = pmp::Point2(200, 400);
	const float outerRadius = 100.0f;
	const float innerRadius = 40.0f;

	const std::vector<pmp::Point2> squareVerticesLarge = {
		pmp::Point2{-50.0f, -50.0f} + center,
		pmp::Point2{50.0f, -50.0f} + center,
		pmp::Point2{50.0f, 50.0f} + center,
		pmp::Point2{-50.0f, 50.0f} + center
	};
	const std::vector<pmp::Point2> triangleVerticesLarge = {
		pmp::Point2{-0.5f, -sqrtf(3.0f) / 6.0f} * (2.0f * outerRadius) + center,
		pmp::Point2{0.5f, -sqrtf(3.0f) / 6.0f} * (2.0f * outerRadius) + center,
		pmp::Point2{0.0f, sqrtf(3.0f) / 3.0f} * (2.0f * outerRadius) + center
	};
	const std::vector<pmp::Point2> squareVerticesSmall = {
		pmp::Point2{-20.0f, -20.0f} + center,
		pmp::Point2{20.0f, -20.0f} + center,
		pmp::Point2{20.0f, 20.0f} + center,
		pmp::Point2{-20.0f, 20.0f} + center
	};
	const std::vector<pmp::Point2> triangleVerticesSmall = {
		pmp::Point2{-0.5f, -sqrtf(3.0f) / 6.0f} * innerRadius + center,
		pmp::Point2{0.5f, -sqrtf(3.0f) / 6.0f} * innerRadius + center,
		pmp::Point2{0.0f, sqrtf(3.0f) / 3.0f} * innerRadius + center
	};
	const unsigned int segments = 30;

	// List of curve pairs to evolve
	const std::vector<std::pair<pmp::ManifoldCurve2D, pmp::ManifoldCurve2D>> curvePairs{
		// Pair 1: Outer chamfered square and inner circle
		{pmp::CurveFactory::sampled_polygon(squareVerticesLarge, segments, true),
		 pmp::CurveFactory::circle(center, 0.5f * innerRadius, segments)},

		// Pair 2: Outer chamfered equilateral triangle and inner circle
		{pmp::CurveFactory::sampled_polygon(triangleVerticesLarge, segments, true),
		pmp::CurveFactory::circle(center, innerRadius, segments)},

		// Pair 3: Outer circle and inner chamfered square
		{pmp::CurveFactory::circle(center, outerRadius, segments),
		pmp::CurveFactory::sampled_polygon(squareVerticesSmall, segments, true)},

		// Pair 4: Outer circle and inner chamfered equilateral triangle
		{pmp::CurveFactory::circle(center, outerRadius, segments),
		pmp::CurveFactory::sampled_polygon(triangleVerticesSmall, segments, true)}
	};

	// Prepare the settings for the evolver
	ManifoldEvolutionSettings strategySettings;
	strategySettings.UseInnerManifolds = true;

	strategySettings.AdvectionInteractWithOtherManifolds = true;
	// use PreComputeAdvectionDiffusionParams?
	strategySettings.OuterManifoldEpsilon = [](double distance)
	{
		return 1.0 * (1.0 - exp(-distance * distance / 1.0));
	};
	strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
	{
		return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
	};
	strategySettings.InnerManifoldEpsilon = [](double distance)
	{
		return 0.001 * TRIVIAL_EPSILON(distance);
	};
	strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
	{
		return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
	};
	strategySettings.LevelOfDetail = 3;
	strategySettings.TangentialVelocityWeight = 0.05;

	strategySettings.TimeStep = 0.05;
	strategySettings.TangentialVelocityWeight = 0.05;
	strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14f;
	strategySettings.RemeshingSettings.UseBackProjection = false;
	strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
	strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = static_cast<float>(M_PI_2);
	strategySettings.FieldSettings.NVoxelsPerMinDimension = 40;
	strategySettings.FieldSettings.FieldIsoLevel = 2.0;

	strategySettings.ExportVariableScalarFieldsDimInfo = true;
	strategySettings.ExportVariableVectorFieldsDimInfo = true;

	// Global settings
	GlobalManifoldEvolutionSettings globalSettings;
	globalSettings.NSteps = 180;
	globalSettings.DoRemeshing = true;
	globalSettings.DetectFeatures = false;
	globalSettings.ExportPerTimeStep = true;
	globalSettings.ExportTargetDistanceFieldAsImage = true;
	globalSettings.OutputPath = dataOutPath;
	globalSettings.ExportResult = false;

	globalSettings.RemeshingResizeFactor = 0.7f;
	globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

	for (unsigned int pairId = 0; const auto & [outerCurve, innerCurve] : curvePairs)
	{
		std::cout << "EquilibriumPairedManifoldTest for pair: " << pairId << "\n";

		globalSettings.ProcedureName = "equilibriumPair" + std::to_string(pairId);

		std::vector innerCurves{ innerCurve };
		innerCurves[0].negate_orientation();
		auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
			strategySettings, outerCurve, innerCurves, nullptr);

		ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

		try
		{
			evolver.Evolve();
		}
		catch (std::invalid_argument& ex)
		{
			std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
		}
		catch (std::runtime_error& ex)
		{
			std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
		}
		catch (...)
		{
			std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
		}

		pairId++;
	}
}

void VisualizeMCF()
{
	std::vector curves{
		pmp::CurveFactory::circle(pmp::Point2(0, 0), 10.0, 40),
		pmp::CurveFactory::sine_deformed_circle(pmp::Point2(0, 0), 10.0, 40, 10.0, 4)
	};

	// remesh the deformed circle
	constexpr float edgeLength = 2.0f;
	constexpr unsigned int iterations = 10;
	pmp::CurveRemeshing remesher(curves[1]);
	pmp::AdaptiveRemeshingSettings settings;
	settings.MinEdgeLength = edgeLength;
	settings.MaxEdgeLength = 1.5f * edgeLength;
	settings.ApproxError = 0.05f * edgeLength;
	settings.NRemeshingIterations = iterations;
	settings.NTangentialSmoothingIters = 6;
	settings.UseProjection = true;
	remesher.adaptive_remeshing(settings);

	constexpr unsigned int NTimeSteps = 30;

	for (unsigned int curveId = 0; const auto& curve : curves)
	{
		std::cout << "Setting up ManifoldEvolutionSettings.\n";

		ManifoldEvolutionSettings strategySettings;
		strategySettings.UseInnerManifolds = true;
		strategySettings.OuterManifoldEpsilon = [](double distance)
		{
			return 0.02 * TRIVIAL_EPSILON(distance);
		};
		strategySettings.TimeStep = 0.01;
		strategySettings.LevelOfDetail = 4;
		strategySettings.TangentialVelocityWeight = 0.05;

		strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.22f;
		strategySettings.RemeshingSettings.UseBackProjection = false;

		strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
		strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

		strategySettings.FieldSettings.NVoxelsPerMinDimension = 40;
		strategySettings.FieldSettings.FieldIsoLevel = 0.01;

		std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

		GlobalManifoldEvolutionSettings globalSettings;
		globalSettings.NSteps = NTimeSteps;
		globalSettings.DoRemeshing = true;
		globalSettings.DetectFeatures = false;
		globalSettings.ExportPerTimeStep = true;
		globalSettings.ExportTargetDistanceFieldAsImage = true;
		globalSettings.ProcedureName = "mcfCurve" + std::to_string(curveId);
		globalSettings.OutputPath = dataOutPath;
		globalSettings.ExportResult = false;

		globalSettings.RemeshingResizeFactor = 0.7f;
		globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

		std::vector<pmp::ManifoldCurve2D> innerCurves{};

		auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
			strategySettings, curve, innerCurves, nullptr);

		std::cout << "Setting up ManifoldEvolver.\n";

		ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

		std::cout << "ManifoldEvolver::Evolve ... ";

		try
		{
			evolver.Evolve();
		}
		catch (std::invalid_argument& ex)
		{
			std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
		}
		catch (std::runtime_error& ex)
		{
			std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
		}
		catch (...)
		{
			std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
		}

		curveId++;
	}
}


void VisualizeMultipleInnerCurves()
{
	const std::vector<std::string> meshForPtCloudNames{
		//"armadillo",
		//"bunny",
		//"maxPlanck",
		"nefertiti"
	};
	const std::map<std::string, double> timeStepSizesForPtClouds{
		{"armadillo", 0.05 },
		{ "bunny", 0.05 },
		{ "maxPlanck", 0.05 },
		{ "nefertiti", 0.05 },
	};
	const std::map<std::string, double> isoLevelOffsetFactors{
		{"armadillo", 0.5 },
		{ "bunny", 0.5 },
		{ "maxPlanck", 0.5 },
		{ "nefertiti", 0.5 },
	};

	const std::map<std::string, pmp::Point> slicingPlaneRefPts{
	{"armadillo", pmp::Point{-0.10348621158928151f, 21.427067319905646f, 9.79369240592005f}},
	{"bunny", pmp::Point{-0.01684039831161499f, 0.11015420407056808f, 0.0012007840834242693f} },
	{"maxPlanck", pmp::Point{30.59686279296875f, -18.105804443359375f, 82.29149055480957f} },
	{"nefertiti", pmp::Point{0.0f, 0.0f, 0.0f} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"armadillo", pmp::vec3{-0.03070969905335075f, 0.12876712096541565f, 0.9911992448253433f}},
		{"bunny", pmp::vec3{0.0f, 0.0f, 1.0f} },
		{"maxPlanck", pmp::vec3{1.0f, 0.0f, 0.0f} },
		{"nefertiti", pmp::vec3{1.0f, 0.0f, 0.0f} }
	};

	const std::map<std::string, Circle2D> outerCircles{
		//{"armadillo", Circle2D{pmp::Point2{0.372234f, 16.6515f}, 121.558f} },
		//{"bunny", Circle2D{pmp::Point2{-0.0155906f, 0.102261f}, 0.142831f} },
		//{"maxPlanck", Circle2D{pmp::Point2{-17.82f, 82.5006f}, 292.263f} },
		{"nefertiti", Circle2D{pmp::Point2{0.178497f, -0.0410004f}, 441.436f} }
	};
	const std::map<std::string, std::vector<Circle2D>> innerCircles{
		//{"armadillo", { Circle2D{pmp::Point2{-3.0f, 52.0f}, 20.0f} }},
		//{"bunny", { Circle2D{pmp::Point2{-0.025f, 0.08f}, 0.025f} } },
		//{"maxPlanck", { Circle2D{pmp::Point2{8.0f, 85.0f}, 50.0f} }},

		{"nefertiti", {
			Circle2D{pmp::Point2{-20.0f, 100.0f}, 55.0f},
			//Circle2D{pmp::Point2{-75.0f, -50.0f}, 25.0f}
			Circle2D{pmp::Point2{-10.0f, -200.0f}, 35.0f}
		}}
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;

	constexpr size_t samplingLevel = 3;
	constexpr size_t nSamplings = 10;
	constexpr size_t minVerts = 9; // Minimum number of vertices to sample

	constexpr unsigned int seed = 5000; // seed for the pt cloud sampling RNG

	//SetRemeshingAdjustmentTimeIndices({}); // no remeshing adjustment
	SetRemeshingAdjustmentTimeIndices({ /*3, 10,*/ 30 /*, 50 , 100, 120, 140, 145*/ });

	constexpr unsigned int NTimeSteps = 180;

	for (const auto& meshName : meshForPtCloudNames)
	{
		// =======================================================================
		//   - - - - - - - - - - - - -   Data   Prep   - - - - - - - - - - - -
		// -----------------------------------------------------------------------

		std::cout << "==================================================================\n";
		std::cout << "Mesh To Pt Cloud: " << meshName << ".obj -> " << meshName << "Pts_" << samplingLevel << ".ply\n";
		std::cout << "------------------------------------------------------------------\n";
		const auto baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", false);
		if (!baseDataOpt.has_value())
		{
			std::cerr << "baseDataOpt == nullopt!\n";
			break;
		}
		std::cout << "meshName.obj" << " imported as BaseMeshGeometryData.\n";
		const auto& baseData = baseDataOpt.value();
		const size_t maxVerts = baseData.Vertices.size(); // Maximum number of vertices available
		size_t nVerts = minVerts + (maxVerts - minVerts) * samplingLevel / (nSamplings - 1);
		nVerts = std::max(minVerts, std::min(nVerts, maxVerts));

		std::cout << "Sampling " << nVerts << "/" << maxVerts << " vertices...\n";

		std::string filename = dataOutPath + meshName + "Pts_" + std::to_string(samplingLevel) + ".ply";
		if (!ExportSampledVerticesToPLY(baseData, nVerts, filename, seed))
		{
			std::cerr << "ExportSampledVerticesToPLY failed!\n";
			break;
		}

		const auto ptCloudName = meshName + "Pts_" + std::to_string(samplingLevel);
		const auto ptCloudOpt = Geometry::ImportPLYPointCloudData(dataOutPath + ptCloudName + ".ply", true);
		if (!ptCloudOpt.has_value())
		{
			std::cerr << "ptCloudOpt == nullopt!\n";
			break;
		}

		const auto& ptCloud = ptCloudOpt.value();
		const pmp::BoundingBox ptCloudBBox(ptCloud);
		const auto center = ptCloudBBox.center();
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		//const float maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float cellSize = minSize / nVoxelsPerMinDimension;

		const double isoLvlOffsetFactor = (isoLevelOffsetFactors.contains(ptCloudName) ? isoLevelOffsetFactors.at(ptCloudName) : defaultOffsetFactor);
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		const double tau = (timeStepSizesForPtClouds.contains(ptCloudName) ? timeStepSizesForPtClouds.at(ptCloudName) : defaultTimeStep); // time step

		// ==========================================================================
		// - - - - - - - - -  New Manifold Evolver (Curve)  - - - - - - - - - - - - 
		// ==========================================================================

		const float distTolerance = 0.01f * minSize;
		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0f, 0.0f, 0.0f });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}

		if (!Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
		{
			std::cerr << "Export2DPointCloudToPLY: internal error during export!\n";
			continue;
		}

		// CustomCurveEvolver
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				return 0.0 * TRIVIAL_EPSILON(distance);
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 1.5 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 4;
			strategySettings.TangentialVelocityWeight = 0.05;

			//strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.22f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "newEvol_Pts2D" + std::to_string(samplingLevel) + "_RepulsionlessTwo";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerCircles.contains(meshName))
			{
				std::cerr << "!outerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}
			if (!innerCircles.contains(meshName))
			{
				std::cerr << "!innerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			const auto outerCircle = outerCircles.at(meshName);
			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);

			std::vector<pmp::ManifoldCurve2D> innerCurves{};
			innerCurves.reserve(innerCircles.at(meshName).size());

			for (const auto& innerCircle : innerCircles.at(meshName))
			{
				const auto nInnerSegments = static_cast<unsigned int>(static_cast<pmp::Scalar>(nSegments) * innerCircle.Radius / outerCircle.Radius * 2);
				auto innerCurve = pmp::CurveFactory::circle(innerCircle.Center, innerCircle.Radius, nInnerSegments);
				innerCurve.negate_orientation();
				innerCurves.push_back(innerCurve);
			}

			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, outerCurve, innerCurves,
				std::make_shared<std::vector<pmp::Point2>>(pts2D));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	}
}

void ExportSlicingPlanes()
{
	const std::vector<std::string> meshForPtCloudNames{
		"armadillo",
		"bunny",
		"maxPlanck",
		"nefertiti"
	};

	const std::map<std::string, pmp::Point> slicingPlaneRefPts{
		{"armadillo", pmp::Point{-0.10348621158928151f, 21.427067319905646f, 9.79369240592005f}},
		{"bunny", pmp::Point{-0.01684039831161499f, 0.11015420407056808f, 0.0012007840834242693f} },
		{"maxPlanck", pmp::Point{30.59686279296875f, -18.105804443359375f, 82.29149055480957f} },
		{"nefertiti", pmp::Point{0.0f, 0.0f, 0.0f} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"armadillo", pmp::vec3{-0.03070969905335075f, 0.12876712096541565f, 0.9911992448253433f}},
		{"bunny", pmp::vec3{0.0f, 0.0f, 1.0f} },
		{"maxPlanck", pmp::vec3{1.0f, 0.0f, 0.0f} },
		{"nefertiti", pmp::vec3{1.0f, 0.0f, 0.0f} }
	};

	for (const auto& meshName : meshForPtCloudNames)
	{
		Geometry::PlaneSettings planeSettings;
		planeSettings.Width = 200.0f;
		planeSettings.Depth = 250.0f;
		planeSettings.nWidthSegments = 1;
		planeSettings.nDepthSegments = 1;
		planeSettings.Origin = slicingPlaneRefPts.at(meshName) - 0.5f * pmp::vec3{ planeSettings.Width, planeSettings.Depth, 0.0f };

		Geometry::PlaneBuilder pb(planeSettings);
		pb.BuildBaseData();
		pb.BuildPMPSurfaceMesh();
		auto pMesh = pb.GetPMPSurfaceMeshResult();

		const pmp::vec3 defaultNormal(0.0f, 0.0f, 1.0f);
		const pmp::vec3 targetNormal = slicingPlaneNormals.at(meshName);

		const pmp::mat4 rotationMat = rotation_matrix(targetNormal, defaultNormal);
		pMesh *= rotationMat;

		const pmp::mat4 translationMat = translation_matrix(slicingPlaneRefPts.at(meshName));
		pMesh *= translationMat;

		pMesh.write(dataOutPath + meshName + "_SlicingPlane.vtk");
	}
}

void AdvectionDrivenInnerOuterCircleTests()
{
	const std::map<std::string, pmp::ManifoldCurve2D> targetCurves{
		{ "Circle", pmp::CurveFactory::circle(pmp::Point2{-3.0f, 52.0f}, 35.0f, 25, 0.0, 2.0 * M_PI)},
		{ "IncompleteCircle", pmp::CurveFactory::circle(pmp::Point2{-3.0f, 52.0f}, 35.0f, 25, M_PI_2, 2.0 * M_PI)},
		{ "SineDeformedCircle", pmp::CurveFactory::sine_deformed_circle(pmp::Point2{-3.0f, 52.0f}, 35.0f, 25, 7.0f, 4.0f, 0.0, 2.0 * M_PI)},
		{ "SineDeformedIncompleteCircle", pmp::CurveFactory::sine_deformed_circle(pmp::Point2{-3.0f, 52.0f}, 35.0f, 25, 7.0f, 4.0f, M_PI_2, 2.0 * M_PI)},
		{ "ChamferedRectangle", pmp::CurveFactory::rectangle(pmp::Point2{-3.0f, 52.0f}, 60.0f, 70.0f, 15, true)},
		{ "IncompleteChamferedRectangle", pmp::CurveFactory::sampled_polygon({
			pmp::Point2{-30.0f, -35.0f} + pmp::Point2{-3.0f, 52.0f},
			pmp::Point2{30.0f, -35.0f} + pmp::Point2{-3.0f, 52.0f},
			pmp::Point2{30.0f, 35.0f} + pmp::Point2{-3.0f, 52.0f},
			pmp::Point2{-30.0f, 35.0f} + pmp::Point2{-3.0f, 52.0f}}, 30, true, false)},
		{ "ChamferedTriangle", pmp::CurveFactory::sampled_polygon({
			pmp::Point2{-0.5f, -sqrtf(3.0f) / 6.0f} *120.0f + pmp::Point2{-3.0f, 52.0f},
			pmp::Point2{0.5f, -sqrtf(3.0f) / 6.0f} *120.0f + pmp::Point2{-3.0f, 52.0f},
			pmp::Point2{0.0f, sqrtf(3.0f) / 3.0f} *120.0f + pmp::Point2{-3.0f, 52.0f}}, 30, true)},
		{ "IncompleteChamferedTriangle", pmp::CurveFactory::sampled_polygon({
			pmp::Point2{-0.5f, -sqrtf(3.0f) / 6.0f} *120.0f + pmp::Point2{-3.0f, 52.0f},
			pmp::Point2{0.5f, -sqrtf(3.0f) / 6.0f} *120.0f + pmp::Point2{-3.0f, 52.0f},
			pmp::Point2{0.0f, sqrtf(3.0f) / 3.0f} *120.0f + pmp::Point2{-3.0f, 52.0f}}, 30, true, false)}
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;
	constexpr unsigned int NTimeSteps = 180;

	const auto innerCircle = Circle2D{ pmp::Point2{-3.0f, 52.0f}, 20.0f };
	const auto outerCircle = Circle2D{ pmp::Point2{-3.0f, 52.0f}, 110.0f };

	for (const auto& [ptCloudName, curve] : targetCurves)
	{
		std::cout << "========================================================\n";
		std::cout << "         inner/outer circle LSW for: " << ptCloudName << " ... \n";
		std::cout << " ------------------------------------------------------ \n";

		const auto& pts2D = curve.positions();
		const pmp::BoundingBox2 bbox{ pts2D };

		const auto bboxSize = bbox.max() - bbox.min();
		const float minSize = std::min(bboxSize[0], bboxSize[1]);
		//const float maxSize = std::max(bboxSize[0], bboxSize[1]);
		const float cellSize = minSize / nVoxelsPerMinDimension;

		const double isoLvlOffsetFactor = defaultOffsetFactor;
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.AdvectionInteractWithOtherManifolds = true;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 2.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				return 0.001 * TRIVIAL_EPSILON(distance);
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 2.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = defaultTimeStep;
			strategySettings.LevelOfDetail = 4;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = "innerOuter" + ptCloudName;
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);

			const auto nInnerSegments = static_cast<unsigned int>(static_cast<pmp::Scalar>(nSegments) * innerCircle.Radius / outerCircle.Radius * 2);
			auto innerCurve = pmp::CurveFactory::circle(innerCircle.Center, innerCircle.Radius, nInnerSegments);
			innerCurve.negate_orientation();
			std::vector innerCurves{ innerCurve };

			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, outerCurve, innerCurves,
				std::make_shared<std::vector<pmp::Point2>>(pts2D));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	}
}

void OuterOnlySimpleShapeTests()
{
	const std::map<std::string, pmp::ManifoldCurve2D> targetCurves{
	{ "Circle", pmp::CurveFactory::circle(pmp::Point2{-3.0f, 52.0f}, 35.0f, 25, 0.0, 2.0 * M_PI)},
	{ "IncompleteCircle", pmp::CurveFactory::circle(pmp::Point2{-3.0f, 52.0f}, 35.0f, 25, M_PI_2, 2.0 * M_PI)},
	{ "SineDeformedCircle", pmp::CurveFactory::sine_deformed_circle(pmp::Point2{-3.0f, 52.0f}, 35.0f, 25, 7.0f, 4.0f, 0.0, 2.0 * M_PI)},
	{ "SineDeformedIncompleteCircle", pmp::CurveFactory::sine_deformed_circle(pmp::Point2{-3.0f, 52.0f}, 35.0f, 25, 7.0f, 4.0f, M_PI_2, 2.0 * M_PI)},
	{ "ChamferedRectangle", pmp::CurveFactory::rectangle(pmp::Point2{-3.0f, 52.0f}, 60.0f, 70.0f, 15, true)},
	{ "IncompleteChamferedRectangle", pmp::CurveFactory::sampled_polygon({
		pmp::Point2{-30.0f, -35.0f} + pmp::Point2{-3.0f, 52.0f},
		pmp::Point2{30.0f, -35.0f} + pmp::Point2{-3.0f, 52.0f},
		pmp::Point2{30.0f, 35.0f} + pmp::Point2{-3.0f, 52.0f},
		pmp::Point2{-30.0f, 35.0f} + pmp::Point2{-3.0f, 52.0f}}, 30, true, false)},
	{ "ChamferedTriangle", pmp::CurveFactory::sampled_polygon({
		pmp::Point2{-0.5f, -sqrtf(3.0f) / 6.0f} *120.0f + pmp::Point2{-3.0f, 52.0f},
		pmp::Point2{0.5f, -sqrtf(3.0f) / 6.0f} *120.0f + pmp::Point2{-3.0f, 52.0f},
		pmp::Point2{0.0f, sqrtf(3.0f) / 3.0f} *120.0f + pmp::Point2{-3.0f, 52.0f}}, 30, true)},
	{ "IncompleteChamferedTriangle", pmp::CurveFactory::sampled_polygon({
		pmp::Point2{-0.5f, -sqrtf(3.0f) / 6.0f} *120.0f + pmp::Point2{-3.0f, 52.0f},
		pmp::Point2{0.5f, -sqrtf(3.0f) / 6.0f} *120.0f + pmp::Point2{-3.0f, 52.0f},
		pmp::Point2{0.0f, sqrtf(3.0f) / 3.0f} *120.0f + pmp::Point2{-3.0f, 52.0f}}, 30, true, false)}
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;
	constexpr unsigned int NTimeSteps = 180;

	const auto outerCircle = Circle2D{ pmp::Point2{-3.0f, 52.0f}, 110.0f };

	for (const auto& [ptCloudName, curve] : targetCurves)
	{
		std::cout << "========================================================\n";
		std::cout << "         inner/outer circle LSW for: " << ptCloudName << " ... \n";
		std::cout << " ------------------------------------------------------ \n";

		const auto& pts2D = curve.positions();
		const pmp::BoundingBox2 bbox{ pts2D };

		const auto bboxSize = bbox.max() - bbox.min();
		const float minSize = std::min(bboxSize[0], bboxSize[1]);
		//const float maxSize = std::max(bboxSize[0], bboxSize[1]);
		const float cellSize = minSize / nVoxelsPerMinDimension;

		const double isoLvlOffsetFactor = defaultOffsetFactor;
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.AdvectionInteractWithOtherManifolds = true;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 2.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = defaultTimeStep;
			strategySettings.LevelOfDetail = 4;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = "outer" + ptCloudName;
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);

			std::vector<pmp::ManifoldCurve2D> innerCurves{};

			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, outerCurve, innerCurves,
				std::make_shared<std::vector<pmp::Point2>>(pts2D));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	}
}

void SimpleMeshesIOLSWTests()
{
	const std::vector<std::string> meshForPtCloudNames{
		"boxWithAHole"
	};
	const std::map<std::string, double> timeStepSizesForPtClouds{
		{ "boxWithAHole", 0.05 }
	};
	const std::map<std::string, double> isoLevelOffsetFactors{
		{"boxWithAHole", 1.0 }
	};

	const std::map<std::string, Sphere3D> outerSpheres{
		{"boxWithAHole", Sphere3D{pmp::Point{0, 0, 0}, 2.0f} },
	};
	const std::map<std::string, std::vector<Sphere3D>> innerSpheres{
		{"boxWithAHole", std::vector{ Sphere3D{pmp::Point{0, 0, 0}, 0.9f}} },
	};

	const std::map<std::string, pmp::Point> slicingPlaneRefPts{
		{"boxWithAHole", pmp::Point{0.0f, 0.0f, 0.0f}},
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"boxWithAHole", pmp::vec3{0.0f, -1.0f, 0.0f}},
	};

	const std::map<std::string, Circle2D> outerCircles{
		{"boxWithAHole", Circle2D{pmp::Point2{0.0f, 0.0f}, 2.0f} },
	};
	const std::map<std::string, std::vector<Circle2D>> innerCircles{
		{"boxWithAHole", std::vector{ Circle2D{pmp::Point2{0.0f, 0.0f}, 0.9f}} }
	};

	const std::map<std::string, std::vector<Sphere3D>> cutSpheres{
		{"boxWithAHole", std::vector{ Sphere3D{pmp::Point{1.0f, 0.0f, 0.0f}, 0.6f} }},
	};

	constexpr unsigned int nVoxelsPerMinDimension = 30;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;

	constexpr size_t samplingLevel = 9;
	constexpr size_t nSamplings = 10;
	constexpr size_t minVerts = 9; // Minimum number of vertices to sample

	constexpr unsigned int seed = 5000; // seed for the pt cloud sampling RNG

	//SetRemeshingAdjustmentTimeIndices({}); // no remeshing adjustment
	SetRemeshingAdjustmentTimeIndices({ /*3, 10,*/ 30 /*, 50 , 100, 120, 140, 145*/ });

	constexpr unsigned int NTimeSteps = 80;

	constexpr bool executeCurveLSW = false;
	constexpr bool executeCurveIOLSW = false;
	constexpr bool executeSurfaceLSW = false;
	constexpr bool executeSurfaceIOLSW = true;

	for (const auto& meshName : meshForPtCloudNames)
	{
		// =======================================================================
		//   - - - - - - - - - - - - -   Data   Prep   - - - - - - - - - - - -
		// -----------------------------------------------------------------------

		std::cout << "==================================================================\n";
		std::cout << "Mesh To Pt Cloud: " << meshName << ".obj -> " << meshName << "Pts_" << samplingLevel << ".ply\n";
		std::cout << "------------------------------------------------------------------\n";
		const auto baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", false);
		if (!baseDataOpt.has_value())
		{
			std::cerr << "baseDataOpt == nullopt!\n";
			break;
		}
		std::cout << "meshName.obj" << " imported as BaseMeshGeometryData.\n";
		const auto& baseData = baseDataOpt.value();
		const size_t maxVerts = baseData.Vertices.size(); // Maximum number of vertices available
		size_t nVerts = minVerts + (maxVerts - minVerts) * samplingLevel / (nSamplings - 1);
		nVerts = std::max(minVerts, std::min(nVerts, maxVerts));

		std::cout << "Sampling " << nVerts << "/" << maxVerts << " vertices...\n";

		auto ptCloud = SamplePointsFromMeshData(baseData, nVerts, seed);
		const auto ptCloudName = meshName + "Pts_" + std::to_string(samplingLevel);

		const pmp::BoundingBox ptCloudBBox(ptCloud);
		const auto center = ptCloudBBox.center();
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float cellSize = minSize / nVoxelsPerMinDimension;
		constexpr float volExpansionFactor = 1.0f;
		const SDF::PointCloudDistanceFieldSettings dfSettings{
			cellSize,
			volExpansionFactor,
			Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::BlurPostprocessingType::None
		};

		if (!cutSpheres.contains(meshName))
		{
			std::cerr << "!cutSpheres.contains(\"" << meshName << "\") ... skipping.\n";
			continue;
		}

		DeletePointsContainedInSpheres(ptCloud, cutSpheres.at(meshName));

		std::cout << "Cutting " << nVerts - ptCloud.size() << " / " << nVerts << " points ... \n";

		std::string filename = dataOutPath + meshName + "Pts_" + std::to_string(samplingLevel) + ".ply";
		Geometry::BaseMeshGeometryData ptData;
		ptData.Vertices = ptCloud;
		if (!ExportPointsToPLY(ptData, filename))
		{
			std::cerr << "ExportPointsToPLY failed!\n";
			break;
		}
		const double isoLvlOffsetFactor = (timeStepSizesForPtClouds.contains(ptCloudName) ? isoLevelOffsetFactors.at(ptCloudName) : defaultOffsetFactor);
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		const double tau = (timeStepSizesForPtClouds.contains(ptCloudName) ? timeStepSizesForPtClouds.at(ptCloudName) : defaultTimeStep); // time step

		// ==========================================================================
		// - - - - - - - - -  New Manifold Evolver (Curve)  - - - - - - - - - - - - 
		// ==========================================================================

		const float distTolerance = 0.01f * minSize;
		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0f, 0.0f, 0.0f });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}

		if (!Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
		{
			std::cerr << "Export2DPointCloudToPLY: internal error during export!\n";
			continue;
		}

		// ------------------------- only outer ---------------------------
		if (executeCurveLSW)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = false;
			strategySettings.AdvectionInteractWithOtherManifolds = false;
			strategySettings.OuterManifoldEpsilon = [](double distance)
				{
					if (distance <= 0.0)
						return 0.0;
					return 1.0 * (1.0 - exp(-distance * distance / 1.0));
				};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
				{
					if (distance <= 0.0)
						return 0.0;
					return 1.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
				};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 4;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "_CurveLSW";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerCircles.contains(meshName))
			{
				std::cerr << "!outerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			const auto outerCircle = outerCircles.at(meshName);
			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);

			std::vector<pmp::ManifoldCurve2D> innerCurves;  // no inner curves to evolve

			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, outerCurve, innerCurves,
				std::make_shared<std::vector<pmp::Point2>>(pts2D));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}

		// --------------------- both inner and outer ---------------------
		if (executeCurveIOLSW)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.AdvectionInteractWithOtherManifolds = true;
			strategySettings.OuterManifoldEpsilon = [](double distance)
				{
					if (distance <= 0.0)
						return 0.0;
					return 1.0 * (1.0 - exp(-distance * distance / 1.0));
				};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
				{
					if (distance <= 0.0)
						return 0.0;
					return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
				};
			strategySettings.InnerManifoldEpsilon = [](double distance)
				{
					if (distance <= 0.0)
						return 0.0;
					return 0.0005 * TRIVIAL_EPSILON(distance);
				};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
				{
					if (distance <= 0.0)
						return 0.0;
					return 0.4 * distance * (negGradDotNormal - 1.3 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
				};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 4;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "_CurveIOLSW";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerCircles.contains(meshName))
			{
				std::cerr << "!outerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			if (!innerCircles.contains(meshName))
			{
				std::cerr << "!innerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			const auto outerCircle = outerCircles.at(meshName);
			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);

			std::vector<pmp::ManifoldCurve2D> innerCurves;
			for (const auto& innerCircle : innerCircles.at(meshName))
			{
				const auto nInnerSegments = static_cast<unsigned int>(static_cast<pmp::Scalar>(nSegments) * (innerCircle.Radius * 2.0) / outerCircle.Radius);
				auto innerCurve = pmp::CurveFactory::circle(innerCircle.Center, innerCircle.Radius, nInnerSegments);
				innerCurve.negate_orientation();
				innerCurves.push_back(innerCurve);
			}

			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, outerCurve, innerCurves,
				std::make_shared<std::vector<pmp::Point2>>(pts2D));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}

		// ==========================================================================
		// - - - - - - - - -  New Manifold Evolver (Surface)  - - - - - - - - - - - - 
		// ==========================================================================

		// ------------------------------ Only outer --------------------------------
		if (executeSurfaceLSW)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = false;
			strategySettings.AdvectionInteractWithOtherManifolds = false;

			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				if (distance <= 0.0)
					return 0.0;
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				if (distance <= 0.0)
					return 0.0;
				return 1.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "_SurfaceLSW";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerSpheres.contains(meshName))
			{
				std::cerr << "!outerSpheres.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}
			const auto outerSphere = outerSpheres.at(meshName);

			// construct outer ico-sphere
			Geometry::IcoSphereBuilder icoBuilder({ strategySettings.LevelOfDetail, outerSphere.Radius });
			icoBuilder.BuildBaseData();
			icoBuilder.BuildPMPSurfaceMesh();
			auto outerSurface = icoBuilder.GetPMPSurfaceMeshResult();
			const pmp::mat4 transfMatrixGeomMove{
				1.0f, 0.0f, 0.0f, outerSphere.Center[0],
				0.0f, 1.0f, 0.0f, outerSphere.Center[1],
				0.0f, 0.0f, 1.0f, outerSphere.Center[2],
				0.0f, 0.0f, 0.0f, 1.0f
			};
			outerSurface *= transfMatrixGeomMove;
			std::vector<pmp::SurfaceMesh> innerSurfaces; // no inner surfaces to evolve

			auto surfaceStrategy = std::make_shared<CustomManifoldSurfaceEvolutionStrategy>(
				strategySettings, MeshLaplacian::Voronoi,
				outerSurface, innerSurfaces,
				std::make_shared<std::vector<pmp::Point>>(ptCloud));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(surfaceStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}

		// -------------------------- Both inner and outer ---------------------------
		if (executeSurfaceIOLSW)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.AdvectionInteractWithOtherManifolds = true;

			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				if (distance <= 0.0)
					return 0.0;
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				if (distance <= 0.0)
					return 0.0;
				return 1.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				if (distance <= 0.0)
					return 0.0;
				return 0.0005 * TRIVIAL_EPSILON(distance);
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				if (distance <= 0.0)
					return 0.0;
				return 0.4 * distance * (negGradDotNormal - 1.3 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			strategySettings.ExportVariableScalarFieldsDimInfo = true;
			//strategySettings.ExportVariableVectorFieldsDimInfo = true;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "_SurfaceIOLSW";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerSpheres.contains(meshName))
			{
				std::cerr << "!outerSpheres.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}
			const auto outerSphere = outerSpheres.at(meshName);

			// construct outer ico-sphere
			Geometry::IcoSphereBuilder icoBuilder({ strategySettings.LevelOfDetail, outerSphere.Radius });
			icoBuilder.BuildBaseData();
			icoBuilder.BuildPMPSurfaceMesh();
			auto outerSurface = icoBuilder.GetPMPSurfaceMeshResult();
			const pmp::mat4 transfMatrixGeomMove{
				1.0f, 0.0f, 0.0f, outerSphere.Center[0],
				0.0f, 1.0f, 0.0f, outerSphere.Center[1],
				0.0f, 0.0f, 1.0f, outerSphere.Center[2],
				0.0f, 0.0f, 0.0f, 1.0f
			};
			outerSurface *= transfMatrixGeomMove;

			if (!innerSpheres.contains(meshName))
			{
				std::cerr << "!innerSpheres.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}
			const auto innerSpheresPerMesh = innerSpheres.at(meshName);
			std::vector<pmp::SurfaceMesh> innerSurfaces;
			innerSurfaces.reserve(innerSpheresPerMesh.size());
			for (const auto& innerSphere : innerSpheresPerMesh)
			{
				// construct innrt ico-sphere
				const auto innerSubdiv = static_cast<unsigned int>(static_cast<pmp::Scalar>(strategySettings.LevelOfDetail) * (innerSphere.Radius * 2.0) / outerSphere.Radius);
				Geometry::IcoSphereBuilder innerIcoBuilder({ innerSubdiv, innerSphere.Radius });
				innerIcoBuilder.BuildBaseData();
				innerIcoBuilder.BuildPMPSurfaceMesh();
				auto innerSurface = innerIcoBuilder.GetPMPSurfaceMeshResult();
				const pmp::mat4 transfMatrixGeomMove{
					1.0f, 0.0f, 0.0f, innerSphere.Center[0],
					0.0f, 1.0f, 0.0f, innerSphere.Center[1],
					0.0f, 0.0f, 1.0f, innerSphere.Center[2],
					0.0f, 0.0f, 0.0f, 1.0f
				};
				innerSurface *= transfMatrixGeomMove;
				innerSurface.negate_orientation();
				innerSurfaces.push_back(innerSurface);
			}

			auto surfaceStrategy = std::make_shared<CustomManifoldSurfaceEvolutionStrategy>(
				strategySettings, MeshLaplacian::Voronoi,
				outerSurface, innerSurfaces,
				std::make_shared<std::vector<pmp::Point>>(ptCloud));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(surfaceStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	}
}

void StandardMeshesIOLSWTests()
{
	const std::vector<std::string> meshForPtCloudNames{
		//"armadillo",
		//"blub",
		"bunny",
		"maxPlanck",
		//"nefertiti",
		//"ogre",
		//"spot"
	};
	const std::map<std::string, double> timeStepSizesForPtClouds{
		{ "armadillo", 0.05 },
		{ "blub", 0.05 },
		{ "bunny", 0.05 },
		{ "maxPlanck", 0.05 },
		{ "nefertiti", 0.05 },
		{ "ogre", 0.05 },
		{ "spot", 0.05 }
	};
	const std::map<std::string, double> isoLevelOffsetFactors{
		{"armadillo", 1.5 },
		{ "blub", 0.5 },
		{ "bunny", 1.5 },
		{ "maxPlanck", 1.5 },
		{ "nefertiti", 0.5 },
		{ "ogre", 0.5 },
		{ "spot", 0.5 }
	};

	const std::map<std::string, pmp::Point> slicingPlaneRefPts{
		{"armadillo", pmp::Point{-0.10348621158928151f, 21.427067319905646f, 9.79369240592005f}},
		{"bunny", pmp::Point{-0.01684039831161499f, 0.11015420407056808f, 0.0012007840834242693f} },
		//{"bunny", pmp::Point{-0.021311134388792542f, 0.11143290480956525f, 0.007428090888817638f} },
		{"maxPlanck", pmp::Point{30.59686279296875f, -18.105804443359375f, 82.29149055480957f} },
		{"nefertiti", pmp::Point{0.0f, 0.0f, 0.0f} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"armadillo", pmp::vec3{-0.03070969905335075f, 0.12876712096541565f, 0.9911992448253433f}},
		{"bunny", pmp::vec3{0.0f, 0.0f, 1.0f} },
		//{"bunny", pmp::vec3{-0.18754182700901353f, -0.029010506066971985f, 0.9818281181855913f} },
		{"maxPlanck", pmp::vec3{1.0f, 0.0f, 0.0f} },
		{"nefertiti", pmp::vec3{1.0f, 0.0f, 0.0f} }
	};

	const std::map<std::string, Circle2D> outerCircles{
		{"armadillo", Circle2D{pmp::Point2{0.372234f, 16.6515f}, 121.558f} },
		{"bunny", Circle2D{pmp::Point2{-0.0155906f, 0.102261f}, 0.142831f} },
		{"maxPlanck", Circle2D{pmp::Point2{-17.82f, 82.5006f}, 292.263f} },
		{"nefertiti", Circle2D{pmp::Point2{0.178497f, -0.0410004f}, 441.436f} }
	};
	const std::map<std::string, std::vector<Circle2D>> innerCircles{
		//{"armadillo", std::vector{ Circle2D{pmp::Point2{-3.0f, 52.0f}, 20.0f}} },
		{"bunny", std::vector{ Circle2D{pmp::Point2{0.0f, 0.08f}, 0.025f}} },
		{"maxPlanck", std::vector{ Circle2D{pmp::Point2{8.0f, 85.0f}, 50.0f}} },
		//{"nefertiti", std::vector{ Circle2D{pmp::Point2{-20.0f, 100.0f}, 55.0f}} }
	};

	const std::map<std::string, Sphere3D> outerSpheres{
		{"armadillo", Sphere3D{pmp::Point{0.0122509f, 21.4183f, -0.000249863f}, 136.963f} },
		{"bunny", Sphere3D{pmp::Point{-0.0168297f, 0.110217f, -0.0015718f}, 0.141622f} },
		{"maxPlanck", Sphere3D{pmp::Point{30.658f, -17.9765f, 82.2885f}, 271.982f} },
		{"nefertiti", Sphere3D{pmp::Point{0.0144997f, -0.00499725f, -0.0215073f}, 392.184f} }
	};
	const std::map<std::string, std::vector<Sphere3D>> innerSpheres{
		//{"armadillo", std::vector{ Sphere3D{pmp::Point{-3.0f, 52.0f}, 20.0f}} },
		{"bunny", std::vector{ Sphere3D{pmp::Point{0.0f, 0.082f, 0.012f}, 0.03f}} },
		{"maxPlanck", std::vector{ Sphere3D{pmp::Point{8.0f, 50.0f, 100.0f}, 50.0f}} },
		//{"nefertiti", std::vector{ Sphere3D{pmp::Point{-20.0f, 100.0f}, 55.0f}} }
	};

	const std::map<std::string, std::vector<Sphere3D>> cutSpheres{
		{"armadillo", {}},
		{"bunny", std::vector{ Sphere3D{pmp::Point{-0.01f, 0.06f, 0.012f}, 0.032f}, Sphere3D{pmp::Point{0.01f, 0.12f, 0.01f}, 0.025f}/**/}},
		{"maxPlanck", std::vector{ Sphere3D{pmp::Point{8.0f, 85.0f, 0.0f}, 50.0f}, Sphere3D{pmp::Point{30.0f, -120.0f, 160.0f}, 100.0f} /**/}},
		{"nefertiti", {}}
	};

	constexpr unsigned int nVoxelsPerMinDimension = 30;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;

	constexpr size_t samplingLevel = 3;
	constexpr size_t nSamplings = 10;
	constexpr size_t minVerts = 9; // Minimum number of vertices to sample

	constexpr unsigned int seed = 5000; // seed for the pt cloud sampling RNG

	//SetRemeshingAdjustmentTimeIndices({}); // no remeshing adjustment
	SetRemeshingAdjustmentTimeIndices({ /*3, 10,*/ 30 /*, 50 , 100, 120, 140, 145*/ });

	constexpr unsigned int NTimeSteps = 180;

	constexpr bool executeCurveLSW = false;
	constexpr bool executeCurveIOLSW = false;
	constexpr bool executeSurfaceLSW = false;
	constexpr bool executeSurfaceIOLSW = true;

	for (const auto& meshName : meshForPtCloudNames)
	{
		// =======================================================================
		//   - - - - - - - - - - - - -   Data   Prep   - - - - - - - - - - - -
		// -----------------------------------------------------------------------

		std::cout << "==================================================================\n";
		std::cout << "Mesh To Pt Cloud: " << meshName << ".obj -> " << meshName << "Pts_" << samplingLevel << ".ply\n";
		std::cout << "------------------------------------------------------------------\n";
		const auto baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", false);
		if (!baseDataOpt.has_value())
		{
			std::cerr << "baseDataOpt == nullopt!\n";
			break;
		}
		std::cout << "meshName.obj" << " imported as BaseMeshGeometryData.\n";
		const auto& baseData = baseDataOpt.value();
		const size_t maxVerts = baseData.Vertices.size(); // Maximum number of vertices available
		size_t nVerts = minVerts + (maxVerts - minVerts) * samplingLevel / (nSamplings - 1);
		nVerts = std::max(minVerts, std::min(nVerts, maxVerts));

		std::cout << "Sampling " << nVerts << "/" << maxVerts << " vertices...\n";

		auto ptCloud = SamplePointsFromMeshData(baseData, nVerts, seed);
		const auto ptCloudName = meshName + "Pts_" + std::to_string(samplingLevel);

		if (!cutSpheres.contains(meshName))
		{
			std::cerr << "!cutSpheres.contains(\"" << meshName << "\") ... skipping.\n";
			continue;
		}

		DeletePointsContainedInSpheres(ptCloud, cutSpheres.at(meshName));

		std::cout << "Cutting " << nVerts - ptCloud.size() << " / " << nVerts << " points ... \n";

		std::string filename = dataOutPath + meshName + "Pts_" + std::to_string(samplingLevel) + ".ply";
		Geometry::BaseMeshGeometryData ptData;
		ptData.Vertices = ptCloud;
		if (!ExportPointsToPLY(ptData, filename))
		{
			std::cerr << "ExportPointsToPLY failed!\n";
			break;
		}

		const pmp::BoundingBox ptCloudBBox(ptCloud);
		const auto center = ptCloudBBox.center();
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float cellSize = minSize / nVoxelsPerMinDimension;
		constexpr float volExpansionFactor = 1.0f;
		const SDF::PointCloudDistanceFieldSettings dfSettings{
			cellSize,
			volExpansionFactor,
			Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::BlurPostprocessingType::None
		};

		const double isoLvlOffsetFactor = (timeStepSizesForPtClouds.contains(ptCloudName) ? isoLevelOffsetFactors.at(ptCloudName) : defaultOffsetFactor);
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		const double tau = (timeStepSizesForPtClouds.contains(ptCloudName) ? timeStepSizesForPtClouds.at(ptCloudName) : defaultTimeStep); // time step

		// ==========================================================================
		// - - - - - - - - -  New Manifold Evolver (Curve)  - - - - - - - - - - - - 
		// ==========================================================================

		const float distTolerance = 0.01f * minSize;
		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0f, 0.0f, 0.0f });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}

		if (!Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
		{
			std::cerr << "Export2DPointCloudToPLY: internal error during export!\n";
			continue;
		}

		// ------------------------- only outer ---------------------------
		if (executeCurveLSW)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = false;
			strategySettings.AdvectionInteractWithOtherManifolds = false;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 4;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "_CurveLSW";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerCircles.contains(meshName))
			{
				std::cerr << "!outerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			const auto outerCircle = outerCircles.at(meshName);
			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);

			std::vector<pmp::ManifoldCurve2D> innerCurves;  // no inner curves to evolve

			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, outerCurve, innerCurves,
				std::make_shared<std::vector<pmp::Point2>>(pts2D));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}

		// --------------------- both inner and outer ---------------------
		if (executeCurveIOLSW)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.AdvectionInteractWithOtherManifolds = true;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				return 0.0005 * TRIVIAL_EPSILON(distance);
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 0.6 * distance * (negGradDotNormal - 1.2 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 4;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			strategySettings.ExportVariableScalarFieldsDimInfo = true;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "_CurveIOLSW";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerCircles.contains(meshName))
			{
				std::cerr << "!outerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			if (!innerCircles.contains(meshName))
			{
				std::cerr << "!innerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			const auto outerCircle = outerCircles.at(meshName);
			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);

			std::vector<pmp::ManifoldCurve2D> innerCurves;
			for (const auto& innerCircle : innerCircles.at(meshName))
			{
				const auto nInnerSegments = static_cast<unsigned int>(static_cast<pmp::Scalar>(nSegments) * (innerCircle.Radius * 2.0) / outerCircle.Radius);
				auto innerCurve = pmp::CurveFactory::circle(innerCircle.Center, innerCircle.Radius, nInnerSegments);
				innerCurve.negate_orientation();
				innerCurves.push_back(innerCurve);
			}

			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, outerCurve, innerCurves,
				std::make_shared<std::vector<pmp::Point2>>(pts2D));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}

		// ==========================================================================
		// - - - - - - - - -  New Manifold Evolver (Surface)  - - - - - - - - - - - - 
		// ==========================================================================

		// ------------------------------ Only outer --------------------------------
		if (executeSurfaceLSW)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = false;
			strategySettings.AdvectionInteractWithOtherManifolds = false;

			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "_SurfaceLSW";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerSpheres.contains(meshName))
			{
				std::cerr << "!outerSpheres.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}
			const auto outerSphere = outerSpheres.at(meshName);

			// construct outer ico-sphere
			Geometry::IcoSphereBuilder icoBuilder({ strategySettings.LevelOfDetail, outerSphere.Radius });
			icoBuilder.BuildBaseData();
			icoBuilder.BuildPMPSurfaceMesh();
			auto outerSurface = icoBuilder.GetPMPSurfaceMeshResult();
			const pmp::mat4 transfMatrixGeomMove{
				1.0f, 0.0f, 0.0f, outerSphere.Center[0],
				0.0f, 1.0f, 0.0f, outerSphere.Center[1],
				0.0f, 0.0f, 1.0f, outerSphere.Center[2],
				0.0f, 0.0f, 0.0f, 1.0f
			};
			outerSurface *= transfMatrixGeomMove;			
			std::vector<pmp::SurfaceMesh> innerSurfaces; // no inner surfaces to evolve

			auto surfaceStrategy = std::make_shared<CustomManifoldSurfaceEvolutionStrategy>(
				strategySettings, MeshLaplacian::Voronoi,
				outerSurface, innerSurfaces,
				std::make_shared<std::vector<pmp::Point>>(ptCloud));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(surfaceStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}

		// -------------------------- Both inner and outer ---------------------------
		if (executeSurfaceIOLSW)
		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.AdvectionInteractWithOtherManifolds = true;

			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				if (distance <= 0.0)
					return 0.0;
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				if (distance <= 0.0)
					return 0.0;
				return 1.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				if (distance <= 0.0)
					return 0.0;
				return 0.0005 * TRIVIAL_EPSILON(distance);
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				if (distance <= 0.0)
					return 0.0;
				return 0.4 * distance * (negGradDotNormal - 1.3 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "_SurfaceIOLSW";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerSpheres.contains(meshName))
			{
				std::cerr << "!outerSpheres.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}
			const auto outerSphere = outerSpheres.at(meshName);

			// construct outer ico-sphere
			Geometry::IcoSphereBuilder icoBuilder({ strategySettings.LevelOfDetail, outerSphere.Radius });
			icoBuilder.BuildBaseData();
			icoBuilder.BuildPMPSurfaceMesh();
			auto outerSurface = icoBuilder.GetPMPSurfaceMeshResult();
			const pmp::mat4 transfMatrixGeomMove{
				1.0f, 0.0f, 0.0f, outerSphere.Center[0],
				0.0f, 1.0f, 0.0f, outerSphere.Center[1],
				0.0f, 0.0f, 1.0f, outerSphere.Center[2],
				0.0f, 0.0f, 0.0f, 1.0f
			};
			outerSurface *= transfMatrixGeomMove;

			if (!innerSpheres.contains(meshName))
			{
				std::cerr << "!innerSpheres.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}
			const auto innerSpheresPerMesh = innerSpheres.at(meshName);
			std::vector<pmp::SurfaceMesh> innerSurfaces;
			innerSurfaces.reserve(innerSpheresPerMesh.size());
			for (const auto& innerSphere : innerSpheresPerMesh)
			{
				// construct innrt ico-sphere
				const auto innerSubdiv = static_cast<unsigned int>(static_cast<pmp::Scalar>(strategySettings.LevelOfDetail) * (innerSphere.Radius * 2.0) / outerSphere.Radius);
				Geometry::IcoSphereBuilder innerIcoBuilder({ innerSubdiv, innerSphere.Radius });
				innerIcoBuilder.BuildBaseData();
				innerIcoBuilder.BuildPMPSurfaceMesh();
				auto innerSurface = innerIcoBuilder.GetPMPSurfaceMeshResult();
				const pmp::mat4 transfMatrixGeomMove{
					1.0f, 0.0f, 0.0f, innerSphere.Center[0],
					0.0f, 1.0f, 0.0f, innerSphere.Center[1],
					0.0f, 0.0f, 1.0f, innerSphere.Center[2],
					0.0f, 0.0f, 0.0f, 1.0f
				};
				innerSurface *= transfMatrixGeomMove;
				innerSurface.negate_orientation();
				innerSurfaces.push_back(innerSurface);
			}

			auto surfaceStrategy = std::make_shared<CustomManifoldSurfaceEvolutionStrategy>(
				strategySettings, MeshLaplacian::Voronoi,
				outerSurface, innerSurfaces,
				std::make_shared<std::vector<pmp::Point>>(ptCloud));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(surfaceStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	}
}

void MedialAxisTests()
{
	pmp::ManifoldCurve2D hyperEllipse = pmp::CurveFactory::hyper_ellipse(pmp::Point2{ 0, 0 }, 300.0, 150.0, 4, 10);
	if (!pmp::write_to_ply(hyperEllipse, dataOutPath + "hyperEllipse.ply"))
		std::cerr << "Error writing hyperEllipse.ply!\n";

	//hyperEllipse.negate_orientation();
	//const auto medialAxis = Geometry::CalculateApproxMedialAxisFromCurve(hyperEllipse);
	//if (medialAxis.has_value())
	//{
	//	if (!Geometry::ExportBaseCurveGeometryDataToPLY(*medialAxis, dataOutPath + "hyperEllipse_medialAxis.ply"))
	//		std::cerr << "Error writing hyperEllipse_medialAxis.ply!\n";
	//}

	// Rohan Sawney's code is amazingly non-universal! Only his examples work.

	{
		const auto medialAxis = Geometry::GetMedialAxisOfSawhneysStupidMATAlgorithm('q');
		if (medialAxis.has_value())
		{
			if (!Geometry::ExportBaseCurveGeometryDataToPLY(*medialAxis, dataOutPath + "sawhneyShapeQ.ply"))
				std::cerr << "Error writing sawhneyShapeQ.ply!\n";
		}
	}

	{
		const auto medialAxis = Geometry::GetMedialAxisOfSawhneysStupidMATAlgorithm('w');
		if (medialAxis.has_value())
		{
			if (!Geometry::ExportBaseCurveGeometryDataToPLY(*medialAxis, dataOutPath + "sawhneyShapeW.ply"))
				std::cerr << "Error writing sawhneyShapeW.ply!\n";
		}
	}

	{
		const auto medialAxis = Geometry::GetMedialAxisOfSawhneysStupidMATAlgorithm('e');
		if (medialAxis.has_value())
		{
			if (!Geometry::ExportBaseCurveGeometryDataToPLY(*medialAxis, dataOutPath + "sawhneyShapeE.ply"))
				std::cerr << "Error writing sawhneyShapeE.ply!\n";
		}
	}

	{
		const auto medialAxis = Geometry::GetMedialAxisOfSawhneysStupidMATAlgorithm('r');
		if (medialAxis.has_value())
		{
			if (!Geometry::ExportBaseCurveGeometryDataToPLY(*medialAxis, dataOutPath + "sawhneyShapeR.ply"))
				std::cerr << "Error writing sawhneyShapeR.ply!\n";
		}
	}

	{
		const auto medialAxis = Geometry::GetMedialAxisOfSawhneysStupidMATAlgorithm('t');
		if (medialAxis.has_value())
		{
			if (!Geometry::ExportBaseCurveGeometryDataToPLY(*medialAxis, dataOutPath + "sawhneyShapeT.ply"))
				std::cerr << "Error writing sawhneyShapeT.ply!\n";
		}
	}
}

void HyperellipseEllipsoidEquilibriumTests()
{
	pmp::ManifoldCurve2D hyperEllipse = pmp::CurveFactory::hyper_ellipse(pmp::Point2{ 0, 0 }, 200.0, 110.0, 4, 40);
	//hyperEllipse.negate_orientation();
	if (!pmp::write_to_ply(hyperEllipse, dataOutPath + "hyperEllipse.ply"))
		std::cerr << "Error writing hyperEllipse.ply!\n";
	// remesh hyper ellipse
	constexpr float edgeLength = 15.0f;
	constexpr unsigned int iterations = 10;
	pmp::CurveRemeshing remesher(hyperEllipse);
	pmp::AdaptiveRemeshingSettings settings;
	settings.MinEdgeLength = edgeLength;
	settings.MaxEdgeLength = 1.5f * edgeLength;
	settings.ApproxError = 0.05f * edgeLength;
	settings.NRemeshingIterations = iterations;
	settings.NTangentialSmoothingIters = 6;
	settings.UseProjection = true;
	remesher.adaptive_remeshing(settings);	

	const std::vector<Circle2D> testInnerCircles{
		//{ pmp::Point2{ 0.0f, 0.0f }, 70.0f },
		//{ pmp::Point2{ 100.0f, 0.0f }, 60.0f },
		//{ pmp::Point2{ 140.0f, 0.0f }, 40.0f }
		//{ pmp::Point2{ 130.0f, 0.0f }, 60.0f },
		{ pmp::Point2{ 130.0f, 40.0f }, 40.0f }
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;
	constexpr unsigned int NTimeSteps = 180;
	const double fieldIsoLevel = defaultOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(5.0);

	for (size_t hyperellipseId = 9; const auto& innerCircle : testInnerCircles)
	{
		std::cout << "Setting up ManifoldEvolutionSettings.\n";

		ManifoldEvolutionSettings strategySettings;
		strategySettings.UseInnerManifolds = true;
		strategySettings.AdvectionInteractWithOtherManifolds = true;
		strategySettings.OuterManifoldEpsilon = [](double distance)
			{
			return 0.00; // * (1.0 - exp(-distance * distance / 1.0));
			};
		strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
			return 0.00; // * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
		strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				return 0.0025 * TRIVIAL_EPSILON(distance);
			};
		strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
		strategySettings.TimeStep = defaultTimeStep;
		strategySettings.LevelOfDetail = 4;
		strategySettings.TangentialVelocityWeight = 0.05;

		strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.22f;
		strategySettings.RemeshingSettings.UseBackProjection = false;

		strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
		strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

		strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
		strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

		std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

		GlobalManifoldEvolutionSettings globalSettings;
		globalSettings.NSteps = NTimeSteps;
		globalSettings.DoRemeshing = true;
		globalSettings.DetectFeatures = false;
		globalSettings.ExportPerTimeStep = true;
		globalSettings.ExportTargetDistanceFieldAsImage = true;
		globalSettings.ProcedureName = "hyperellipse" + std::to_string(hyperellipseId);
		globalSettings.OutputPath = dataOutPath;
		globalSettings.ExportResult = false;

		globalSettings.RemeshingResizeFactor = 0.7f;
		globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

		const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
		const auto outerCurve = hyperEllipse;

		std::vector<pmp::ManifoldCurve2D> innerCurves{
			pmp::CurveFactory::circle(innerCircle.Center, innerCircle.Radius, nSegments)
		};
		innerCurves[0].negate_orientation();

		auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
			strategySettings, outerCurve, innerCurves, nullptr);

		std::cout << "Setting up ManifoldEvolver.\n";

		ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

		std::cout << "ManifoldEvolver::Evolve ... ";

		try
		{
			evolver.Evolve();
		}
		catch (std::invalid_argument& ex)
		{
			std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
		}
		catch (std::runtime_error& ex)
		{
			std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
		}
		catch (...)
		{
			std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
		}

		hyperellipseId++;
	}
}

void JunkCan2DTests()
{
	const std::vector<std::string> meshForPtCloudNames{
		//"canStraight",
		"canStraightMissingBottom",
		//"crushedCan",
		"crushedCanMissingBottom"
	};
	const std::map<std::string, double> timeStepSizesForPtClouds{
		{ "canStraight", 0.05 },
		{ "canStraightMissingBottom", 0.05 },
		{ "crushedCan", 0.05 },
		{ "crushedCanMissingBottom", 0.05 }
	};
	const std::map<std::string, double> isoLevelOffsetFactors{
		{ "canStraight", 0.5 },
		{ "canStraightMissingBottom", 0.5 },
		{ "crushedCan", 0.5 },
		{ "crushedCanMissingBottom", 0.5 }
	};

	const std::map<std::string, pmp::Point> slicingPlaneRefPts{
		{"canStraight", pmp::Point{0.07009050250053406f, 0.3348689912818372f, -0.018547505140304565f}},
		{"canStraightMissingBottom", pmp::Point{0.07009050250053406f, 0.3348689912818372f, -0.018547505140304565f} },
		{"crushedCan", pmp::Point{0.0007844995707273483f, 0.030990499537438154f, -0.006461501121520996} },
		{"crushedCanMissingBottom", pmp::Point{0.0007844995707273483f, 0.030990499537438154f, -0.006461501121520996} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"canStraight", pmp::vec3{1.0f, 0.0f, 0.0f}},
		{"canStraightMissingBottom", pmp::vec3{1.0f, 0.0f, 0.0f} },
		{"crushedCan", pmp::vec3{1.0f, 0.0f, 0.0f} },
		{"crushedCanMissingBottom", pmp::vec3{1.0f, 0.0f, 0.0f} }
	};

	const std::map<std::string, Circle2D> outerCircles{
		{"canStraight", Circle2D{pmp::Point2{0.25f, 0.0f}, 0.9f} },
		{"canStraightMissingBottom", Circle2D{pmp::Point2{0.25f, 0.0f}, 0.9f} },
		{"crushedCan", Circle2D{pmp::Point2{0.0205f, 0.0f}, 0.07f} },
		{"crushedCanMissingBottom", Circle2D{pmp::Point2{0.0205f, 0.0f}, 0.07f} }
	};
	const std::map<std::string, std::vector<Circle2D>> innerCircles{
		{"canStraight", std::vector{ Circle2D{pmp::Point2{0.35f, -0.25f}, 0.16f} } },
		{"canStraightMissingBottom", std::vector{ Circle2D{pmp::Point2{0.35f, -0.25f}, 0.16f}, Circle2D{pmp::Point2{0.25f, 0.35f}, 0.16f}} },
		{"crushedCan", std::vector{ Circle2D{pmp::Point2{0.025f, -0.035f}, 0.02f}} },
		{"crushedCanMissingBottom", std::vector{ Circle2D{pmp::Point2{0.030f, -0.035f}, 0.02f}, Circle2D{pmp::Point2{0.025f, 0.035f}, 0.02f} } }
	};

	constexpr unsigned int nVoxelsPerMinDimension = 30;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;

	constexpr size_t samplingLevel = 7;
	constexpr size_t nSamplings = 10;
	constexpr size_t minVerts = 9; // Minimum number of vertices to sample

	constexpr unsigned int seed = 5000; // seed for the pt cloud sampling RNG

	//SetRemeshingAdjustmentTimeIndices({}); // no remeshing adjustment
	SetRemeshingAdjustmentTimeIndices({ /*3, 10,*/ 30 /*, 50 , 100, 120, 140, 145*/ });

	constexpr unsigned int NTimeSteps = 180;

	constexpr bool executeNewSurfaceCustomEvolver = false;
	constexpr bool executeNewCurveCustomEvolver = true;

	for (const auto& meshName : meshForPtCloudNames)
	{
		// =======================================================================
		//   - - - - - - - - - - - - -   Data   Prep   - - - - - - - - - - - -
		// -----------------------------------------------------------------------

		std::cout << "==================================================================\n";
		std::cout << "Mesh To Pt Cloud: " << meshName << ".obj -> " << meshName << "Pts_" << samplingLevel << ".ply\n";
		std::cout << "------------------------------------------------------------------\n";
		const auto baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", false);
		if (!baseDataOpt.has_value())
		{
			std::cerr << "baseDataOpt == nullopt!\n";
			break;
		}
		std::cout << "meshName.obj" << " imported as BaseMeshGeometryData.\n";
		const auto& baseData = baseDataOpt.value();
		const size_t maxVerts = baseData.Vertices.size(); // Maximum number of vertices available
		size_t nVerts = minVerts + (maxVerts - minVerts) * samplingLevel / (nSamplings - 1);
		nVerts = std::max(minVerts, std::min(nVerts, maxVerts));

		std::cout << "Sampling " << nVerts << "/" << maxVerts << " vertices...\n";

		auto ptCloud = SamplePointsFromMeshData(baseData, nVerts, seed);
		const auto ptCloudName = meshName + "Pts_" + std::to_string(samplingLevel);

		std::cout << "Cutting " << nVerts - ptCloud.size() << " / " << nVerts << " points ... \n";

		std::string filename = dataOutPath + meshName + "Pts_" + std::to_string(samplingLevel) + ".ply";
		Geometry::BaseMeshGeometryData ptData;
		ptData.Vertices = ptCloud;
		if (!ExportPointsToPLY(ptData, filename))
		{
			std::cerr << "ExportPointsToPLY failed!\n";
			break;
		}

		const pmp::BoundingBox ptCloudBBox(ptCloud);
		const auto center = ptCloudBBox.center();
		const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
		const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const float cellSize = minSize / nVoxelsPerMinDimension;
		constexpr float volExpansionFactor = 1.0f;
		const SDF::PointCloudDistanceFieldSettings dfSettings{
			cellSize,
			volExpansionFactor,
			Geometry::DEFAULT_SCALAR_GRID_INIT_VAL,
			SDF::BlurPostprocessingType::None
		};

		const double isoLvlOffsetFactor = (timeStepSizesForPtClouds.contains(ptCloudName) ? isoLevelOffsetFactors.at(ptCloudName) : defaultOffsetFactor);
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		const double tau = (timeStepSizesForPtClouds.contains(ptCloudName) ? timeStepSizesForPtClouds.at(ptCloudName) : defaultTimeStep); // time step

		// ==========================================================================
		// - - - - - - - - -  New Manifold Evolver (Curve)  - - - - - - - - - - - - 
		// ==========================================================================

		const float distTolerance = 0.05f * minSize;
		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0f, 0.0f, 0.0f });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}

		if (!Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
		{
			std::cerr << "Export2DPointCloudToPLY: internal error during export!\n";
			continue;
		}


		{
			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.AdvectionInteractWithOtherManifolds = true;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				return 0.0005 * TRIVIAL_EPSILON(distance);
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 0.6 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 4;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14f;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = meshName + "_CurveIOLSW";
			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7f;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			if (!outerCircles.contains(meshName))
			{
				std::cerr << "!outerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			if (!innerCircles.contains(meshName))
			{
				std::cerr << "!innerCircles.contains(\"" << meshName << "\") ... skipping.\n";
				continue;
			}

			const auto outerCircle = outerCircles.at(meshName);
			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);

			std::vector<pmp::ManifoldCurve2D> innerCurves;
			for (const auto& innerCircle : innerCircles.at(meshName))
			{
				const auto nInnerSegments = static_cast<unsigned int>(static_cast<pmp::Scalar>(nSegments) * (innerCircle.Radius * 2.0) / outerCircle.Radius);
				auto innerCurve = pmp::CurveFactory::circle(innerCircle.Center, innerCircle.Radius, nInnerSegments);
				innerCurve.negate_orientation();
				innerCurves.push_back(innerCurve);
			}

			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, outerCurve, innerCurves,
				std::make_shared<std::vector<pmp::Point2>>(pts2D));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			try
			{
				evolver.Evolve();
			}
			catch (std::invalid_argument& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (std::runtime_error& ex)
			{
				std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
			}
			catch (...)
			{
				std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
			}
		}
	}
}

void InscribedCircleCalculatorVisualization()
{
	// Data (incomplete circle)
	pmp::ManifoldCurve2D targetCurve = pmp::CurveFactory::circle(pmp::Point2(0.0f, 0.0f), 0.75f, 16);
	auto targetPts = targetCurve.positions();
	targetPts.erase(targetPts.begin());
	targetPts.erase(targetPts.begin());
	targetPts.erase(targetPts.begin());
	InscribedCircleInputData inputData;
	inputData.Points = targetPts;

	const auto pointBBox = pmp::BoundingBox2(inputData.Points);
	const auto pointBBoxSize = pointBBox.max() - pointBBox.min();
	const float minSize = std::min(pointBBoxSize[0], pointBBoxSize[1]);
	const float cellSize = minSize / 20.0f;
	const SDF::PointCloudDistanceField2DSettings sdfSettings{
		cellSize,
		0.5f,
		DBL_MAX
	};
	inputData.DistanceField = std::make_shared<Geometry::ScalarGrid2D>(SDF::PlanarPointCloudDistanceFieldGenerator::Generate(inputData.Points, sdfSettings));
	// export png and gdim2d
	ExportScalarGridDimInfo2D(dataOutPath + "incompleteCircleDF.gdim2d", *inputData.DistanceField);
	constexpr double colorMapPlotScaleFactor = 1.0; // scale the distance field color map down to show more detail
	ExportScalarGrid2DToPNG(dataOutPath + "incompleteCircleDF.png", *inputData.DistanceField,
		Geometry::BilinearInterpolateScalarValue, 
		//Geometry::GetNearestNeighborScalarValue2D,
		10, 10, RAINBOW_TO_WHITE_MAP * colorMapPlotScaleFactor);

	//const auto [dfMin, dfMax] = std::ranges::minmax_element(inputData.DistanceField->Values());
	//constexpr double colorLegendScaleFactor = 2.0; // stretch color map keys (ratios) to fit the gradient to the full legend bar
	//constexpr double dfValueClampFactor = colorMapPlotScaleFactor / colorLegendScaleFactor; // value cutoff for the chosen color legend.
	//constexpr bool verticalLegend = true;
	//constexpr unsigned int resolutionFactor = 4;
	//constexpr unsigned int legendPxHeight = (verticalLegend ? 400 : 100) * resolutionFactor;
	//constexpr unsigned int legendPxWidth = (verticalLegend ? 100 : 600) * resolutionFactor;
	//ExportColorScaleToPNG(
	//	dataOutPath + "incompleteCircleDF_Scale.png",
	//	(*dfMin),
	//	dfValueClampFactor * (*dfMax),
	//	RAINBOW_TO_WHITE_MAP * colorLegendScaleFactor,
	//	legendPxHeight, legendPxWidth);

	// Distance field (brute force)
	std::cout << "DistanceFieldInscribedCircleCalculator::Calculate ... ";
	DistanceFieldInscribedCircleCalculator dfCalculator;
	const auto dfCircles = dfCalculator.Calculate(inputData);
	std::cout << "done.\n";

	// Hierarchical
	std::cout << "HierarchicalDistanceFieldInscribedCircleCalculator::Calculate ... ";
	HierarchicalDistanceFieldInscribedCircleCalculator hCalculator;
	std::ofstream hCalculatorFile(dataOutPath + "incompleteCircleDF_HierarchicalLog.json");
	if (!hCalculatorFile)
		throw std::runtime_error("Failed to open file \"incompleteCircleDF_HierarchicalLog.txt\" for writing.");
	hCalculator.SetOutputStream(hCalculatorFile);
	const auto hCircles = hCalculator.Calculate(inputData);
	hCalculatorFile.close();
	std::cout << "done.\n";

	// Particle swarm
	std::cout << "ParticleSwarmDistanceFieldInscribedCircleCalculator::Calculate ... ";
	ParticleSwarmDistanceFieldInscribedCircleCalculator psCalculator;
	std::ofstream psCalculatorFile(dataOutPath + "incompleteCircleDF_ParticleSwarmLog.json");
	if (!psCalculatorFile)
		throw std::runtime_error("Failed to open file \"incompleteCircleDF_ParticleSwarmLog.txt\" for writing.");
	psCalculator.SetOutputStream(psCalculatorFile);
	const auto psCircles = psCalculator.Calculate(inputData);
	psCalculatorFile.close();
	std::cout << "done.\n";

	// =============================================================
	// Ellipsoid

	{
		Geometry::IcoSphereBuilder icoBuilder({ 2, 1.0f });
		icoBuilder.BuildBaseData();
		icoBuilder.BuildPMPSurfaceMesh();
		auto sphereMesh = icoBuilder.GetPMPSurfaceMeshResult();
		constexpr float a = 1.0f;
		constexpr float b = 1.5f;
		constexpr float c = 2.0f;
		const auto scalingMat = scaling_matrix(pmp::vec3{ a, b, c });
		sphereMesh *= scalingMat;

		const auto points = sphereMesh.positions();
		const auto pointBBox3D = pmp::BoundingBox(points);
		const auto pointBBox3DSize = pointBBox3D.max() - pointBBox3D.min();
		const float minSize3D = std::min({ pointBBox3DSize[0], pointBBox3DSize[1], pointBBox3DSize[2] });
		const float cellSize3D = minSize3D / 20.0f;

		const SDF::PointCloudDistanceFieldSettings ellipsoidDfSettings{
			cellSize3D,
			0.5f,
			DBL_MAX
		};
		const auto dfEllipsoid = std::make_shared<Geometry::ScalarGrid>(SDF::PointCloudDistanceFieldGenerator::Generate(points, ellipsoidDfSettings));

		ExportToVTI(dataOutPath + "\\EllipsoidSampling_DF", *dfEllipsoid);
	}
	// =============================================================
}

void StandardMeshesExportWithNormals()
{
	const std::vector<std::string> meshForPtCloudNames{
		"bunnyWNormals",
		"maxPlanckWNormals",
	};

	const std::map<std::string, std::vector<Sphere3D>> cutSpheres{
		{"bunnyWNormals", std::vector{ Sphere3D{pmp::Point{-0.01f, 0.06f, 0.012f}, 0.032f}, Sphere3D{pmp::Point{0.01f, 0.12f, 0.01f}, 0.025f}/**/}},
		{"maxPlanckWNormals", std::vector{ Sphere3D{pmp::Point{8.0f, 85.0f, 0.0f}, 50.0f}, Sphere3D{pmp::Point{30.0f, -120.0f, 160.0f}, 100.0f} /**/}},
	};

	constexpr size_t samplingLevel = 3;
	constexpr size_t nSamplings = 10;
	constexpr size_t minVerts = 9; // Minimum number of vertices to sample

	constexpr unsigned int seed = 5000; // seed for the pt cloud sampling RNG

	for (const auto& meshName : meshForPtCloudNames)
	{
		// =======================================================================
		//   - - - - - - - - - - - - -   Data   Prep   - - - - - - - - - - - -
		// -----------------------------------------------------------------------

		std::cout << "==================================================================\n";
		std::cout << "Mesh To Pt Cloud: " << meshName << ".obj -> " << meshName << "Pts_" << samplingLevel << ".ply\n";
		std::cout << "------------------------------------------------------------------\n";
		const auto baseDataOpt = Geometry::ImportOBJMeshGeometryData(dataDirPath + meshName + ".obj", false);
		if (!baseDataOpt.has_value())
		{
			std::cerr << "baseDataOpt == nullopt!\n";
			break;
		}
		std::cout << "meshName.obj" << " imported as BaseMeshGeometryData.\n";
		const auto& baseData = baseDataOpt.value();
		const size_t maxVerts = baseData.Vertices.size(); // Maximum number of vertices available
		size_t nVerts = minVerts + (maxVerts - minVerts) * samplingLevel / (nSamplings - 1);
		nVerts = std::max(minVerts, std::min(nVerts, maxVerts));

		std::cout << "Sampling " << nVerts << "/" << maxVerts << " vertices...\n";

		auto orientedPtCloud = SamplePointsWithNormalsFromMeshData(baseData, nVerts, seed);
		const auto ptCloudName = meshName + "PtsWNormals_" + std::to_string(samplingLevel);

		if (!cutSpheres.contains(meshName))
		{
			std::cerr << "!cutSpheres.contains(\"" << meshName << "\") ... skipping.\n";
			continue;
		}

		DeletePointsWithNormalsContainedInSpheres(orientedPtCloud, cutSpheres.at(meshName));

		std::cout << "Cutting " << nVerts - orientedPtCloud.size() << " / " << nVerts << " points ... \n";

		std::string filename = dataOutPath + meshName + "Pts_" + std::to_string(samplingLevel) + ".ply";
		Geometry::BaseMeshGeometryData ptData;
        ptData.Vertices.resize(orientedPtCloud.size());
        ptData.VertexNormals.resize(orientedPtCloud.size());
        // Extract vertices and vertex normals from orientedPtCloud
        std::ranges::transform(orientedPtCloud, ptData.Vertices.begin(), [](const auto& pair) { return pair.first; });
        std::ranges::transform(orientedPtCloud, ptData.VertexNormals.begin(), [](const auto& pair) { return pair.second; });

		if (!ExportPointsToPLY(ptData, filename))
		{
			std::cerr << "ExportPointsToPLY failed!\n";
			break;
		}
	}
}

void ExportBallPivotingBoundaryEdges()
{
	const std::vector<std::string> meshNames{
		"bunnyPts_3_BallPivoting",
		"maxPlanckPts_3_BallPivoting10",
	};

	for (const auto& meshName : meshNames)
	{
		std::cout << "==================================================================\n";
		std::cout << "Mesh: " << meshName << ".obj\n";
		std::cout << "------------------------------------------------------------------\n";
		
		pmp::SurfaceMesh mesh;
		//mesh.read(dataDirPath + meshName + ".obj");
		// TODO: implement for BaseMeshGeometryData because this data is non-manifold

		//if (!Geometry::ExportBoundaryEdgesToPLY(mesh, dataOutPath + meshName + "_boundaries.ply"))
		//{
		//	std::cerr << "Error while exporting " << (dataOutPath + meshName + "_boundaries.ply") << "!\n";
		//	continue;
		//}

		std::cout << "meshName.obj" << " imported as BaseMeshGeometryData.\n";
	}
}

void TestProblematicMedialAxisPtClouds()
{
	const std::vector unevenCrossPolyPts{
		pmp::Point2{39.507142f, 14.544772f},
		pmp::Point2{46.104261f, 5.542495f},
		pmp::Point2{61.36143f, 4.906308f},
		pmp::Point2{68.282948f, 13.11281f},
		pmp::Point2{66.153916f, 31.426095f},
		pmp::Point2{69.924933f, 39.754365f},
		pmp::Point2{111.082270f, 39.723965f},
		pmp::Point2{117.795930f, 46.528945f},
		pmp::Point2{117.765430f, 66.419005f},
		pmp::Point2{113.358230f, 72.215125f},
		pmp::Point2{89.030514f, 71.743755f},
		pmp::Point2{82.788235f, 77.337235f},
		pmp::Point2{87.565954f, 122.613040f},
		pmp::Point2{80.332640f, 129.222760f},
		pmp::Point2{68.372952f, 128.366110f},
		pmp::Point2{61.451434f, 122.392570f},
		pmp::Point2{57.089489f, 85.394595f},
		pmp::Point2{36.929786f, 84.297475f},
		pmp::Point2{10.835265f, 83.544745f},
		pmp::Point2{3.558908f, 76.519305f},
		pmp::Point2{3.558908f, 57.450225f},
		pmp::Point2{10.584357f, 47.664785f},
		pmp::Point2{32.413427f, 48.166595f},
		pmp::Point2{40.865615f, 41.985195f}
	};
	pmp::ManifoldCurve2D unevenCrossCurve = pmp::CurveFactory::sampled_polygon(unevenCrossPolyPts, 100, false);
	if (!pmp::write_to_ply(unevenCrossCurve, dataOutPath + "unevenCrossCurve.ply"))
		std::cerr << "Error writing unevenCrossCurve.ply!\n";

	const auto unevenCrossPts = unevenCrossCurve.positions();

	const std::vector<Circle2D> testInnerCircles{
		{ pmp::Point2{ 60.322582f, 63.604839f }, 6.6532254f },
		{ pmp::Point2{ 87.733871f, 55.975807f }, 9.08871 },
		{ pmp::Point2{ 62.3629f, 61.741936f }, 19.161289f },
		{ pmp::Point2{ 73.274193f, 103.87903f }, 8.5161285f },
		{ pmp::Point2{ 24.749998f, 67.330643f }, 10.733871f },
		{ pmp::Point2{ 53.935482f, 24.129032f }, 7.6290321f },
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;
	constexpr unsigned int NTimeSteps = 400;
	const double fieldIsoLevel = defaultOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(5.0);

	for (size_t unevenCrossId = 0; const auto & innerCircle : testInnerCircles)
	{
		std::cout << "Setting up ManifoldEvolutionSettings.\n";

		ManifoldEvolutionSettings strategySettings;
		strategySettings.UseInnerManifolds = true;
		//strategySettings.AdvectionInteractWithOtherManifolds = true;
		//strategySettings.OuterManifoldEpsilon = [](double distance)
		//{
		//	return 0.00; // * (1.0 - exp(-distance * distance / 1.0));
		//};
		//strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
		//{
		//	return 0.00; // * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		//};
		strategySettings.InnerManifoldEpsilon = [](double distance)
		{
			return 0.0025 * TRIVIAL_EPSILON(distance);
		};
		strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
		{
			return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		};
		strategySettings.TimeStep = defaultTimeStep;
		strategySettings.LevelOfDetail = 4;
		strategySettings.TangentialVelocityWeight = 0.05;

		strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.22f;
		strategySettings.RemeshingSettings.UseBackProjection = false;

		strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2f;
		strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0f * static_cast<float>(M_PI_2);

		strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
		strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;
		strategySettings.FieldSettings.FieldExpansionFactor = 0.5f;

		std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

		GlobalManifoldEvolutionSettings globalSettings;
		globalSettings.NSteps = NTimeSteps;
		globalSettings.DoRemeshing = true;
		globalSettings.DetectFeatures = false;
		globalSettings.ExportPerTimeStep = true;
		globalSettings.ExportTargetDistanceFieldAsImage = true;
		globalSettings.ProcedureName = "unevenCross" + std::to_string(unevenCrossId);
		globalSettings.OutputPath = dataOutPath;
		globalSettings.ExportResult = false;

		globalSettings.RemeshingResizeFactor = 0.7f;
		globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

		const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;

		std::vector<pmp::ManifoldCurve2D> innerCurves{
			pmp::CurveFactory::circle(innerCircle.Center, innerCircle.Radius, nSegments)
		};
		innerCurves[0].negate_orientation();

		auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
			strategySettings, std::nullopt, innerCurves, 
			std::make_shared<std::vector<pmp::Point2>>(unevenCrossPts));

		std::cout << "Setting up ManifoldEvolver.\n";

		ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

		std::cout << "ManifoldEvolver::Evolve ... ";

		try
		{
			evolver.Evolve();
		}
		catch (std::invalid_argument& ex)
		{
			std::cerr << "> > > > > > > > > > > > > > std::invalid_argument: " << ex.what() << " Continue... < < < < < \n";
		}
		catch (std::runtime_error& ex)
		{
			std::cerr << "> > > > > > > > > > > > > > std::runtime_error: " << ex.what() << " Continue... < < < < < \n";
		}
		catch (...)
		{
			std::cerr << "> > > > > > > > > > > > > > ManifoldEvolver::Evolve has thrown an exception! Continue... < < < < < \n";
		}

		unevenCrossId++;
	}
}
