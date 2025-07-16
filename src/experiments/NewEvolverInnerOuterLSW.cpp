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
#include "geometry/GeometryIOUtils.h"

#include "sdf/SDF.h"

#include "utils/TimingUtils.h"
#include "utils/NumericalUtils.h"

#include "core/ConversionUtils.h"
#include "core/EvolverUtilsCommon.h"
#include "core/InscribedManifold.h"
#include "core/ManifoldEvolver.h"
#include "core/SurfaceEvolver.h"
#include "core/CharacteristicsBuilder.h"

#include "IOEnvironment.h"
#include "PolygonalDatasets.h"

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
	const pmp::Scalar minSize = std::min({ meshBBoxSize[0], meshBBoxSize[1], meshBBoxSize[2] });
	const pmp::Scalar cellSize = minSize / 10.0;

	const SDF::DistanceFieldSettings sdfSettings{
		cellSize,
		1.0,
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
		1.0,
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
			1.0,
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
		const pmp::Scalar minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
		const pmp::Scalar cellSize = minSize / 10.0;

		const SDF::DistanceField2DSettings sdfSettings{
			cellSize,
			1.0,
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
		const pmp::Scalar minSize = std::min(curveBBoxSize[0], curveBBoxSize[1]);
		const pmp::Scalar cellSize = minSize / 10.0;

		const SDF::DistanceField2DSettings sdfSettings{
			cellSize,
			1.0,
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
		const pmp::Scalar minSize = std::min(pointBBoxSize[0], pointBBoxSize[1]);
		const pmp::Scalar cellSize = minSize / 10.0;

		const SDF::PointCloudDistanceField2DSettings sdfSettings{
			cellSize,
			1.0,
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
		auto ellipse = pmp::CurveFactory::circle(pmp::Point2(0, 0), 1.0, 32);
		ellipse *= scaling_matrix_2d(pmp::vec2(2.0, 1.0));

		const auto pointBBox = pmp::BoundingBox2(ellipse.positions());
		const auto pointBBoxSize = pointBBox.max() - pointBBox.min();
		const pmp::Scalar minSize = std::min(pointBBoxSize[0], pointBBoxSize[1]);
		const pmp::Scalar cellSize = minSize / 20.0;

		const SDF::PointCloudDistanceField2DSettings sdfSettings{
			cellSize,
			1.0,
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

	if (!Geometry::ExportManifoldCurve2DToPLY(closedArc, dataOutPath + "closedArc.ply"))
	{
		std::cerr << "ExportManifoldCurve2DToPLY: internal error!\n";
	}

	const auto openArc = pmp::CurveFactory::circle(pmp::Point2(0, 0), 1.0, 16, 0, M_PI);

	if (!Geometry::ExportManifoldCurve2DToPLY(openArc, dataOutPath + "openArc.ply"))
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
		{"armadillo", pmp::Point{-0.10348621158928151, 21.427067319905646, 9.79369240592005}},
		{"bunny", pmp::Point{-0.01684039831161499, 0.11015420407056808, 0.0012007840834242693} },
		{"maxPlanck", pmp::Point{30.59686279296875, -18.105804443359375, 82.29149055480957} },
		{"nefertiti", pmp::Point{0.0, 0.0, 0.0} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"armadillo", pmp::vec3{-0.03070969905335075, 0.12876712096541565, 0.9911992448253433}},
		{"bunny", pmp::vec3{0.0, 0.0, 1.0} },
		{"maxPlanck", pmp::vec3{1.0, 0.0, 0.0} },
		{"nefertiti", pmp::vec3{1.0, 0.0, 0.0} }
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
		const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		//const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar distTolerance = 0.01 * minSize;

		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0, 0.0, 0.0 });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(pts3D, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}

		if (!Geometry::Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
		{
			std::cerr << "Export2DPointCloudToPLY: internal error during export!\n";
			continue;
		}
	}
}


void CurveReorientTests()
{
	auto circleCurve = pmp::CurveFactory::circle(pmp::Point2{}, 1.0, 25);
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

void DeletePointsContainedInCircles(std::vector<pmp::Point2>& pointsToFilter, const std::vector<Geometry::Circle2D>& cutCircles)
{
	auto isPointInCircle = [](const pmp::Point2& point, const Geometry::Circle2D& circle) {
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

void DeletePointsContainedInSpheres(std::vector<pmp::Point>& pointsToFilter, const std::vector<Geometry::Sphere3D>& cutSpheres)
{
	auto isPointInSphere = [](const pmp::Point& point, const Geometry::Sphere3D& sphere) {
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

void DeletePointsWithNormalsContainedInSpheres(std::vector<std::pair<pmp::Point, pmp::vec3>>& pointsToFilter, const std::vector<Geometry::Sphere3D>& cutSpheres)
{
	auto isPointInSphere = [](const pmp::Point& point, const Geometry::Sphere3D& sphere) {
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
		{"armadillo", pmp::Point{-0.10348621158928151, 21.427067319905646, 9.79369240592005}},
		{"bunny", pmp::Point{-0.01684039831161499, 0.11015420407056808, 0.0012007840834242693} },
		{"maxPlanck", pmp::Point{30.59686279296875, -18.105804443359375, 82.29149055480957} },
		{"nefertiti", pmp::Point{0.0, 0.0, 0.0} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"armadillo", pmp::vec3{-0.03070969905335075, 0.12876712096541565, 0.9911992448253433}},
		{"bunny", pmp::vec3{0.0, 0.0, 1.0} },
		{"maxPlanck", pmp::vec3{1.0, 0.0, 0.0} },
		{"nefertiti", pmp::vec3{1.0, 0.0, 0.0} }
	};

	const std::map<std::string, Geometry::Circle2D> outerCircles{
		{"armadillo", Geometry::Circle2D{pmp::Point2{0.372234, 16.6515}, 121.558} },
		{"bunny", Geometry::Circle2D{pmp::Point2{-0.0155906, 0.102261}, 0.142831} },
		{"maxPlanck", Geometry::Circle2D{pmp::Point2{-17.82, 82.5006}, 292.263} },
		{"nefertiti", Geometry::Circle2D{pmp::Point2{0.178497, -0.0410004}, 441.436} }
	};

	const std::map<std::string, Geometry::Sphere3D> outerSpheres{
		{"armadillo", Geometry::Sphere3D{pmp::Point{0.0122509, 21.4183, -0.000249863}, 136.963} },
		{"bunny", Geometry::Sphere3D{pmp::Point{-0.0168297, 0.110217, -0.0015718}, 0.141622} },
		{"maxPlanck", Geometry::Sphere3D{pmp::Point{30.658, -17.9765, 82.2885}, 271.982} },
		{"nefertiti", Geometry::Sphere3D{pmp::Point{0.0144997, -0.00499725, -0.0215073}, 392.184} }
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
		const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
		constexpr pmp::Scalar volExpansionFactor = 1.0;
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
			topoParams.MinEdgeMultiplier = 0.14;
			topoParams.UseBackProjection = false;
			topoParams.PrincipalCurvatureFactor = 3.2;
			topoParams.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);
			topoParams.EdgeLengthDecayFactor = 0.7;
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
				0.05,
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

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			std::cout << "Setting up ManifoldSurfaceEvolutionStrategy.\n";

			// Set up the evolution strategy
			auto surfaceStrategy = std::make_shared<ManifoldSurfaceEvolutionStrategy>(
				strategySettings, MeshLaplacian::Voronoi,
				std::make_shared<std::vector<pmp::Point>>(ptCloud));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver newEvolver(globalSettings, std::move(surfaceStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			newEvolver.Evolve();
		}

		// ==========================================================================
		// - - - - - - - - -  New Manifold Evolver (Curve)  - - - - - - - - - - - - 
		// ==========================================================================

		const pmp::Scalar distTolerance = 0.01 * minSize;
		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0, 0.0, 0.0 });
		auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}
		const std::vector cutCircles{ 
			Geometry::Circle2D{pmp::Point2{-0.02, 0.04}, 0.015 },
			Geometry::Circle2D{pmp::Point2{0.03, 0.04}, 0.015 } };
		DeletePointsContainedInCircles(pts2D, cutCircles);

		if (!Geometry::Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
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

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			auto curveStrategy = std::make_shared<ManifoldCurveEvolutionStrategy>(
				strategySettings, std::make_shared<std::vector<pmp::Point2>>(pts2D));

			std::cout << "Setting up ManifoldEvolver.\n";

			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			evolver.Evolve();
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

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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
				1.0, 0.0, 0.0, outerSphere.Center[0],
				0.0, 1.0, 0.0, outerSphere.Center[1],
				0.0, 0.0, 1.0, outerSphere.Center[2],
				0.0, 0.0, 0.0, 1.0
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

			evolver.Evolve();
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

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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

			evolver.Evolve();
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
		{"armadillo", pmp::Point{-0.10348621158928151, 21.427067319905646, 9.79369240592005}},
		{"bunny", pmp::Point{-0.01684039831161499, 0.11015420407056808, 0.0012007840834242693} },
		{"maxPlanck", pmp::Point{30.59686279296875, -18.105804443359375, 82.29149055480957} },
		{"nefertiti", pmp::Point{0.0, 0.0, 0.0} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"armadillo", pmp::vec3{-0.03070969905335075, 0.12876712096541565, 0.9911992448253433}},
		{"bunny", pmp::vec3{0.0, 0.0, 1.0} },
		{"maxPlanck", pmp::vec3{1.0, 0.0, 0.0} },
		{"nefertiti", pmp::vec3{1.0, 0.0, 0.0} }
	};

	const std::map<std::string, Geometry::Circle2D> outerCircles{
		{"armadillo", Geometry::Circle2D{pmp::Point2{0.372234, 16.6515}, 121.558} },
		{"bunny", Geometry::Circle2D{pmp::Point2{-0.0155906, 0.102261}, 0.142831} },
		{"maxPlanck", Geometry::Circle2D{pmp::Point2{-17.82, 82.5006}, 292.263} },
		{"nefertiti", Geometry::Circle2D{pmp::Point2{0.178497, -0.0410004}, 441.436} }
	};
	const std::map<std::string, Geometry::Circle2D> innerCircles{
		{"armadillo", Geometry::Circle2D{pmp::Point2{-3.0, 52.0}, 20.0}},
		{"bunny", Geometry::Circle2D{pmp::Point2{-0.025, 0.08}, 0.025}},
		{"maxPlanck", Geometry::Circle2D{pmp::Point2{8.0, 85.0}, 50.0}},
		{"nefertiti", Geometry::Circle2D{pmp::Point2{-20.0, 100.0}, 55.0}}
	};

	const std::map<std::string, Geometry::Sphere3D> outerSpheres{
		{"armadillo", Geometry::Sphere3D{pmp::Point{0.0122509, 21.4183, -0.000249863}, 136.963} },
		{"bunny", Geometry::Sphere3D{pmp::Point{-0.0168297, 0.110217, -0.0015718}, 0.141622} },
		{"maxPlanck", Geometry::Sphere3D{pmp::Point{30.658, -17.9765, 82.2885}, 271.982} },
		{"nefertiti", Geometry::Sphere3D{pmp::Point{0.0144997, -0.00499725, -0.0215073}, 392.184} }
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
		const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		// const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;

		const double isoLvlOffsetFactor = (isoLevelOffsetFactors.contains(ptCloudName) ? isoLevelOffsetFactors.at(ptCloudName) : defaultOffsetFactor);
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		const double tau = (timeStepSizesForPtClouds.contains(ptCloudName) ? timeStepSizesForPtClouds.at(ptCloudName) : defaultTimeStep); // time step

		// ==========================================================================
		// - - - - - - - - -  New Manifold Evolver (Curve)  - - - - - - - - - - - - 
		// ==========================================================================

		const pmp::Scalar distTolerance = 0.01 * minSize;
		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0, 0.0, 0.0 });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}

		if (!Geometry::Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
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

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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
				1.0, 0.0, 0.0, outerSphere.Center[0],
				0.0, 1.0, 0.0, outerSphere.Center[1],
				0.0, 0.0, 1.0, outerSphere.Center[2],
				0.0, 0.0, 0.0, 1.0
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

			evolver.Evolve();
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

			//strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.22;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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

			evolver.Evolve();
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
		{"armadillo", pmp::Point{-0.10348621158928151, 21.427067319905646, 9.79369240592005}},
		{"bunny", pmp::Point{-0.01684039831161499, 0.11015420407056808, 0.0012007840834242693} },
		{"maxPlanck", pmp::Point{30.59686279296875, -18.105804443359375, 82.29149055480957} },
		{"nefertiti", pmp::Point{0.0, 0.0, 0.0} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"armadillo", pmp::vec3{-0.03070969905335075, 0.12876712096541565, 0.9911992448253433}},
		{"bunny", pmp::vec3{0.0, 0.0, 1.0} },
		{"maxPlanck", pmp::vec3{1.0, 0.0, 0.0} },
		{"nefertiti", pmp::vec3{1.0, 0.0, 0.0} }
	};

	const std::map<std::string, Geometry::Circle2D> outerCircles{
		{"armadillo", Geometry::Circle2D{pmp::Point2{0.372234, 16.6515}, 121.558} },
		{"bunny", Geometry::Circle2D{pmp::Point2{-0.0155906, 0.102261}, 0.142831} },
		{"maxPlanck", Geometry::Circle2D{pmp::Point2{-17.82, 82.5006}, 292.263} },
		{"nefertiti", Geometry::Circle2D{pmp::Point2{0.178497, -0.0410004}, 441.436} }
	};
	const std::map<std::string, Geometry::Circle2D> innerCircles{
		{"armadillo", Geometry::Circle2D{pmp::Point2{-3.0, 52.0}, 20.0}},
		{"bunny", Geometry::Circle2D{pmp::Point2{-0.025, 0.08}, 0.025}},
		{"maxPlanck", Geometry::Circle2D{pmp::Point2{8.0, 85.0}, 50.0}},
		{"nefertiti", Geometry::Circle2D{pmp::Point2{-20.0, 100.0}, 55.0}}
	};

	const std::map<std::string, Geometry::Sphere3D> outerSpheres{
		{"armadillo", Geometry::Sphere3D{pmp::Point{0.0122509, 21.4183, -0.000249863}, 136.963} },
		{"bunny", Geometry::Sphere3D{pmp::Point{-0.0168297, 0.110217, -0.0015718}, 0.141622} },
		{"maxPlanck", Geometry::Sphere3D{pmp::Point{30.658, -17.9765, 82.2885}, 271.982} },
		{"nefertiti", Geometry::Sphere3D{pmp::Point{0.0144997, -0.00499725, -0.0215073}, 392.184} }
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
		const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;

		const double isoLvlOffsetFactor = (timeStepSizesForPtClouds.contains(ptCloudName) ? isoLevelOffsetFactors.at(ptCloudName) : defaultOffsetFactor);
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		const double tau = (timeStepSizesForPtClouds.contains(ptCloudName) ? timeStepSizesForPtClouds.at(ptCloudName) : defaultTimeStep); // time step

		// ==========================================================================
		// - - - - - - - - -  New Manifold Evolver (Curve)  - - - - - - - - - - - - 
		// ==========================================================================

		const pmp::Scalar distTolerance = 0.01 * minSize;
		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0, 0.0, 0.0 });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}

		if (!Geometry::Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
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

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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
				1.0, 0.0, 0.0, outerSphere.Center[0],
				0.0, 1.0, 0.0, outerSphere.Center[1],
				0.0, 0.0, 1.0, outerSphere.Center[2],
				0.0, 0.0, 0.0, 1.0
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

			evolver.Evolve();
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

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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

			const auto& outerCircle = outerCircles.at(meshName);
			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);

			const auto& innerCircle = innerCircles.at(meshName);
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

			evolver.Evolve();
		}
	}
}


void OutwardEvolvingInnerCircleTest()
{
	const Geometry::Circle2D innerTestCircle{ pmp::Point2{-3.0, 52.0}, 20.0 };

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
	const pmp::Scalar minSize = std::min(bboxSize[0], bboxSize[1]);
	const pmp::Scalar maxSize = std::max(bboxSize[0], bboxSize[1]);
	const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;

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
		strategySettings.LevelOfDetail = 4;
		strategySettings.TangentialVelocityWeight = 0.05;

		strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.22;
		strategySettings.RemeshingSettings.UseBackProjection = false;

		strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
		strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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
		globalSettings.ProcedureName = "singleInnerCircleTestPhaseTwo";
		globalSettings.OutputPath = dataOutPath;
		globalSettings.ExportResult = false;

		globalSettings.RemeshingResizeFactor = 0.7;
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

		evolver.Evolve();
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
		pmp::Point2{0.0, 0.0},
		pmp::Point2{1.0, 0.0},
		pmp::Point2{1.0, 1.0},
		pmp::Point2{0.0, 1.0}
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
		pmp::Point2{0.0, 0.0},
		pmp::Point2{1.0, 0.0},
		pmp::Point2{1.5, 0.5},
		pmp::Point2{2.0, 0.0}
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
		pmp::Point2{0.0, 1.0},
		pmp::Point2{0.866, 0.5},
		pmp::Point2{0.866, -0.5},
		pmp::Point2{0.0, -1.0},
		pmp::Point2{-0.866, -0.5},
		pmp::Point2{-0.866, 0.5}
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
		{ "circle", pmp::CurveFactory::circle(pmp::Point2{-3.0, 52.0}, 35.0, 25, 0.0, 2.0 * M_PI)},
		{ "incompleteCircle", pmp::CurveFactory::circle(pmp::Point2{-3.0, 52.0}, 35.0, 25, M_PI_2, 2.0 * M_PI)},
		{ "sineDeformedCircle", pmp::CurveFactory::sine_deformed_circle(pmp::Point2{-3.0, 52.0}, 35.0, 25, 7.0, 4.0, 0.0, 2.0 * M_PI)},
		{ "sineDeformedIncompleteCircle", pmp::CurveFactory::sine_deformed_circle(pmp::Point2{-3.0, 52.0}, 35.0, 25, 7.0, 4.0, M_PI_2, 2.0 * M_PI)},
		{ "chamferedRectangle", pmp::CurveFactory::rectangle(pmp::Point2{-3.0, 52.0}, 60.0, 70.0, 15, true)},
		{ "incompleteChamferedRectangle", pmp::CurveFactory::sampled_polygon({
			pmp::Point2{-30.0, -35.0} + pmp::Point2{-3.0, 52.0},
			pmp::Point2{30.0, -35.0} + pmp::Point2{-3.0, 52.0},
			pmp::Point2{30.0, 35.0} + pmp::Point2{-3.0, 52.0},
			pmp::Point2{-30.0, 35.0} + pmp::Point2{-3.0, 52.0}}, 30, true, false)},
		{ "chamferedTriangle", pmp::CurveFactory::sampled_polygon({
			pmp::Point2{-0.5, -(pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)6.0} *120.0 + pmp::Point2{-3.0, 52.0},
			pmp::Point2{0.5, -(pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)6.0} *120.0 + pmp::Point2{-3.0, 52.0},
			pmp::Point2{0.0, (pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)3.0} *120.0 + pmp::Point2{-3.0, 52.0}}, 30, true)},
		{ "incompleteChamferedTriangle", pmp::CurveFactory::sampled_polygon({
			pmp::Point2{-0.5, -(pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)6.0} *120.0 + pmp::Point2{-3.0, 52.0},
			pmp::Point2{0.5, -(pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)6.0} *120.0 + pmp::Point2{-3.0, 52.0},
			pmp::Point2{0.0, (pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)3.0} *120.0 + pmp::Point2{-3.0, 52.0}}, 30, true, false)}
	};
	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;
	constexpr unsigned int NTimeSteps = 180;
	const auto innerTestCircle = Geometry::Circle2D{ pmp::Point2{ -3.0, 52.0 }, 20.0 };

	for (const auto& [ptCloudName, curve] : targetCurves)
	{
		std::cout << "========================================================\n";
		std::cout << "         inner circle LSW for: " << ptCloudName << " ... \n";
		std::cout << " ------------------------------------------------------ \n";
		const auto& pts2D = curve.positions();
		const pmp::BoundingBox2 bbox{ pts2D };

		const auto bboxSize = bbox.max() - bbox.min();
		const pmp::Scalar minSize = std::min(bboxSize[0], bboxSize[1]);
		//const pmp::Scalar maxSize = std::max(bboxSize[0], bboxSize[1]);
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;

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

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.22;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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

			evolver.Evolve();
		}
	}
}


void ConcentricCirclesTests()
{
	// Define the inner and outer circle pairs directly
	const std::vector<std::pair<Geometry::Circle2D, Geometry::Circle2D>> circlePairs{
		{Geometry::Circle2D{pmp::Point2{-3.0, 52.0}, 100.0}, Geometry::Circle2D{pmp::Point2{-3.0, 52.0}, 121.558}},
		{Geometry::Circle2D{pmp::Point2{-25.0, 8.0}, 0.055}, Geometry::Circle2D{pmp::Point2{-25.0, 8.0}, 0.142831}},
		{Geometry::Circle2D{pmp::Point2{8.0, 85.0}, 50.0}, Geometry::Circle2D{pmp::Point2{8.0, 85.0}, 292.263}},
		{Geometry::Circle2D{pmp::Point2{-20.0, 90.0}, 55.0}, Geometry::Circle2D{pmp::Point2{-20.0, 90.0}, 441.436}}
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
		const pmp::Scalar minSize = std::min(bboxSize[0], bboxSize[1]);
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;

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
				return -1.0 * distance * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = defaultTimeStep;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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

			evolver.Evolve();
		}

		curveId++;
	}
}


void NonConcentricCirclesTest()
{
	// Define the inner and outer circle pairs directly
	const std::vector<std::pair<Geometry::Circle2D, Geometry::Circle2D>> circlePairs{
		{Geometry::Circle2D{pmp::Point2{-3.0, 52.0}, 20.0}, Geometry::Circle2D{pmp::Point2{0.372234, 16.6515}, 121.558}},
		{Geometry::Circle2D{pmp::Point2{-0.025, 0.08}, 0.025}, Geometry::Circle2D{pmp::Point2{-0.0155906, 0.102261}, 0.142831}},
		{Geometry::Circle2D{pmp::Point2{8.0, 85.0}, 50.0}, Geometry::Circle2D{pmp::Point2{-17.82, 82.5006}, 292.263}},
		{Geometry::Circle2D{pmp::Point2{-20.0, 90.0}, 55.0}, Geometry::Circle2D{pmp::Point2{0.178497, -0.0410004}, 441.436}}
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;
	constexpr unsigned int NTimeSteps = 1500;

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
		const pmp::Scalar minSize = std::min(bboxSize[0], bboxSize[1]);
		//const pmp::Scalar maxSize = std::max(bboxSize[0], bboxSize[1]);
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;

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
				//return 1.0 * (1.0 - exp(-distance * distance / 1.0));
				return 0.0;
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				return 0.0025 * (1.0 - exp(-distance * distance / 1.0)); //TRIVIAL_EPSILON(distance);
				//return 0.001 * TRIVIAL_EPSILON(distance);
				// return 0.0;
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				return -1.5 * distance * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = defaultTimeStep;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;
			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.22;

			strategySettings.RemeshingSettings.UseBackProjection = true;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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

			evolver.Evolve();
		}

		curveId++;
	}
}

namespace
{
	void RemeshWithDefaultSettings(pmp::ManifoldCurve2D& curve, const std::shared_ptr<pmp::EvolvingArcLengthCalculator>& calc = nullptr, const pmp::Scalar& factor = 1.0)
	{
		// DISCLAIMER: Some curves from svg paths might have duplicate points
		const pmp::Scalar edgeLength = 2.0 * factor;
		constexpr unsigned int iterations = 10;
		pmp::CurveRemeshing remesher(curve, calc);
		pmp::AdaptiveRemeshingSettings settings;
		settings.MinEdgeLength = edgeLength;
		settings.MaxEdgeLength = 1.5 * edgeLength * factor;
		settings.ApproxError = 0.05 * edgeLength * factor;
		settings.NRemeshingIterations = iterations;
		settings.NTangentialSmoothingIters = 6;
		settings.UseProjection = false;
		remesher.adaptive_remeshing(settings);
	}

	[[nodiscard]] pmp::ManifoldCurve2D GetCurveRotatedAboutCenterPoint(const pmp::ManifoldCurve2D& curve, const pmp::Scalar& angle, const bool& remesh = false)
	{
		if (curve.n_vertices() == 0)
		{
			return curve;
		}

		const auto curveBBox = curve.bounds();
		const auto center = curveBBox.center();
		const auto rotationMatrix = pmp::rotation_matrix_2d(center, angle);
		pmp::ManifoldCurve2D resultCurve{ curve };
		resultCurve *= rotationMatrix;
		if (remesh)
			RemeshWithDefaultSettings(resultCurve, nullptr, 2.0);
		return resultCurve;
	}
} // anonymous namespace

void EquilibriumPairedManifoldTests()
{
	const auto center = pmp::Point2(200, 400);
	const pmp::Scalar outerRadius = 100.0;
	const pmp::Scalar innerRadius = 40.0;

	const std::vector<pmp::Point2> squareVerticesLarge = {
		pmp::Point2{-50.0, -50.0} + center,
		pmp::Point2{50.0, -50.0} + center,
		pmp::Point2{50.0, 50.0} + center,
		pmp::Point2{-50.0, 50.0} + center
	};
	const std::vector<pmp::Point2> triangleVerticesLarge = {
		pmp::Point2{-0.5, -(pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)6.0} *(2.0 * outerRadius) + center,
		pmp::Point2{0.5, -(pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)6.0} *(2.0 * outerRadius) + center,
		pmp::Point2{0.0, (pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)3.0} *(2.0 * outerRadius) + center
	};
	const std::vector<pmp::Point2> squareVerticesSmall = {
		pmp::Point2{-20.0, -20.0} + center,
		pmp::Point2{20.0, -20.0} + center,
		pmp::Point2{20.0, 20.0} + center,
		pmp::Point2{-20.0, 20.0} + center
	};
	const std::vector<pmp::Point2> triangleVerticesSmall = {
		pmp::Point2{-0.5, -(pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)6.0} *innerRadius + center,
		pmp::Point2{0.5, -(pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)6.0} *innerRadius + center,
		pmp::Point2{0.0, (pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)3.0} *innerRadius + center
	};
	const unsigned int segments = 40;

	//const pmp::Scalar angle = 46.0; // in degrees
	const pmp::Scalar angle = 0.0;
	const double minDistancePercentageEpsilon = 0.0;
	const double minDistancePercentageEta = 0.05;

	constexpr bool logCtrlFunctionValues{ false };
	constexpr bool logEpsilon{ false };
	constexpr bool logErrorValues{ false };

	// List of curve pairs to evolve
	const std::vector<std::pair<pmp::ManifoldCurve2D, pmp::ManifoldCurve2D>> curvePairs{
		// Pair 1: Outer chamfered square and inner circle
		//{pmp::CurveFactory::sampled_polygon(squareVerticesLarge, segments, true),
		// pmp::CurveFactory::circle(center, 0.5 * innerRadius, segments)},

		// Pair 2: Outer chamfered equilateral triangle and inner circle
		//{pmp::CurveFactory::sampled_polygon(triangleVerticesLarge, segments, true),
		//pmp::CurveFactory::circle(center, innerRadius, segments)},

		//// Pair 3: Outer circle and inner chamfered square
		//{pmp::CurveFactory::circle(center, outerRadius, segments),
		//pmp::CurveFactory::sampled_polygon(squareVerticesSmall, segments, true)},

		// Pair 4: Outer circle and inner chamfered equilateral triangle
		{GetCurveRotatedAboutCenterPoint(pmp::CurveFactory::circle(center, outerRadius, segments), angle, false),
		GetCurveRotatedAboutCenterPoint(pmp::CurveFactory::sampled_polygon(triangleVerticesSmall, segments / 2, true), angle, false)}
	};

	// Prepare the settings for the evolver
	ManifoldEvolutionSettings strategySettings;
	strategySettings.UseInnerManifolds = true;

	strategySettings.AdvectionInteractWithOtherManifolds = true;
	// use PreComputeAdvectionDiffusionParams?
	strategySettings.OuterManifoldEpsilon.Bind(minDistancePercentageEpsilon, [](double distance)
		{
			//return 0.0;
			return 0.5 * (1.0 - exp(-distance * distance / 1.0));
		});
	strategySettings.OuterManifoldEta.Bind(minDistancePercentageEta, [](double distance, double negGradDotNormal)
		{
			return -1.0 * distance * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			//return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			//return 1.0 * (1.0 - exp(-distance * distance / 0.5)) * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		});
	strategySettings.InnerManifoldEpsilon.Bind(minDistancePercentageEpsilon, [](double distance)
		{
			//return 0.0;
			return 0.001 * TRIVIAL_EPSILON(distance);
		});
	strategySettings.InnerManifoldEta.Bind(minDistancePercentageEta, [](double distance, double negGradDotNormal)
		{
			return 1.0 * distance * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			//return -1.0 * (1.0 - exp(-distance * distance / 0.5)) * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		});

	strategySettings.TimeStep = 0.05;
	strategySettings.TangentialVelocityWeight = 0.05;
	strategySettings.RemeshingSettings.MinEdgeMultiplier = 1.0;
	strategySettings.RemeshingSettings.UseBackProjection = true;
	strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
	strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = static_cast<pmp::Scalar>(M_PI_2);
	strategySettings.FieldSettings.NVoxelsPerMinDimension = 40;
	//strategySettings.FieldSettings.FieldIsoLevel = 2.0;
	strategySettings.UseLinearGridInterpolation = false;

	strategySettings.ExportVariableScalarFieldsDimInfo = true;
	strategySettings.ExportVariableVectorFieldsDimInfo = true;

	if (logCtrlFunctionValues)
	{
		if (logEpsilon)
		{
			strategySettings.DiagSettings.LogOuterManifoldEpsilon = true;
			strategySettings.DiagSettings.LogInnerManifoldsEpsilon = true;
		}
		else
		{
			strategySettings.DiagSettings.LogOuterManifoldEta = true;
			strategySettings.DiagSettings.LogInnerManifoldsEta = true;
		}
	}
	if (logErrorValues)
	{
		//strategySettings.DiagSettings.LogOuterManifoldErrors = true;
		//strategySettings.DiagSettings.LogInnerManifoldsErrors = true;

		//strategySettings.DiagSettings.LogOuterManifoldXErrors = true;
		//strategySettings.DiagSettings.LogInnerManifoldsXErrors = true;

		strategySettings.DiagSettings.LogOuterManifoldYErrors = true;
		strategySettings.DiagSettings.LogInnerManifoldsYErrors = true;
	}

	// Global settings
	GlobalManifoldEvolutionSettings globalSettings;
	globalSettings.NSteps = 2500;
	//globalSettings.NSteps = 1378;
	//globalSettings.NSteps = 1000;
	//globalSettings.NSteps = 500;
	//globalSettings.NSteps = 20;
	globalSettings.DoRemeshing = true;
	globalSettings.DetectFeatures = false;
	globalSettings.ExportPerTimeStep = true;
	globalSettings.ExportTargetDistanceFieldAsImage = true;
	globalSettings.OutputPath = dataOutPath;
	globalSettings.ExportResult = false;

	globalSettings.RemeshingResizeFactor = 0.7;
	globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

	for (unsigned int pairId = 3; const auto & [outerCurve, innerCurve] : curvePairs)
	{
		std::cout << "EquilibriumPairedManifoldTest for pair: " << pairId << "\n";

		globalSettings.ProcedureName = "equilibriumPair" + std::to_string(pairId);

		const auto innerCurveMinDim = Geometry::GetCurveBoundsMinDimension(innerCurve);
		strategySettings.FieldSettings.FieldIsoLevel = 4.0 * innerCurveMinDim / strategySettings.FieldSettings.NVoxelsPerMinDimension;

		// DEBUG translation
		const pmp::Scalar smallXYShift = strategySettings.FieldSettings.FieldIsoLevel * 0.1;
		const auto smallTranslation = translation_matrix(pmp::vec2{ smallXYShift , (pmp::Scalar)-0.66 * smallXYShift });
		pmp::ManifoldCurve2D outerCurveCopy{ outerCurve };
		outerCurveCopy *= smallTranslation;

		std::vector innerCurves{ innerCurve };
		//innerCurves[0].negate_orientation();
		//auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
		//	strategySettings, outerCurve, innerCurves, nullptr);
		auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
			strategySettings, outerCurveCopy, innerCurves, nullptr);

		ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

		evolver.Evolve();

		pairId++;
	}
}

void EquilibriumPairedConcaveManifoldTests()
{
	//const pmp::Scalar angle = 46.0; // in degrees
	const pmp::Scalar angle = 0.0;

	// Pair 0
	const auto path00Vertices = ParsePolygonalSVGPath(svgPathPair00);
	auto path00Curve = pmp::CurveFactory::sampled_polygon(path00Vertices, 140, false);
	RemeshWithDefaultSettings(path00Curve);
	if (!pmp::write_to_ply(path00Curve, dataOutPath + "path00Curve.ply"))
		std::cerr << "Error writing path00Curve.ply!\n";

	const auto path01Vertices = ParsePolygonalSVGPath(svgPathPair01);
	auto path01Curve = pmp::CurveFactory::sampled_polygon(path01Vertices, 80, false);
	RemeshWithDefaultSettings(path01Curve);
	if (!pmp::write_to_ply(path01Curve, dataOutPath + "path01Curve.ply"))
		std::cerr << "Error writing path01Curve.ply!\n";

	// Pair 1
	const auto path10Vertices = ParsePolygonalSVGPath(svgPathPair10);
	auto path10Curve = pmp::CurveFactory::sampled_polygon(path10Vertices, 140, false);
	RemeshWithDefaultSettings(path10Curve);
	if (!pmp::write_to_ply(path10Curve, dataOutPath + "path10Curve.ply"))
		std::cerr << "Error writing path10Curve.ply!\n";

	const auto path11Vertices = ParsePolygonalSVGPath(svgPathPair11);
	auto path11Curve = pmp::CurveFactory::sampled_polygon(path11Vertices, 80, false);
	RemeshWithDefaultSettings(path11Curve);
	if (!pmp::write_to_ply(path11Curve, dataOutPath + "path11Curve.ply"))
		std::cerr << "Error writing path11Curve.ply!\n";

	// Pair 2
	const auto path20Vertices = ParsePolygonalSVGPath(svgPathPair20);
	auto path20Curve = pmp::CurveFactory::sampled_polygon(path20Vertices, 100, true);
	RemeshWithDefaultSettings(path20Curve);
	if (!pmp::write_to_ply(path20Curve, dataOutPath + "path20Curve.ply"))
		std::cerr << "Error writing path20Curve.ply!\n";

	const auto path21Vertices = ParsePolygonalSVGPath(svgPathPair21);
	auto path21Curve = pmp::CurveFactory::sampled_polygon(path21Vertices, 80, false);
	RemeshWithDefaultSettings(path21Curve);
	if (!pmp::write_to_ply(path21Curve, dataOutPath + "path21Curve.ply"))
		std::cerr << "Error writing path21Curve.ply!\n";

	// Pair 3
	const auto path30Vertices = ParsePolygonalSVGPath(svgPathPair30);
	auto path30Curve = pmp::CurveFactory::sampled_polygon(path30Vertices, 100, true);
	RemeshWithDefaultSettings(path30Curve);
	if (!pmp::write_to_ply(path30Curve, dataOutPath + "path30Curve.ply"))
		std::cerr << "Error writing path30Curve.ply!\n";

	const auto path31Vertices = ParsePolygonalSVGPath(svgPathPair31);
	auto path31Curve = pmp::CurveFactory::sampled_polygon(path31Vertices, 100, false);
	RemeshWithDefaultSettings(path31Curve);
	if (!pmp::write_to_ply(path31Curve, dataOutPath + "path31Curve.ply"))
		std::cerr << "Error writing path31Curve.ply!\n";

	// Pair 4
	const auto path40Vertices = ParsePolygonalSVGPath(svgPathPair40);
	auto path40Curve = pmp::CurveFactory::sampled_polygon(path40Vertices, 100, true);
	RemeshWithDefaultSettings(path40Curve);
	if (!pmp::write_to_ply(path40Curve, dataOutPath + "path40Curve.ply"))
		std::cerr << "Error writing path40Curve.ply!\n";

	const auto path41Vertices = ParsePolygonalSVGPath(svgPathPair41);
	auto path41Curve = pmp::CurveFactory::sampled_polygon(path41Vertices, 100, false);
	RemeshWithDefaultSettings(path41Curve);
	if (!pmp::write_to_ply(path41Curve, dataOutPath + "path41Curve.ply"))
		std::cerr << "Error writing path41Curve.ply!\n";

	// Pair 5
	const auto path50Vertices = ParsePolygonalSVGPath(svgPathPair40);
	auto path50Curve = pmp::CurveFactory::sampled_polygon(path50Vertices, 100, true);
	RemeshWithDefaultSettings(path50Curve);
	if (!pmp::write_to_ply(path50Curve, dataOutPath + "path50Curve.ply"))
		std::cerr << "Error writing path50Curve.ply!\n";

	auto path51Curve = pmp::CurveFactory::circle(pmp::Point2{ 98.689514, 54.600803 }, 9.625, 80);
	//auto path51Curve = pmp::CurveFactory::circle(pmp::Point2{ 60.244923, 60.766129 }, 14.346979, 100);
	if (!pmp::write_to_ply(path51Curve, dataOutPath + "path51Curve.ply"))
		std::cerr << "Error writing path51Curve.ply!\n";

	// Pair 6
	const auto path60Vertices = ParsePolygonalSVGPath(svgPathPair40);
	auto path60Curve = pmp::CurveFactory::sampled_polygon(path60Vertices, 100, true);
	RemeshWithDefaultSettings(path60Curve);
	if (!pmp::write_to_ply(path60Curve, dataOutPath + "path60Curve.ply"))
		std::cerr << "Error writing path60Curve.ply!\n";

	const auto path61Vertices = ParsePolygonalSVGPath(svgPathPair61);
	auto path61Curve = pmp::CurveFactory::sampled_polygon(path61Vertices, 100, false);
	RemeshWithDefaultSettings(path61Curve);
	if (!pmp::write_to_ply(path61Curve, dataOutPath + "path61Curve.ply"))
		std::cerr << "Error writing path61Curve.ply!\n";

	// Pair 7
	const auto path70Vertices = ParsePolygonalSVGPath(svgPathPair40);
	auto path70Curve = pmp::CurveFactory::sampled_polygon(path70Vertices, 100, true);
	RemeshWithDefaultSettings(path70Curve);
	if (!pmp::write_to_ply(path70Curve, dataOutPath + "path70Curve.ply"))
		std::cerr << "Error writing path70Curve.ply!\n";

	const auto path71Vertices = ParsePolygonalSVGPath(svgPathPair71);
	auto path71Curve = pmp::CurveFactory::sampled_polygon(path71Vertices, 100, false);
	RemeshWithDefaultSettings(path71Curve);
	if (!pmp::write_to_ply(path71Curve, dataOutPath + "path71Curve.ply"))
		std::cerr << "Error writing path71Curve.ply!\n";

	// Pair 8
	const auto path80Vertices = ParsePolygonalSVGPath(svgPathPair40);
	auto path80Curve = pmp::CurveFactory::sampled_polygon(path80Vertices, 100, true);
	RemeshWithDefaultSettings(path80Curve);
	if (!pmp::write_to_ply(path80Curve, dataOutPath + "path80Curve.ply"))
		std::cerr << "Error writing path80Curve.ply!\n";

	const auto path81Vertices = ParsePolygonalSVGPath(svgPathPair81);
	auto path81Curve = pmp::CurveFactory::sampled_polygon(path81Vertices, 80, true);
	RemeshWithDefaultSettings(path81Curve);
	if (!pmp::write_to_ply(path81Curve, dataOutPath + "path81Curve.ply"))
		std::cerr << "Error writing path81Curve.ply!\n";

	// Pair 9
	auto path90Curve = pmp::CurveFactory::circle(pmp::Point2{ 60.411285, 68.040321 }, 56.508064, 180);
	RemeshWithDefaultSettings(path90Curve);
	if (!pmp::write_to_ply(path90Curve, dataOutPath + "path90Curve.ply"))
		std::cerr << "Error writing path90Curve.ply!\n";

	const auto path91Vertices = ParsePolygonalSVGPath(svgPathPair31);
	auto path91Curve = pmp::CurveFactory::sampled_polygon(path91Vertices, 50, true);
	RemeshWithDefaultSettings(path91Curve);
	if (!pmp::write_to_ply(path91Curve, dataOutPath + "path91Curve.ply"))
		std::cerr << "Error writing path91Curve.ply!\n";

	// Pair 10
	auto path100Curve = pmp::CurveFactory::circle(pmp::Point2{ 61.653221, 60.588707 }, 56.508064, 180);
	RemeshWithDefaultSettings(path100Curve);
	if (!pmp::write_to_ply(path100Curve, dataOutPath + "path100Curve.ply"))
		std::cerr << "Error writing path90Curve.ply!\n";

	const auto path101Vertices = ParsePolygonalSVGPath(svgPathPair31);
	auto path101Curve = pmp::CurveFactory::sampled_polygon(path101Vertices, 50, true);
	RemeshWithDefaultSettings(path101Curve);
	if (!pmp::write_to_ply(path101Curve, dataOutPath + "path101Curve.ply"))
		std::cerr << "Error writing path101Curve.ply!\n";

	// Pair 11
	const auto path110Vertices = ParsePolygonalSVGPath(svgPathPair40);
	auto path110Curve = pmp::CurveFactory::sampled_polygon(path110Vertices, 100, true);
	RemeshWithDefaultSettings(path110Curve);
	if (!pmp::write_to_ply(path110Curve, dataOutPath + "path110Curve.ply"))
		std::cerr << "Error writing path110Curve.ply!\n";

	const auto path111Vertices = ParsePolygonalSVGPath(svgPathPair111);
	auto path111Curve = pmp::CurveFactory::sampled_polygon(path111Vertices, 80, false);
	RemeshWithDefaultSettings(path111Curve);
	if (!pmp::write_to_ply(path111Curve, dataOutPath + "path111Curve.ply"))
		std::cerr << "Error writing path111Curve.ply!\n";

	// Pair 12
	const auto path120Vertices = ParsePolygonalSVGPath(svgPathPair40);
	auto path120Curve = pmp::CurveFactory::sampled_polygon(path120Vertices, 100, true);
	RemeshWithDefaultSettings(path120Curve);
	if (!pmp::write_to_ply(path120Curve, dataOutPath + "path120Curve.ply"))
		std::cerr << "Error writing path120Curve.ply!\n";

	const auto path121Vertices = ParsePolygonalSVGPath(svgPathPair121);
	auto path121Curve = pmp::CurveFactory::sampled_polygon(path121Vertices, 80, false);
	RemeshWithDefaultSettings(path121Curve);
	if (!pmp::write_to_ply(path121Curve, dataOutPath + "path121Curve.ply"))
		std::cerr << "Error writing path121Curve.ply!\n";

	// Pair 13
	auto path130Curve = pmp::CurveFactory::circle(pmp::Point2{ 60.411285, 68.040321 }, 56.508064, 100);
	RemeshWithDefaultSettings(path130Curve);
	if (!pmp::write_to_ply(path130Curve, dataOutPath + "path130Curve.ply"))
		std::cerr << "Error writing path130Curve.ply!\n";

	auto path131Curve = pmp::CurveFactory::circle(pmp::Point2{ 59.879032, 51.983871 }, 14.104838, 80);
	RemeshWithDefaultSettings(path131Curve);
	if (!pmp::write_to_ply(path131Curve, dataOutPath + "path131Curve.ply"))
		std::cerr << "Error writing path131Curve.ply!\n";

	// Pair 14
	auto path140Curve = pmp::CurveFactory::circle(pmp::Point2{ 60.411285, 68.040321 }, 56.508064, 100);
	RemeshWithDefaultSettings(path140Curve);
	if (!pmp::write_to_ply(path140Curve, dataOutPath + "path140Curve.ply"))
		std::cerr << "Error writing path140Curve.ply!\n";

	auto path141Curve = pmp::CurveFactory::circle(pmp::Point2{ 53.669357, 34.419353 }, 13.217741, 80);
	RemeshWithDefaultSettings(path141Curve);
	if (!pmp::write_to_ply(path141Curve, dataOutPath + "path141Curve.ply"))
		std::cerr << "Error writing path141Curve.ply!\n";

	// Pair 15
	auto path150Curve = pmp::CurveFactory::circle(pmp::Point2{ 64.497322, 62.582409 }, 92.768562, 250);
	RemeshWithDefaultSettings(path150Curve);
	if (!pmp::write_to_ply(path150Curve, dataOutPath + "path150Curve.ply"))
		std::cerr << "Error writing path150Curve.ply!\n";

	const auto path151Vertices = ParsePolygonalSVGPath(svgPathPair00);
	auto path151Curve = pmp::CurveFactory::sampled_polygon(path151Vertices, 140, false);
	RemeshWithDefaultSettings(path151Curve);
	if (!pmp::write_to_ply(path151Curve, dataOutPath + "path151Curve.ply"))
		std::cerr << "Error writing path151Curve.ply!\n";

	// Path 16
	const auto path160Vertices = ParsePolygonalSVGPath(svgPathPair160);
	auto path160Curve = pmp::CurveFactory::sampled_polygon(path160Vertices, 120, false);
	RemeshWithDefaultSettings(path160Curve);
	if (!pmp::write_to_ply(path160Curve, dataOutPath + "path160Curve.ply"))
		std::cerr << "Error writing path160Curve.ply!\n";

	const auto path161Vertices = ParsePolygonalSVGPath(svgPathPair161);
	auto path161Curve = pmp::CurveFactory::sampled_polygon(path161Vertices, 100, false);
	RemeshWithDefaultSettings(path161Curve);
	if (!pmp::write_to_ply(path161Curve, dataOutPath + "path161Curve.ply"))
		std::cerr << "Error writing path161Curve.ply!\n";

	// List of curve pairs to evolve
	const std::vector<std::pair<pmp::ManifoldCurve2D, pmp::ManifoldCurve2D>> curvePairs{
		//{path00Curve, path01Curve},
		//{path10Curve, path11Curve},
		//{path20Curve, path21Curve},
		//{path30Curve, path31Curve},
		//{path40Curve, path41Curve},
		//{path50Curve, path51Curve},
		//{path60Curve, path61Curve},
		//{path70Curve, path71Curve},
		//{path80Curve, path81Curve},
		//{path90Curve, path91Curve},
		//{path100Curve, path101Curve},
		//{path110Curve, path111Curve},
		//{path120Curve, path121Curve},
		//{path130Curve, path131Curve},
		//{path140Curve, path141Curve},
		//{path150Curve, path151Curve},
		{path160Curve, path161Curve},
	};

	// Prepare the settings for the evolver
	ManifoldEvolutionSettings strategySettings;
	strategySettings.UseInnerManifolds = true;

	strategySettings.AdvectionInteractWithOtherManifolds = true;
	// use PreComputeAdvectionDiffusionParams?
	strategySettings.OuterManifoldEpsilon = [](double distance)
	{
		return 0.0;

		//return 1.0 * (1.0 - exp(-distance * distance / 1.0));
	};
	strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
	{
		return 0.0;

		//return -1.0 * distance * (std::abs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));

		//return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		//return 1.0 * (1.0 - exp(-distance * distance / 0.5)) * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
	};
	strategySettings.InnerManifoldEpsilon = [](double distance)
	{
		return 0.0;

		//return 0.001 * TRIVIAL_EPSILON(distance);
	};
	strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
	{
		//return 0.0;

		return 1.0 * distance * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		//return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		// return 1.0 * distance * (std::fabs(-negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		// return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		//return 1.0 * distance * (std::fabs(negGradDotNormal) - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		//return -1.0 * (1.0 - exp(-distance * distance / 0.5)) * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
	};
	strategySettings.LevelOfDetail = 4;

	strategySettings.TimeStep = 0.05;
	strategySettings.TangentialVelocityWeight = 0.05;
	strategySettings.RemeshingSettings.MinEdgeMultiplier = 1.0;
	strategySettings.RemeshingSettings.UseBackProjection = true;
	strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
	strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = static_cast<pmp::Scalar>(M_PI_2);
	strategySettings.FieldSettings.NVoxelsPerMinDimension = 50;

	strategySettings.ExportVariableScalarFieldsDimInfo = true;
	strategySettings.ExportVariableVectorFieldsDimInfo = true;

	//strategySettings.DiagSettings.LogOuterManifoldEpsilon = true;
	//strategySettings.DiagSettings.LogInnerManifoldsEpsilon = true;
	//strategySettings.DiagSettings.LogOuterManifoldEta = true;
	strategySettings.DiagSettings.LogInnerManifoldsEta = true;

	// Global settings
	GlobalManifoldEvolutionSettings globalSettings;
	globalSettings.NSteps = 2500;
	//globalSettings.NSteps = 1000;
	//globalSettings.NSteps = 500;
	//globalSettings.NSteps = 20;
	globalSettings.DoRemeshing = true;
	globalSettings.DetectFeatures = false;
	globalSettings.ExportPerTimeStep = true;
	globalSettings.ExportTargetDistanceFieldAsImage = true;
	globalSettings.OutputPath = dataOutPath;
	globalSettings.ExportResult = false;

	globalSettings.RemeshingResizeFactor = 0.7;
	globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

	for (unsigned int pairId = 19; const auto & [outerCurve, innerCurve] : curvePairs)
	{
		std::cout << "EquilibriumPairedConcaveManifoldTests for pair: " << pairId << "\n";

		globalSettings.ProcedureName = "equilibriumConcavePair" + std::to_string(pairId);

		const auto innerCurveMinDim = Geometry::GetCurveBoundsMinDimension(innerCurve);
		strategySettings.FieldSettings.FieldIsoLevel = 4.0 * innerCurveMinDim / strategySettings.FieldSettings.NVoxelsPerMinDimension;

		std::vector innerCurves{ GetCurveRotatedAboutCenterPoint(innerCurve, angle) };
		//innerCurves[0].negate_orientation();
		auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
			strategySettings, outerCurve, innerCurves, nullptr);

		ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

		evolver.Evolve();

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
	constexpr pmp::Scalar edgeLength = 2.0;
	constexpr unsigned int iterations = 10;
	pmp::CurveRemeshing remesher(curves[1]);
	pmp::AdaptiveRemeshingSettings settings;
	settings.MinEdgeLength = edgeLength;
	settings.MaxEdgeLength = 1.5 * edgeLength;
	settings.ApproxError = 0.05 * edgeLength;
	settings.NRemeshingIterations = iterations;
	settings.NTangentialSmoothingIters = 6;
	settings.UseProjection = true;
	remesher.adaptive_remeshing(settings);

	constexpr unsigned int NTimeSteps = 30;

	for (unsigned int curveId = 0; const auto & curve : curves)
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

		strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.22;
		strategySettings.RemeshingSettings.UseBackProjection = false;

		strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
		strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

		globalSettings.RemeshingResizeFactor = 0.7;
		globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

		std::vector<pmp::ManifoldCurve2D> innerCurves{};

		auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
			strategySettings, curve, innerCurves, nullptr);

		std::cout << "Setting up ManifoldEvolver.\n";

		ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

		std::cout << "ManifoldEvolver::Evolve ... ";

		evolver.Evolve();

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
	{"armadillo", pmp::Point{-0.10348621158928151, 21.427067319905646, 9.79369240592005}},
	{"bunny", pmp::Point{-0.01684039831161499, 0.11015420407056808, 0.0012007840834242693} },
	{"maxPlanck", pmp::Point{30.59686279296875, -18.105804443359375, 82.29149055480957} },
	{"nefertiti", pmp::Point{0.0, 0.0, 0.0} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"armadillo", pmp::vec3{-0.03070969905335075, 0.12876712096541565, 0.9911992448253433}},
		{"bunny", pmp::vec3{0.0, 0.0, 1.0} },
		{"maxPlanck", pmp::vec3{1.0, 0.0, 0.0} },
		{"nefertiti", pmp::vec3{1.0, 0.0, 0.0} }
	};

	const std::map<std::string, Geometry::Circle2D> outerCircles{
		//{"armadillo", Geometry::Circle2D{pmp::Point2{0.372234, 16.6515}, 121.558} },
		//{"bunny", Geometry::Circle2D{pmp::Point2{-0.0155906, 0.102261}, 0.142831} },
		//{"maxPlanck", Geometry::Circle2D{pmp::Point2{-17.82, 82.5006}, 292.263} },
		{"nefertiti", Geometry::Circle2D{pmp::Point2{0.178497, -0.0410004}, 441.436} }
	};
	const std::map<std::string, std::vector<Geometry::Circle2D>> innerCircles{
		//{"armadillo", { Geometry::Circle2D{pmp::Point2{-3.0, 52.0}, 20.0} }},
		//{"bunny", { Geometry::Circle2D{pmp::Point2{-0.025, 0.08}, 0.025} } },
		//{"maxPlanck", { Geometry::Circle2D{pmp::Point2{8.0, 85.0}, 50.0} }},

		{"nefertiti", {
			Geometry::Circle2D{pmp::Point2{-20.0, 100.0}, 55.0},
			//Geometry::Circle2D{pmp::Point2{-75.0, -50.0}, 25.0}
			Geometry::Circle2D{pmp::Point2{-10.0, -200.0}, 35.0}
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
		const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		//const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;

		const double isoLvlOffsetFactor = (isoLevelOffsetFactors.contains(ptCloudName) ? isoLevelOffsetFactors.at(ptCloudName) : defaultOffsetFactor);
		const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

		const double tau = (timeStepSizesForPtClouds.contains(ptCloudName) ? timeStepSizesForPtClouds.at(ptCloudName) : defaultTimeStep); // time step

		// ==========================================================================
		// - - - - - - - - -  New Manifold Evolver (Curve)  - - - - - - - - - - - - 
		// ==========================================================================

		const pmp::Scalar distTolerance = 0.01 * minSize;
		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0, 0.0, 0.0 });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}

		if (!Geometry::Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
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

			//strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.22;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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
		{"armadillo", pmp::Point{-0.10348621158928151, 21.427067319905646, 9.79369240592005}},
		{"bunny", pmp::Point{-0.01684039831161499, 0.11015420407056808, 0.0012007840834242693} },
		{"maxPlanck", pmp::Point{30.59686279296875, -18.105804443359375, 82.29149055480957} },
		{"nefertiti", pmp::Point{0.0, 0.0, 0.0} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"armadillo", pmp::vec3{-0.03070969905335075, 0.12876712096541565, 0.9911992448253433}},
		{"bunny", pmp::vec3{0.0, 0.0, 1.0} },
		{"maxPlanck", pmp::vec3{1.0, 0.0, 0.0} },
		{"nefertiti", pmp::vec3{1.0, 0.0, 0.0} }
	};

	for (const auto& meshName : meshForPtCloudNames)
	{
		Geometry::PlaneSettings planeSettings;
		planeSettings.Width = 200.0;
		planeSettings.Depth = 250.0;
		planeSettings.nWidthSegments = 1;
		planeSettings.nDepthSegments = 1;
		planeSettings.Origin = slicingPlaneRefPts.at(meshName) - 0.5 * pmp::vec3{ planeSettings.Width, planeSettings.Depth, 0.0 };

		Geometry::PlaneBuilder pb(planeSettings);
		pb.BuildBaseData();
		pb.BuildPMPSurfaceMesh();
		auto pMesh = pb.GetPMPSurfaceMeshResult();

		const pmp::vec3 defaultNormal(0.0, 0.0, 1.0);
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
		{ "Circle", pmp::CurveFactory::circle(pmp::Point2{-3.0, 52.0}, 35.0, 25, 0.0, 2.0 * M_PI)},
		{ "IncompleteCircle", pmp::CurveFactory::circle(pmp::Point2{-3.0, 52.0}, 35.0, 25, M_PI_2, 2.0 * M_PI)},
		{ "SineDeformedCircle", pmp::CurveFactory::sine_deformed_circle(pmp::Point2{-3.0, 52.0}, 35.0, 25, 7.0, 4.0, 0.0, 2.0 * M_PI)},
		{ "SineDeformedIncompleteCircle", pmp::CurveFactory::sine_deformed_circle(pmp::Point2{-3.0, 52.0}, 35.0, 25, 7.0, 4.0, M_PI_2, 2.0 * M_PI)},
		{ "ChamferedRectangle", pmp::CurveFactory::rectangle(pmp::Point2{-3.0, 52.0}, 60.0, 70.0, 15, true)},
		{ "IncompleteChamferedRectangle", pmp::CurveFactory::sampled_polygon({
			pmp::Point2{-30.0, -35.0} + pmp::Point2{-3.0, 52.0},
			pmp::Point2{30.0, -35.0} + pmp::Point2{-3.0, 52.0},
			pmp::Point2{30.0, 35.0} + pmp::Point2{-3.0, 52.0},
			pmp::Point2{-30.0, 35.0} + pmp::Point2{-3.0, 52.0}}, 30, true, false)},
		{ "ChamferedTriangle", pmp::CurveFactory::sampled_polygon({
			pmp::Point2{-0.5, (pmp::Scalar)-sqrtf(3.0) / (pmp::Scalar)6.0} *120.0 + pmp::Point2{-3.0, 52.0},
			pmp::Point2{0.5, (pmp::Scalar)-sqrtf(3.0) / (pmp::Scalar)6.0} *120.0 + pmp::Point2{-3.0, 52.0},
			pmp::Point2{0.0, (pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)3.0} *120.0 + pmp::Point2{-3.0, 52.0}}, 30, true)},
		{ "IncompleteChamferedTriangle", pmp::CurveFactory::sampled_polygon({
			pmp::Point2{-0.5, (pmp::Scalar)-sqrtf(3.0) / (pmp::Scalar)6.0} *120.0 + pmp::Point2{-3.0, 52.0},
			pmp::Point2{0.5, (pmp::Scalar)-sqrtf(3.0) / (pmp::Scalar)6.0} *120.0 + pmp::Point2{-3.0, 52.0},
			pmp::Point2{0.0, (pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)3.0} *120.0 + pmp::Point2{-3.0, 52.0}}, 30, true, false)}
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;
	constexpr unsigned int NTimeSteps = 180;

	const auto innerCircle = Geometry::Circle2D{ pmp::Point2{ -3.0, 52.0 }, 20.0 };
	const auto outerCircle = Geometry::Circle2D{ pmp::Point2{ -3.0, 52.0 }, 110.0 };

	for (const auto& [ptCloudName, curve] : targetCurves)
	{
		std::cout << "========================================================\n";
		std::cout << "         inner/outer circle LSW for: " << ptCloudName << " ... \n";
		std::cout << " ------------------------------------------------------ \n";

		const auto& pts2D = curve.positions();
		const pmp::BoundingBox2 bbox{ pts2D };

		const auto bboxSize = bbox.max() - bbox.min();
		const pmp::Scalar minSize = std::min(bboxSize[0], bboxSize[1]);
		//const pmp::Scalar maxSize = std::max(bboxSize[0], bboxSize[1]);
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;

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

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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
	{ "Circle", pmp::CurveFactory::circle(pmp::Point2{-3.0, 52.0}, 35.0, 25, 0.0, 2.0 * M_PI)},
	{ "IncompleteCircle", pmp::CurveFactory::circle(pmp::Point2{-3.0, 52.0}, 35.0, 25, M_PI_2, 2.0 * M_PI)},
	{ "SineDeformedCircle", pmp::CurveFactory::sine_deformed_circle(pmp::Point2{-3.0, 52.0}, 35.0, 25, 7.0, 4.0, 0.0, 2.0 * M_PI)},
	{ "SineDeformedIncompleteCircle", pmp::CurveFactory::sine_deformed_circle(pmp::Point2{-3.0, 52.0}, 35.0, 25, 7.0, 4.0, M_PI_2, 2.0 * M_PI)},
	{ "ChamferedRectangle", pmp::CurveFactory::rectangle(pmp::Point2{-3.0, 52.0}, 60.0, 70.0, 15, true)},
	{ "IncompleteChamferedRectangle", pmp::CurveFactory::sampled_polygon({
		pmp::Point2{-30.0, -35.0} + pmp::Point2{-3.0, 52.0},
		pmp::Point2{30.0, -35.0} + pmp::Point2{-3.0, 52.0},
		pmp::Point2{30.0, 35.0} + pmp::Point2{-3.0, 52.0},
		pmp::Point2{-30.0, 35.0} + pmp::Point2{-3.0, 52.0}}, 30, true, false)},
	{ "ChamferedTriangle", pmp::CurveFactory::sampled_polygon({
		pmp::Point2{-0.5, (pmp::Scalar)-sqrtf(3.0) / (pmp::Scalar)6.0} *120.0 + pmp::Point2{-3.0, 52.0},
		pmp::Point2{0.5, (pmp::Scalar)-sqrtf(3.0) / (pmp::Scalar)6.0} *120.0 + pmp::Point2{-3.0, 52.0},
		pmp::Point2{0.0, (pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)3.0} *120.0 + pmp::Point2{-3.0, 52.0}}, 30, true)},
	{ "IncompleteChamferedTriangle", pmp::CurveFactory::sampled_polygon({
		pmp::Point2{-0.5, (pmp::Scalar)-sqrtf(3.0) / (pmp::Scalar)6.0} *120.0 + pmp::Point2{-3.0, 52.0},
		pmp::Point2{0.5, (pmp::Scalar)-sqrtf(3.0) / (pmp::Scalar)6.0} *120.0 + pmp::Point2{-3.0, 52.0},
		pmp::Point2{0.0, (pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)3.0} *120.0 + pmp::Point2{-3.0, 52.0}}, 30, true, false)}
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;
	constexpr unsigned int NTimeSteps = 180;

	const auto outerCircle = Geometry::Circle2D{ pmp::Point2{ -3.0, 52.0 }, 110.0 };

	for (const auto& [ptCloudName, curve] : targetCurves)
	{
		std::cout << "========================================================\n";
		std::cout << "         inner/outer circle LSW for: " << ptCloudName << " ... \n";
		std::cout << " ------------------------------------------------------ \n";

		const auto& pts2D = curve.positions();
		const pmp::BoundingBox2 bbox{ pts2D };

		const auto bboxSize = bbox.max() - bbox.min();
		const pmp::Scalar minSize = std::min(bboxSize[0], bboxSize[1]);
		//const pmp::Scalar maxSize = std::max(bboxSize[0], bboxSize[1]);
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;

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

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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

	const std::map<std::string, Geometry::Sphere3D> outerSpheres{
		{"boxWithAHole", Geometry::Sphere3D{pmp::Point{0, 0, 0}, 2.0} },
	};
	const std::map<std::string, std::vector<Geometry::Sphere3D>> innerSpheres{
		{"boxWithAHole", std::vector{ Geometry::Sphere3D{pmp::Point{0, 0, 0}, 0.9}} },
	};

	const std::map<std::string, pmp::Point> slicingPlaneRefPts{
		{"boxWithAHole", pmp::Point{0.0, 0.0, 0.0}},
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"boxWithAHole", pmp::vec3{0.0, -1.0, 0.0}},
	};

	const std::map<std::string, Geometry::Circle2D> outerCircles{
		{"boxWithAHole", Geometry::Circle2D{pmp::Point2{0.0, 0.0}, 2.0} },
	};
	const std::map<std::string, std::vector<Geometry::Circle2D>> innerCircles{
		{"boxWithAHole", std::vector{ Geometry::Circle2D{pmp::Point2{0.0, 0.0}, 0.9}} }
	};

	const std::map<std::string, std::vector<Geometry::Sphere3D>> cutSpheres{
		{"boxWithAHole", std::vector{ Geometry::Sphere3D{pmp::Point{1.0, 0.0, 0.0}, 0.6} }},
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
		const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
		constexpr pmp::Scalar volExpansionFactor = 1.0;
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

		const pmp::Scalar distTolerance = 0.01 * minSize;
		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0, 0.0, 0.0 });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}

		if (!Geometry::Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
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

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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
				return 0.4 * distance * (negGradDotNormal - 1.3 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 4;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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
				1.0, 0.0, 0.0, outerSphere.Center[0],
				0.0, 1.0, 0.0, outerSphere.Center[1],
				0.0, 0.0, 1.0, outerSphere.Center[2],
				0.0, 0.0, 0.0, 1.0
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
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				return 0.0005 * TRIVIAL_EPSILON(distance);
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 0.4 * distance * (negGradDotNormal - 1.3 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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
				1.0, 0.0, 0.0, outerSphere.Center[0],
				0.0, 1.0, 0.0, outerSphere.Center[1],
				0.0, 0.0, 1.0, outerSphere.Center[2],
				0.0, 0.0, 0.0, 1.0
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
					1.0, 0.0, 0.0, innerSphere.Center[0],
					0.0, 1.0, 0.0, innerSphere.Center[1],
					0.0, 0.0, 1.0, innerSphere.Center[2],
					0.0, 0.0, 0.0, 1.0
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
		{"armadillo", pmp::Point{-0.10348621158928151, 21.427067319905646, 9.79369240592005}},
		{"bunny", pmp::Point{-0.01684039831161499, 0.11015420407056808, 0.0012007840834242693} },
		//{"bunny", pmp::Point{-0.021311134388792542, 0.11143290480956525, 0.007428090888817638} },
		{"maxPlanck", pmp::Point{30.59686279296875, -18.105804443359375, 82.29149055480957} },
		{"nefertiti", pmp::Point{0.0, 0.0, 0.0} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"armadillo", pmp::vec3{-0.03070969905335075, 0.12876712096541565, 0.9911992448253433}},
		{"bunny", pmp::vec3{0.0, 0.0, 1.0} },
		//{"bunny", pmp::vec3{-0.18754182700901353, -0.029010506066971985, 0.9818281181855913} },
		{"maxPlanck", pmp::vec3{1.0, 0.0, 0.0} },
		{"nefertiti", pmp::vec3{1.0, 0.0, 0.0} }
	};

	const std::map<std::string, Geometry::Circle2D> outerCircles{
		{"armadillo", Geometry::Circle2D{pmp::Point2{0.372234, 16.6515}, 121.558} },
		{"bunny", Geometry::Circle2D{pmp::Point2{-0.0155906, 0.102261}, 0.142831} },
		{"maxPlanck", Geometry::Circle2D{pmp::Point2{-17.82, 82.5006}, 292.263} },
		{"nefertiti", Geometry::Circle2D{pmp::Point2{0.178497, -0.0410004}, 441.436} }
	};
	const std::map<std::string, std::vector<Geometry::Circle2D>> innerCircles{
		//{"armadillo", std::vector{ Geometry::Circle2D{pmp::Point2{-3.0, 52.0}, 20.0}} },
		{"bunny", std::vector{ Geometry::Circle2D{pmp::Point2{0.0, 0.08}, 0.025}} },
		{"maxPlanck", std::vector{ Geometry::Circle2D{pmp::Point2{8.0, 85.0}, 50.0}} },
		//{"nefertiti", std::vector{ Geometry::Circle2D{pmp::Point2{-20.0, 100.0}, 55.0}} }
	};

	const std::map<std::string, Geometry::Sphere3D> outerSpheres{
		{"armadillo", Geometry::Sphere3D{pmp::Point{0.0122509, 21.4183, -0.000249863}, 136.963} },
		{"bunny", Geometry::Sphere3D{pmp::Point{-0.0168297, 0.110217, -0.0015718}, 0.141622} },
		{"maxPlanck", Geometry::Sphere3D{pmp::Point{30.658, -17.9765, 82.2885}, 271.982} },
		{"nefertiti", Geometry::Sphere3D{pmp::Point{0.0144997, -0.00499725, -0.0215073}, 392.184} }
	};
	const std::map<std::string, std::vector<Geometry::Sphere3D>> innerSpheres{
		//{"armadillo", std::vector{ Geometry::Sphere3D{pmp::Point{-3.0, 52.0}, 20.0}} },
		{"bunny", std::vector{ Geometry::Sphere3D{pmp::Point{0.0, 0.082, 0.012}, 0.03}} },
		{"maxPlanck", std::vector{ Geometry::Sphere3D{pmp::Point{8.0, 50.0, 100.0}, 50.0}} },
		//{"nefertiti", std::vector{ Geometry::Sphere3D{pmp::Point{-20.0, 100.0}, 55.0}} }
	};

	const std::map<std::string, std::vector<Geometry::Sphere3D>> cutSpheres{
		{"armadillo", {}},
		{"bunny", std::vector{ Geometry::Sphere3D{pmp::Point{-0.01, 0.06, 0.012}, 0.032}, Geometry::Sphere3D{pmp::Point{0.01, 0.12, 0.01}, 0.025}/**/}},
		{"maxPlanck", std::vector{ Geometry::Sphere3D{pmp::Point{8.0, 85.0, 0.0}, 50.0}, Geometry::Sphere3D{pmp::Point{30.0, -120.0, 160.0}, 100.0} /**/}},
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
		const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
		constexpr pmp::Scalar volExpansionFactor = 1.0;
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

		const float distTolerance = 0.01 * minSize;
		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0, 0.0, 0.0 });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}

		if (!Geometry::Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
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

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<float>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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
				1.0, 0.0, 0.0, outerSphere.Center[0],
				0.0, 1.0, 0.0, outerSphere.Center[1],
				0.0, 0.0, 1.0, outerSphere.Center[2],
				0.0, 0.0, 0.0, 1.0
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
				return 1.0 * (1.0 - exp(-distance * distance / 1.0));
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (negGradDotNormal - 2.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.InnerManifoldEpsilon = [](double distance)
			{
				return 0.0005 * TRIVIAL_EPSILON(distance);
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 0.4 * distance * (negGradDotNormal - 1.3 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = tau;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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
				1.0, 0.0, 0.0, outerSphere.Center[0],
				0.0, 1.0, 0.0, outerSphere.Center[1],
				0.0, 0.0, 1.0, outerSphere.Center[2],
				0.0, 0.0, 0.0, 1.0
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
					1.0, 0.0, 0.0, innerSphere.Center[0],
					0.0, 1.0, 0.0, innerSphere.Center[1],
					0.0, 0.0, 1.0, innerSphere.Center[2],
					0.0, 0.0, 0.0, 1.0
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
			if (!ExportBaseCurveGeometryDataToPLY(*medialAxis, dataOutPath + "sawhneyShapeQ.ply"))
				std::cerr << "Error writing sawhneyShapeQ.ply!\n";
		}
	}

	{
		const auto medialAxis = Geometry::GetMedialAxisOfSawhneysStupidMATAlgorithm('w');
		if (medialAxis.has_value())
		{
			if (!ExportBaseCurveGeometryDataToPLY(*medialAxis, dataOutPath + "sawhneyShapeW.ply"))
				std::cerr << "Error writing sawhneyShapeW.ply!\n";
		}
	}

	{
		const auto medialAxis = Geometry::GetMedialAxisOfSawhneysStupidMATAlgorithm('e');
		if (medialAxis.has_value())
		{
			if (!ExportBaseCurveGeometryDataToPLY(*medialAxis, dataOutPath + "sawhneyShapeE.ply"))
				std::cerr << "Error writing sawhneyShapeE.ply!\n";
		}
	}

	{
		const auto medialAxis = Geometry::GetMedialAxisOfSawhneysStupidMATAlgorithm('r');
		if (medialAxis.has_value())
		{
			if (!ExportBaseCurveGeometryDataToPLY(*medialAxis, dataOutPath + "sawhneyShapeR.ply"))
				std::cerr << "Error writing sawhneyShapeR.ply!\n";
		}
	}

	{
		const auto medialAxis = Geometry::GetMedialAxisOfSawhneysStupidMATAlgorithm('t');
		if (medialAxis.has_value())
		{
			if (!ExportBaseCurveGeometryDataToPLY(*medialAxis, dataOutPath + "sawhneyShapeT.ply"))
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
	constexpr pmp::Scalar edgeLength = 15.0;
	constexpr unsigned int iterations = 10;
	pmp::CurveRemeshing remesher(hyperEllipse);
	pmp::AdaptiveRemeshingSettings settings;
	settings.MinEdgeLength = edgeLength;
	settings.MaxEdgeLength = 1.5 * edgeLength;
	settings.ApproxError = 0.05 * edgeLength;
	settings.NRemeshingIterations = iterations;
	settings.NTangentialSmoothingIters = 6;
	settings.UseProjection = true;
	remesher.adaptive_remeshing(settings);

	const std::vector<Geometry::Circle2D> testInnerCircles{
		{ pmp::Point2{ 0.0, 0.0 }, 70.0 },
		{ pmp::Point2{ 100.0, 0.0 }, 60.0 },
		{ pmp::Point2{ 140.0, 0.0 }, 40.0 },
		{ pmp::Point2{ 130.0, 0.0 }, 60.0 },
		{ pmp::Point2{ 130.0, 40.0 }, 40.0 }
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.5;
	constexpr unsigned int NTimeSteps = 180;
	const double fieldIsoLevel = defaultOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(5.0);

	for (size_t hyperellipseId = 9; const auto & innerCircle : testInnerCircles)
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
			return -1.0 * distance * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		};
		strategySettings.TimeStep = defaultTimeStep;
		strategySettings.LevelOfDetail = 4;
		strategySettings.TangentialVelocityWeight = 0.05;

		strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.22;
		strategySettings.RemeshingSettings.UseBackProjection = false;

		strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
		strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

		globalSettings.RemeshingResizeFactor = 0.7;
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

		evolver.Evolve();

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
		{"canStraight", pmp::Point{0.07009050250053406, 0.3348689912818372, -0.018547505140304565}},
		{"canStraightMissingBottom", pmp::Point{0.07009050250053406, 0.3348689912818372, -0.018547505140304565} },
		{"crushedCan", pmp::Point{0.0007844995707273483, 0.030990499537438154, -0.006461501121520996} },
		{"crushedCanMissingBottom", pmp::Point{0.0007844995707273483, 0.030990499537438154, -0.006461501121520996} }
	};

	const std::map<std::string, pmp::vec3> slicingPlaneNormals{
		{"canStraight", pmp::vec3{1.0, 0.0, 0.0}},
		{"canStraightMissingBottom", pmp::vec3{1.0, 0.0, 0.0} },
		{"crushedCan", pmp::vec3{1.0, 0.0, 0.0} },
		{"crushedCanMissingBottom", pmp::vec3{1.0, 0.0, 0.0} }
	};

	const std::map<std::string, Geometry::Circle2D> outerCircles{
		{"canStraight", Geometry::Circle2D{pmp::Point2{0.25, 0.0}, 0.9} },
		{"canStraightMissingBottom", Geometry::Circle2D{pmp::Point2{0.25, 0.0}, 0.9} },
		{"crushedCan", Geometry::Circle2D{pmp::Point2{0.0205, 0.0}, 0.07} },
		{"crushedCanMissingBottom", Geometry::Circle2D{pmp::Point2{0.0205, 0.0}, 0.07} }
	};
	const std::map<std::string, std::vector<Geometry::Circle2D>> innerCircles{
		{"canStraight", std::vector{ Geometry::Circle2D{pmp::Point2{0.35, -0.25}, 0.16} } },
		{"canStraightMissingBottom", std::vector{ Geometry::Circle2D{pmp::Point2{0.35, -0.25}, 0.16}, Geometry::Circle2D{pmp::Point2{0.25, 0.35}, 0.16}} },
		{"crushedCan", std::vector{ Geometry::Circle2D{pmp::Point2{0.025, -0.035}, 0.02}} },
		{"crushedCanMissingBottom", std::vector{ Geometry::Circle2D{pmp::Point2{0.030, -0.035}, 0.02}, Geometry::Circle2D{pmp::Point2{0.025, 0.035}, 0.02} } }
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
		const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;
		constexpr pmp::Scalar volExpansionFactor = 1.0;
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

		const pmp::Scalar distTolerance = 0.05 * minSize;
		const auto planeRefPt = (slicingPlaneRefPts.contains(meshName) ? slicingPlaneRefPts.at(meshName) : center);
		const auto planeNormal = (slicingPlaneNormals.contains(meshName) ? slicingPlaneNormals.at(meshName) : pmp::vec3{ -1.0, 0.0, 0.0 });
		const auto pts2D = Geometry::GetSliceOfThePointCloud(ptCloud, planeRefPt, planeNormal, distTolerance);
		if (pts2D.empty())
		{
			std::cerr << "GetSliceOfThePointCloud sampled no 2D points during slicing for mesh " << meshName << "!\n";
			continue;
		}

		if (!Geometry::Export2DPointCloudToPLY(pts2D, dataOutPath + meshName + "_Pts_2D.ply"))
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

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

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

			globalSettings.RemeshingResizeFactor = 0.7;
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
	using namespace Geometry;
	// Data (incomplete circle)
	pmp::ManifoldCurve2D targetCurve = pmp::CurveFactory::circle(pmp::Point2(0.0, 0.0), 0.75, 16);
	auto targetPts = targetCurve.positions();
	targetPts.erase(targetPts.begin());
	targetPts.erase(targetPts.begin());
	targetPts.erase(targetPts.begin());
	InscribedCircleInputData inputData;
	inputData.Points = targetPts;

	const auto pointBBox = pmp::BoundingBox2(inputData.Points);
	const auto pointBBoxSize = pointBBox.max() - pointBBox.min();
	const pmp::Scalar minSize = std::min(pointBBoxSize[0], pointBBoxSize[1]);
	const pmp::Scalar cellSize = minSize / 20.0;
	const SDF::PointCloudDistanceField2DSettings sdfSettings{
		cellSize,
		0.5,
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
		Geometry::IcoSphereBuilder icoBuilder({ 2, 1.0 });
		icoBuilder.BuildBaseData();
		icoBuilder.BuildPMPSurfaceMesh();
		auto sphereMesh = icoBuilder.GetPMPSurfaceMeshResult();
		constexpr pmp::Scalar a = 1.0;
		constexpr pmp::Scalar b = 1.5;
		constexpr pmp::Scalar c = 2.0;
		const auto scalingMat = scaling_matrix(pmp::vec3{ a, b, c });
		sphereMesh *= scalingMat;

		const auto points = sphereMesh.positions();
		const auto pointBBox3D = pmp::BoundingBox(points);
		const auto pointBBox3DSize = pointBBox3D.max() - pointBBox3D.min();
		const pmp::Scalar minSize3D = std::min({ pointBBox3DSize[0], pointBBox3DSize[1], pointBBox3DSize[2] });
		const pmp::Scalar cellSize3D = minSize3D / 20.0;

		const SDF::PointCloudDistanceFieldSettings ellipsoidDfSettings{
			cellSize3D,
			0.5,
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

	const std::map<std::string, std::vector<Geometry::Sphere3D>> cutSpheres{
		{"bunnyWNormals", std::vector{ Geometry::Sphere3D{pmp::Point{-0.01, 0.06, 0.012}, 0.032}, Geometry::Sphere3D{pmp::Point{0.01, 0.12, 0.01}, 0.025}/**/}},
		{"maxPlanckWNormals", std::vector{ Geometry::Sphere3D{pmp::Point{8.0, 85.0, 0.0}, 50.0}, Geometry::Sphere3D{pmp::Point{30.0, -120.0, 160.0}, 100.0} /**/}},
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
	pmp::ManifoldCurve2D unevenCrossCurve = pmp::CurveFactory::sampled_polygon(unevenCrossPolyPts, 100, false);
	if (!pmp::write_to_ply(unevenCrossCurve, dataOutPath + "unevenCrossCurve.ply"))
		std::cerr << "Error writing unevenCrossCurve.ply!\n";

	const auto unevenCrossPts = unevenCrossCurve.positions();

	const std::vector<Geometry::Circle2D> testInnerCircles{
		//{ pmp::Point2{ 60.322582, 63.604839 }, 6.6532254 },
		//{ pmp::Point2{ 87.733871, 55.975807 }, 9.08871 },
		//{ pmp::Point2{ 62.3629, 61.741936 }, 19.161289 },
		//{ pmp::Point2{ 73.274193, 103.87903 }, 8.5161285 },
		//{ pmp::Point2{ 24.749998, 67.330643 }, 10.733871 },
		{ pmp::Point2{ 53.935482, 30.129032 }, 7.6290321 },
		//{ pmp::Point2{ 51.0, 72.0 }, 5.0 },

	};

	constexpr unsigned int nVoxelsPerMinDimension = 80;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 1.0;
	constexpr unsigned int NTimeSteps = 1500;
	const double fieldIsoLevel = defaultOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(5.0);

	for (size_t unevenCrossId = 10; const auto & innerCircle : testInnerCircles)
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
			//return 0.0025 * TRIVIAL_EPSILON(distance);
			return 0.0;
		};
		strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
		{
			// return 1.5 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			return -1.5 * distance * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		};
		strategySettings.TimeStep = defaultTimeStep;
		strategySettings.LevelOfDetail = 4;
		strategySettings.TangentialVelocityWeight = 0.05;

		strategySettings.RemeshingSettings.MinEdgeMultiplier = 1.0;
		strategySettings.RemeshingSettings.UseBackProjection = true;

		strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
		strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

		strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
		strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;
		strategySettings.FieldSettings.FieldExpansionFactor = 0.5;

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

		globalSettings.RemeshingResizeFactor = 0.7;
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

		evolver.Evolve();

		unevenCrossId++;
	}
}

void TestDFDivergence2D()
{
	using namespace Geometry;
	//const std::vector unevenCrossPolyPts{
	//	pmp::Point2{39.507142, 14.544772},
	//	pmp::Point2{46.104261, 5.542495},
	//	pmp::Point2{61.36143, 4.906308},
	//	pmp::Point2{68.282948, 13.11281},
	//	pmp::Point2{66.153916, 31.426095},
	//	pmp::Point2{69.924933, 39.754365},
	//	pmp::Point2{111.082270, 39.723965},
	//	pmp::Point2{117.795930, 46.528945},
	//	pmp::Point2{117.765430, 66.419005},
	//	pmp::Point2{113.358230, 72.215125},
	//	pmp::Point2{89.030514, 71.743755},
	//	pmp::Point2{82.788235, 77.337235},
	//	pmp::Point2{87.565954, 122.613040},
	//	pmp::Point2{80.332640, 129.222760},
	//	pmp::Point2{68.372952, 128.366110},
	//	pmp::Point2{61.451434, 122.392570},
	//	pmp::Point2{57.089489, 85.394595},
	//	pmp::Point2{36.929786, 84.297475},
	//	pmp::Point2{10.835265, 83.544745},
	//	pmp::Point2{3.558908, 76.519305},
	//	pmp::Point2{3.558908, 57.450225},
	//	pmp::Point2{10.584357, 47.664785},
	//	pmp::Point2{32.413427, 48.166595},
	//	pmp::Point2{40.865615, 41.985195}
	//};

	const auto unevenCrossPolyPts = ParsePolygonalSVGPath(svgPathPair160);
	pmp::ManifoldCurve2D unevenCrossCurve = pmp::CurveFactory::sampled_polygon(unevenCrossPolyPts, 100, false);
	if (!pmp::write_to_ply(unevenCrossCurve, dataOutPath + "unevenCrossCurve.ply"))
		std::cerr << "Error writing unevenCrossCurve.ply!\n";
	const auto unevenCrossPts = unevenCrossCurve.positions();

	constexpr unsigned int nVoxelsPerMinDimension = 350;

	const pmp::BoundingBox2 ptCloudBBox(unevenCrossPts);
	const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
	const pmp::Scalar minSize = std::min(ptCloudBBoxSize[0], ptCloudBBoxSize[1]);
	const pmp::Scalar maxSize = std::max(ptCloudBBoxSize[0], ptCloudBBoxSize[1]);
	const pmp::Scalar cellSize = minSize / static_cast<pmp::Scalar>(nVoxelsPerMinDimension);
	const SDF::PointCloudDistanceField2DSettings dfSettings{
		cellSize,
		0.6,
		DBL_MAX
	};
	const auto df = SDF::PlanarPointCloudDistanceFieldGenerator::Generate(unevenCrossPts, dfSettings);

	//const SDF::DistanceField2DSettings dfSettings{
	//	cellSize,
	//	0.6,
	//	DBL_MAX
	//};
	//const Geometry::ManifoldCurve2DAdapter unevenCrossCurveAdapter(std::make_shared<pmp::ManifoldCurve2D>(unevenCrossCurve));
	//const auto df = SDF::PlanarDistanceFieldGenerator::Generate(unevenCrossCurveAdapter, dfSettings);

	auto blurredDf = df;
	ApplyWideGaussianBlur2D(blurredDf);
	const unsigned int nPx = 2;
	ExportScalarGrid2DToPNG(dataOutPath + "unevenCrossDF.png", blurredDf,
		Geometry::BilinearInterpolateScalarValue,
		//Geometry::GetNearestNeighborScalarValue2D,
		nPx, nPx, RAINBOW_TO_WHITE_MAP);

	//const auto negGradDF = ComputeNormalizedNegativeGradient(blurredDf);
	const auto gradDF = ComputeNormalizedGradient(blurredDf);

	// Calculate stream lines
	Geometry::StreamLineSettings slSettings{};
	slSettings.NSamplePts = 1000;
	slSettings.StepSize = 0.5 * cellSize;
	slSettings.NSteps = 100;
	const auto streamLines = CalculateStreamLines(gradDF, slSettings);
	ExportPolyLinesToPLY(streamLines, dataOutPath + "unevenCrossNegGradDF_Streamlines.ply");

	// Calculate euler stream lines
	//Geometry::EulerStreamLineSettings eslSettings{};
	//eslSettings.NSamplePts = 1000;
	//eslSettings.StepSize = 0.5 * cellSize;
	//eslSettings.NSteps = 100;
	//const auto characteristics = CalculateStreamLinesEuler(negGradDF, eslSettings);
	//ExportPolyLinesToPLY(characteristics, dataOutPath + "unevenCrossNegGradDF_StreamlinesEuler.ply");

	// Calculate characteristics
	CharacteristicsBuilderSettings chBuilderSettings;
	chBuilderSettings.DFSettings = dfSettings;
	chBuilderSettings.ConstructInwardCharacteristics = true;
	chBuilderSettings.ConstructOutwardCharacteristics = true;
	chBuilderSettings.DivFieldThresholdFactor = -0.15;
	PlanarPointCloudCharacteristicsBuilder charBuilder{ unevenCrossPts, chBuilderSettings };

	//PlanarManifoldCurveCharacteristicsBuilder charBuilder{ unevenCrossCurve, chBuilderSettings };
	const auto characteristics = charBuilder.Build();
	ExportPolyLinesToPLY(characteristics, dataOutPath + "unevenCrossNegGradDF_Characteristics.ply");

	//auto negGradDF = ComputeGradient(blurredDf);
	//NegateGrid(negGradDF);
	const auto divNegGradDF = ComputeDivergenceField(gradDF);
	const auto [dMin, dMax] = std::ranges::minmax_element(divNegGradDF.Values());
	std::cout << "divNegGradDF value range: [" << *dMin << ", " << *dMax << "]\n";
	const auto colorMap = AdjustColorMapForZeroMidpoint(SIGN_TEMP_MAP, *dMin, *dMax);
	const double rangeVal = std::min(std::abs(*dMin), *dMax);
	ExportScalarGridDimInfo2D(dataOutPath + "unevenCrossDivNegGradDF.gdim2d", divNegGradDF);
	ExportScalarGrid2DToPNG(dataOutPath + "unevenCrossDivNegGradDF.png", divNegGradDF,
		Geometry::BilinearInterpolateScalarValue,
		//Geometry::GetNearestNeighborScalarValue2D,
		nPx, nPx, colorMap);
	constexpr bool verticalLegend = true;
	constexpr unsigned int resolutionFactor = 4;
	constexpr unsigned int legendPxHeight = (verticalLegend ? 400 : 100) * resolutionFactor;
	constexpr unsigned int legendPxWidth = (verticalLegend ? 100 : 600) * resolutionFactor;
	ExportColorScaleToPNG(
		dataOutPath + "unevenCrossDivNegGradDF_Scale.png",
		-rangeVal,
		rangeVal,
		colorMap,
		legendPxHeight, legendPxWidth);
}

void TestArcLengthCalculation()
{
	auto testCircleCurve = pmp::CurveFactory::circle(pmp::Point2{ 53.669357, 34.419353 }, 13.217741, 10);

	const auto arcLengthCalc = std::make_shared<pmp::EvolvingArcLengthCalculator>(testCircleCurve);

	std::cout << "TestArcLengthCalculation: test 1: circle with 10 segments\n";
	const auto arcLengths1 = arcLengthCalc->CalculateArcLengths();
	Geometry::PrintCurveValuesInTopologicalOrder(testCircleCurve, arcLengths1, std::cout);

	std::cout << "TestArcLengthCalculation: test 2: circle with 10 segments after remeshing\n";
	RemeshWithDefaultSettings(testCircleCurve, arcLengthCalc);
	const auto arcLengths2 = arcLengthCalc->CalculateArcLengths();
	Geometry::PrintCurveValuesInTopologicalOrder(testCircleCurve, arcLengths2, std::cout);

	std::cout << "TestArcLengthCalculation: test 3: circle with 20 segments after remeshing with factor 0.7\n";
	RemeshWithDefaultSettings(testCircleCurve, arcLengthCalc, 0.7);
	const auto arcLengths3 = arcLengthCalc->CalculateArcLengths();
	Geometry::PrintCurveValuesInTopologicalOrder(testCircleCurve, arcLengths3, std::cout);

	std::cout << "TestArcLengthCalculation: test 4: after remeshing with factor 0.5\n";
	RemeshWithDefaultSettings(testCircleCurve, arcLengthCalc, 0.5);
	const auto arcLengths4 = arcLengthCalc->CalculateArcLengths();
	Geometry::PrintCurveValuesInTopologicalOrder(testCircleCurve, arcLengths4, std::cout);
}

void TestCurve2DRotation()
{
	const auto center = pmp::Point2(200, 400);
	const pmp::Scalar innerRadius = 40.0;
	const std::vector<pmp::Point2> triangleVerticesSmall = {
		pmp::Point2{-0.5, (pmp::Scalar)-sqrtf(3.0) / (pmp::Scalar)6.0} *innerRadius + center,
		pmp::Point2{0.5, (pmp::Scalar)-sqrtf(3.0) / (pmp::Scalar)6.0} *innerRadius + center,
		pmp::Point2{0.0, (pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)3.0} *innerRadius + center
	};
	const unsigned int segments = 30;

	const pmp::Scalar angle = 45 * 1.2; // in degrees
	auto curve = pmp::CurveFactory::sampled_polygon(triangleVerticesSmall, segments, true);
	if (!pmp::write_to_ply(curve, dataOutPath + "curve_BeforeRot.ply"))
		std::cerr << "Error writing curve_BeforeRot.ply!\n";

	curve = GetCurveRotatedAboutCenterPoint(curve, angle);
	if (!pmp::write_to_ply(curve, dataOutPath + "curve_AfterRot.ply"))
		std::cerr << "Error writing curve_AfterRot.ply!\n";
}

void TestSmoothingAdvectionEquilibrium()
{
	// Define the inner and outer circle pairs directly
	const std::vector<std::pair<Geometry::Circle2D, Geometry::Circle2D>> circlePairs{
		{Geometry::Circle2D{pmp::Point2{-3.0, 52.0}, 100.0}, Geometry::Circle2D{pmp::Point2{-3.0, 52.0}, 121.558}},
		//{Geometry::Circle2D{pmp::Point2{-25.0, 8.0}, 0.055}, Geometry::Circle2D{pmp::Point2{-25.0, 8.0}, 0.142831}},
		//{Geometry::Circle2D{pmp::Point2{8.0, 85.0}, 50.0}, Geometry::Circle2D{pmp::Point2{8.0, 85.0}, 292.263}},
		//{Geometry::Circle2D{pmp::Point2{-20.0, 90.0}, 55.0}, Geometry::Circle2D{pmp::Point2{-20.0, 90.0}, 441.436}}
	};

	constexpr unsigned int nVoxelsPerMinDimension = 40;
	constexpr double defaultTimeStep = 0.05;
	constexpr double defaultOffsetFactor = 0.25; // 1.5;
	constexpr unsigned int NTimeSteps = 1800;

	constexpr double innerEpsilonStart = 0.01;
	constexpr double innerEpsilonEnd = 0.001;
	//constexpr double innerEpsilonExpFactor = 0.5;
	constexpr double innerEpsilonStep = -0.001;

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
		const pmp::Scalar minSize = std::min(bboxSize[0], bboxSize[1]);
		const pmp::Scalar cellSize = minSize / nVoxelsPerMinDimension;

		//for (unsigned int epsilonId = 0; double epsConstantValue : Utils::GetExponentialValueRange(innerEpsilonStart, innerEpsilonEnd, innerEpsilonExpFactor))
		for (unsigned int epsilonId = 0; double epsConstantValue : Utils::GetLinearValueRange(innerEpsilonStart, innerEpsilonEnd, innerEpsilonStep))
		{
			const double isoLvlOffsetFactor = defaultOffsetFactor;
			const double fieldIsoLevel = isoLvlOffsetFactor * sqrt(3.0) / 2.0 * static_cast<double>(cellSize);

			std::cout << "Setting up ManifoldEvolutionSettings.\n";

			ManifoldEvolutionSettings strategySettings;
			strategySettings.UseInnerManifolds = true;
			strategySettings.AdvectionInteractWithOtherManifolds = true;
			strategySettings.OuterManifoldEpsilon = [](double distance)
			{
				//return 1.0 * (1.0 - exp(-distance * distance / 1.0));
				return 0.0;
			};
			strategySettings.OuterManifoldEta = [](double distance, double negGradDotNormal)
			{
				//return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
				return 0.0;
			};
			strategySettings.InnerManifoldEpsilon = [&epsConstantValue](double distance)
			{
				return epsConstantValue * TRIVIAL_EPSILON(distance);
			};
			strategySettings.InnerManifoldEta = [](double distance, double negGradDotNormal)
			{
				return 1.0 * distance * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
			};
			strategySettings.TimeStep = defaultTimeStep;
			strategySettings.LevelOfDetail = 3;
			strategySettings.TangentialVelocityWeight = 0.05;

			strategySettings.RemeshingSettings.MinEdgeMultiplier = 0.14;
			strategySettings.RemeshingSettings.UseBackProjection = false;

			strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
			strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = 1.0 * static_cast<pmp::Scalar>(M_PI_2);

			strategySettings.FieldSettings.NVoxelsPerMinDimension = nVoxelsPerMinDimension;
			strategySettings.FieldSettings.FieldIsoLevel = fieldIsoLevel;

			strategySettings.ExportVariableScalarFieldsDimInfo = true;
			strategySettings.ExportVariableVectorFieldsDimInfo = true;

			//std::cout << "Setting up GlobalManifoldEvolutionSettings.\n";

			GlobalManifoldEvolutionSettings globalSettings;
			globalSettings.NSteps = NTimeSteps;
			globalSettings.DoRemeshing = true;
			globalSettings.DetectFeatures = false;
			globalSettings.ExportPerTimeStep = true;
			globalSettings.ExportTargetDistanceFieldAsImage = true;
			globalSettings.ProcedureName = "concentricCircles" + std::to_string(curveId) + "_eps" + std::to_string(epsilonId);

			std::cout << globalSettings.ProcedureName << " with epsilon: " << epsConstantValue << "\n";

			globalSettings.OutputPath = dataOutPath;
			globalSettings.ExportResult = false;

			globalSettings.RemeshingResizeFactor = 0.7;
			globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

			const auto nSegments = static_cast<unsigned int>(pow(2, strategySettings.LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
			auto outerCurve = pmp::CurveFactory::circle(outerCircle.Center, outerCircle.Radius, nSegments);

			const auto nInnerSegments = static_cast<unsigned int>(static_cast<pmp::Scalar>(nSegments) * innerCircle.Radius / outerCircle.Radius * 2);
			auto innerCurve = pmp::CurveFactory::circle(innerCircle.Center, innerCircle.Radius, nInnerSegments);
			//innerCurve.negate_orientation();
			std::vector innerCurves{ innerCurve };

			auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
				strategySettings, outerCurve, innerCurves, nullptr);

			std::cout << "Setting up ManifoldEvolver.\n";


			ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

			std::cout << "ManifoldEvolver::Evolve ... ";

			evolver.Evolve();

			epsilonId++;
		}

		curveId++;
	}
}

void TestImageToDistanceField()
{
	using namespace Geometry;
	// const std::string imgName = "room";
	const std::string imgName = "shape";
	SDF::ImageDistanceField2DSettings dfSettings;
	dfSettings.CellSize = (pmp::Scalar)1.0;
	dfSettings.ImageScaleFactor = (pmp::Scalar)4.0;
	auto imgDf = SDF::ImageDistanceFieldGenerator::Generate(dataDirPath + imgName + ".png", dfSettings);

	ExportScalarGridDimInfo2D(dataOutPath + imgName + ".gdim2d", imgDf);
	constexpr double colorMapPlotScaleFactor = 1.0; // scale the distance field color map down to show more detail
	ExportScalarGrid2DToPNG(dataOutPath + imgName + "_df.png", imgDf,
		Geometry::BilinearInterpolateScalarValue,
		//Geometry::GetNearestNeighborScalarValue2D,
		10, 10, RAINBOW_TO_WHITE_MAP * colorMapPlotScaleFactor);

	auto blurredDf = imgDf;
	ApplyWideGaussianBlur2D(blurredDf);
	const unsigned int nPx = 2;
	const auto gradDF = ComputeNormalizedGradient(blurredDf);
	const auto divNegGradDF = ComputeDivergenceField(gradDF);
	const auto [dMin, dMax] = std::ranges::minmax_element(divNegGradDF.Values());
	std::cout << "divNegGradDF value range: [" << *dMin << ", " << *dMax << "]\n";
	const auto colorMap = AdjustColorMapForZeroMidpoint(SIGN_TEMP_MAP, *dMin, *dMax);
	const double rangeVal = std::min(std::abs(*dMin), *dMax);
	ExportScalarGridDimInfo2D(dataOutPath + imgName + "_divGradDF.gdim2d", divNegGradDF);
	ExportScalarGrid2DToPNG(dataOutPath + imgName + "_divGradDF.png", divNegGradDF,
		Geometry::BilinearInterpolateScalarValue,
		//Geometry::GetNearestNeighborScalarValue2D,
		nPx, nPx, colorMap);
	constexpr bool verticalLegend = true;
	constexpr unsigned int resolutionFactor = 4;
	constexpr unsigned int legendPxHeight = (verticalLegend ? 400 : 100) * resolutionFactor;
	constexpr unsigned int legendPxWidth = (verticalLegend ? 100 : 600) * resolutionFactor;
	ExportColorScaleToPNG(
		dataOutPath + imgName + "_divGradDF_Scale.png",
		-rangeVal,
		rangeVal,
		colorMap,
		legendPxHeight, legendPxWidth);
}

void TestImageSegmentation()
{
	const std::vector<std::string> imageNames{
		"room",
		//"U64",
	};

	constexpr size_t nCurvePts = 50;
	constexpr pmp::Scalar imageScale = 2.0;
	constexpr size_t nPixels = 64 * imageScale;
	std::map<std::string, std::pair<pmp::ManifoldCurve2D, pmp::ManifoldCurve2D>> circles{
		{"room", 
		{pmp::CurveFactory::circle(pmp::Point2{ nPixels * 0.49, nPixels * 0.5 }, 1.7 * nPixels / 4.0, nCurvePts),
		 pmp::CurveFactory::circle(pmp::Point2{ nPixels * 0.49, nPixels * 0.47 }, 0.15 * nPixels / 4.0, nCurvePts)}},
		 {"U64",
		{pmp::CurveFactory::circle(pmp::Point2{ nPixels * 0.49, nPixels * 0.5 }, 1.8 * nPixels / 4.0, nCurvePts),
		 pmp::CurveFactory::circle(pmp::Point2{ nPixels * 0.49, nPixels * 0.63 }, 0.2 * nPixels / 4.0, nCurvePts)}}
	};

	constexpr double minDistancePercentageEpsilon = 0.01;
	constexpr double minDistancePercentageEta = 0.01;

	// Prepare the settings for the evolver
	ManifoldEvolutionSettings strategySettings;
	strategySettings.UseInnerManifolds = true;

	strategySettings.AdvectionInteractWithOtherManifolds = true;
	// use PreComputeAdvectionDiffusionParams?
	strategySettings.OuterManifoldEpsilon.Bind(minDistancePercentageEpsilon, [](double distance)
	{
		//return 0.0;

		return 0.2 * (1.0 - exp(-distance * distance / 1.0));
	});
	strategySettings.OuterManifoldEta.Bind(minDistancePercentageEta, [](double distance, double negGradDotNormal)
	{
		//return 0.0;

		return -1.0 * distance * (std::abs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));

		//return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		//return 1.0 * (1.0 - exp(-distance * distance / 0.5)) * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
	});
	strategySettings.InnerManifoldEpsilon.Bind(minDistancePercentageEpsilon, [](double distance)
	{
		//return 0.0;

		return 0.001 * TRIVIAL_EPSILON(distance);
	});
	strategySettings.InnerManifoldEta.Bind(minDistancePercentageEta, [](double distance, double negGradDotNormal)
	{
		//return 0.0;

		return 0.8 * distance * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		//return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		// return 1.0 * distance * (std::fabs(-negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		// return 1.0 * distance * (negGradDotNormal - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		//return 1.0 * distance * (std::fabs(negGradDotNormal) - 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
		//return -1.0 * (1.0 - exp(-distance * distance / 0.5)) * (std::fabs(negGradDotNormal) + 1.0 * sqrt(1.0 - negGradDotNormal * negGradDotNormal));
	});
	strategySettings.LevelOfDetail = 4;

	strategySettings.TimeStep = 0.05;
	strategySettings.TangentialVelocityWeight = 0.05;
	strategySettings.RemeshingSettings.MinEdgeMultiplier = 1.0;
	strategySettings.RemeshingSettings.UseBackProjection = true;
	strategySettings.FeatureSettings.PrincipalCurvatureFactor = 3.2;
	strategySettings.FeatureSettings.CriticalMeanCurvatureAngle = static_cast<pmp::Scalar>(M_PI_2);
	strategySettings.FieldSettings.NVoxelsPerMinDimension = 40;

	//strategySettings.DistanceSelection = DistanceSelectionType::QuadricBlend;
	//strategySettings.DistanceBlendingRadius = 5.0;
	strategySettings.NormalActivation.On = true;
	strategySettings.NormalActivation.TargetDFCriticalRadius = 10.0;
	strategySettings.NormalActivation.ManifoldCriticalRadius = 15.0;
	strategySettings.NormalActivation.NPointsFromCriticalBound = 4;

	strategySettings.ExportVariableScalarFieldsDimInfo = true;
	strategySettings.ExportVariableVectorFieldsDimInfo = true;

	//strategySettings.DiagSettings.LogOuterManifoldEpsilon = true;
	//strategySettings.DiagSettings.LogInnerManifoldsEpsilon = true;
	//strategySettings.DiagSettings.LogOuterManifoldEta = true;
	//strategySettings.DiagSettings.LogInnerManifoldsEta = true;

	// Global settings
	GlobalManifoldEvolutionSettings globalSettings;
	globalSettings.NSteps = 2500;
	//globalSettings.NSteps = 218;
	//globalSettings.NSteps = 1000;
	//globalSettings.NSteps = 500;
	//globalSettings.NSteps = 20;
	globalSettings.DoRemeshing = true;
	globalSettings.DetectFeatures = false;
	globalSettings.ExportPerTimeStep = true;
	globalSettings.ExportTargetDistanceFieldAsImage = true;
	globalSettings.OutputPath = dataOutPath;
	globalSettings.ExportResult = false;

	globalSettings.RemeshingResizeFactor = 0.7;
	globalSettings.RemeshingResizeTimeIds = GetRemeshingAdjustmentTimeIndices();

	for (const auto& imgName : imageNames)
	{
		std::cout << "==================================================================\n";
		std::cout << "Image: " << imgName << ".png\n";
		std::cout << "------------------------------------------------------------------\n";

		globalSettings.ProcedureName = imgName + "Segment";

		SDF::ImageDistanceField2DSettings dfSettings;
		dfSettings.CellSize = (pmp::Scalar)1.0;
		dfSettings.ImageScaleFactor = imageScale;
		auto imgDf = SDF::ImageDistanceFieldGenerator::Generate(dataDirPath + imgName + ".png", dfSettings);

		//ExportScalarGridDimInfo2D(dataOutPath + imgName + ".gdim2d", imgDf);
		//constexpr double colorMapPlotScaleFactor = 1.0; // scale the distance field color map down to show more detail
		//ExportScalarGrid2DToPNG(dataOutPath + imgName + "_df.png", imgDf,
		//	Geometry::BilinearInterpolateScalarValue,
		//	//Geometry::GetNearestNeighborScalarValue2D,
		//	10, 10, RAINBOW_TO_WHITE_MAP * colorMapPlotScaleFactor);

		if (!circles.contains(imgName))
			continue;

		auto [outerCurve, innerCurve] = circles.at(imgName);
		RemeshWithDefaultSettings(outerCurve, nullptr, 1.0);
		RemeshWithDefaultSettings(innerCurve, nullptr, 1.0);
		std::vector innerCurves{ innerCurve };

		strategySettings.FieldSettings.FieldIsoLevel = imgDf.CellSize() * imageScale * 0.1;

		auto curveStrategy = std::make_shared<CustomManifoldCurveEvolutionStrategy>(
			strategySettings, outerCurve, innerCurves, nullptr, std::make_shared<Geometry::ScalarGrid2D>(imgDf));

		ManifoldEvolver evolver(globalSettings, std::move(curveStrategy));

		evolver.Evolve();

	}
}

void TestPointCloudGaps()
{
	const std::map<std::string, pmp::ManifoldCurve2D> targetCurves{
	{ "IncompleteCircle", pmp::CurveFactory::circle(pmp::Point2{-3.0, 52.0}, 35.0, 25, M_PI_2, 2.0 * M_PI)},
	{ "SineDeformedIncompleteCircle", pmp::CurveFactory::sine_deformed_circle(pmp::Point2{-3.0, 52.0}, 35.0, 25, 7.0, 4.0, M_PI_2, 2.0 * M_PI)},
	{ "IncompleteChamferedRectangle", pmp::CurveFactory::sampled_polygon({
		pmp::Point2{-30.0, -35.0} + pmp::Point2{-3.0, 52.0},
		pmp::Point2{30.0, -35.0} + pmp::Point2{-3.0, 52.0},
		pmp::Point2{30.0, 35.0} + pmp::Point2{-3.0, 52.0},
		pmp::Point2{-30.0, 35.0} + pmp::Point2{-3.0, 52.0}}, 30, true, false)},
	{ "IncompleteChamferedTriangle", pmp::CurveFactory::sampled_polygon({
		pmp::Point2{-0.5, (pmp::Scalar)-sqrtf(3.0) / (pmp::Scalar)6.0} *120.0 + pmp::Point2{-3.0, 52.0},
		pmp::Point2{0.5, (pmp::Scalar)-sqrtf(3.0) / (pmp::Scalar)6.0} *120.0 + pmp::Point2{-3.0, 52.0},
		pmp::Point2{0.0, (pmp::Scalar)sqrtf(3.0) / (pmp::Scalar)3.0} *120.0 + pmp::Point2{-3.0, 52.0}}, 30, true, false)}
	};

	for (const auto& [ptCloudName, curve] : targetCurves)
	{
		std::cout << "========================================================\n";
		std::cout << "   Point cloud gap detection 2D for: " << ptCloudName << " ... \n";
		std::cout << " ------------------------------------------------------ \n";

		const auto& pts2D = curve.positions();
		const pmp::BoundingBox2 bbox{ pts2D };

		const auto bdPtIds = Geometry::GetBoundaryPointsOfPointCloudGaps2D(pts2D);
	}
}

void TestNormalActivation()
{
	const auto pathNormalActivationTarget00Vertices = ParsePolygonalSVGPath(svgPathNormalActivationTarget00);
	auto pathNormalActivationTarget00Curve = pmp::CurveFactory::sampled_polygon(pathNormalActivationTarget00Vertices, 50, false, false);
	RemeshWithDefaultSettings(pathNormalActivationTarget00Curve);
	if (!pmp::write_to_ply(pathNormalActivationTarget00Curve, dataOutPath + "pathNormalActivationTarget00Curve.ply"))
		std::cerr << "Error writing pathNormalActivationTarget00Curve.ply!\n";

	const auto pathNormalActivationTarget01Vertices = ParsePolygonalSVGPath(svgPathNormalActivationTarget01);
	auto pathNormalActivationTarget01Curve = pmp::CurveFactory::sampled_polygon(pathNormalActivationTarget01Vertices, 50, false, false);
	RemeshWithDefaultSettings(pathNormalActivationTarget01Curve);
	if (!pmp::write_to_ply(pathNormalActivationTarget01Curve, dataOutPath + "pathNormalActivationTarget01Curve.ply"))
		std::cerr << "Error writing pathNormalActivationTarget01Curve.ply!\n";

	const auto pathNormalActivationInner0Vertices = ParsePolygonalSVGPath(svgPathNormalActivationInner0);
	auto pathNormalActivationInner0Curve = pmp::CurveFactory::sampled_polygon(pathNormalActivationInner0Vertices, 60, false);
	RemeshWithDefaultSettings(pathNormalActivationInner0Curve);
	if (!pmp::write_to_ply(pathNormalActivationInner0Curve, dataOutPath + "pathNormalActivationInner0Curve.ply"))
		std::cerr << "Error writing pathNormalActivationInner0Curve.ply!\n";

	const auto pathNormalActivationOuter0Vertices = ParsePolygonalSVGPath(svgPathNormalActivationOuter0);
	auto pathNormalActivationOuter0Curve = pmp::CurveFactory::sampled_polygon(pathNormalActivationOuter0Vertices, 100, false);
	RemeshWithDefaultSettings(pathNormalActivationOuter0Curve);
	if (!pmp::write_to_ply(pathNormalActivationOuter0Curve, dataOutPath + "pathNormalActivationOuter0Curve.ply"))
		std::cerr << "Error writing pathNormalActivationOuter0Curve.ply!\n";

	std::vector<pmp::Point2> targetPtCloud;
	targetPtCloud.reserve(pathNormalActivationTarget00Curve.n_vertices() + pathNormalActivationTarget01Curve.n_vertices());
	const auto& pts0 = pathNormalActivationTarget00Curve.positions();
	const auto& pts1 = pathNormalActivationTarget01Curve.positions();
	targetPtCloud.insert(targetPtCloud.end(), pts0.begin(), pts0.end());
	targetPtCloud.insert(targetPtCloud.end(), pts1.begin(), pts1.end());

	constexpr unsigned int nVoxelsPerMinDimension = 40;

	const pmp::BoundingBox2 ptCloudBBox(targetPtCloud);
	const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
	const pmp::Scalar minSize = std::min(ptCloudBBoxSize[0], ptCloudBBoxSize[1]);
	const pmp::Scalar maxSize = std::max(ptCloudBBoxSize[0], ptCloudBBoxSize[1]);
	const pmp::Scalar cellSize = minSize / static_cast<pmp::Scalar>(nVoxelsPerMinDimension);
	const SDF::PointCloudDistanceField2DSettings dfSettings{
		cellSize,
		0.6,
		DBL_MAX
	};
	const auto targetDf = std::make_shared<Geometry::ScalarGrid2D>(SDF::PlanarPointCloudDistanceFieldGenerator::Generate(targetPtCloud, dfSettings));

	const SDF::DistanceField2DSettings curveDFSettings{
		cellSize,
		0.6,
		DBL_MAX,
		SDF::KDTreeSplitType::Center,
		SDF::SignComputation2D::None,
		SDF::PreprocessingType2D::Quadtree
	};

	NormalActivationSettings naSettings;
	naSettings.On = true;
	naSettings.TargetDFCriticalRadius = 20.0;
	naSettings.ManifoldCriticalRadius = 25.0;
	naSettings.NPointsFromCriticalBound = 1;

	const Geometry::ManifoldCurve2DAdapter outerCurveAdapter(std::make_shared<pmp::ManifoldCurve2D>(pathNormalActivationOuter0Curve));
	const auto outerCurveDf = std::make_shared<Geometry::ScalarGrid2D>(SDF::PlanarDistanceFieldGenerator::Generate(outerCurveAdapter, curveDFSettings));

	const Geometry::ManifoldCurve2DAdapter innerCurveAdapter(std::make_shared<pmp::ManifoldCurve2D>(pathNormalActivationInner0Curve));
	const auto innerCurveDf = std::make_shared<Geometry::ScalarGrid2D>(SDF::PlanarDistanceFieldGenerator::Generate(innerCurveAdapter, curveDFSettings, outerCurveDf->Box()));

	const auto vOuterGap = Geometry::GetVerticesWithinMinDistance(pathNormalActivationOuter0Curve, 
		{ innerCurveDf }, naSettings.ManifoldCriticalRadius, "v:gap_activated", Geometry::BilinearInterpolateScalarValue);
	const auto vInnerGap = Geometry::GetVerticesWithinMinDistance(pathNormalActivationInner0Curve,
		{ outerCurveDf }, naSettings.ManifoldCriticalRadius, "v:gap_activated", Geometry::BilinearInterpolateScalarValue);

	std::vector<pmp::Point2> outerCurveGapActivated;
	for (const auto v : pathNormalActivationOuter0Curve.vertices())
	{
		if (!vOuterGap[v])
			continue;

		outerCurveGapActivated.push_back(pathNormalActivationOuter0Curve.position(v));
	}

	std::vector<pmp::Point2> innerCurveGapActivated;
	for (const auto v : pathNormalActivationInner0Curve.vertices())
	{
		if (!vInnerGap[v])
			continue;

		innerCurveGapActivated.push_back(pathNormalActivationInner0Curve.position(v));
	}

	if (!Geometry::Export2DPointCloudToPLY(outerCurveGapActivated, dataOutPath + "outerCurveGapActivated.ply"))
	{
		std::cerr << "Export2DPointCloudToPLY: internal error during export!\n";
	}

	if (!Geometry::Export2DPointCloudToPLY(innerCurveGapActivated, dataOutPath + "innerCurveGapActivated.ply"))
	{
		std::cerr << "Export2DPointCloudToPLY: internal error during export!\n";
	}

	const auto [outerNextBoundary, outerPrevBoundary] = GetNearestGapBoundaryVertices(pathNormalActivationOuter0Curve,
		targetDf, { innerCurveDf }, Geometry::BilinearInterpolateScalarValue, naSettings);
	if (!outerNextBoundary || !outerPrevBoundary)
		return;

	const auto [innerNextBoundary, innerPrevBoundary] = GetNearestGapBoundaryVertices(pathNormalActivationInner0Curve,
		targetDf, { outerCurveDf }, Geometry::BilinearInterpolateScalarValue, naSettings);
	if (!innerNextBoundary || !innerPrevBoundary)
		return;

	VertexValueLogger<pmp::ManifoldCurve2D> curveLogger;
	curveLogger.AddManifold(&pathNormalActivationOuter0Curve);
	curveLogger.AddManifold(&pathNormalActivationInner0Curve);

	std::unordered_map<pmp::ManifoldCurve2D*, std::shared_ptr<pmp::EvolvingArcLengthCalculator>> arcLengthCalculators{
		{&pathNormalActivationOuter0Curve, std::make_shared<pmp::EvolvingArcLengthCalculator>(pathNormalActivationOuter0Curve) },
		{&pathNormalActivationInner0Curve, std::make_shared<pmp::EvolvingArcLengthCalculator>(pathNormalActivationInner0Curve) }
	};

	curveLogger.Init(dataOutPath + "activationBoundary_log.json");
	curveLogger.StartNewTimeStep(1);

	const auto outerArcLengths = arcLengthCalculators[&pathNormalActivationOuter0Curve]->CalculateArcLengths();
	if (!outerArcLengths.empty())
	{
		for (const auto v : pathNormalActivationOuter0Curve.vertices())
		{
			curveLogger.LogValue(&pathNormalActivationOuter0Curve, "arcLength", v.idx(), outerArcLengths[v.idx()]);

			if ((*outerNextBoundary)[v].is_valid())
			{
				const auto nextBdArcLength = outerArcLengths[(*outerNextBoundary)[v].idx()];
				curveLogger.LogValue(&pathNormalActivationOuter0Curve, "nextBoundaryArcLength", v.idx(), nextBdArcLength);
			}
			else
			{
				curveLogger.LogValue(&pathNormalActivationOuter0Curve, "nextBoundaryArcLength", v.idx(), -1.0);
			}

			if ((*outerPrevBoundary)[v].is_valid())
			{
				const auto prevBdArcLength = outerArcLengths[(*outerPrevBoundary)[v].idx()];
				curveLogger.LogValue(&pathNormalActivationOuter0Curve, "prevBoundaryArcLength", v.idx(), prevBdArcLength);
			}
			else
			{
				curveLogger.LogValue(&pathNormalActivationOuter0Curve, "prevBoundaryArcLength", v.idx(), -1.0);
			}
		}
	}
	const auto innerArcLengths = arcLengthCalculators[&pathNormalActivationInner0Curve]->CalculateArcLengths();
	if (!innerArcLengths.empty())
	{
		for (const auto v : pathNormalActivationInner0Curve.vertices())
		{
			curveLogger.LogValue(&pathNormalActivationInner0Curve, "arcLength", v.idx(), innerArcLengths[v.idx()]);

			if ((*innerNextBoundary)[v].is_valid())
			{
				const auto nextBdArcLength = innerArcLengths[(*innerNextBoundary)[v].idx()];
				curveLogger.LogValue(&pathNormalActivationInner0Curve, "nextBoundaryArcLength", v.idx(), nextBdArcLength);
			}
			else
			{
				curveLogger.LogValue(&pathNormalActivationInner0Curve, "nextBoundaryArcLength", v.idx(), -1.0);
			}

			if ((*innerPrevBoundary)[v].is_valid())
			{
				const auto prevBdArcLength = innerArcLengths[(*innerPrevBoundary)[v].idx()];
				curveLogger.LogValue(&pathNormalActivationInner0Curve, "prevBoundaryArcLength", v.idx(), prevBdArcLength);
			}
			else
			{
				curveLogger.LogValue(&pathNormalActivationInner0Curve, "prevBoundaryArcLength", v.idx(), -1.0);
			}
		}
	}

	curveLogger.Save(false);
}

void TestGapSpecificBehaviorForRealData()
{
	const std::string& imgName = "room";
	constexpr pmp::Scalar imageScale = 2.0;

	SDF::ImageDistanceField2DSettings dfSettings;
	dfSettings.CellSize = (pmp::Scalar)1.0;
	dfSettings.ImageScaleFactor = imageScale;
	const auto targetDf = std::make_shared<Geometry::ScalarGrid2D>(SDF::ImageDistanceFieldGenerator::Generate(dataDirPath + imgName + ".png", dfSettings));
	if (!targetDf)
	{
		std::cerr << "SDF::ImageDistanceFieldGenerator::Generate: error!\n";
		return;
	}

	Geometry::ExportScalarGridDimInfo2D(dataOutPath + imgName + ".gdim2d", *targetDf);
	constexpr double colorMapPlotScaleFactor = 1.0; // scale the distance field color map down to show more detail
	ExportScalarGrid2DToPNG(dataOutPath + imgName + "_df.png", *targetDf,
		Geometry::BilinearInterpolateScalarValue,
		//Geometry::GetNearestNeighborScalarValue2D,
		10, 10, Geometry::RAINBOW_TO_WHITE_MAP * colorMapPlotScaleFactor);

	const std::vector<std::pair<size_t, size_t>> timeStepPairs{
		{75, 76},
		{117, 118}
	};

	const pmp::Scalar cellSize = targetDf->CellSize();
	const SDF::DistanceField2DSettings curveDFSettings{
		cellSize,
		1.0,
		DBL_MAX,
		SDF::KDTreeSplitType::Center,
		SDF::SignComputation2D::None,
		SDF::PreprocessingType2D::Quadtree
	};

	NormalActivationSettings naSettings;
	naSettings.On = true;
	naSettings.TargetDFCriticalRadius = 10.0;
	naSettings.ManifoldCriticalRadius = 15.0;
	naSettings.NPointsFromCriticalBound = 4;

	for (const auto [ts0, ts1] : timeStepPairs)
	{
		const auto coverTimeStep = [&](const size_t& ts) {
			pmp::ManifoldCurve2D outerCurve;
			if (!pmp::read_from_ply(outerCurve, dataOutPath + "roomSegment_Outer_Evol_" + std::to_string(ts) + ".ply"))
			{
				std::cerr << "pmp::read_from_ply: internal error!\n";
				return false;
			}

			pmp::ManifoldCurve2D innerCurve;
			if (!pmp::read_from_ply(innerCurve, dataOutPath + "roomSegment_Inner0_Evol_" + std::to_string(ts) + ".ply"))
			{
				std::cerr << "pmp::read_from_ply: internal error!\n";
				return false;
			}

			const Geometry::ManifoldCurve2DAdapter outerCurveAdapter(std::make_shared<pmp::ManifoldCurve2D>(outerCurve));
			const auto outerCurveDf = std::make_shared<Geometry::ScalarGrid2D>(SDF::PlanarDistanceFieldGenerator::Generate(outerCurveAdapter, curveDFSettings));

			const Geometry::ManifoldCurve2DAdapter innerCurveAdapter(std::make_shared<pmp::ManifoldCurve2D>(innerCurve));
			const auto innerCurveDf = std::make_shared<Geometry::ScalarGrid2D>(SDF::PlanarDistanceFieldGenerator::Generate(innerCurveAdapter, curveDFSettings, outerCurveDf->Box()));

			const auto vOuterGap = Geometry::GetVerticesWithinMinDistance(outerCurve,
				{ innerCurveDf }, naSettings.ManifoldCriticalRadius, "v:gap_activated", Geometry::BilinearInterpolateScalarValue);
			const auto vInnerGap = Geometry::GetVerticesWithinMinDistance(innerCurve,
				{ outerCurveDf }, naSettings.ManifoldCriticalRadius, "v:gap_activated", Geometry::BilinearInterpolateScalarValue);

			std::vector<pmp::Point2> outerCurveGapActivated;
			for (const auto v : outerCurve.vertices())
			{
				if (!vOuterGap[v])
					continue;

				outerCurveGapActivated.push_back(outerCurve.position(v));
			}

			std::vector<pmp::Point2> innerCurveGapActivated;
			for (const auto v : innerCurve.vertices())
			{
				if (!vInnerGap[v])
					continue;

				innerCurveGapActivated.push_back(innerCurve.position(v));
			}

			if (!Geometry::Export2DPointCloudToPLY(outerCurveGapActivated, dataOutPath + "outerCurveGapActivated_step" + std::to_string(ts) + ".ply"))
			{
				std::cerr << "Export2DPointCloudToPLY: internal error during export!\n";
				return false;
			}

			if (!Geometry::Export2DPointCloudToPLY(innerCurveGapActivated, dataOutPath + "innerCurveGapActivated_step" + std::to_string(ts) + ".ply"))
			{
				std::cerr << "Export2DPointCloudToPLY: internal error during export!\n";
				return false;
			}

			const auto [outerNextBoundary, outerPrevBoundary] = GetNearestGapBoundaryVertices(outerCurve,
				targetDf, { innerCurveDf }, Geometry::BilinearInterpolateScalarValue, naSettings);
			if (!outerNextBoundary || !outerPrevBoundary)
				return false;

			const auto [innerNextBoundary, innerPrevBoundary] = GetNearestGapBoundaryVertices(innerCurve,
				targetDf, { outerCurveDf }, Geometry::BilinearInterpolateScalarValue, naSettings);
			if (!innerNextBoundary || !innerPrevBoundary)
				return false;

			VertexValueLogger<pmp::ManifoldCurve2D> curveLogger;
			curveLogger.AddManifold(&outerCurve);
			curveLogger.AddManifold(&innerCurve);

			std::unordered_map<pmp::ManifoldCurve2D*, std::shared_ptr<pmp::EvolvingArcLengthCalculator>> arcLengthCalculators{
				{&outerCurve, std::make_shared<pmp::EvolvingArcLengthCalculator>(outerCurve) },
				{&innerCurve, std::make_shared<pmp::EvolvingArcLengthCalculator>(innerCurve) }
			};

			curveLogger.Init(dataOutPath + "gapBehavior_step" + std::to_string(ts) + "_log.json");
			curveLogger.StartNewTimeStep(1);

			const auto outerArcLengths = arcLengthCalculators[&outerCurve]->CalculateArcLengths();
			if (!outerArcLengths.empty())
			{
				for (const auto v : outerCurve.vertices())
				{
					curveLogger.LogValue(&outerCurve, "arcLength", v.idx(), outerArcLengths[v.idx()]);

					if ((*outerNextBoundary)[v].is_valid())
					{
						const auto nextBdArcLength = outerArcLengths[(*outerNextBoundary)[v].idx()];
						curveLogger.LogValue(&outerCurve, "nextBoundaryArcLength", v.idx(), nextBdArcLength);
					}
					else
					{
						curveLogger.LogValue(&outerCurve, "nextBoundaryArcLength", v.idx(), -1.0);
					}

					if ((*outerPrevBoundary)[v].is_valid())
					{
						const auto prevBdArcLength = outerArcLengths[(*outerPrevBoundary)[v].idx()];
						curveLogger.LogValue(&outerCurve, "prevBoundaryArcLength", v.idx(), prevBdArcLength);
					}
					else
					{
						curveLogger.LogValue(&outerCurve, "prevBoundaryArcLength", v.idx(), -1.0);
					}
				}
			}
			const auto innerArcLengths = arcLengthCalculators[&innerCurve]->CalculateArcLengths();
			if (!innerArcLengths.empty())
			{
				for (const auto v : innerCurve.vertices())
				{
					curveLogger.LogValue(&innerCurve, "arcLength", v.idx(), innerArcLengths[v.idx()]);

					if ((*innerNextBoundary)[v].is_valid())
					{
						const auto nextBdArcLength = innerArcLengths[(*innerNextBoundary)[v].idx()];
						curveLogger.LogValue(&innerCurve, "nextBoundaryArcLength", v.idx(), nextBdArcLength);
					}
					else
					{
						curveLogger.LogValue(&innerCurve, "nextBoundaryArcLength", v.idx(), -1.0);
					}

					if ((*innerPrevBoundary)[v].is_valid())
					{
						const auto prevBdArcLength = innerArcLengths[(*innerPrevBoundary)[v].idx()];
						curveLogger.LogValue(&innerCurve, "prevBoundaryArcLength", v.idx(), prevBdArcLength);
					}
					else
					{
						curveLogger.LogValue(&innerCurve, "prevBoundaryArcLength", v.idx(), -1.0);
					}
				}
			}

			curveLogger.Save(false);
			return true;
		};

		if (!coverTimeStep(ts0))
		{
			std::cerr << "coverTimeStep: internal error for step " << ts0 << "\n";
			continue;
		}
		if (!coverTimeStep(ts1))
		{
			std::cerr << "coverTimeStep: internal error for step " << ts1 << "\n";
			continue;
		}
	}
}
