﻿/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
#pragma once

#ifdef WIN32
#include <windows.h>
#include <Psapi.h>
#endif
#include <cstdio>

#include "Poisson/MyTime.h"
#include "Poisson/MarchingCubes.h"
#include "Poisson/Octree.h"
#include "Poisson/SparseMatrix.h"
#include "Poisson/CmdLineParser.h"
#include "Poisson/PPolynomial.h"

#include "Poisson/MultiGridOctreeData.h"

#include "VCGAdapter.h"
#include "vcg/complex/algorithms/update/position.h"

inline void DumpOutput(const char* format, ...)
{
	char buf[4096];
	va_list marker;
	va_start(marker, format);

	vsprintf(buf, format, marker);
	va_end(marker);

	//std::cerr << buf;
}

inline void DumpOutput2(std::vector< char* >&, const char* format, ...)
{
	char buf[4096];
	va_list marker;
	va_start(marker, format);

	vsprintf(buf, format, marker);
	va_end(marker);
	//std::cerr << buf;
}

#if defined( _WIN32 ) || defined( _WIN64 )
inline double PeakMemoryUsageMB(void)
{
	HANDLE h = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	return GetProcessMemoryInfo(h, &pmc, sizeof(pmc)) ? ((double)pmc.PeakWorkingSetSize) / (1 << 20) : 0;
}

inline double to_seconds(const FILETIME& ft)
{
	const double low_to_sec = 100e-9; // 100 nanoseconds
	const double high_to_sec = low_to_sec * 4294967296.0;
	return ft.dwLowDateTime * low_to_sec + ft.dwHighDateTime * high_to_sec;
}
#endif // _WIN32 || _WIN64

template< class Real >
struct OctreeProfiler
{
	Octree< Real >& tree;
	double t;

	OctreeProfiler(Octree< Real >& t) : tree(t) { ; }
	void start(void) { t = Time(), tree.resetLocalMemoryUsage(); }
	void print(const char* header) const
	{
		tree.memoryUsage();
//#if defined( _WIN32 ) || defined( _WIN64 )
//		if (header) printf("%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n", header, Time() - t, tree.localMemoryUsage(), tree.maxMemoryUsage(), PeakMemoryUsageMB());
//		else         printf("%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n", Time() - t, tree.localMemoryUsage(), tree.maxMemoryUsage(), PeakMemoryUsageMB());
//#else // !_WIN32 && !_WIN64
//		if (header) printf("%s %9.1f (s), %9.1f (MB) / %9.1f (MB)\n", header, Time() - t, tree.localMemoryUsage(), tree.maxMemoryUsage());
//		else         printf("%9.1f (s), %9.1f (MB) / %9.1f (MB)\n", Time() - t, tree.localMemoryUsage(), tree.maxMemoryUsage());
//#endif // _WIN32 || _WIN64
	}
	void dumpOutput(const char* header) const
	{
		tree.memoryUsage();
//#if defined( _WIN32 ) || defined( _WIN64 )
//		if (header) DumpOutput("%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n", header, Time() - t, tree.localMemoryUsage(), tree.maxMemoryUsage(), PeakMemoryUsageMB());
//		else         DumpOutput("%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n", Time() - t, tree.localMemoryUsage(), tree.maxMemoryUsage(), PeakMemoryUsageMB());
//#else // !_WIN32 && !_WIN64
//		if (header) DumpOutput("%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n", header, Time() - t, tree.localMemoryUsage(), tree.maxMemoryUsage());
//		else         DumpOutput("%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n", Time() - t, tree.localMemoryUsage(), tree.maxMemoryUsage());
//#endif // _WIN32 || _WIN64
	}
	void dumpOutput2(std::vector< char* >& comments, const char* header) const
	{
		tree.memoryUsage();
//#if defined( _WIN32 ) || defined( _WIN64 )
//		if (header) DumpOutput2(comments, "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n", header, Time() - t, tree.localMemoryUsage(), tree.maxMemoryUsage(), PeakMemoryUsageMB());
//		else         DumpOutput2(comments, "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n", Time() - t, tree.localMemoryUsage(), tree.maxMemoryUsage(), PeakMemoryUsageMB());
//#else // !_WIN32 && !_WIN64
//		if (header) DumpOutput2(comments, "%s %9.1f (s), %9.1f (MB) / %9.1f (MB)\n", header, Time() - t, tree.localMemoryUsage(), tree.maxMemoryUsage());
//		else         DumpOutput2(comments, "%9.1f (s), %9.1f (MB) / %9.1f (MB)\n", Time() - t, tree.localMemoryUsage(), tree.maxMemoryUsage());
//#endif // _WIN32 || _WIN64
	}
};

template <class Real>
class PoissonParam
{
public:
	int MaxDepthVal;
	int MaxSolveDepthVal;
	int KernelDepthVal;
	int MinDepthVal;
	int FullDepthVal;
	Real SamplesPerNodeVal;
	Real ScaleVal;
	bool ConfidenceFlag;
	bool CleanFlag;
	bool DensityFlag;
	Real PointWeightVal;
	int AdaptiveExponentVal;
	int BoundaryTypeVal;
	bool CompleteFlag;
	bool NonManifoldFlag;
	bool ShowResidualFlag;
	int CGDepthVal;
	int ItersVal;
	Real CSSolverAccuracyVal;

	bool VerboseFlag;
	int ThreadsVal;
	bool LinearFitFlag;
	float LowResIterMultiplierVal;
	float ColorVal;

	PoissonParam() {
		MaxDepthVal = 8;
		MaxSolveDepthVal = -1;
		KernelDepthVal = -1;
		MinDepthVal = 0;
		FullDepthVal = 5;
		SamplesPerNodeVal = 1.5f;
		ScaleVal = 1.1f;
		ConfidenceFlag = false;
		CleanFlag = false;
		DensityFlag = false;
		PointWeightVal = 4.0f;
		AdaptiveExponentVal = 1;
		BoundaryTypeVal = 1;
		CompleteFlag = false;
		NonManifoldFlag = false;
		ShowResidualFlag = false;
		CGDepthVal = 0;
		ItersVal = 8;
		CSSolverAccuracyVal = 1e-3f;

		VerboseFlag = false;
		ThreadsVal = omp_get_num_procs();
		LinearFitFlag = false;
		LowResIterMultiplierVal = 1.f;
		ColorVal = 16.0f;
	}
};

template< class Real >
class MeshModelPointStream : public OrientedPointStreamWithData< Real, Point3m >
{
	VCG_Mesh& _m;
	size_t _curPos;
public:
	MeshModelPointStream(VCG_Mesh& m) :_m(m), _curPos(0)
	{
		vcg::tri::RequireCompactness(m);
	}

	~MeshModelPointStream(void) {}

	void reset(void) { _curPos = 0; }

	bool nextPoint(OrientedPoint3D< Real >& pt, Point3m& d)
	{
		if (_curPos >= (unsigned int)_m.vn)
			return false;
		Point3m& nn = _m.vert[_curPos].N();
		Point3m tp = _m.vert[_curPos].P();
		Point4m np = Point4m(nn[0], nn[1], nn[2], 0);

		pt.p[0] = tp[0];
		pt.p[1] = tp[1];
		pt.p[2] = tp[2];
		pt.n[0] = np[0];
		pt.n[1] = np[1];
		pt.n[2] = np[2];

		d[0] = Real(_m.vert[_curPos].C()[0]);
		d[1] = Real(_m.vert[_curPos].C()[1]);
		d[2] = Real(_m.vert[_curPos].C()[2]);

		++_curPos;
		return true;
	}
};

template< class Real>
XForm4x4<Real> GetPointStreamScale(vcg::Box3<Real>& bb, float expFact)
{
	//printf("bbox %f %f %f - %f %f %f ", bb.min[0], bb.min[1], bb.min[2], bb.max[0], bb.max[1], bb.max[2]);
	Real scale = bb.Dim()[bb.MaxDim()] * expFact;
	Point3m center = bb.Center();
	for (int i = 0; i < 3; i++)
		center[i] -= scale / 2;
	XForm4x4< Real > tXForm = XForm4x4< Real >::Identity(), sXForm = XForm4x4< Real >::Identity();
	for (int i = 0; i < 3; i++) {
		sXForm(i, i) = (Real)(1. / scale);
		tXForm(3, i) = -center[i];
	}
	return sXForm * tXForm;
}

template< class Real, int Degree, BoundaryType BType, class Vertex >
int _Execute(
	OrientedPointStream< Real >* pointStream,
	Box3m bb, VCG_Mesh& pm,
	PoissonParam<Real>& pp,
	vcg::CallBackPos* cb)
{
	typedef typename Octree< Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;
	typedef typename Octree< Real >::template InterpolationInfo< false > InterpolationInfo;
	typedef OrientedPointStreamWithData< Real, Point3D< Real > > PointStreamWithData;
	typedef TransformedOrientedPointStreamWithData< Real, Point3D< Real > > XPointStreamWithData;
	Reset< Real >();
	std::vector< char* > comments;

	XForm4x4< Real > xForm = GetPointStreamScale(bb, pp.ScaleVal);
	XForm4x4< Real > iXForm = xForm.inverse();
	//DumpOutput2(comments, "Running Screened Poisson Reconstruction (Version 9.0)\n");
	double startTime = Time();

	OctNode< TreeNodeData >::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);
	Octree< Real > tree;
	OctreeProfiler< Real > profiler(tree);
	tree.threads = pp.ThreadsVal;
	if (pp.MaxSolveDepthVal < 0) pp.MaxSolveDepthVal = pp.MaxDepthVal;

	//std::cerr << "Using " << pp.ThreadsVal << " threads\n";
	//	int kernelDepth = KernelDepth.set ? KernelDepth.value : Depth.value-2;
	if (pp.KernelDepthVal < 0) pp.KernelDepthVal = pp.MaxDepthVal - 2;
	if (pp.KernelDepthVal > pp.MaxDepthVal) {
		printf("kernelDepth cannot be greateer Depth.value\n");
		return false;
	}

	int pointCount;

	Real pointWeightSum;
	std::vector< typename Octree< Real >::PointSample >* samples = new std::vector< typename Octree< Real >::PointSample >();
	std::vector< ProjectiveData< Point3D< Real >, Real > >* sampleData = NULL;
	DensityEstimator* density = NULL;
	SparseNodeData< Point3D< Real >, NORMAL_DEGREE >* normalInfo = NULL;
	Real targetValue = (Real)0.5;

	// Read in the samples (and color data)
	{
		profiler.start();
		//		PointStream* pointStream;

		//		char* ext = GetFileExtension( In.value );
		//		if( Color.set && Color.value>0 )
		//		{
		//			sampleData = new std::vector< ProjectiveData< Point3D< Real > , Real > >();
		//			if     ( !strcasecmp( ext , "bnpts" ) ) pointStream = new BinaryOrientedPointStreamWithData< Real , Point3D< Real > , float , Point3D< unsigned char > >( In.value );
		//			else if( !strcasecmp( ext , "ply"   ) ) pointStream = new    PLYOrientedPointStreamWithData< Real , Point3D< Real > >( In.value , ColorInfo< Real >::PlyProperties , 6 , ColorInfo< Real >::ValidPlyProperties );
		//			else                                    pointStream = new  ASCIIOrientedPointStreamWithData< Real , Point3D< Real > >( In.value , ColorInfo< Real >::ReadASCII );
		//		}
		//		else
		//		{
		//			if     ( !strcasecmp( ext , "bnpts" ) ) pointStream = new BinaryOrientedPointStream< Real , float >( In.value );
		//			else if( !strcasecmp( ext , "ply"   ) ) pointStream = new    PLYOrientedPointStream< Real >( In.value );
		//			else                                    pointStream = new  ASCIIOrientedPointStream< Real >( In.value );
		//		}
		//		delete[] ext;
		sampleData = new std::vector< ProjectiveData< Point3D< Real >, Real > >();
		XPointStreamWithData _pointStream(xForm, (PointStreamWithData&)*pointStream);
		pointCount = tree.template init< Point3D< Real > >(_pointStream, pp.MaxDepthVal, pp.ConfidenceFlag, *samples, sampleData);

#pragma omp parallel for num_threads( pp.ThreadsVal )
		for (int i = 0; i < (int)samples->size(); i++)
			(*samples)[i].sample.data.n *= (Real)-1;

		//DumpOutput("Input Points / Samples: %d / %d\n", pointCount, samples->size());
		profiler.dumpOutput2(comments, "# Read input into tree:");
	}

	DenseNodeData< Real, Degree > solution;

	{
		DenseNodeData< Real, Degree > constraints;
		InterpolationInfo* iInfo = NULL;
		int solveDepth = pp.MaxSolveDepthVal;

		tree.resetNodeIndices();

		// Get the kernel density estimator [If discarding, compute anew. Otherwise, compute once.]
		{
			profiler.start();
			density = tree.template setDensityEstimator< WEIGHT_DEGREE >(*samples, pp.KernelDepthVal, pp.SamplesPerNodeVal);
			profiler.dumpOutput2(comments, "#   Got kernel density:");
		}

		// Transform the Hermite samples into a vector field [If discarding, compute anew. Otherwise, compute once.]
		{
			profiler.start();
			normalInfo = new SparseNodeData< Point3D< Real >, NORMAL_DEGREE >();
			*normalInfo = tree.template setNormalField< NORMAL_DEGREE >(*samples, *density, pointWeightSum, BType == BOUNDARY_NEUMANN);
			profiler.dumpOutput2(comments, "#     Got normal field:");
		}

		if (!pp.DensityFlag) {
			delete density;
			density = NULL;
		}

		// Trim the tree and prepare for multigrid
		{
			profiler.start();
			std::vector< int > indexMap;

			constexpr int MAX_DEGREE = NORMAL_DEGREE > Degree ? NORMAL_DEGREE : Degree;
			tree.template inalizeForBroodedMultigrid< MAX_DEGREE, Degree, BType >(pp.FullDepthVal, typename Octree< Real >::template HasNormalDataFunctor< NORMAL_DEGREE >(*normalInfo), &indexMap);

			if (normalInfo) normalInfo->remapIndices(indexMap);
			if (density) density->remapIndices(indexMap);
			profiler.dumpOutput2(comments, "#       Finalized tree:");
		}

		// Add the FEM constraints
		{
			profiler.start();
			constraints = tree.template initDenseNodeData< Degree >();
			tree.template addFEMConstraints< Degree, BType, NORMAL_DEGREE, BType >(FEMVFConstraintFunctor< NORMAL_DEGREE, BType, Degree, BType >(1., 0.), *normalInfo, constraints, solveDepth);
			profiler.dumpOutput2(comments, "#  Set FEM constraints:");
		}

		// Free up the normal info [If we don't need it for subseequent iterations.]
		delete normalInfo;
		normalInfo = NULL;

		// Add the interpolation constraints
		if (pp.PointWeightVal > 0) {
			profiler.start();
			iInfo = new InterpolationInfo(tree, *samples, targetValue, pp.AdaptiveExponentVal, (Real)pp.PointWeightVal * pointWeightSum, (Real)0);
			tree.template addInterpolationConstraints< Degree, BType >(*iInfo, constraints, solveDepth);
			profiler.dumpOutput2(comments, "#Set point constraints:");
		}

		//DumpOutput("Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n", (int)tree.leaves(), (int)tree.nodes(), (int)tree.ghostNodes());
		//DumpOutput( "Memory Usage: %.3f MB\n" , float( MemoryInfo::Usage())/(1<<20) );

		// Solve the linear system
		{
			profiler.start();
			typename Octree< Real >::SolverInfo solverInfo;
			solverInfo.cgDepth = pp.CGDepthVal, solverInfo.iters = pp.ItersVal, solverInfo.cgAccuracy = pp.CSSolverAccuracyVal, solverInfo.verbose = pp.VerboseFlag, solverInfo.showResidual = pp.ShowResidualFlag, solverInfo.lowResIterMultiplier = std::max< double >(1., pp.LowResIterMultiplierVal);
			solution = tree.template solveSystem< Degree, BType >(FEMSystemFunctor< Degree, BType >(0, 1., 0), iInfo, constraints, solveDepth, solverInfo);
			profiler.dumpOutput2(comments, "# Linear system solved:");
			if (iInfo) delete iInfo, iInfo = NULL;
		}
	}

	CoredFileMeshData< Vertex > mesh;

	{
		profiler.start();
		double valueSum = 0, weightSum = 0;
		typename Octree< Real >::template MultiThreadedEvaluator< Degree, BType > evaluator(&tree, solution, pp.ThreadsVal);
#pragma omp parallel for num_threads( pp.ThreadsVal ) reduction( + : valueSum , weightSum )
		for (int j = 0; j < (int)samples->size(); j++) {
			ProjectiveData< OrientedPoint3D< Real >, Real >& sample = (*samples)[j].sample;
			Real w = sample.weight;
			if (w > 0) weightSum += w, valueSum += evaluator.value(sample.data.p / sample.weight, omp_get_thread_num(), (*samples)[j].node) * w;
		}
		Real isoValue = (Real)(valueSum / weightSum);
		//		if( samples ) delete samples , samples = NULL;
		profiler.dumpOutput("Got average:");
		//DumpOutput("Iso-Value: %e\n", isoValue);

		profiler.start();
		SparseNodeData< ProjectiveData< Point3D< Real >, Real >, DATA_DEGREE >* colorData = NULL;
		if (sampleData) {
			colorData = new SparseNodeData< ProjectiveData< Point3D< Real >, Real >, DATA_DEGREE >();
			*colorData = tree.template setDataField< DATA_DEGREE, false >(*samples, *sampleData, (DensityEstimator*)NULL);
			delete sampleData, sampleData = NULL;
			for (const OctNode< TreeNodeData >* n = tree.tree().nextNode(); n; n = tree.tree().nextNode(n)) {
				ProjectiveData< Point3D< Real >, Real >* clr = (*colorData)(n);
				if (clr)
					(*clr) *= (Real)pow(pp.ColorVal, tree.depth(n));
			}
		}
		tree.template getMCIsoSurface< Degree, BType, WEIGHT_DEGREE, DATA_DEGREE >(density, colorData, solution, isoValue, mesh, !pp.LinearFitFlag, !pp.NonManifoldFlag, false /*PolygonMesh.set*/);
		//DumpOutput("Vertices / Polygons: %d / %d\n", mesh.outOfCorePointCount() + mesh.inCorePoints.size(), mesh.polygonCount());
		profiler.dumpOutput2(comments, "#        Got triangles:");
	}

	//        FreePointer( solution );

	if (cb) cb(90, "Creating Mesh");
	mesh.resetIterator();
	//int vm = mesh.outOfCorePointCount()+mesh.inCorePoints.size();
	for (auto pt = mesh.inCorePoints.begin(); pt != mesh.inCorePoints.end(); ++pt) {
		Point3D<Real> pp = iXForm * pt->point;
		vcg::tri::Allocator<VCG_Mesh>::AddVertex(pm, Point3m(pp[0], pp[1], pp[2]));
		pm.vert.back().Q() = pt->value;
		pm.vert.back().C()[0] = pt->color[0];
		pm.vert.back().C()[1] = pt->color[1];
		pm.vert.back().C()[2] = pt->color[2];
	}
	for (int ii = 0; ii < mesh.outOfCorePointCount(); ii++) {
		Vertex pt;
		mesh.nextOutOfCorePoint(pt);
		Point3D<Real> pp = iXForm * pt.point;
		vcg::tri::Allocator<VCG_Mesh>::AddVertex(pm, Point3m(pp[0], pp[1], pp[2]));
		pm.vert.back().Q() = pt.value;
		pm.vert.back().C()[0] = pt.color[0];
		pm.vert.back().C()[1] = pt.color[1];
		pm.vert.back().C()[2] = pt.color[2];
	}

	std::vector< CoredVertexIndex > polygon;
	while (mesh.nextPolygon(polygon)) {
		assert(polygon.size() == 3);
		int indV[3];
		for (int i = 0; i<int(polygon.size()); i++) {
			if (polygon[i].inCore)
				indV[i] = polygon[i].idx;
			else
				indV[i] = polygon[i].idx + int(mesh.inCorePoints.size());
		}
		vcg::tri::Allocator<VCG_Mesh>::AddFace(pm, &pm.vert[indV[0]], &pm.vert[indV[1]], &pm.vert[indV[2]]);
	}
	if (cb) cb(100, "Done");

	//if( colorData ) delete colorData , colorData = NULL;


	if (density) delete density, density = NULL;
	//DumpOutput2(comments, "#          Total Solve: %9.1f (s), %9.1f (MB)\n", Time() - startTime, tree.maxMemoryUsage());

	return 1;
}

template <class MeshType>
void PoissonClean(MeshType& m, bool scaleNormal, bool cleanFlag)
{
	vcg::tri::UpdateBounding<MeshType>::Box(m);
	vcg::tri::UpdateNormal<MeshType>::NormalizePerVertex(m);

	if (cleanFlag) {
		for (auto vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
			if (vcg::SquaredNorm(vi->N()) < std::numeric_limits<float>::min() * 10.0)
				vcg::tri::Allocator<MeshType>::DeleteVertex(m, *vi);
		}

		for (auto fi = m.face.begin(); fi != m.face.end(); ++fi)
			if (fi->V(0)->IsD() || fi->V(1)->IsD() || fi->V(2)->IsD())
				vcg::tri::Allocator<MeshType>::DeleteFace(m, *fi);
	}

	vcg::tri::Allocator<MeshType>::CompactEveryVector(m);
	if (scaleNormal) {
		for (auto vi = m.vert.begin(); vi != m.vert.end(); ++vi)
			vi->N() *= vi->Q();
	}
}

bool HasGoodNormal(VCG_Mesh& m)
{
	for (auto vi = m.vert.begin(); vi != m.vert.end(); ++vi)
		if (vcg::SquaredNorm(vi->N()) < std::numeric_limits<float>::min() * 10.0)
			return false;

	return true;
}
