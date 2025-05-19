#pragma once

#include "Poisson/Octree.h"
#include "Poisson/MultiGridOctreeData.h"
#include "Poisson/MarchingCubes.h"
#include "Poisson/SparseMatrix.h"
#include "Poisson/PointStream.h"
#include "Poisson/Vector.h"

using Real = float;  // or double

// 1) Wrap your points+normals in their Stream interface
struct VectorPointStream : OrientedPointStream<Real> {
    const std::vector<Point3<Real>>& P;
    const std::vector<Point3<Real>>& N;
    size_t cur = 0;
    VectorPointStream(const auto& pts, const auto& nms) :P(pts), N(nms) {
        assert(P.size() == N.size());
    }
    void reset() override { cur = 0; }
    bool nextPoint(OrientedPoint3D<Real>& o) override {
        if (cur >= P.size()) return false;
        o.point = P[cur];
        o.normal = N[cur];
        ++cur;
        return true;
    }
};

// 2) Bundle parameters into the same struct as PoissonParam<Real>
struct MyPoissonParams {
    int    depth = 8;
    int    fullDepth = 5;
    int    cgDepth = 0;
    int    threads = omp_get_max_threads();
    Real   scale = 1.1f;
    Real   samplesPerNode = 1.5f;
    Real   pointWeight = 4.0f;
    int    iters = 8;
    bool   confidence = false;
    bool   preClean = false;
};

// 3) Your user‐facing API:
struct Mesh {
    std::vector<Point3<Real>> vertices;
    std::vector<std::array<int, 3>> faces;
};

inline Mesh ComputePoisson(
    const std::vector<Point3<Real>>& pts,
    const std::vector<Point3<Real>>& nms,
    const MyPoissonParams& up)
{
    // a) Prepare stream & bounding box
    VectorPointStream stream(pts, nms);
    // You can compute bounds however you like; here we use the code's helper:
    vcg::Box3<Real> bb; bb.Set(pts[0]);
    for (auto& p : pts) bb.Add(p);

    // b) Fill PoissonParam<Real>
    PoissonParam<Real> pp;
    pp.MaxDepthVal = up.depth;
    pp.FullDepthVal = up.fullDepth;
    pp.CGDepthVal = up.cgDepth;
    pp.ThreadsVal = up.threads;
    pp.ScaleVal = up.scale;
    pp.SamplesPerNodeVal = up.samplesPerNode;
    pp.PointWeightVal = up.pointWeight;
    pp.ItersVal = up.iters;
    pp.ConfidenceFlag = up.confidence;
    pp.CleanFlag = up.preClean;
    pp.DensityFlag = true;

    // c) Call the core _Execute<...>, writing into an in-memory mesh
    CoredFileMeshData< typename Octree<Real>::Vertex > meshData;
    ::Reset<Real>();         // clear any static multigrid state
    stream.reset();
    _Execute<Real, 2, BOUNDARY_NEUMANN, sizeof(typename Octree<Real>::Vertex)>(
        &stream, bb, meshData, pp, /*callback*/ nullptr
        );

    // d) Extract vertices + triangles
    Mesh out;
    meshData.resetIterator();
    // in-core points
    for (auto& v : meshData.inCorePoints) {
        out.vertices.push_back(v.point);
    }
    // out-of-core points
    for (int i = 0; i < meshData.outOfCorePointCount(); ++i) {
        typename Octree<Real>::Vertex v;
        meshData.nextOutOfCorePoint(v);
        out.vertices.push_back(v.point);
    }
    // faces
    std::vector<CoredVertexIndex> poly;
    while (meshData.nextPolygon(poly)) {
        std::array<int, 3> f;
        for (int i = 0; i < 3; ++i) {
            f[i] = poly[i].inCore
                ? poly[i].idx
                : poly[i].idx + (int)meshData.inCorePoints.size();
        }
        out.faces.push_back(f);
    }

    return out;
}