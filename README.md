![ThreeMeshReps](https://github.com/MCInversion/MCInversionPhDCollection/blob/main/images/ThreeMeshReps.png)

# Introduction

Welcome to **MCInversionPhDCollection**! This application is a platform for my PhD research on meshes and other 3D data representations.

See more in Section ["**Functionality**"](https://github.com/MCInversion/MCInversionPhDCollection/blob/main/README.md#functionality) or at: https://mshgrid.com/2023/02/28/my-research/

This work supports a range of applications from engineering, where models are used for Finite Element Method (FEM) simulations, to purely visual applications in virtual environments. By leveraging advancements in data storage and processing, our methods push the boundaries of what's possible in digital representations of physical objects.

This repo is a fork of the PMP Library: https://www.pmp-library.org.

Additional dependencies include:

- Eigen
- Glew
- Glfw
- GoogleTest
- Imgui
- [Nanoflann](https://github.com/jlblancoc/nanoflann)
- Rply
- [VCGLib](https://github.com/cnr-isti-vclab/vcglib)
- [Quickhull](https://github.com/akuukka/quickhull)

## Get Started

Clone the repository:

```sh
git clone --recursive https://github.com/MCInversion/MCInversionPhDCollection.git
```

Configure and build:

```sh
cd MCInversionPhDCollection && mkdir build && cd build && cmake .. && make
```

## Functionality 

This repository explores advanced techniques in 3D data storage and display, focusing primarily on polygonal meshes, especially triangular meshes for tessellating surfaces of solid objects. Meshes, although finite in information, provide a vast array of configurations creating the need for a variety of processing algorithms. Our research delves into extracting volumetric data from meshes and reconstructing optimized meshes, enhancing them with the latest techniques in decimation, remeshing, and specific combinatorial adjustments.
Key features of our project include:

![AllSDFs](https://github.com/MCInversion/MCInversionPhDCollection/blob/main/images/AllSDFs.png)
**Signed Distance Fields Computation**: From triangle soup meshes and point cloud data.

![TitlePicBunnyPtCloud](https://github.com/MCInversion/MCInversionPhDCollection/blob/main/images/TitlePicBunnyPtCloud.png)
**Lagrangian Shrink-Wrapping**: An innovative approach to mesh optimization that adapts to the underlying geometry dynamically.

![IcoSphereWith2HolesSubdiv](https://github.com/MCInversion/MCInversionPhDCollection/blob/main/images/IcoSphereWith2HolesSubdiv.jpg)
**Counting Formulas for Subdivision Surfaces**: These formulas help predict the complexity and behavior of surfaces as they are subdivided, aiding in efficient mesh refinement.

![LargeMeshBufferingNoTransparency](https://github.com/MCInversion/MCInversionPhDCollection/blob/main/images/LargeMeshBufferingNoTransparency.png)
**Massive 3D File Simplification**: Utilizing a multithreaded raw data sampling strategy, this feature significantly reduces the complexity of 3D models without compromising essential details, facilitating faster processing and lower memory requirements.


Here's the list of tests carried out in Main.cpp:

- SDF Tests
- Sphere Test
- Evolver Tests
- Old Armadillo LSW Test
- Isosurface Evolver Tests
- Sheet Evolver Test
- Brain Evolver Tests
- Subdivision Tests 1
- Subdivision Tests 2
- Subdivision Tests 3
- Subdivision Test 4
- Subdiv Tests Boundary
- Subdiv Tests Multi Torus
- Subdiv Preallocation Tests
- New Icosphere Tests
- Icosphere Performance Tests
- Catmull Clark Counting
- Remeshing Tests
- Mobius Strip Voxelization
- Metaball Test
- Imported Obj Metrics Eval
- MMap Import Test
- MMap OBJ Chunk Marking Test
- Simple Bunny OBJ Sampling Demo
- PDaniel Pt Cloud PLY Export
- Pt Cloud To DF
- PDaniel Pt Cloud Comparison Test
- Repulsive Surf Result Evaluation
- Histogram Result Evaluation
- Old Result Jacobian Metric Eval
- Hausdorff Distance Measurements Per Time Step
- Direct Higher Genus Pt Cloud Sampling
- Higher Genus Pt Cloud LSW
- TriTri Intersection Tests
- Mesh Self Intersection Tests
- Hurtado Meshes Isosurface Evolver Tests
- Hurtado Trex Icosphere LSW
- Import VTI Debug Tests
- Convex Hull Tests
- Convex Hull Remeshing Tests
- Convex Hull Evolver Tests
- IcoSphere Evolver Tests
- BPATest
- Incremental Mesh Builder Tests
- 2GB Apollon Mesh Builder Test
- Nanoflann Distance Tests
- Apollon LSW Saliency Eval
- Incremental Mesh Builder Hausdorff Eval
- Apollon Artec Eva LSW Hausdorff Eval
