![CoverBunnyEvol](https://github.com/MCInversion/ImplicitSurfaceWrap/blob/main/images/BunnyEvolCoverPic.png)

### Article Abstract:

Fairing methods, frequently used for smoothing noisy features of surfaces, evolve a surface towards the simplest shape. The inverse process - shaping a simple surface into a more complex object - requires a field defined in the ambient space driving the surface towards a target shape. Sharp protruding features combined with deep chasms in the target may give rise to severe reduction of triangle quality as the surface stretches to fit into concave regions. Our key contribution is a combination of adaptive remeshing and curvature-based feature detection to mitigate these issues while also maintaining the stability of a semi-implicit formulation of the method. We analyze the results of this approach using triangle quality metrics.

![ISWArchitecture](https://github.com/MCInversion/ImplicitSurfaceWrap/blob/main/images/ShrinkWrapMainUML.png)

# Introduction

Welcome to **ImplicitSurfaceWrap**! This application is a platform for the research of [Lagrangian surface evolution](chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/http://www.math.sk/mikula/mrss_SISC.pdf) on triangular meshes. Besides the core functionality (`SurfaceEvolver`) the codebase provides a distance field computation functionality in `DistanceFieldGenerator` class.

This repo is a fork of the PMP Library: https://www.pmp-library.org.

## Get Started

Clone the repository:

```sh
git clone --recursive https://github.com/MCInversion/ImplicitSurfaceWrap.git
```

Configure and build:

```sh
cd ImplicitSurfaceWrap && mkdir build && cd build && cmake .. && make
```

# Features

At the moment, there is no argument parsing, so all input is read from `./data` folder, and the application outputs into `./output`. The `main()` function is defined in `ShrinkWrapMain.cpp` with the vector of input *.obj mesh names (without the extension):

```
const std::vector<std::string> meshNames{
    "armadillo", // Stanford
    "BentChair", // custom
    "blub", // Keenan Crane
    "bunny", // Stanford
    "maxPlanck", // Kobbelt
    "nefertiti", // Cosmo Wenman
    "ogre", // Keenan Crane
    "spot" // Keenan Crane
};
```

Feel free to import a geometry of your choosing. It should be noted that the model is tested on the above meshes only. 

### Distance field computation

![SDFPic](https://github.com/MCInversion/ImplicitSurfaceWrap/blob/main/images/SDFsSixMeshes.jpg)

To test distance field, turn boolean flage `performSDFTests` on for your `meshNames`. Then feel free to edit parameters:

```
const SDF::DistanceFieldSettings sdfSettings{
		 cellSize,
		 1.0f,
		 0.2,
		 SDF::KDTreeSplitType::Center,
		 SDF::SignComputation::VoxelFloodFill,
		 SDF::BlurPostprocessingType::None,
		 SDF::PreprocessingType::Octree
};
```

according to your needs. `SDF::DistanceFieldGenerator::Generate` is called afterwards.

### Surface Evolver

![EvolverResults](https://github.com/MCInversion/ImplicitSurfaceWrap/blob/main/images/MeshAnalysisResultsSizing_LowRes.jpg)

The same holds for `SurfaceEvolver` functionality which is turned on by `performEvolverTests` flag. In the for loop for `meshNames` we also pre-compute the distance field using `SDF::DistanceFieldGenerator::Generate`. `SurfaceEvolver` provides a wide range of parameters defined in `SurfaceEvolver.h` (with descriptions). 

There is a variant of `SurfaceEvolver` class, called `SphereTest` which verifies the rate of convergence of mean curvature flow with respect to exact solution with radius `r(t) = sqrt(r0 * r0 - 4 * t)`. The test is applicable to remeshed surface, but be careful with the error and sizing settings for adaptive remeshing because it might crash.

### WIP

`BrainSurfaceEvolver` does not work yet. It just evolves a sphere and detects no brain. We also prepare `IsosurfaceEvolver` for future research.

## License

PMP is provided under a simple and flexible MIT-style [license](https://github.com/pmp-library/pmp-library/blob/master/LICENSE.txt) allowing for both open-source and commercial usage.
