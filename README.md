![ThreeMeshReps](https://github.com/MCInversion/MCInversionPhDCollection/blob/main/images/ThreeMeshReps.png)

# Introduction

Welcome to **MCInversionPhDCollection*! This application is a platform for my PhD research on meshes and other 3D data representations.

This repository explores advanced techniques in 3D data storage and display, focusing primarily on polygonal meshes, especially triangular meshes for tessellating surfaces of solid objects. Meshes, although finite in information, provide a vast array of configurations necessitating diverse algorithms for processing. Our research delves into extracting volumetric data from meshes and reconstructing optimized meshes, enhancing them with the latest advancements in 3D modeling technologies.

Key features of our project include:

![TitlePicBunnyPtCloud](https://github.com/MCInversion/MCInversionPhDCollection/blob/main/images/TitlePicBunnyPtCloud.png)
**Lagrangian Shrink-Wrapping**: An innovative approach to mesh optimization that adapts to the underlying geometry dynamically.

![IcoSphereWith2HolesSubdiv](https://github.com/MCInversion/MCInversionPhDCollection/blob/main/images/IcoSphereWith2HolesSubdiv.jpg)
**Counting Formulas for Subdivision Surfaces**: These formulas help predict the complexity and behavior of surfaces as they are subdivided, aiding in efficient mesh refinement.

![LargeMeshBufferingNoTransparency](https://github.com/MCInversion/MCInversionPhDCollection/blob/main/images/LargeMeshBufferingNoTransparency.png)
**Massive 3D File Simplification**: Utilizing a multithreaded raw data sampling strategy, this feature significantly reduces the complexity of 3D models without compromising essential details, facilitating faster processing and lower memory requirements.

This work supports a range of applications from engineering, where models are used for Finite Element Method (FEM) simulations, to purely visual applications in virtual environments. By leveraging advancements in data storage and processing, our methods push the boundaries of what's possible in digital representations of physical objects.

This repo is a fork of the PMP Library: https://www.pmp-library.org.

Additional dependencies include:

- Eigen
- Glew
- Glfw
- GoogleTest
- Imgui
- Nanoflann
- Rply
- [VCGLib](https://github.com/cnr-isti-vclab/vcglib)

## Get Started

Clone the repository:

```sh
git clone --recursive https://github.com/MCInversion/MCInversionPhDCollection.git
```

Configure and build:

```sh
cd ImplicitSurfaceWrap && mkdir build && cd build && cmake .. && make
```

# Features

- 

## License

PMP is provided under a simple and flexible MIT-style [license](https://github.com/pmp-library/pmp-library/blob/master/LICENSE.txt) allowing for both open-source and commercial usage.
