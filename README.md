![CoverBunnyEvol](https://github.com/MCInversion/ImplicitSurfaceWrap/blob/main/images/BunnyEvolCoverPic.png)

### Article Abstract:

Fairing methods, frequently used for smoothing noisy features of surfaces, evolve a surface towards the simplest shape. The inverse process - shaping a simple surface into a more complex object - requires a field defined in the ambient space driving the surface towards a target shape. Sharp protruding features combined with deep chasms in the target may give rise to severe reduction of triangle quality as the surface stretches to fit into concave regions. Our key contribution is a combination of adaptive remeshing and curvature-based feature detection to mitigate these issues while also maintaining the stability of a semi-implicit formulation of the method. We analyze the results of this approach using triangle quality metrics.

![ISWArchitecture](https://github.com/MCInversion/ImplicitSurfaceWrap/blob/main/images/ShrinkWrapMainUML.png)

# Introduction

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

## License

PMP is provided under a simple and flexible MIT-style [license](https://github.com/pmp-library/pmp-library/blob/master/LICENSE.txt) allowing for both open-source and commercial usage.
