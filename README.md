# unstructured-data-interpolation (UDI)

UDI provides efficient interpolation of unstructured 2D data points using Delaunay triangulation. This project aims to replicate the functionality of SciPy's griddata interpolation in a standalone C library.

- Linear interpolation using barycentric coordinates
- Delaunay triangulation using Qhull (and soon [Delaunator](https://github.com/mapbox/delaunator/tree/main))
- Support for 2D point interpolation

⚠️ **Work in Progress**
UDI is currently functional and can perform basic 2D interpolation tasks, but it's actively under development. While the core functionality works, you may encounter:
- Limited error handling
- Ongoing API refinements
- Documentation updates

To run:
```bash
make
make test
```

An example of how to use UDI:
```bash
./bin/example.exe
```

## Installation prerequisites

[QHull](http://www.qhull.org/) is necessary to compile udi.

[qhull github]: https://github.com/qhull/qhull

### Qhull steps


- Clone [qhull](https://github.com/qhull/qhull) within the `lib` directory.
- `cd` into `qhull` directory and run `make` and `make test`.
- `export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH`
- Now you should be able to run `make` in the project root directory and then run `./udi.exe`.

## Delaunator

The Delaunator component is based on two main sources:

- A C translation of the JavaScript [Delaunator](https://github.com/mapbox/delaunator/tree/main) library
- Jonathan Shewchuk's robust geometric [predicates](https://www.cs.cmu.edu/afs/cs/project/quake/public/code/). The link I used to access code: https://www.cs.cmu.edu/afs/cs/project/quake/public/code/.


Initially, this project considered using Delaunator as an alternative meshing approach for interpolation tasks. However, after implementation, I decided to maintain only the triangulation functionality while deferring the full Delaunay-based interpolation to Qhull. While Delaunator's implementation remains in the codebase, I need to study its Delaunay properties more thoroughly to understand why its triangulations differ from Qhull's results.

## Barycentric Coordinates notes

UDI uses [barycentric coordinates](https://en.wikipedia.org/wiki/Barycentric_coordinate_system) for interpolation within triangles. The coordinates (λ₁, λ₂, λ₃) for a point (x,y) in a triangle with vertices (x₁,y₁), (x₂,y₂), (x₃,y₃) are calculated as:

λ₁ = ((y₂-y₃)(x-x₃) + (x₃-x₂)(y-y₃)) / D
λ₂ = ((y₃-y₁)(x-x₃) + (x₁-x₃)(y-y₃)) / D
λ₃ = 1 - λ₁ - λ₂

where:
D = (y₂-y₃)(x₁-x₃) + (x₃-x₂)(y₁-y₃)

#### Properties
- λ₁ + λ₂ + λ₃ = 1
- Point is inside triangle if all λᵢ > 0
- Point is on edge if one λᵢ = 0
- Point is outside triangle if any λᵢ < 0

#### Interpolation
For a value v at point (x,y):

v = λ₁v₁ + λ₂v₂ + λ₃v₃

where v₁, v₂, v₃ are the values at the triangle vertices.

This interpolation provides a continuous, linear variation of values across the triangle, making it ideal for unstructured data interpolation.


## Future work
- Improve documentation and examples
- Implement different grid interpolation procedures
- Add support for 3D convex hulls
- Add support for [Fortran](https://fortran-lang.org/)

