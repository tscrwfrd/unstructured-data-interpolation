# unstructured-data-interpolation (UDI)

UDI provides efficient interpolation of unstructured 2D and 3D data points using Delaunay triangulation. This project aims to replicate the functionality of SciPy's `griddata` interpolation in a standalone C library.

- Linear interpolation using barycentric coordinates
- Delaunay triangulation via [Qhull](http://www.qhull.org/)
- Unified API supporting both 2D (triangles) and 3D (tetrahedra) interpolation

⚠️ **Work in Progress**
UDI is currently functional and can perform basic 2D and 3D interpolation tasks, but it's actively under development. While the core functionality works, you may encounter:
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
./bin/examples.exe
```

## Installation prerequisites

[QHull](http://www.qhull.org/) development headers and library are required. Install via your system package manager:

```bash
# Fedora/RHEL
sudo dnf install qhull-devel

# Debian/Ubuntu
sudo apt install libqhull-dev
```

Tests require the [cmocka](https://cmocka.org/) library:

```bash
# Fedora/RHEL
sudo dnf install libcmocka-devel

# Debian/Ubuntu
sudo apt install libcmocka-dev
```

## API

Both public functions take a `dims` parameter (2 or 3) to select 2D or 3D mode. All point arrays use flat interleaved layout: `[x₀, y₀, x₁, y₁, ...]` for 2D or `[x₀, y₀, z₀, x₁, y₁, z₁, ...]` for 3D.

```c
/* Interpolate values at query points (ipoints → ivalues).
 * Points outside the convex hull receive fill_value.
 * Requires num_pts >= dims+2, inum_pts >= 1. */
int griddata(double* points, double* values, int num_pts,
             double* ipoints, double* ivalues, int inum_pts,
             double fill_value, int dims);

/* Return the Delaunay simplex mesh (triangles for dims=2, tetrahedra for dims=3).
 * Output simplex_list is a flat array with stride dims+1. Caller must free it. */
int griddata_triangles(double* points, int num_pts,
                       int** simplex_list, int* num, int dims);
```

Error codes: `QHULL_GENERAL_ERROR (-30)`, `INTERP_MIN_ERROR (-20)`, `INT_ALLOCATION_ERROR (-10)`.

### Example

```c
// 2D interpolation
double points[]  = {0.0,0.0, 2.0,0.0, 0.0,2.0, 2.0,2.0};
double values[]  = {1.0, 2.0, 3.0, 4.0};
double ipoints[] = {1.0, 1.0};
double ivalues[1];
griddata(points, values, 4, ipoints, ivalues, 1, -9999.0, 2);

// 3D interpolation
double pts3[]  = {0,0,0, 2,0,0, 0,2,0, 0,0,2, 2,2,2};
double vals3[] = {1.0, 1.0, 1.0, 1.0, 1.0};
double ipt3[]  = {0.5, 0.5, 0.5};
double ival3[1];
griddata(pts3, vals3, 5, ipt3, ival3, 1, -9999.0, 3);
```

## Delaunator

The Delaunator component is based on two main sources:

- A C translation of the JavaScript [Delaunator](https://github.com/mapbox/delaunator/tree/main) library
- Jonathan Shewchuk's robust geometric [predicates](https://www.cs.cmu.edu/afs/cs/project/quake/public/code/)

Delaunator is a 2D-only algorithm and is not used for interpolation — it exists in the codebase for comparison with Qhull's triangulations. Its triangulations currently differ from Qhull's results and require further study.

## Barycentric Coordinates

UDI uses [barycentric coordinates](https://en.wikipedia.org/wiki/Barycentric_coordinate_system) for interpolation within simplices.

### 2D (triangles)

For a point (x,y) in a triangle with vertices (x₁,y₁), (x₂,y₂), (x₃,y₃):

```
λ₁ = ((y₂-y₃)(x-x₃) + (x₃-x₂)(y-y₃)) / D
λ₂ = ((y₃-y₁)(x-x₃) + (x₁-x₃)(y-y₃)) / D
λ₃ = 1 - λ₁ - λ₂

D  = (y₂-y₃)(x₁-x₃) + (x₃-x₂)(y₁-y₃)
```

### 3D (tetrahedra)

For a point p in a tetrahedron with vertices v₀, v₁, v₂, v₃, substituting λ₄ = 1 - λ₁ - λ₂ - λ₃ yields the 3×3 linear system:

```
T · [λ₁, λ₂, λ₃]ᵀ = p - v₃
```

where T has columns (v₀-v₃), (v₁-v₃), (v₂-v₃). Solved via Cramer's rule; no external linear algebra library required.

### Properties (both dimensions)
- All λᵢ sum to 1
- Point is inside the simplex if all λᵢ ≥ 0
- Interpolated value: v = λ₁v₁ + λ₂v₂ + ... + λₙvₙ

## Future work
- Improve documentation and examples
- Implement different interpolation methods (nearest-neighbor, cubic)
- Add support for [Fortran](https://fortran-lang.org/)

