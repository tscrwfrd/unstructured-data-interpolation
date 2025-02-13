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
- Jonathan Shewchuk's robust geometric [predicates](https://www.cs.cmu.edu/afs/cs/project/quake/public/code/). The link I used to access code: https://www.cs.cmu.edu/afs/cs/project/quake/public/code/


## Future work
- Improve documentation and examples
- Implement different grid interpolation procedures
- Add support for 3D convex hulls
- Add support for [Fortran](https://fortran-lang.org/)

