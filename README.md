# unstructured-data-interpolation

This is a replica of Scipy's griddata implementation.

The unstructured-data-interpolation (udi) aims to provide 
the same functionality of Scipy's interpolation functionality. Current 
devlopement efforts is only on linux and is tested with gcc version 14.2.1. **THIS IS A WORK IN PROGRESS.**

## Installation prerequisites

[QHull](http://www.qhull.org/) is necessary to compile udi.

[qhull github]: https://github.com/qhull/qhull

### Qhull steps

- Clone [qhull](https://github.com/qhull/qhull) within the `lib` directory.
- `cd` into `qhull` directory and run `make` and `make test`.
- `export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH`
- Now you should be able to run `make` in the project root directory and then run `./udi.exe`.

## Delaunator

Delaunator is mostly a direct translation from the [delaunator](https://github.com/mapbox/delaunator/tree/main) javascript library.  I had to heavily lift code from Jonathon Shewchuk's [predicate](https://www.cs.cmu.edu/afs/cs/project/quake/public/code/) routines from his robust [Predicates](https://www.cs.cmu.edu/~quake/robust.html) page.

The link I used to access code: https://www.cs.cmu.edu/afs/cs/project/quake/public/code/



## Future work
- Different grid interpolation procedures
- 3D convex hulls and interpolation procedures
- possible optimization applications

