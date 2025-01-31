# unstructured-data-interpolation

This is a replica of Scipy's griddata implementation.

The scatttered-data-interpolation (SDI) aims to provide 
the same functionality of Scipy's interpolation functionality. Current 
devlopement efforts is only on linux. **THIS IS A WORK IN PROGRESS.**

## Installation prerequisites

The following is necessary to use SDI

- [GNU Make][gmake] 

- [GCC compiler][gcc]

- [QHull][qhull]


[gmake]: https://www.gnu.org/software/make/
[gcc]: https://gcc.gnu.org/
[qhull]: http://www.qhull.org/

For QHull, run basic generic makefile with default settings.  Once compiled, SDI needs to know how to include "qhull_ra.h"
and link .so files libqhullstatic_r, libqhullstatic, and qhull_r.

### QHULL steps

- clone [qhull](https://github.com/qhull/qhull) within the root directory.
- `cd` into `qhull` directory and run `make` and `make test`.
- `export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH`
- Now you should be able to run `make` in the root directory and then run `./sdi.exe`.

## Delaunator

Delaunator is mostly a direct translation from the [delaunator](https://github.com/mapbox/delaunator/tree/main) javascript library.  I had to heavily lift code from Jonathon Shewchuk's [predicate](https://www.cs.cmu.edu/afs/cs/project/quake/public/code/) routines from his robust [Predicates](https://www.cs.cmu.edu/~quake/robust.html) page.

The link I used to access code: https://www.cs.cmu.edu/afs/cs/project/quake/public/code/



## Future work
- Different grid interpolation procedures
- 3D convex hulls and interpolation procedures
- possible optimization applications
- Unit test cases
- Add [delaunator](https://github.com/mapbox/delaunator) translation

