#include "../include/interpolation.h"

#include <stdbool.h>

#include "../lib/qhull/src/libqhull_r/qhull_ra.h"

/**
 * Helper print function...
 */
static void print_triangles_facets(qhT *qh, double *pts) {
  facetT *facet = qh->facet_list;

  while (facet->id != 0) {
    if (!facet->upperdelaunay) {
      vertexT *v1 = facet->vertices->e[0].p;
      vertexT *v2 = facet->vertices->e[1].p;
      vertexT *v3 = facet->vertices->e[2].p;
      int h1 = qh_pointid(qh, v1->point);
      int h2 = qh_pointid(qh, v2->point);
      int h3 = qh_pointid(qh, v3->point);
      printf("%d -> %d -> %d \n", h1, h2, h3);
      printf(" >> (%f, %f) -> (%f, %f) -> (%f, %f) \n", pts[h1 * 2],
             pts[h1 * 2 + 1], pts[h2 * 2], pts[h2 * 2 + 1], pts[h3 * 2],
             pts[h3 * 2 + 1]);
    }
    facet = facet->next;
  }
}

/**
 * Creates an array of triangle indices from a Qhull triangulation.
 *
 * This function extracts triangle indices from a Qhull triangulation,
 * considering only non-upperdelaunay facets. Each triangle is represented
 * by three consecutive indices in the output array.
 *
 * @param qh       Pointer to Qhull's state data structure
 * @param indices  Double pointer to store the allocated array of indices.
 *                 The function will allocate memory for this array.
 *                 Caller is responsible for freeing this memory.
 * @param num      Pointer to store the number of non-upperdelaunay facets.
 *                 The actual size of indices array will be num * 3.
 *
 * @return         0 on success, INT_ALLOCATION_ERROR if memory allocation fails
 *
 * @note The indices array will be allocated with size num * 3 * sizeof(int)
 * @note Only processes non-upperdelaunay facets (lower hull triangles)
 * @note Caller must free the allocated memory in *indices when no longer needed
 */
int create_triangles_indices(qhT *qh, int **indices, int *num) {
  qh_findgood_all(qh, qh->facet_list);
  *num = qh->num_good;

  *indices = malloc(sizeof(double *) * (*num) * 3);
  if (*indices == NULL) {
    fprintf(stderr,
            "ERROR: allocation error while creating triangle indices.\n");
    return INT_ALLOCATION_ERROR;
  }

  int count = 0;
  facetT *facet = qh->facet_list;
  while (facet->id != 0) {
    if (!facet->upperdelaunay) {
      vertexT *v1 = facet->vertices->e[0].p;
      vertexT *v2 = facet->vertices->e[1].p;
      vertexT *v3 = facet->vertices->e[2].p;
      (*indices)[count] = qh_pointid(qh, v1->point);
      (*indices)[count + 1] = qh_pointid(qh, v2->point);
      (*indices)[count + 2] = qh_pointid(qh, v3->point);
      count += 3;
    }
    facet = facet->next;
  }
  return 0;
}

/**
 * Performs linear interpolation on 2D points using barycentric coordinates.
 *
 * This function interpolates values for a set of 2D points using triangular
 * facets from a Delaunay triangulation. For each input point, it:
 * 1. Finds which triangle (facet) contains the point
 * 2. Calculates barycentric coordinates for that point
 * 3. Interpolates the value using the triangle's vertex values
 *
 * @param qh          Pointer to Qhull's state data structure
 * @param ipoints     Array of 2D interpolation points [x1,y1,x2,y2,...,xn,yn]
 * @param ipval       Output array to store interpolated values for each point
 * @param inum_pts    Number of points to interpolate
 * @param pval        Array of known values at the triangulation vertices
 * @param fill_value  Value to use for points outside the convex hull
 *
 * @note Points outside the convex hull will:
 *       - Be assigned the fill_value
 *       - Generate a warning message to stdout
 * @warning Input arrays must be properly allocated:
 *          - ipoints must be size inum_pts * 2
 *          - ipval must be size inum_pts
 *          - pval must match tohe size of original points in triangulation
 * Example usage:
 * @code
 *     double points[] = {1.0, 2.0, 3.0, 4.0};  // Two points to interpolate
 *     double values[2];                         // Output interpolated values
 *     linear_interp2d_facet(qh, points, values, 2, original_values, -999.0);
 * @endcode
 */
static void linear_interp2d_facet(qhT *qh, double *ipoints, double *ipval,
                                  int inum_pts, double *pval,
                                  double fill_value) {
  // facet reference
  facetT *facet;

  // vertices references
  vertexT *vertex, **vertexp;

  // coordinates
  coordT single_point[2];

  // coords always in 2D
  int ncoords = inum_pts * 2;

  // loop through all coordinates!
  for (int i = 0; i < ncoords; i += 2) {
    single_point[0] = ipoints[i];
    single_point[1] = ipoints[i + 1];
    const int idx = i >> 1;

    // facet list reference
    facet = qh->facet_list;

    // flag to mark points inside convex hull
    bool inside_hull = false;

    // default value
    ipval[idx] = -9999.0;

    while (facet->id != 0) {
      if (!facet->upperdelaunay) {
        double tript_x[3];
        double tript_y[3];

        int count = 0;

        FOREACHvertex_(facet->vertices) {
          tript_x[count] = vertex->point[0];
          tript_y[count] = vertex->point[1];
          ++count;
        } // end collecting triangle points

        // More efficient with this procedure:
        // https://math.stackexchange.com/questions/51326/determining-if-an-arbitrary-point-lies-inside-a-triangle-defined-by-three-points

        double x = single_point[0];
        double y = single_point[1];

        vertex = facet->vertices->e[0].p;
        double x1 = tript_x[0];
        double y1 = tript_y[0];
        double f1 = pval[qh_pointid(qh, vertex->point)];

        vertex = facet->vertices->e[1].p;
        double x2 = tript_x[1];
        double y2 = tript_y[1];
        double f2 = pval[qh_pointid(qh, vertex->point)];

        vertex = facet->vertices->e[2].p;
        double x3 = tript_x[2];
        double y3 = tript_y[2];
        double f3 = pval[qh_pointid(qh, vertex->point)];

        /*
           Converting Cartesian to barycentric coordinates
           https://en.wikipedia.org/wiki/Barycentric_coordinate_system

           denomr = ((y2-y3)*(x1-x3) + (x3-x2)*(y1-y3))
           L1 = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3)) / denomr
           L2 = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3)) / denomr
           L3 = 1.0 - L1 - L2
         */

        double denom =
            y2 * x1 - y2 * x3 - y3 * x1 + x3 * y1 - x2 * y1 + x2 * y3;

        double L1 = y2 * x - y2 * x3 - y3 * x + x3 * y - x2 * y + x2 * y3;
        L1 = L1 / denom;

        double L2 = y3 * x - y1 * x + y1 * x3 + x1 * y - x1 * y3 - x3 * y;
        L2 = L2 / denom;

        double L3 = 1 - L1 - L2;

        if (L1 > 0.0 && L2 > 0.0 && L3 > 0.0) {
          inside_hull = true;

          // apply barycentric coordinates
          double f = L1 * f1 + L2 * f2 + L3 * f3;
          ipval[idx] = f;
          break;
        }

      } // end !upperdelauney

      facet = facet->next;
    } // end looking through facets

    // if false, current point is not inside any triangle
    // and not inside convex hull.
    if (!inside_hull) {
      ipval[idx] = fill_value;
      fprintf(stdout, " -- WARNING:: point (%f, %f) outside convex hull \n",
              single_point[0], single_point[1]);
    }

  } // end for
}

/**
 * Interpolates values for unstructured 2D points using Delaunay triangulation.
 *
 * This function performs interpolation by:
 * 1. Creating a Delaunay triangulation of the input points
 * 2. Using linear interpolation within the triangles
 * 3. Applying barycentric coordinates for the interpolation
 *
 * @param points    Array of known point coordinates [x1,y1,x2,y2,...,xn,yn]
 * @param values    Array of known values at each point location [v1,v2,...,vn]
 * @param num_pts   Number of known points (must be >= 4)
 * @param ipoints   Array of points to interpolate [x1,y1,x2,y2,...,xn,yn]
 * @param ivalues   Output array for interpolated values
 * @param inum_pts  Number of points to interpolate (must be >= 1)
 * @param fill_value Value to use for points outside the convex hull
 *
 * @return          0 on success, error code on failure:
 *                  INTERP_MIN_ERROR if insufficient points
 *                  QHULL_GENERAL_ERROR if Qhull fails
 *
 * @warning This function does not check for:
 *          - Collinear points
 *          - Degenerate triangles
 *          These conditions will produce undefined behavior
 *
 * @note Array sizes must be:
 *       - points:  num_pts * 2 elements
 *       - values:  num_pts elements
 *       - ipoints: inum_pts * 2 elements
 *       - ivalues: inum_pts elements
 *
 * Example usage:
 * @code
 *     double points[] = {0.0, 1.0, 0.1, 1.1};  // 4 points
 *     double values[] = {0, 1, 1, 2};          // values at points
 *     double ipoints[] = {0.5, 0.5};           // point to interpolate
 *     double ivalues[1];                       // interpolated result
 *     int result = griddata(points, values, 4, ipoints, ivalues, 1, -999.0);
 * @endcode
 *
 * @note Uses Qhull library for Delaunay triangulation with options:
 *       "QVn QJ d" for:
 *       - Quiet output
 *       - Verify results
 *       - Jiggle input to avoid precision problems
 *       - Delaunay triangulation
 */
int griddata(double *points, double *values, int num_pts, double *ipoints,
             double *ivalues, int inum_pts, double fill_value) {
  const int DIMS2D = 2;

  if (num_pts < 4) {
    fprintf(stdout, " -- ERROR: Qhull needs a minimum of four points.\n");
    return INTERP_MIN_ERROR;
  } else if (inum_pts < 1) {
    fprintf(stdout, " -- ERROR: implementation needs a minimum ");
    fprintf(stdout, "of one interpolated point location.\n");
    return INTERP_MIN_ERROR;
  }

  //   d - Delaunay triangulation by lifting points to a paraboloid
  //   J  - slightly joggles the point location
  //   char noptions[CMDOPTS] = "QVn d";
  char *noptions = "QVn QJ d";

  // new instance for qhull
  qhT qh_qh;
  qhT *qh = &qh_qh;

  // True if qh_freeqhull should 'free(array)'
  boolT ismalloc = False;

  // initiate qhull -
  // legacy call to arm variables that are generally used for command line
  qh_init_A(qh, stdin, stdout, stderr, 0, NULL);
  int exitcode = setjmp(qh->errexit);
  if (exitcode < 0)
    return QHULL_GENERAL_ERROR;

  // not sure what this means yet
  qh->NOerrexit = False;

  // commandline option for:
  qh_initflags(qh, noptions);

  // true for delaunay and will allow read-in points
  qh->PROJECTdelaunay = True;

  coordT *allpts = points;

  // second initialization with points
  qh_init_B(qh, allpts, num_pts, DIMS2D, ismalloc);

  // convex hull
  qh_qhull(qh);

  // triangluate hull points
  qh_triangulate(qh);

  // interpolate with known point values
  linear_interp2d_facet(qh, ipoints, ivalues, inum_pts, values, fill_value);

  if (exitcode < 0)
    return -1;
  else
    return 0;
}

/**
 * Produces a triangle mesh from an array of 2-D data points.
 *
 * This function does not check for conditions where points form lines
 * or skinny triangle elements - which will produce undefined behaviors.
 *
 * @param points x and y coordinates of the points with known values
 * @param value Given values for each point location. Must be half the
 *                   size of the length of points.
 * @param num_pts The number of known values, or use following:
 *                int num_pts  = sizeof(values)/sizeof(values[0]); should be
 *                at least three points
 * @param ipoints Array of x and y coordinate of the point(s) whose
 *                unknown value to interpolate.
 * @param ivalues The value array for qhull interpolation results
 * @param inum_pts The number of unknown values. or use following:
 *                 int inum_pts = sizeof(ipoints)/(2 * sizeof(ipoints[0]));
 * @param fill_Value Default fill value to use if point location can't
 *                   interpolated.
 */
int griddata_triangles(double *points, int num_pts, int **indices, int *num) {
  const int DIMS2D = 2;

  if (num_pts < 4) {
    fprintf(stdout, " -- ERROR: Qhull needs a minimum of four points.\n");
    return INTERP_MIN_ERROR;
  }

  //   d - Delaunay triangulation by lifting points to a paraboloid
  //   J  - slightly joggles the point location
  //   char noptions[CMDOPTS] = "QVn d";
  char *noptions = "QVn QJ d";

  // new instance for qhull
  qhT qh_qh;
  qhT *qh = &qh_qh;

  // True if qh_freeqhull should 'free(array)'
  boolT ismalloc = False;

  // initiate qhull -
  // legacy call to arm variables that are generally used for command line
  qh_init_A(qh, stdin, stdout, stderr, 0, NULL);
  int exitcode = setjmp(qh->errexit);
  if (exitcode < 0)
    return QHULL_GENERAL_ERROR;

  // not sure what this means yet
  qh->NOerrexit = False;

  // commandline option for:
  qh_initflags(qh, noptions);

  // true for delaunay and will allow read-in points
  qh->PROJECTdelaunay = True;

  coordT *allpts = points;

  // second initialization with points
  qh_init_B(qh, allpts, num_pts, DIMS2D, ismalloc);

  // convex hull
  qh_qhull(qh);

  // triangluate hull points
  qh_triangulate(qh);

  exitcode = create_triangles_indices(qh, indices, num);

  if (exitcode < 0)
    return exitcode;
  else
    return 0;
}
