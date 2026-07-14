#include "../include/interpolation.h"

#include <math.h>
#include <stdbool.h>

#include <libqhull_r/qhull_ra.h>

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
int create_simplex_indices(qhT *qh, int **indices, int *num, int dims) {
  qh_findgood_all(qh, qh->facet_list);
  *num = qh->num_good;

  int verts_per_simplex = dims + 1;
  *indices = malloc(sizeof(int) * (*num) * verts_per_simplex);
  if (*indices == NULL) {
    fprintf(stderr,
            "ERROR: allocation error while creating simplex indices.\n");
    return INT_ALLOCATION_ERROR;
  }

  int count = 0;
  for (facetT *facet = qh->facet_list; facet->id != 0; facet = facet->next) {
    if (!facet->upperdelaunay) {
      for (int v = 0; v < verts_per_simplex; v++) {
        vertexT *vtx = facet->vertices->e[v].p;
        (*indices)[count + v] = qh_pointid(qh, vtx->point);
      }
      count += verts_per_simplex;
    }
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
static void linear_interp2d(qhT *qh, const double *ipoints, double *ipval,
                            int inum_pts, const double *pval,
                            double fill_value) {
  // Process each point
  for (int i = 0; i < inum_pts; i++) {
    const double x = ipoints[i * 2];
    const double y = ipoints[i * 2 + 1];
    bool found = false;

    // Default to fill value
    ipval[i] = fill_value;

    // Check each facet
    for (facetT *facet = qh->facet_list; facet && facet->id != 0;
         facet = facet->next) {
      if (facet->upperdelaunay)
        continue;

      // Get triangle vertices
      vertexT *v0 = facet->vertices->e[0].p;
      vertexT *v1 = facet->vertices->e[1].p;
      vertexT *v2 = facet->vertices->e[2].p;

      // Get vertex coordinates
      const double x1 = v0->point[0], y1 = v0->point[1];
      const double x2 = v1->point[0], y2 = v1->point[1];
      const double x3 = v2->point[0], y3 = v2->point[1];

      // Calculate barycentric coordinates
      const double denom = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
      const double L1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denom;
      const double L2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denom;
      const double L3 = 1.0 - L1 - L2;

      // Check if point is inside triangle
      if (L1 >= 0.0 && L2 >= 0.0 && L3 >= 0.0) {
        // Get values at vertices
        const double f1 = pval[qh_pointid(qh, v0->point)];
        const double f2 = pval[qh_pointid(qh, v1->point)];
        const double f3 = pval[qh_pointid(qh, v2->point)];

        // Interpolate value
        ipval[i] = L1 * f1 + L2 * f2 + L3 * f3;
        found = true;
        break;
      }
    }

    if (!found) {
      fprintf(stdout, " -- WARNING:: point (%f, %f) outside convex hull\n", x,
              y);
    }
  }
}

/* Interpolates values at 3D query points using barycentric coordinates within
 * Qhull tetrahedra. Points outside the convex hull receive fill_value.
 *
 * For a query point p in tetrahedron (v0,v1,v2,v3), substituting
 * L4 = 1 - L1 - L2 - L3 yields the 3x3 system T*[L1,L2,L3]^T = p - v3,
 * where T has columns (v0-v3), (v1-v3), (v2-v3). Solved via Cramer's rule.
 * Vertex orientation is absorbed into the sign of det, so containment
 * (all Li >= 0) is correct regardless of Qhull's winding order. */
static void linear_interp3d(qhT *qh, const double *ipoints, double *ipval,
                             int inum_pts, const double *pval,
                             double fill_value) {
  for (int i = 0; i < inum_pts; i++) {
    const double px = ipoints[i * 3];
    const double py = ipoints[i * 3 + 1];
    const double pz = ipoints[i * 3 + 2];
    bool found = false;

    ipval[i] = fill_value;

    for (facetT *facet = qh->facet_list; facet && facet->id != 0;
         facet = facet->next) {
      if (facet->upperdelaunay)
        continue;

      vertexT *vt0 = facet->vertices->e[0].p;
      vertexT *vt1 = facet->vertices->e[1].p;
      vertexT *vt2 = facet->vertices->e[2].p;
      vertexT *vt3 = facet->vertices->e[3].p;

      const double x0 = vt0->point[0], y0 = vt0->point[1], z0 = vt0->point[2];
      const double x1 = vt1->point[0], y1 = vt1->point[1], z1 = vt1->point[2];
      const double x2 = vt2->point[0], y2 = vt2->point[1], z2 = vt2->point[2];
      const double x3 = vt3->point[0], y3 = vt3->point[1], z3 = vt3->point[2];

      /* Columns of T = (v0-v3), (v1-v3), (v2-v3) */
      const double a = x0 - x3, b = x1 - x3, c = x2 - x3;
      const double d = y0 - y3, e = y1 - y3, f = y2 - y3;
      const double g = z0 - z3, h = z1 - z3, k = z2 - z3;

      const double det = a * (e * k - f * h) - b * (d * k - f * g) + c * (d * h - e * g);
      if (fabs(det) < 1e-15)
        continue;

      const double rx = px - x3, ry = py - y3, rz = pz - z3;

      const double L1 = (rx * (e * k - f * h) - b * (ry * k - rz * h) + c * (ry * h - rz * e)) / det;
      const double L2 = (a * (ry * k - rz * h) - rx * (d * k - f * g) + c * (d * rz - ry * g)) / det;
      const double L3 = (a * (e * rz - ry * h) - b * (d * rz - ry * g) + rx * (d * h - e * g)) / det;
      const double L4 = 1.0 - L1 - L2 - L3;

      if (L1 >= 0.0 && L2 >= 0.0 && L3 >= 0.0 && L4 >= 0.0) {
        const double f0 = pval[qh_pointid(qh, vt0->point)];
        const double f1 = pval[qh_pointid(qh, vt1->point)];
        const double f2 = pval[qh_pointid(qh, vt2->point)];
        const double f3 = pval[qh_pointid(qh, vt3->point)];
        ipval[i] = L1 * f0 + L2 * f1 + L3 * f2 + L4 * f3;
        found = true;
        break;
      }
    }

    if (!found) {
      fprintf(stdout, " -- WARNING:: point (%f, %f, %f) outside convex hull\n",
              px, py, pz);
    }
  }
}

/**
 * Helper function to initialize and setup Qhull for Delaunay triangulation.
 *
 * @param points Input points array (flat interleaved, length num_pts * dims)
 * @param num_pts Number of points (must be >= dims + 2)
 * @param dims Spatial dimension (2 or 3)
 * @return Initialized qhull structure or NULL on error
 */
static qhT *init_qhull(double *points, int num_pts, int dims) {
  if (num_pts < dims + 2) {
    fprintf(stdout, " -- ERROR: Qhull needs a minimum of %d points for %dD triangulation.\n",
            dims + 2, dims);
    return NULL;
  }

  char *noptions = "QVn QJ d";
  qhT *qh = (qhT *)malloc(sizeof(qhT));
  if (!qh)
    return NULL;

  // Initialize qhull
  boolT ismalloc = False;
  qh_init_A(qh, stdin, stdout, stderr, 0, NULL);

  if (setjmp(qh->errexit) < 0) {
    free(qh);
    return NULL;
  }

  qh->NOerrexit = False;
  qh_initflags(qh, noptions);
  qh->PROJECTdelaunay = True;

  // Initialize with points
  qh_init_B(qh, points, num_pts, dims, ismalloc);

  // Create triangulation
  qh_qhull(qh);
  qh_triangulate(qh);

  return qh;
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
             double *ivalues, int inum_pts, double fill_value, int dims) {
  if (inum_pts < 1) {
    fprintf(stdout, " -- ERROR: implementation needs a minimum of one "
                    "interpolated point.\n");
    return INTERP_MIN_ERROR;
  }
  if (dims != 2 && dims != 3) {
    fprintf(stderr, " -- ERROR: dims must be 2 or 3.\n");
    return QHULL_GENERAL_ERROR;
  }

  qhT *qh = init_qhull(points, num_pts, dims);
  if (!qh)
    return QHULL_GENERAL_ERROR;

  if (dims == 2)
    linear_interp2d(qh, ipoints, ivalues, inum_pts, values, fill_value);
  else
    linear_interp3d(qh, ipoints, ivalues, inum_pts, values, fill_value);

  qh_freeqhull(qh, False);
  free(qh);
  return 0;
}

int griddata_triangles(double *points, int num_pts, int **simplex_list,
                       int *num, int dims) {
  if (dims != 2 && dims != 3) {
    fprintf(stderr, " -- ERROR: dims must be 2 or 3.\n");
    return QHULL_GENERAL_ERROR;
  }

  qhT *qh = init_qhull(points, num_pts, dims);
  if (!qh)
    return QHULL_GENERAL_ERROR;

  int result = create_simplex_indices(qh, simplex_list, num, dims);

  qh_freeqhull(qh, False);
  free(qh);
  return result;
}
