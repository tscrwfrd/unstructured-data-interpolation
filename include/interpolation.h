#ifndef _INTERPOLATE_H
#define _INTERPOLATE_H

#define TRIANGLE_INDICES_ERROR -40
#define QHULL_GENERAL_ERROR -30
#define INTERP_MIN_ERROR -20
#define INT_ALLOCATION_ERROR -10

/* dims must be 2 or 3. Points are flat interleaved [x0,y0,...] for 2D or
 * [x0,y0,z0,...] for 3D. fill_value is used for query points outside the
 * convex hull. */
int griddata(double* points, double* values, int num_pts, double* ipoints,
             double* ivalues, int inum_pts, double fill_value, int dims);

/* Output simplex_list is a flat array of vertex indices with stride dims+1
 * (triangles for 2D, tetrahedra for 3D). Caller must free simplex_list. */
int griddata_triangles(double* points, int num_pts, int** simplex_list,
                       int* num, int dims);

#endif
