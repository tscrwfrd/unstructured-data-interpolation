#ifndef _INTERPOLATE_H
#define _INTERPOLATE_H

#define TRIANGLE_INDICES_ERROR -40
#define QHULL_GENERAL_ERROR -30
#define INTERP_MIN_ERROR -20
#define INT_ALLOCATION_ERROR -10

int griddata(double* points, double* values, int num_pts, double* ipoints,
             double* ivalues, int inum_pts, double fill_value);

int griddata_triangles(double* points, int num_pts, int** tri_list, int* num);

#endif
