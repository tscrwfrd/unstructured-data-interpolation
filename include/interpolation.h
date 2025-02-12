#ifndef _INTERPOLATE_H
#define _INTERPOLATE_H

#define INTERP_MIN_ERROR -30

int griddata(double* points, double* values, int num_pts, double* ipoints,
             double* ivalues, int inum_pts, double fill_value);

int griddata_triangles(double* points, int num_pts);

#endif
