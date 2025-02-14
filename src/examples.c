#include <stdio.h>
#include <stdlib.h>

#include "../../qhull/src/libqhull_r/libqhull_r.h"
#include "../include/delaunator.h"
#include "../include/interpolation.h"

//   EXAMPLE INPUT POINTS AND FACETS
// +++++++++++++++++++++++++++++++++++++++++++++++++
//
//           *(0,6)
//
//
//                           * (3,4)
//
//
//             * (1,3)
//
//                                * (4,2)
//
//                    *(2,1)
//
//                             X (3, 0) <<-- point location
//     *
//   (-1,-1)
//
//                           * (3,-2)
//
// +++++++++++++++++++++++++++++++++++++++++++++++++
//
//  FACETS: (1, 3) -> (0, 6) -> (-1,-1)
//  FACETS: (3, 4) -> (1, 3) -> ( 0, 6)
//  FACETS: (2, 1) -> (3,-2) -> ( 4, 2)
//  FACETS: (2, 1) -> (3,-2) -> (-1,-1)
//  FACETS: (2, 1) -> (1, 3) -> (-1,-1)
//  FACETS: (2, 1) -> (3, 4) -> ( 1, 3)
//  FACETS: (2, 1) -> (3, 4) -> ( 4, 2)
//  ======================================
int interpolation_qhull(void);
int interpolation_delaunator(void);
int triangle_mesh_qhull(void);
int triangle_mesh_delaunator(void);

int main() {
  if (interpolation_qhull() != 0) {
    fprintf(stderr, "Interpolation failed.\n");
  }
  if (interpolation_delaunator() != 0) {
    fprintf(stderr, "Interpolation failed.\n");
  }
  if (triangle_mesh_qhull() != 0) {
    fprintf(stderr, "Triangle meshing (qhull) failed.\n");
  }
  if (triangle_mesh_delaunator() != 0) {
    fprintf(stderr, "Triangle meshing (delaunator) failed.\n");
  }

  return 0;
}

int interpolation_qhull(void) {

  int result = 0;

  const int NUM_KNOWN_POINTS = 7;
  const int NUM_INTERP_POINTS = 2;
  const double FILL_VALUE = 9.5;
  // given point data
  double points[] = {0.0, 6.0, -1.0, -1.0, 1.0, 3.0, 2.0,
                     1.0, 3.0, 4.0,  4.0,  2.0, 3.0, -2.0};
  double point_values[] = {3.4, 4.4, 5.4, 6.4, 7.4, 8.4, 9.4};
  // points to interpolate
  double interp_points[] = {3.0, 0.0, 1.0, 40.0};
  double interp_values[] = {0.0, 0.0};

  result = griddata(points, point_values, NUM_KNOWN_POINTS, interp_points,
                    interp_values, NUM_INTERP_POINTS, FILL_VALUE);
  if (result != 0) {
    fprintf(stderr, "Interpolation failed\n");
  }

  for (int i = 0; i < NUM_INTERP_POINTS; i++) {
    printf("Interpolated value %d: %lf\n", i + 1, interp_values[i]);
  }

  return result;
}

int interpolation_delaunator(void) {

  int result = 0;

  const int NUM_KNOWN_POINTS = 7;
  const int NUM_INTERP_POINTS = 2;
  const double FILL_VALUE = 9.5;
  // given point data
  double points[] = {0.0, 6.0, -1.0, -1.0, 1.0, 3.0, 2.0,
                     1.0, 3.0, 4.0,  4.0,  2.0, 3.0, -2.0};
  double point_values[] = {3.4, 4.4, 5.4, 6.4, 7.4, 8.4, 9.4};
  // points to interpolate
  double interp_points[] = {3.0, 0.0, 1.0, 40.0};
  double interp_values[] = {0.0, 0.0};

  result = griddata(points, point_values, NUM_KNOWN_POINTS, interp_points,
                    interp_values, NUM_INTERP_POINTS, FILL_VALUE);
  if (result != 0) {
    fprintf(stderr, "Interpolation failed\n");
  }

  for (int i = 0; i < NUM_INTERP_POINTS; i++) {
    printf("Interpolated value %d: %lf\n", i + 1, interp_values[i]);
  }

  return result;
}

int triangle_mesh_qhull(void) {
  int result = 0;
  const int NUM_KNOWN_POINTS = 7;
  // given point data
  double points[] = {0.0, 6.0, -1.0, -1.0, 1.0, 3.0, 2.0,
                     1.0, 3.0, 4.0,  4.0,  2.0, 3.0, -2.0};

  int *triangle_list = NULL;
  int num_triangles = -1;
  griddata_triangles(points, NUM_KNOWN_POINTS, &triangle_list, &num_triangles);

  printf(" >>>>>>  QHULL  <<<<<<<\n");
  printf("Number of triangles ==> %d\n", num_triangles);
  int tnum = 0;
  for (int i = 0; i < num_triangles * 3; i += 3) {
    printf("Triangle: %d\n", tnum);
    printf("   %d -> %d -> %d \n", triangle_list[i], triangle_list[i + 1],
           triangle_list[i + 2]);
    tnum += 1;
  }

  free(triangle_list);

  return result;
}

int triangle_mesh_delaunator(void) {
  int result = 0;
  const int NUM_KNOWN_POINTS = 7;
  double points[] = {0.0, 6.0, -1.0, -1.0, 1.0, 3.0, 2.0,
                     1.0, 3.0, 4.0,  4.0,  2.0, 3.0, -2.0};

  Delaunator *d = delaunator_create(points, NUM_KNOWN_POINTS*2);
  update(d);

  printf(" >>>>>>  DELAUNATOR  <<<<<<<\n");
  printf("Number of triangles ==> %d\n", d->triangles_len / 3);
  int tnum = 0;
  for (int i = 0; i < d->triangles_len; i += 3) {
    printf("Triangle: %d\n", tnum);
    printf("   %d -> %d -> %d \n", d->triangles[i], d->triangles[i + 1],
           d->triangles[i + 2]);
    tnum += 1;
  }

  delaunator_free(d);
  return result;
}

