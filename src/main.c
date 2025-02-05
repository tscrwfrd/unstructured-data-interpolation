#include <stdio.h>
#include <stdlib.h>

#include "../../qhull/src/libqhull_r/libqhull_r.h"
#include "../include/delaunator.h"
#include "../include/interpolation.h"

//   EXAMPLE 1 INPUT POINTS AND FACETS
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

// void example1();
// void example2();
void example3();
double *create_points1();
double *create_points2();
double *create_func_values1();
double *create_func_values2();
double *create_points_test();
double *create_test_values();
void _print_all_facets(qhT *qh);

int main() {
  // example1();
  // example2();
  example3();
}

void example1() {
  double *pts = create_points1();
  double *pval = create_func_values1();
  double fill_value = 9.5;

  double *ipts = malloc(sizeof(double *) * 6);
  ipts[0] = 3.0;
  ipts[1] = 0.0;
  ipts[2] = 1.0;
  ipts[3] = 4.2;
  ipts[4] = 1.0;
  ipts[5] = 40.2;

  double *ipval = malloc(sizeof(double *) * 3);
  ipval[0] = 0.0;
  ipval[1] = 0.0;
  ipval[2] = 0.0;

  griddata(pts, pval, 7, ipts, ipval, 3, fill_value);

  printf("interpolated value 1: %lf \n", ipval[0]);
  printf("interpolated value 2: %lf \n", ipval[1]);
  printf("interpolated value 3: %lf \n", ipval[2]);

  free(pts);
  free(pval);
  free(ipts);
  free(ipval);
}

void example2() {
  double *pts = create_points2();
  double *pval = create_func_values2();
  double fill_value = 0.838;

  double *ipts = malloc(sizeof(double *) * 2);
  ipts[0] = 237.7;
  ipts[1] = 47.6;

  double *ipval = malloc(sizeof(double *) * 1);
  ipval[0] = 0.0;

  griddata(pts, pval, 7, ipts, ipval, 1, fill_value);

  printf("interpolated value: %lf \n", ipval[0]);

  free(pts);
  free(pval);
  free(ipts);
  free(ipval);
}

void example3() {
  double *pts = create_points_test();
  double *pval = create_test_values();
  double fill_value = -9999.0;

  double *ipts = malloc(sizeof(double *) * 4);
  ipts[0] = 3.0;
  ipts[1] = 0.0;
  ipts[2] = 1.0;
  ipts[3] = 7.2;
  double *ipval = calloc(3, sizeof(double *));

  griddata(pts, pval, 5, ipts, ipval, 2, fill_value);

  double *test_points = create_points_test();
  Delaunator *d = delaunator_create(test_points, 10);
  update(d);

  for (int i = 0; i < d->triangles_len; i += 3) {
    printf("triangle: %d, %d, %d\n", d->triangles[i], d->triangles[i + 1],
           d->triangles[i + 2]);
  }
  printf(" <> Triangle length: %d \n", d->triangles_len);

  delaunator_free(d);

  printf(" <---->  <---->   \n");
  printf("interpolated value 1: %f \n", ipval[0]);
  printf("interpolated value 2: %f \n", ipval[1]);

  free(pts);
  free(pval);
  free(ipts);
  free(ipval);
}

double *create_points_test() {
  double *T;

  T = malloc(sizeof(double *) * 10);

  T[0] = 0.0;
  T[1] = 0.0;

  T[2] = 2.0;
  T[3] = 4.0;

  T[4] = 8.0;
  T[5] = -2.0;

  T[6] = 5.0;
  T[7] = 7.0;

  T[8] = 10.0;
  T[9] = 4.0;

  return T;
}

double *create_test_values() {
  double *f_val = malloc(sizeof(double *) * 5);

  f_val[0] = 110.4;
  f_val[1] = 150.4;
  f_val[2] = 111.4;
  f_val[3] = 112.4;
  f_val[4] = 115.4;

  return f_val;
}

/**
 * Test function that creates 2D points.
 */
// +++++++++++++++++++++++++++++++++++++++++++++++++
//
//           *(0,6)[0]
//
//
//                           * (3,4)[4]
//
//
//             * (1,3)[2]
//
//                                * (4,2)[5]
//
//                    *(2,1)[3]
//
//
//     *
//   (-1,-1)[1]
//
//                           * (3,-2)[6]
//
// +++++++++++++++++++++++++++++++++++++++++++++++++
//  FACETS: (1, 3) -> (0, 6) -> (-1,-1)
//  FACETS: (3, 4) -> (1, 3) -> ( 0, 6)
//  FACETS: (2, 1) -> (3,-2) -> ( 4, 2)
//  FACETS: (2, 1) -> (3,-2) -> (-1,-1)
//  FACETS: (2, 1) -> (1, 3) -> (-1,-1)
//  FACETS: (2, 1) -> (3, 4) -> ( 1, 3)
//  FACETS: (2, 1) -> (3, 4) -> ( 4, 2)
double *create_points1() {
  double *T;
  int num_pts = 14;

  T = malloc(sizeof(double *) * num_pts);

  T[0] = 0.0;
  T[1] = 6.0;

  T[2] = -1.0;
  T[3] = -1.0;

  T[4] = 1.0;
  T[5] = 3.0;

  T[6] = 2.0;
  T[7] = 1.0;

  T[8] = 3.0;
  T[9] = 4.0;

  T[10] = 4.0;
  T[11] = 2.0;

  T[12] = 3.0;
  T[13] = -2.0;

  return T;
}

/**
 * Test function that creates 2D points.
 */
double *create_points2() {
  double *T;
  int num_pts = 14;

  T = malloc(sizeof(double *) * num_pts);

  T[0] = 262.3;
  T[1] = 30.4;

  T[2] = 254.7;
  T[3] = 40.0;

  T[4] = 273.5;
  T[5] = 30.5;

  T[6] = 164.7;
  T[7] = -14.600000000000001;

  T[8] = 363.5;
  T[9] = -14.600000000000001;

  T[10] = 363.5;
  T[11] = 85.0;

  T[12] = 164.7;
  T[13] = 85.0;

  return T;
}

/**
 * Test function that creates function values
 *
 *
 */
double *create_func_values1() {
  double *f_val = malloc(sizeof(double *) * 7);

  f_val[0] = 3.4;
  f_val[1] = 4.4;
  f_val[2] = 5.4;
  f_val[3] = 6.4;
  f_val[4] = 7.4;
  f_val[5] = 8.4;
  f_val[6] = 9.4;

  return f_val;
}

/**
 * Test function that creates function values
 *
 *
 */
double *create_func_values2() {
  double *f_val = malloc(sizeof(double *) * 7);

  f_val[0] = 0.8166254934793893;
  f_val[1] = 0.972422825781967;
  f_val[2] = 0.7249863954040258;
  f_val[3] = 0.8380115715551274;
  f_val[4] = 0.8380115715551274;
  f_val[5] = 0.8380115715551274;
  f_val[6] = 0.8380115715551274;

  return f_val;
}
