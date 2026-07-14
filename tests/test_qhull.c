#include <stdio.h>
#include <stdlib.h>

#include <cmocka.h>

#include "../include/interpolation.h"

/**
 * Interpolating values from a diamond shape.
 */
static void test_equilateral_triangle(void **state) {
  (void)state;
  double points[] = {0.0, 0.0, 1.0, 2.0, 2.0, 0.0, 1.0, -2.0};
  double point_values[] = {1.0, 1.0, 1.0, 1.0};
  double fill_value = -9999.0;
  double ipoints_loc[] = {1.0, 1.0, 1.0, -1.0};
  double ipoint_values[] = {0.0, 0.0};

  double expected_values[] = {1.0, 1.0};

  int value = griddata(points, point_values, 4, ipoints_loc, ipoint_values, 2,
                       fill_value, 2);
  assert_int_equal(value, 0);
  assert_float_equal(ipoint_values[0], expected_values[0], 1.0e-6);
  assert_float_equal(ipoint_values[1], expected_values[1], 1.0e-6);
}

/**
 * Interpolating values from a square.
 */
static void test_qhull_square(void **state) {
  (void)state;
  double points[] = {0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, 0.0};
  double point_values[] = {1.0, 1.0, 1.0, 1.0};
  double fill_value = -9999.0;
  double ipoints_loc[] = {0.5, 0.75, 0.5, 0.75};
  double ipoint_values[] = {0.0, 0.0};

  double expected_values[] = {1.0, 1.0};

  int value = griddata(points, point_values, 4, ipoints_loc, ipoint_values, 2,
                       fill_value, 2);
  assert_int_equal(value, 0);
  assert_float_equal(ipoint_values[0], expected_values[0], 1.0e-6);
}

/**
 * Get triangular mesh results.
 */
static void test_qhull_triangle_list(void **state) {
  (void)state;
  double points[] = {0.0, 0.0, 2.0, 4.0, 8.0, -2.0, 5.0, 7.0, 10.0, 4.0};
  int *triangle_list = NULL;
  int num = -1;

  int value = griddata_triangles(points, 5, &triangle_list, &num, 2);
  assert_int_equal(value, 0);
  assert_int_equal(num, 3);

  assert_int_equal(triangle_list[0], 1);
  assert_int_equal(triangle_list[1], 2);
  assert_int_equal(triangle_list[2], 0);
  assert_int_equal(triangle_list[3], 1);
  assert_int_equal(triangle_list[4], 2);
  assert_int_equal(triangle_list[5], 4);
  assert_int_equal(triangle_list[6], 1);
  assert_int_equal(triangle_list[7], 3);
  assert_int_equal(triangle_list[8], 4);
}


/* Constant field in 3D: all values equal, centroid of a tetrahedron must
 * interpolate to that same value. */
static void test_griddata_3d_constant_field(void **state) {
  (void)state;
  double points[] = {
    0.0, 0.0, 0.0,
    2.0, 0.0, 0.0,
    0.0, 2.0, 0.0,
    0.0, 0.0, 2.0,
    2.0, 2.0, 2.0,
  };
  double values[] = {3.0, 3.0, 3.0, 3.0, 3.0};
  double ipoints[] = {0.5, 0.5, 0.5};
  double ivalues[] = {0.0};

  int rc = griddata(points, values, 5, ipoints, ivalues, 1, -9999.0, 3);
  assert_int_equal(rc, 0);
  assert_float_equal(ivalues[0], 3.0, 1.0e-6);
}

/* 3D linear field: value = x + y + z. Centroid of the tet (0,0,0)(2,0,0)(0,2,0)(0,0,2)
 * is (0.5,0.5,0.5), expected interpolated value = 1.5. */
static void test_griddata_3d_linear_field(void **state) {
  (void)state;
  double points[] = {
    0.0, 0.0, 0.0,
    2.0, 0.0, 0.0,
    0.0, 2.0, 0.0,
    0.0, 0.0, 2.0,
    2.0, 2.0, 2.0,
  };
  double values[] = {0.0, 2.0, 2.0, 2.0, 6.0};
  double ipoints[] = {0.5, 0.5, 0.5};
  double ivalues[] = {0.0};

  int rc = griddata(points, values, 5, ipoints, ivalues, 1, -9999.0, 3);
  assert_int_equal(rc, 0);
  assert_float_equal(ivalues[0], 1.5, 1.0e-6);
}

/* griddata_triangles in 3D: 5 points in general position produce 2 tetrahedra,
 * so num == 2 and the flat output array has 2*4 = 8 valid index entries. */
static void test_griddata_triangles_3d(void **state) {
  (void)state;
  double points[] = {
    0.0, 0.0, 0.0,
    2.0, 0.0, 0.0,
    0.0, 2.0, 0.0,
    0.0, 0.0, 2.0,
    2.0, 2.0, 2.0,
  };
  int *simplex_list = NULL;
  int num = -1;

  int rc = griddata_triangles(points, 5, &simplex_list, &num, 3);
  assert_int_equal(rc, 0);
  assert_int_equal(num, 2);

  for (int i = 0; i < num * 4; i++) {
    assert_in_range(simplex_list[i], 0, 4);
  }
  free(simplex_list);
}

/* Requesting 3D with only 4 points (need >= 5) must return QHULL_GENERAL_ERROR. */
static void test_griddata_3d_min_points_error(void **state) {
  (void)state;
  double points[] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  double values[] = {1.0, 1.0, 1.0, 1.0};
  double ipoints[] = {0.25, 0.25, 0.25};
  double ivalues[] = {0.0};

  int rc = griddata(points, values, 4, ipoints, ivalues, 1, -9999.0, 3);
  assert_int_equal(rc, QHULL_GENERAL_ERROR);
}


int main(void) {
  const struct CMUnitTest test[] = {
    cmocka_unit_test(test_equilateral_triangle),
    cmocka_unit_test(test_qhull_square),
    cmocka_unit_test(test_qhull_triangle_list),
    cmocka_unit_test(test_griddata_3d_constant_field),
    cmocka_unit_test(test_griddata_3d_linear_field),
    cmocka_unit_test(test_griddata_triangles_3d),
    cmocka_unit_test(test_griddata_3d_min_points_error),
  };

  return cmocka_run_group_tests(test, NULL, NULL);
};

