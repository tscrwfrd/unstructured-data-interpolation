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
                       fill_value);
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
                       fill_value);
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

  int value = griddata_triangles(points, 5, &triangle_list, &num);
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

#ifdef TESTING
int main(void) {
  const struct CMUnitTest test[] = {cmocka_unit_test(test_equilateral_triangle),
                                    cmocka_unit_test(test_qhull_square),
                                    cmocka_unit_test(test_qhull_triangle_list)};

  return cmocka_run_group_tests(test, NULL, NULL);
};
#endif
