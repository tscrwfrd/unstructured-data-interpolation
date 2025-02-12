#include <stdio.h>
#include <stdlib.h>

#include <cmocka.h>

#include "../include/interpolation.h"

/**
 * Interpolating values from a diamond shape.
 */
static void test_trianglulation(void **state) {
  (void)state;
  double points[] = {0.0, 0.0, 1.0, 2.0, 2.0, 0.0, 1.0, -2.0};
  double point_values[] = {1.0, 1.0, 1.0, 1.0};
  double fill_value = -9999.0;
  double ipoints_loc[] = {1.0, 1.0, 1.0, -1.0};
  double ipoint_values[] = {0.0, 0.0};

  double expected_values[] = {1.0, 1.0};

  int value = griddata(points, point_values, 4, ipoints_loc, ipoint_values, 2, fill_value);
  assert_int_equal(value, 0);
  assert_float_equal(ipoint_values[0], expected_values[0], 1.0e-6);
  assert_float_equal(ipoint_values[1], expected_values[1], 1.0e-6);
}

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

  int value = griddata(points, point_values, 4, ipoints_loc, ipoint_values, 2, fill_value);
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

  int value = griddata(points, point_values, 4, ipoints_loc, ipoint_values, 2, fill_value);
  assert_int_equal(value, 0);
  assert_float_equal(ipoint_values[0], expected_values[0], 1.0e-6);

}

int main(void) {
  const struct CMUnitTest test[] = {
    cmocka_unit_test(test_equilateral_triangle),
    cmocka_unit_test(test_qhull_square)
  };

  return cmocka_run_group_tests(test, NULL, NULL);
};
