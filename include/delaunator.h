#ifndef DELAUNATOR_H
#define DELAUNATOR_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Constants
#define EPSILON pow(2, -52)
#define MAX_POINTS 10000
#define EDGE_STACK_SIZE 512
#define SPLITTER 134217729
#define INEXACT                     
/* #define INEXACT volatile */

// Macros
#define Absolute(a)  ((a) >= 0.0 ? (a) : -(a))

#define Fast_Two_Sum_Tail(a, b, x, y) \
  bvirt = x - a; \
  y = b - bvirt

#define Fast_Two_Sum(a, b, x, y) \
  x = (double) (a + b); \
  Fast_Two_Sum_Tail(a, b, x, y)

#define Split(a, ahi, alo) \
    c = (double)(SPLITTER * a); \
    abig = (double)(c - a); \
    ahi = c - abig; \
    alo = a - ahi

#define Two_Product_Tail(a, b, x, y) \
    Split(a, ahi, alo); \
    Split(b, bhi, blo); \
    err1 = x - (ahi * bhi); \
    err2 = err1 - (alo * bhi); \
    err3 = err2 - (ahi * blo); \
    y = (alo * blo) - err3

#define Two_Diff_Tail(a, b, x, y) \
  bvirt = (double) (a - x); \
  avirt = x + bvirt; \
  bround = bvirt - b; \
  around = a - avirt; \
  y = around + bround

#define Two_Sum_Tail(a, b, x, y) \
  bvirt = (double) (x - a); \
  avirt = x - bvirt; \
  bround = b - bvirt; \
  around = a - avirt; \
  y = around + bround

#define Two_Sum(a, b, x, y) \
  x = (double) (a + b); \
  Two_Sum_Tail(a, b, x, y)

#define Two_Diff(a, b, x, y) \
  x = (double) (a - b); \
  Two_Diff_Tail(a, b, x, y)

#define Two_One_Sum(a1, a0, b, x2, x1, x0) \
  Two_Sum(a0, b , _i, x0); \
  Two_Sum(a1, _i, x2, x1)

#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
  Two_Diff(a0, b , _i, x0); \
  Two_Sum( a1, _i, x2, x1)

#define Two_Product(a, b, x, y) \
    x = (double)(a * b); \
    Two_Product_Tail(a, b, x, y)

#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
  Two_Diff(a0, b , _i, x0); \
  Two_Sum( a1, _i, x2, x1)

#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Diff(a1, a0, b0, _j, _0, x0); \
  Two_One_Diff(_j, _0, b1, x3, x2, x1)

typedef struct {
    int ncoords;
    int hullstart;
    double* coords;        // Array of coordinates [x0, y0, x1, y1, ...]
    unsigned* triangles;   // Array of triangle indices
    int* halfedges;        // Array of halfedge indices
    unsigned* hull;        // Convex hull indices
    unsigned* hullPrev;    // Previous hull vertex
    unsigned* hullNext;    // Next hull vertex
    unsigned* hullTri;     // Triangle associated with hull vertex
    int* hullHash;         // Hash table for hull vertices
    unsigned* ids;         // Point indices
    double* dists;         // Point distances
    int trianglesLen;      // Number of triangles
    int hullSize;          // Size of convex hull
    int hashSize;          // Size of hash table
    double cx;             // Center x coordinate
    double cy;             // Center y coordinate
    unsigned* edge_stack;  // 
} Delaunator;

// Function declarations
Delaunator* delaunator_create(double* coords, int n_points);
void update(Delaunator* del);
void delaunator_destroy(Delaunator* d);

static double orient2d(double pa0, double pa1, double pb0, double pb1, double pc0, double pc1);
static double orient2dadapt(double pa0, double pa1, double pb0, double pb1, double pc0, double pc1, double detsum);
static int in_circle(const double ax, const double ay, const double bx, const double by, 
                    const double cx, const double cy, const double px, const double py);
static void quicksort(int* ids, double* dists, int left, int right);
static double circumradius(double ax, double ay, double bx, double by, double cx, double cy);
static double* circumcenter(double ax, double ay, double bx, double by, double cx, double cy);
static double dist(double ax, double ay, double bx, double by);
static double pseudo_angle(double dx, double dy);
static void swap(int *arr, int i, int j);
static int hash_key(double x, double y, Delaunator* d);
static int add_triangle(int i0, int i1, int i2, int a, int b, int c, Delaunator* d);
static int legalize(int a, Delaunator* d);
static void link(int a, int b, Delaunator* d);
static double estimate(int elen, double* e);
static double* circumcenter(double ax, double ay, double bx, double by, double cx, double cy);
#endif // DELAUNATOR_H
