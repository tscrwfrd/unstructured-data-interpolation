#ifndef DELAUNATOR_H
#define DELAUNATOR_H

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// Constants
#define EPSILON pow(2, -52)
#define MAX_POINTS 10000
#define EDGE_STACK_SIZE 512
#define SPLITTER 134217729
#define INEXACT

// Macros based on Shewchuk code
// https://www.cs.cmu.edu/afs/cs/project/quake/public/code/
#define Absolute(a) ((a) >= 0.0 ? (a) : -(a))

#define Fast_Two_Sum_Tail(a, b, x, y) \
  bvirt = x - a;                      \
  y = b - bvirt

#define Fast_Two_Sum(a, b, x, y) \
  x = (double)(a + b);           \
  Fast_Two_Sum_Tail(a, b, x, y)

#define Split(a, ahi, alo)    \
  c = (double)(SPLITTER * a); \
  abig = (double)(c - a);     \
  ahi = c - abig;             \
  alo = a - ahi

#define Two_Product_Tail(a, b, x, y) \
  Split(a, ahi, alo);                \
  Split(b, bhi, blo);                \
  err1 = x - (ahi * bhi);            \
  err2 = err1 - (alo * bhi);         \
  err3 = err2 - (ahi * blo);         \
  y = (alo * blo) - err3

#define Two_Diff_Tail(a, b, x, y) \
  bvirt = (double)(a - x);        \
  avirt = x + bvirt;              \
  bround = bvirt - b;             \
  around = a - avirt;             \
  y = around + bround

#define Two_Sum_Tail(a, b, x, y) \
  bvirt = (double)(x - a);       \
  avirt = x - bvirt;             \
  bround = b - bvirt;            \
  around = a - avirt;            \
  y = around + bround

#define Two_Sum(a, b, x, y) \
  x = (double)(a + b);      \
  Two_Sum_Tail(a, b, x, y)

#define Two_Diff(a, b, x, y) \
  x = (double)(a - b);       \
  Two_Diff_Tail(a, b, x, y)

#define Two_One_Sum(a1, a0, b, x2, x1, x0) \
  Two_Sum(a0, b, _i, x0);                  \
  Two_Sum(a1, _i, x2, x1)

#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
  Two_Diff(a0, b, _i, x0);                  \
  Two_Sum(a1, _i, x2, x1)

#define Two_Product(a, b, x, y) \
  x = (double)(a * b);          \
  Two_Product_Tail(a, b, x, y)

#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
  Two_Diff(a0, b, _i, x0);                  \
  Two_Sum(a1, _i, x2, x1)

#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Diff(a1, a0, b0, _j, _0, x0);              \
  Two_One_Diff(_j, _0, b1, x3, x2, x1)

typedef struct {
  int num_coords;
  int hull_start;
  double *coords;        // Array of coordinates [x0, y0, x1, y1, ...]
  unsigned *triangles;   // Array of triangle indices
  int *half_edges;       // Array of halfedge indices
  unsigned *hull;        // Convex hull indices
  unsigned *hull_prev;   // Previous hull vertex
  unsigned *hull_next;   // Next hull vertex
  unsigned *hull_tri;    // Triangle associated with hull vertex
  int *hull_hash;        // Hash table for hull vertices
  unsigned *ids;         // Point indices
  double *dists;         // Point distances
  int triangles_len;     // Number of triangles
  int hull_size;         // Size of convex hull
  int hash_size;         // Size of hash table
  double cx;             // Center x coordinate
  double cy;             // Center y coordinate
  unsigned *edge_stack;  //
} Delaunator;

// Function declarations
Delaunator *delaunator_create(double *coords, int n_points);
void update(Delaunator *del);
void delaunator_free(Delaunator *d);

#endif  // DELAUNATOR_H
