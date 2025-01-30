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

// Macros
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

#define Two_Product(a, b, x, y) \
    x = (double)(a * b); \
    Two_Product_Tail(a, b, x, y)

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
static double orient2adapt(double pa0, double pa1, double pb0, double pb1, double pc0, double pc1, double detsum);
static int in_circle(const double ax, const double ay, const double bx, const double by, 
                    const double cx, const double cy, const double px, const double py);
static void quicksort(int* ids, double* dists, int left, int right);
static double circumradius(double ax, double ay, double bx, double by, double cx, double cy);
static double* circumcenter(double ax, double ay, double bx, double by, double cx, double cy);
static double dist(double ax, double ay, double bx, double by);
static double pseudo_angle(double dx, double dy);
static void swap(double *arr, int i, int j);

#endif // DELAUNATOR_H
