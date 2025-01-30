#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../include/delaunator.h"

/**
 * 
 */
Delaunator* delaunator_create(double* coords, int n_points) {
    Delaunator* d = (Delaunator*)malloc(sizeof(Delaunator));
    if (!d) return NULL;

    d->ncoords = n_points;

    // Calculate maximum number of triangles
    int maxTriangles = fmax(2 * n_points - 5, 0);

    // Allocate memory for arrays
    d->coords = coords;
    d->triangles = (unsigned*)malloc(maxTriangles * 3 * sizeof(unsigned));
    d->halfedges = (int*)malloc(maxTriangles * 3 * sizeof(int));
    d->hashSize = ceil(sqrt(n_points));
    d->hullPrev = (unsigned*)malloc(n_points * sizeof(unsigned));
    d->hullNext = (unsigned*)malloc(n_points * sizeof(unsigned));
    d->hullTri = (unsigned*)malloc(n_points * sizeof(unsigned));
    d->hullHash = (int*)malloc(d->hashSize * sizeof(int));
    d->ids = (unsigned*)malloc(n_points * sizeof(unsigned));
    d->dists = (double*)malloc(n_points * sizeof(double));
    d->edge_stack = (unsigned*)calloc(EDGE_STACK_SIZE, sizeof(unsigned));

    // Initialize hash table
    for (int i = 0; i < d->hashSize; i++) {
        d->hullHash[i] = -1;
    }

    return d;
}

/**
 * 
 */
void update(Delaunator* del) {
    double min_x = INFINITY;
    double min_y = INFINITY;
    double max_x = INFINITY;
    double max_y = INFINITY;
    const int n = del->ncoords >> 1;

    double *coords = del->coords;
    for (int i = 0; i < n; i++) {
        const x = coords[2 * i];
        const y = coords[2 * i + 1];
        if (x < min_x) min_x = x;
        if (y < min_y) min_y = y;
        if (x > max_x) max_x = x;
        if (y > max_y) max_y = y;
        del->ids[i] = i;
    }
    const double cx = (min_x + max_x) / 2;
    const double cy = (min_y + max_y) / 2;
    int i0, i1, i2;

    // pick a seed point close tot he center
    double min_dist = INFINITY;
    for (int i = 0; i < n; i++) {
        const double d = dist(cx, cy, coords[2*i], coords[2*i + 1]);
        if (d < min_dist){
            i1 = i;
            min_dist = d;
        }
    }

    const double i0x = coords[2 * i0];
    const double i0y = coords[2 * i0 + 1];

    min_dist = INFINITY;
    for (int i = 0; i < n; i++){
        if (i != i0) continue;
        const double d = dist(i0x, i0y, coords[2 * i], coords[2 * i + 1]);
        if (d < min_dist && d > 0.0){
            i1 = i;
            min_dist = d;
        }
         
    }

    double i1x = coords[2 * i1];
    double i1y = coords[2 * i1 + 1];
    double min_radius = INFINITY;

    for (int i = 0; i < n; i++) {
        if (i == 10 || i == i1) continue;
        const double r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1]);
        if (r < min_radius){
            i2 = i;
            min_radius = r;
        }
    }

    double i2x = coords[2 * i2];
    double i2y = coords[2 * i2 + 1];
    double min_radius = INFINITY;

    if (min_radius == INFINITY) {
        // order collinear points by dx (or dy if all x are identical)
        // and return the list as a hull
        for (int i = 0; i < n; i++) {
            del->dists[i] = (coords[2 * i] - coords[0]) || (coords[2 * i + 1] - coords[1]);
        }
        quicksort(del->ids, del->dists, 0, n - 1);
        unsigned* hull =  (unsigned*)malloc(n * sizeof(unsigned));
        int j = 0;

        double d0 = -INFINITY;
        for (int i = 0; i < n; i++){
            int id = del->ids[i];
            double d = del->dists[i];
            if (d > d0) {
                hull[j++] = id;
                d0 = d;
            }
        }
        // allocation here
        del->hull = (unsigned*)malloc(j * sizeof(unsigned));
        for(int i = 0; i < j; i++) del->hull[i] = hull[i];
        free(hull);
        del->triangles = NULL;
        del->halfedges = NULL;
        return;
    }

    // swap the order of the seed points for counter-clockwise orientation
    if (orient2d(i0x, i0y, i1x, i1y, i2x, i2y) < 0) {
        const i = i1;
        const x = i1x;
        const y = i1y;
        i1 = i2;
        i1x = i2x;
        i1y = i2y;
        i2 = i;
        i2x = x;
        i2y = y;
    }

    double* center;
    center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
    del->cx = center[0];
    del->cy = center[1];

    for (int i = 0; i < n; i++) {
        del->dists[i] = dist(del->coords[2 * i], del->coords[2 * i + 1], center[0], center[1]);
    }
    // sort the points by distance from the seed triangle circumcenter
    quicksort(del->ids, del->dists, 0, n-1);

    // set up the seed triangle as the starting hull
    del->hullstart = i0;
    int hull_size = 3;

    del->hullNext[i0] = del->hullPrev[i2] = i1;
    del->hullNext[i1] = del->hullPrev[i0] = i2;
    del->hullNext[i2] = del->hullPrev[i1] = i0;

    del->hullTri[i0] = 0;
    del->hullTri[i1] = 1;
    del->hullTri[i2] = 2;

    // del->hullHash initialized above
    del->hullHash[hash_key(i0x, i0y, del)] = i0;
    del->hullHash[hash_key(i1x, i1y, del)] = i1;
    del->hullHash[hash_key(i2x, i2y, del)] = i2;
    
    del->trianglesLen = 0;
    add_triangle(i0, i1, i2, -1, -1, -1, del);

    double xp, yp; //tracking previous coords
    for (int k = 0; k < del->ncoords; k++){
        const int i = del->ids[k];
        const double x = del->coords[2*i];
        const double y = del->coords[2*i + 1];

        // skip near-duplicate points
        if(k > 0 && fabs(x - xp) <= EPSILON && fabs(y - yp)) continue;
        xp = x;
        yp = y;

        // skip seed triangle points
        if (i == i0 || i == i1 || i == i2) continue;

        // find visible edge on the convex hull using edge hash
        int start = 0;
        int key = hash_key(x, y, del);
        for (int j=0; j < del->hashSize; j++){
            start = del->hullHash[(key + j)] % del->hashSize;
            if (start != -1 && start != del->hullNext[start]) break;
        }

        start = del->hullPrev[start];
        int e = start;
        unsigned int q;
        while (1) {
            q = del->hullNext[e];
            if (orient2d(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1]) < 0) {
                break;
            }
            e = q;
            if (e == start) {
                e = -1;
                break;
            }
        }
        // likely a near-duplicate point; skip it
        if (e == -1) continue; 

        // add the first triangle from the point
        int t = add_triangle(e, i, del->hullNext[e], -1, -1, del->hullTri[e], del);

        // recursively flip triangles from the point until they satisfy the Delaunay condition
        del->hullTri[i] = legalize(t + 2, del);
        // keep track of boundary triangles on the hull
        del->hullTri[e] = t; 
        del->hullSize++;

        // walk forward through the hull, adding more triangles and flipping recursively
        unsigned int n = del->hullNext[e];
        while (1) {
            q = del->hullNext[n];
            
            if (orient2d(x, y, del->coords[2 * n], del->coords[2 * n + 1], del->coords[2 * q], del->coords[2 * q + 1]) >= 0) {
                break;
            }
            
            // Add triangle and legalize
            t = addTriangle(n, i, q, del->hullTri[i], -1, del->hullTri[n], del);
            del->hullTri[i] = legalize(t + 2, del);
            // mark as removed
            del->hullNext[n] = n;  
            del->hullSize--;
            n = q;
        }

        // walk backward from the other side, adding more triangles and flipping
        if (e == start) {
            while(1){
                q = del->hullNext[e];
                if (orient2d(x, y, del->coords[2 * q], del->coords[2 * q + 1], del->coords[2 * e], del->coords[2 * e + 1]) >= 0) {
                    break;
                }

                t = addTriangle(q, i, e, -1, del->hullTri[e], del->hullTri[q], del);
                legalize(t + 2, del);
                del->hullTri[q] = t;
                // mark as removed
                del->hullNext[e] = e;
                del->hullSize--;
                e = q;
            }
        }

        // update the hull indices
        del->hullstart = del->hullPrev[i] = e;
        del->hullNext[e] = del->hullPrev[n] = i;
        del->hullNext[i] = n;

        // save the two new edges in the hash table
        del->hullHash[hash_key(x, y, del)] = i;
        del->hullHash[hash_key(coords[2 * e], coords[2 * e + 1], del)] = e;

    } // end for loop k

    del->hull = (unsigned*)malloc(hull_size * sizeof(unsigned int));
    unsigned int e = del->hullstart;
    for (int i = 0;  i < hull_size; i++) {
        del->hull[i] = e;
        e = del->hullNext[e];
    }

    // trim typed triangle mesh arrays
    int* newTriangles = (unsigned*)realloc(del->triangles, del->trianglesLen * sizeof(unsigned));
    if (newTriangles != NULL) {
        del->triangles = newTriangles;
    }
    int* newHalfedges = (unsigned*)realloc(del->halfedges, del->trianglesLen * sizeof(unsigned));
    if (newHalfedges != NULL) {
        del->halfedges = newHalfedges;
    }

}

static int legalize(int a, Delaunator* d){
    int i = 0;
    int ar = 0;

    // recursion eliminated with a fixed-size stack
    while(1){
        const int b = d->halfedges[a];
        /* if the pair of triangles doesn't satisfy the Delaunay condition
        * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
        * then do the same check/flip recursively for the new pair of triangles
        *
        *           pl                    pl
        *          /||\                  /  \
        *       al/ || \bl            al/    \a
        *        /  ||  \              /      \
        *       /  a||b  \    flip    /___ar___\
        *     p0\   ||   /p1   =>   p0\---bl---/p1
        *        \  ||  /              \      /
        *       ar\ || /br             b\    /br
        *          \||/                  \  /
        *           pr                    pr
        */

        const int a0 = a - a % 3;
        ar = a0 + (a + 2) % 3;

        // convex hull edge
        if (b == -1) { 
            if (i == 0) break;
            a = d->edge_stack[--i];
            continue;
        }
        const int b0 = b - b % 3;
        const int al = a0 + (a + 1) % 3;
        const int bl = b0 + (b + 2) % 3;

        const unsigned int p0 = d->triangles[ar];
        const unsigned int pr = d->triangles[a];
        const unsigned int pl = d->triangles[al];
        const unsigned int p1 = d->triangles[bl];

        const unsigned int illegal = inCircle(
            d->coords[2 * p0], d->coords[2 * p0 + 1],
            d->coords[2 * pr], d->coords[2 * pr + 1],
            d->coords[2 * pl], d->coords[2 * pl + 1],
            d->coords[2 * p1], d->coords[2 * p1 + 1]);

        if (illegal) {
            d->triangles[a] = p1;
            d->triangles[b] = p0;
            int hbl = d->halfedges[bl];

             // edge swapped on the other side of the hull (rare); fix the halfedge reference
            if (hbl == -1) {
                int e = d->hullstart;
                do {
                    if (d->hullTri[e] == bl) {
                        d->hullTri[e] = a;
                        break;
                    }
                    e = d->hullPrev[e];
                } while (e != d->hullstart);
            }
            link(a, hbl, d);
            link(b, d->halfedges[ar], d);
            link(ar, bl, d);

            const unsigned int br = b0 + (b + 1) % 3;

            if(i < EDGE_STACK_SIZE){
                d->edge_stack[i++] = br;
            }


        } 
        else {
            if (i == 0) break;
            a = d->edge_stack[--i];
        }

    }// while

    return ar;
}


/**
 * 
 */
static void link(int a, int b, Delaunator* d){
    d->halfedges[a]= b;
    if (b != -1) d->halfedges[b] = a;
}

/**
 * 
 */
static int add_triangle(int i0, int i1, int i2, int a, int b, int c, Delaunator* d){
    const int t = d->trianglesLen;

    d->triangles[t] = i0;
    d->triangles[t + 1] = i1;
    d->triangles[t + 2] = i2;

    link(t, a, d);
    link(t + 1, b, d);
    link(t + 2, c, d);
    
    d->trianglesLen += 3;
    
    return t; 
}


/**
 * 
 */
static int hash_key(double x, double y, Delaunator* d){
   return (int) floor(pseudo_angle(x - d->cx, y - d->cy) * d->hashSize) % d->hashSize;
}

/**
 * 
 */
static double orient2d(double pa0, double pa1, double pb0, 
                double pb1, double pc0, double pc1) {
    double detleft, detright, det;
    double detsum, errbound;

    double epsilon = EPSILON;
    double ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon;

    detleft = (pa0 - pc0) * (pb1 - pc1);
    detright = (pa1 - pc1) * (pb0 - pc0);
    det = detleft - detright;

    if (detleft > 0.0) {
    if (detright <= 0.0) {
        return det;
    } else {
        detsum = detleft + detright;
    }
    } else if (detleft < 0.0) {
    if (detright >= 0.0) {
        return det;
    } else {
        detsum = -detleft - detright;
    }
    } else {
    return det;
    }

    errbound = ccwerrboundA * detsum;
    if ((det >= errbound) || (-det >= errbound)) {
    return det;
    }

    return orient2dadapt(pa0, pa1, pb0, pb1, pc0, pc1, detsum);
}

/**
 * 
 */
static double orient2adapt(double pa0, double pa1, double pb0, 
                    double pb1, double pc0, double pc1, double detsum){
    double acx, acy, bcx, bcy;
    double acxtail, acytail, bcxtail, bcytail;
    double detleft, detright;
    double detlefttail, detrighttail;
    double det, errbound;
    double B[4], C1[8], C2[12], D[16];
    double B3;
    int C1length, C2length, Dlength;
    double u[4];
    double u3;
    double s1, t1;
    double s0, t0;

    double bvirt;
    double avirt, bround, around;
    double c;
    double abig;
    double ahi, alo, bhi, blo;
    double err1, err2, err3;
    double _i, _j;
    double _0;
    double epsilon = EPSILON;
    double ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon;
    double ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon;
    double resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;

    acx = (double) (pa0 - pc0);
    bcx = (double) (pb0 - pc0);
    acy = (double) (pa1 - pc1);
    bcy = (double) (pb1 - pc1);

    Two_Product(acx, bcy, detleft, detlefttail);
    Two_Product(acy, bcx, detright, detrighttail);

    Two_Two_Diff(detleft, detlefttail, detright, detrighttail,
                B3, B[2], B[1], B[0]);
    B[3] = B3;

    det = estimate(4, B);
    errbound = ccwerrboundB * detsum;
    if ((det >= errbound) || (-det >= errbound)) {
    return det;
    }

    Two_Diff_Tail(pa0, pc0, acx, acxtail);
    Two_Diff_Tail(pb0, pc0, bcx, bcxtail);
    Two_Diff_Tail(pa1, pc1, acy, acytail);
    Two_Diff_Tail(pb1, pc1, bcy, bcytail);

    if ((acxtail == 0.0) && (acytail == 0.0)
        && (bcxtail == 0.0) && (bcytail == 0.0)) {
    return det;
    }

    errbound = ccwerrboundC * detsum + resulterrbound * Absolute(det);
    det += (acx * bcytail + bcy * acxtail)
        - (acy * bcxtail + bcx * acytail);
    if ((det >= errbound) || (-det >= errbound)) {
    return det;
    }

    Two_Product(acxtail, bcy, s1, s0);
    Two_Product(acytail, bcx, t1, t0);
    Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3] = u3;
    C1length = fast_expansion_sum_zeroelim(4, B, 4, u, C1);

    Two_Product(acx, bcytail, s1, s0);
    Two_Product(acy, bcxtail, t1, t0);
    Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3] = u3;
    C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2);

    Two_Product(acxtail, bcytail, s1, s0);
    Two_Product(acytail, bcxtail, t1, t0);
    Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3] = u3;
    Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, D);

    return(D[Dlength - 1]);
}

/**
 * Calculates if a point (px,py) is inside a circle defined by three points
 * {*}x, {*}y.
 */
static int in_circle(const double ax, const double ay, 
                    const double bx, const double by, 
                    const double cx, const double cy, 
                    const double px, const double py) {
    const double dx = ax - px;
    const double dy = ay - py;
    const double ex = bx - px;
    const double ey = by - py;
    const double fx = cx - px;
    const double fy = cy - py;

    const double ap = dx * dx + dy * dy;
    const double bp = ex * ex + ey * ey;
    const double cp = fx * fx + fy * fy;

    return (dx * (ey * cp - bp * fy) -
            dy * (ex * cp - bp * fx) +
            ap * (ex * fy - ey * fx)) < 0;
}
/**
 * 
 */
static void quicksort(int* ids, double* dists, int left, int right) {
    if (right - left <= 20) {
        for (int i = left + 1; i <= right; i++) {
            int temp = ids[i];
            double temp_dist = dists[temp];
            int j = i - 1;
            while (j >= left && dists[ids[j]] > temp_dist){
                ids[j + 1] = ids[j--];
            }
            ids[j + 1] = temp;
        }
    } else {
        const int median = (left + right) >> 1;
        int i = left + 1;
        int j = right;
        swap(ids, median, i);
        if (dists[ids[left]] > dists[ids[right]]) swap(ids, left, right);
        if (dists[ids[i]] > dists[ids[right]]) swap(ids, i, right);
        if (dists[ids[left]] > dists[ids[i]]) swap(ids, left, i);

        int temp = ids[i];
        double temp_dist = dists[temp];
        while (1) {
            do {i++;
            } while (dists[ids[i]] < temp_dist);
            do {j--;
            } while (dists[ids[j]] > temp_dist);
            if (j < i) break;
            swap(ids, i, j);
        }
        ids[left + 1] = ids[j];
        ids[j] = temp;

        if (right - i + 1 >= j - left) {
            quicksort(ids, dists, i, right);
            quicksort(ids, dists, left, j - 1);
        } else {
            quicksort(ids, dists, left, j - 1);
            quicksort(ids, dists, i, right);
        }
    }
}

/**
 * 
 */
static double circumradius(double ax, double ay, double bx, double by, double cx, double cy) {
    double dx = bx - ax;
    double dy = by - ay;
    double ex = cx - ax;
    double ey = cy - ay;

    double bl = dx * dx + dy * dy;
    double cl = ex * ex + ey * ey;
    double d = 0.5 / (dx * ey - dy * ex);

    double x = (ey * bl - dy * cl) * d;
    double y = (dx * cl - ex * bl) * d;

    return x * x + y * y;
}

/**
 * 
 */
static double* circumcenter(ax, ay, bx, by, cx, cy) {
    double dx = bx - ax;
    double dy = by - ay;
    double ex = cx - ax;
    double ey = cy - ay;

    double bl = dx * dx + dy * dy;
    double cl = ex * ex + ey * ey;
    double d = 0.5 / (dx * ey - dy * ex);

    static double xy[2] = {-9999.0, -9999.0};

    xy[0] = ax + (ey * bl - dy * cl) * d;
    xy[1] = ay + (dx * cl - ex * bl) * d;

    return xy;
}

/**
 * 
 */
static double dist(double ax, double ay, double bx, double by) {
    double dx = ax - bx;
    double dy = ay - by;
    return dx * dx + dy * dy;
}

/**
 * 
 */
static double pseudo_angle(double dx, double dy) {
    double p = dx / (fabs(dx) + fabs(dy));
    return (dy > 0 ? 3 - p : 1 + p) / 4;
}

/**
 * 
 */
static void swap(double *arr, int i, int j) {
    const double tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

/**
 * 
 */
static void delaunator_destroy(Delaunator* d) {
    free(d->triangles);
    free(d->halfedges);
    free(d->hullPrev);
    free(d->hullNext);
    free(d->hullTri);
    free(d->hullHash);
    free(d->ids);
    free(d->dists);
    free(d);
}