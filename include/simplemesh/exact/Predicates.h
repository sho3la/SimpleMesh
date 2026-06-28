// ============================================================================
//  simplemesh/exact/Predicates.h - exact geometric predicates (filtered)
// ----------------------------------------------------------------------------
//  orient2d / orient3d / in_circle are sign-of-determinant predicates. The
//  determinant is a geometric invariant: any *exact* evaluation of the same
//  formula yields the same sign. Each predicate is written ONCE as a templated
//  determinant formula and evaluated:
//
//      1. fast path: in interval arithmetic (Interval). If the sign is
//         DETERMINED, return it. (the filter never reports a wrong sign.)
//      2. exact path: re-evaluate the SAME formula in exact arithmetic
//         (ExpansionNt) and return its sign.
//
//  Writing the formula once (a function template over the number type) is what
//  guarantees the filtered and exact paths agree by construction.
//
//  Sign convention (standard robust-predicate convention):
//   * orient2d(a,b,c) > 0  <=>  c is left of a->b  (a,b,c CCW)
//   * orient3d(a,b,c,d) > 0 <=> d is below the plane (a,b,c) seen CCW, i.e. the
//                               tetrahedron (a,b,c,d) has positive signed volume
//                               with the convention det(b-a,c-a,d-a).
//   * in_circle(a,b,c,d) > 0 <=> d is strictly INSIDE the circumcircle of the
//                                CCW triangle (a,b,c).
//
//  No symbolic perturbation (SOS) here: ZERO is meaningful (collinear / coplanar
//  / cocircular), which is exactly what the surface-intersection exact path
//  needs. SOS can be layered on later if a call site requires it.
//
//  Compile strict-FP (-ffp-contract=off / /fp:strict) - the interval filter and
//  the expansion transforms both depend on it.
// ============================================================================
#pragma once

#include "Interval.h"
#include "Expansion.h"

namespace sm {
namespace exact {

// --- determinant formulas, written once, evaluated in T = Interval or ExpansionNt

template <class T>
inline T orient2d_det(T ax, T ay, T bx, T by, T cx, T cy) {
    // | bx-ax  by-ay |
    // | cx-ax  cy-ay |
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

template <class T>
inline T orient3d_det(T ax, T ay, T az, T bx, T by, T bz,
                      T cx, T cy, T cz, T dx, T dy, T dz) {
    // det( b-a, c-a, d-a )
    T bax = bx - ax, bay = by - ay, baz = bz - az;
    T cax = cx - ax, cay = cy - ay, caz = cz - az;
    T dax = dx - ax, day = dy - ay, daz = dz - az;
    T m1 = cay * daz - caz * day;
    T m2 = cax * daz - caz * dax;
    T m3 = cax * day - cay * dax;
    return bax * m1 - bay * m2 + baz * m3;
}

template <class T>
inline T in_circle_det(T ax, T ay, T bx, T by, T cx, T cy, T dx, T dy) {
    // | ax-dx  ay-dy  (ax-dx)^2+(ay-dy)^2 |
    // | bx-dx  by-dy  (bx-dx)^2+(by-dy)^2 |
    // | cx-dx  cy-dy  (cx-dx)^2+(cy-dy)^2 |
    T adx = ax - dx, ady = ay - dy;
    T bdx = bx - dx, bdy = by - dy;
    T cdx = cx - dx, cdy = cy - dy;
    T alift = adx * adx + ady * ady;
    T blift = bdx * bdx + bdy * bdy;
    T clift = cdx * cdx + cdy * cdy;
    // Expand along the lift column.
    T bc = bdx * cdy - cdx * bdy;
    T ca = cdx * ady - adx * cdy;
    T ab = adx * bdy - bdx * ady;
    return alift * bc + blift * ca + clift * ab;
}

// --- filtered drivers -------------------------------------------------------

inline int orient2d(double ax, double ay, double bx, double by,
                    double cx, double cy) {
    Interval d = orient2d_det<Interval>(ax, ay, bx, by, cx, cy);
    Sign2 s = d.sign();
    if (sign_is_determined(s)) return convert_sign(s);
    ExpansionNt e = orient2d_det<ExpansionNt>(ax, ay, bx, by, cx, cy);
    return e.sign();
}

inline int orient3d(double ax, double ay, double az,
                    double bx, double by, double bz,
                    double cx, double cy, double cz,
                    double dx, double dy, double dz) {
    Interval d = orient3d_det<Interval>(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz);
    Sign2 s = d.sign();
    if (sign_is_determined(s)) return convert_sign(s);
    ExpansionNt e = orient3d_det<ExpansionNt>(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz);
    return e.sign();
}

inline int in_circle(double ax, double ay, double bx, double by,
                     double cx, double cy, double dx, double dy) {
    Interval d = in_circle_det<Interval>(ax, ay, bx, by, cx, cy, dx, dy);
    Sign2 s = d.sign();
    if (sign_is_determined(s)) return convert_sign(s);
    ExpansionNt e = in_circle_det<ExpansionNt>(ax, ay, bx, by, cx, cy, dx, dy);
    return e.sign();
}

// Convenience pointer overloads (data() interfaces).
inline int orient2d(const double* a, const double* b, const double* c) {
    return orient2d(a[0], a[1], b[0], b[1], c[0], c[1]);
}
inline int orient3d(const double* a, const double* b, const double* c, const double* d) {
    return orient3d(a[0], a[1], a[2], b[0], b[1], b[2],
                    c[0], c[1], c[2], d[0], d[1], d[2]);
}
inline int in_circle(const double* a, const double* b, const double* c, const double* d) {
    return in_circle(a[0], a[1], b[0], b[1], c[0], c[1], d[0], d[1]);
}

// ---- predicates the symbolic triangle intersection needs -------------------

// sign( (b-a) . (c-a) ) in 3D. >0 means the angle at a is acute.
template <class T>
inline T dot3d_det(T ax, T ay, T az, T bx, T by, T bz, T cx, T cy, T cz) {
    return (bx - ax) * (cx - ax) + (by - ay) * (cy - ay) + (bz - az) * (cz - az);
}
inline int dot_3d(const double* a, const double* b, const double* c) {
    Interval d = dot3d_det<Interval>(a[0],a[1],a[2], b[0],b[1],b[2], c[0],c[1],c[2]);
    Sign2 s = d.sign();
    if (sign_is_determined(s)) return convert_sign(s);
    ExpansionNt e = dot3d_det<ExpansionNt>(a[0],a[1],a[2], b[0],b[1],b[2], c[0],c[1],c[2]);
    return e.sign();
}

// True iff p1,p2,p3 are collinear (degenerate triangle). Exact test:
// the plane through p1,p2,p3 is undefined iff it "contains" 4 independent
// reference points (orient3d == 0 against each).
inline bool aligned_3d(const double* p1, const double* p2, const double* p3) {
    static const double q000[3] = {0,0,0}, q001[3] = {0,0,1},
                        q010[3] = {0,1,0}, q100[3] = {1,0,0};
    return orient3d(p1,p2,p3,q000) == 0 && orient3d(p1,p2,p3,q001) == 0 &&
           orient3d(p1,p2,p3,q010) == 0 && orient3d(p1,p2,p3,q100) == 0;
}

// Axis (0/1/2) of the largest-magnitude component of the triangle normal
// N = (p2-p1) x (p3-p1); the coordinate to drop for a non-degenerate 2D
// projection. Comparison of |Ni| done exactly via squares.
inline int triangle_normal_axis(const double* p1, const double* p2, const double* p3) {
    ExpansionNt ux = ExpansionNt(p2[0])-ExpansionNt(p1[0]);
    ExpansionNt uy = ExpansionNt(p2[1])-ExpansionNt(p1[1]);
    ExpansionNt uz = ExpansionNt(p2[2])-ExpansionNt(p1[2]);
    ExpansionNt vx = ExpansionNt(p3[0])-ExpansionNt(p1[0]);
    ExpansionNt vy = ExpansionNt(p3[1])-ExpansionNt(p1[1]);
    ExpansionNt vz = ExpansionNt(p3[2])-ExpansionNt(p1[2]);
    ExpansionNt nx = uy*vz - uz*vy;
    ExpansionNt ny = uz*vx - ux*vz;
    ExpansionNt nz = ux*vy - uy*vx;
    ExpansionNt nx2 = nx*nx, ny2 = ny*ny, nz2 = nz*nz;
    // axis = argmax(nx2, ny2, nz2)
    if (compare(nx2, ny2) >= 0) {
        return (compare(nx2, nz2) >= 0) ? 0 : 2;
    } else {
        return (compare(ny2, nz2) >= 0) ? 1 : 2;
    }
}

} // namespace exact
} // namespace sm
