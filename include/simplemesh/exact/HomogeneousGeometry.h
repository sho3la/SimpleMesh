// ============================================================================
//  simplemesh/exact/HomogeneousGeometry.h - exact homogeneous points & predicates
// ----------------------------------------------------------------------------
//  Exact homogeneous points and the projective predicates the surface
//  intersection uses. An intersection point (e.g. where a segment meets a
//  triangle) has rational coordinates; storing them as a numerator/denominator
//  pair AVOIDS DIVISION and keeps everything exact. A homogeneous point
//  (x,y,w) denotes the affine point (x/w, y/w); (x,y,z,w) denotes (x/w,y/w,z/w).
//
//  Predicates on homogeneous points reduce to integer-style determinants of the
//  homogeneous coordinates, corrected by the signs of the w (denominator) terms:
//
//    orient2d(p0,p1,p2) = sign det3x3(xi,yi,wi) * sign(w0)*sign(w1)*sign(w2)
//    orient3d(p0,p1,p2,p3): U=p1-p0, V=p2-p0, W=p3-p0 (homogeneous differences),
//                           = sign det3x3(U,V,W) * sign(U.w)*sign(V.w)*sign(W.w)
//
//  Each is filtered: evaluate the determinant in Interval; if its sign is
//  determined, combine with the exact w-signs; otherwise fall back to the exact
//  ExpansionNt determinant. The determinant and the homogeneous subtraction are
//  each written ONCE as templates over the number type, so filter and exact paths
//  agree by construction.
//
//  Compile strict-FP (-ffp-contract=off / /fp:strict).
// ============================================================================
#pragma once

#include "Predicates.h"   // pulls in Interval + Expansion
#include "Interval.h"
#include "Expansion.h"

namespace sm {
namespace exact {

// ---- generic 2x2 / 3x3 determinants ----------------------------------------
template <class T> inline T det2x2(const T& a11, const T& a12,
                                   const T& a21, const T& a22) {
    return a11 * a22 - a12 * a21;
}
template <class T> inline T det3x3(
    const T& a11, const T& a12, const T& a13,
    const T& a21, const T& a22, const T& a23,
    const T& a31, const T& a32, const T& a33
) {
    return a11 * det2x2(a22, a23, a32, a33)
         - a21 * det2x2(a12, a13, a32, a33)
         + a31 * det2x2(a12, a13, a22, a23);
}

// ---- homogeneous point types ----------------------------------------------
template <class T> struct Vec2h {
    T x, y, w;
    Vec2h() = default;
    Vec2h(T x_, T y_, T w_) : x(std::move(x_)), y(std::move(y_)), w(std::move(w_)) {}
};
template <class T> struct Vec3h {
    T x, y, z, w;
    Vec3h() = default;
    Vec3h(T x_, T y_, T z_, T w_)
        : x(std::move(x_)), y(std::move(y_)), z(std::move(z_)), w(std::move(w_)) {}
};

// homogeneous difference p1 - p2
template <class T> inline Vec2h<T> operator-(const Vec2h<T>& p1, const Vec2h<T>& p2) {
    return Vec2h<T>(det2x2(p1.x, p1.w, p2.x, p2.w),
                    det2x2(p1.y, p1.w, p2.y, p2.w),
                    p1.w * p2.w);
}
template <class T> inline Vec3h<T> operator-(const Vec3h<T>& p1, const Vec3h<T>& p2) {
    return Vec3h<T>(det2x2(p1.x, p1.w, p2.x, p2.w),
                    det2x2(p1.y, p1.w, p2.y, p2.w),
                    det2x2(p1.z, p1.w, p2.z, p2.w),
                    p1.w * p2.w);
}

// Exact / interval / filter aliases.
using Vec2HE = Vec2h<ExpansionNt>;
using Vec3HE = Vec3h<ExpansionNt>;
using Vec2HI = Vec2h<Interval>;
using Vec3HI = Vec3h<Interval>;

// Build a homogeneous point from an affine (Cartesian) double point (w = 1).
inline Vec2HE make_vec2he(double x, double y) { return Vec2HE(x, y, 1.0); }
inline Vec3HE make_vec3he(double x, double y, double z) { return Vec3HE(x, y, z, 1.0); }

// Interval views of an exact homogeneous point (for the filter).
inline Vec2HI to_interval(const Vec2HE& p) {
    return Vec2HI(Interval(p.x.rep()), Interval(p.y.rep()), Interval(p.w.rep()));
}
inline Vec3HI to_interval(const Vec3HE& p) {
    return Vec3HI(Interval(p.x.rep()), Interval(p.y.rep()),
                  Interval(p.z.rep()), Interval(p.w.rep()));
}

// ---- predicates on homogeneous points --------------------------------------

inline int orient2d(const Vec2HE& p0, const Vec2HE& p1, const Vec2HE& p2) {
    int sw = p0.w.sign() * p1.w.sign() * p2.w.sign();  // exact w-signs
    {   // filter
        Vec2HI a = to_interval(p0), b = to_interval(p1), c = to_interval(p2);
        Interval Delta = det3x3<Interval>(a.x, a.y, a.w, b.x, b.y, b.w, c.x, c.y, c.w);
        Sign2 s = Delta.sign();
        if (sign_is_determined(s)) return convert_sign(s) * sw;
    }
    ExpansionNt Delta = det3x3<ExpansionNt>(p0.x, p0.y, p0.w,
                                            p1.x, p1.y, p1.w,
                                            p2.x, p2.y, p2.w);
    return Delta.sign() * sw;
}

// Affine squared length of a homogeneous point: (x/w)^2 + (y/w)^2, as a double.
// These "given lengths" are passed to incircle to avoid expansion overflow.
inline double vec2h_length(const Vec2HE& p) {
    double x = p.x.estimate(), y = p.y.estimate(), w = p.w.estimate();
    return (x*x + y*y) / (w*w);
}

// Homogeneous in-circle (with given lengths, WITHOUT the SOS
// tie-break): >0 iff p3 is strictly inside the circumcircle of (p0,p1,p2) for a
// CCW triangle. li = given squared length of pi. On a cocircular tie returns 0
// (a valid - if not canonical - Delaunay decision; SOS uniqueness deferred).
template <class T>
inline T incircle_det_h(const Vec2h<T>& P1, const Vec2h<T>& P2, const Vec2h<T>& P3,
                        const T& L1, const T& L2, const T& L3) {
    T M1 = det2x2(P2.x, P2.y, P3.x, P3.y);
    T M2 = det2x2(P1.x, P1.y, P3.x, P3.y);
    T M3 = det2x2(P1.x, P1.y, P2.x, P2.y);
    return L1*P1.w*M1 - L2*P2.w*M2 + L3*P3.w*M3;
}

inline int incircle_2d(const Vec2HE& p0, const Vec2HE& p1,
                       const Vec2HE& p2, const Vec2HE& p3,
                       double l0, double l1, double l2, double l3) {
    {   // filter
        Vec2HI P3i = to_interval(p3);
        Vec2HI P1 = to_interval(p0) - P3i;
        Vec2HI P2 = to_interval(p1) - P3i;
        Vec2HI P3 = to_interval(p2) - P3i;
        Sign2 s1 = P1.w.sign(), s2 = P2.w.sign(), s3 = P3.w.sign();
        if (sign_is_non_zero(s1) && sign_is_non_zero(s2) && sign_is_non_zero(s3)) {
            Interval D = incircle_det_h<Interval>(
                P1, P2, P3, Interval(l0)-Interval(l3),
                Interval(l1)-Interval(l3), Interval(l2)-Interval(l3));
            Sign2 s = D.sign();
            if (sign_is_non_zero(s))
                return convert_sign(s)*convert_sign(s1)*convert_sign(s2)*convert_sign(s3);
        }
    }
    Vec2HE P1 = p0 - p3, P2 = p1 - p3, P3 = p2 - p3;
    ExpansionNt L1 = ExpansionNt(l0)-ExpansionNt(l3);
    ExpansionNt L2 = ExpansionNt(l1)-ExpansionNt(l3);
    ExpansionNt L3 = ExpansionNt(l2)-ExpansionNt(l3);
    ExpansionNt D = incircle_det_h<ExpansionNt>(P1, P2, P3, L1, L2, L3);
    return D.sign() * P1.w.sign() * P2.w.sign() * P3.w.sign();
}

// Homogeneous mix: the point at rational parameter t = tnum/tden along p1->p2
// with NO division. General (differing-w) form.
inline Vec2HE mix(const ExpansionNt& tnum, const ExpansionNt& tden,
                  const Vec2HE& p1, const Vec2HE& p2) {
    ExpansionNt sn = p2.w * (tden - tnum);
    ExpansionNt tn = p1.w * tnum;
    return Vec2HE(sn*p1.x + tn*p2.x, sn*p1.y + tn*p2.y, tden*p1.w*p2.w);
}

inline int orient3d(const Vec3HE& p0, const Vec3HE& p1,
                    const Vec3HE& p2, const Vec3HE& p3) {
    {   // filter: homogeneous differences and det3x3 in interval arithmetic
        Vec3HI P0 = to_interval(p0);
        Vec3HI U = to_interval(p1) - P0;
        Vec3HI V = to_interval(p2) - P0;
        Vec3HI W = to_interval(p3) - P0;
        Sign2 su = U.w.sign(), sv = V.w.sign(), sw = W.w.sign();
        if (sign_is_non_zero(su) && sign_is_non_zero(sv) && sign_is_non_zero(sw)) {
            Interval Delta = det3x3<Interval>(U.x, U.y, U.z, V.x, V.y, V.z, W.x, W.y, W.z);
            Sign2 s = Delta.sign();
            if (sign_is_non_zero(s)) {
                return convert_sign(s) * convert_sign(su) *
                       convert_sign(sv) * convert_sign(sw);
            }
        }
    }
    Vec3HE U = p1 - p0, V = p2 - p0, W = p3 - p0;
    ExpansionNt Delta = det3x3<ExpansionNt>(U.x, U.y, U.z, V.x, V.y, V.z, W.x, W.y, W.z);
    return Delta.sign() * U.w.sign() * V.w.sign() * W.w.sign();
}

} // namespace exact
} // namespace sm
