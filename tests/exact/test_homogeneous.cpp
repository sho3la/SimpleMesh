// ============================================================================
//  test_homogeneous.cpp - test for exact homogeneous predicates.
// ----------------------------------------------------------------------------
//  Validates orient2d/orient3d on homogeneous (numerator/denominator) points:
//    * vs an exact 128-bit integer projective-determinant oracle on integer x,y,(z),w  (zero
//      misses), with random POSITIVE and NEGATIVE denominators,
//    * cross-check: when all w==1 the homogeneous predicate equals the Day-3
//      Cartesian predicate (pins the sign convention),
//    * a real rational case: the exact intersection point of two segments,
//      represented with w != 1 and no division, gives the geometrically correct
//      orientation.
// ============================================================================
#include "simplemesh/exact/HomogeneousGeometry.h"
#include "Int128.h"
#include "simplemesh/exact/Predicates.h"

#include <cstdint>
#include <iostream>
#include <random>
#include <string>

using namespace sm::exact;

static int failures = 0;
#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } \
        else         { std::cout << "ok:   " << msg << "\n"; }             \
    } while (0)

static int isign(i128 v) { return (v > 0) - (v < 0); }

// oracle: orient2d on homogeneous integer points = sign det3x3 * prod sign(w).
static int oracle_orient2d_h(i128 x0,i128 y0,i128 w0, i128 x1,i128 y1,i128 w1,
                             i128 x2,i128 y2,i128 w2) {
    i128 det = x0*(y1*w2 - w1*y2) - x1*(y0*w2 - w0*y2) + x2*(y0*w1 - w0*y1);
    return isign(det) * isign(w0) * isign(w1) * isign(w2);
}
// oracle: orient3d, replicating the homogeneous-difference + det3x3 path.
static int oracle_orient3d_h(
    i128 x0,i128 y0,i128 z0,i128 w0, i128 x1,i128 y1,i128 z1,i128 w1,
    i128 x2,i128 y2,i128 z2,i128 w2, i128 x3,i128 y3,i128 z3,i128 w3) {
    // U = p1 - p0, etc.  (det2x2(a.x,a.w,b.x,b.w) = a.x*b.w - a.w*b.x)
    i128 Ux=x1*w0-w1*x0, Uy=y1*w0-w1*y0, Uz=z1*w0-w1*z0, Uw=w1*w0;
    i128 Vx=x2*w0-w2*x0, Vy=y2*w0-w2*y0, Vz=z2*w0-w2*z0, Vw=w2*w0;
    i128 Wx=x3*w0-w3*x0, Wy=y3*w0-w3*y0, Wz=z3*w0-w3*z0, Ww=w3*w0;
    i128 det = Ux*(Vy*Wz-Vz*Wy) - Vx*(Uy*Wz-Uz*Wy) + Wx*(Uy*Vz-Uz*Vy);
    return isign(det) * isign(Uw) * isign(Vw) * isign(Ww);
}

int main() {
    std::mt19937_64 rng(20260628);

    // ===== orient2d homogeneous vs oracle (with +/- denominators) ============
    {
        std::uniform_int_distribution<int64_t> P(-(1<<14), (1<<14));
        std::uniform_int_distribution<int64_t> Wd(1, (1<<12));
        std::uniform_int_distribution<int> Sg(0, 1);
        auto w = [&]{ int64_t v = Wd(rng); return Sg(rng) ? v : -v; };  // never 0
        int mism = 0;
        for (int i = 0; i < 200000; ++i) {
            int64_t x0=P(rng),y0=P(rng),w0=w(), x1=P(rng),y1=P(rng),w1=w(),
                    x2=P(rng),y2=P(rng),w2=w();
            Vec2HE p0{double(x0),double(y0),double(w0)};
            Vec2HE p1{double(x1),double(y1),double(w1)};
            Vec2HE p2{double(x2),double(y2),double(w2)};
            int got = orient2d(p0,p1,p2);
            int exp = oracle_orient2d_h(x0,y0,w0,x1,y1,w1,x2,y2,w2);
            if (got != exp) ++mism;
        }
        CHECK(mism == 0, "orient2d(homogeneous) matches oracle on 200k configs (+/- w)");
    }

    // ===== orient3d homogeneous vs oracle ====================================
    {
        std::uniform_int_distribution<int64_t> P(-(1<<8), (1<<8));
        std::uniform_int_distribution<int64_t> Wd(1, (1<<6));
        std::uniform_int_distribution<int> Sg(0, 1);
        auto w = [&]{ int64_t v = Wd(rng); return Sg(rng) ? v : -v; };
        int mism = 0;
        for (int i = 0; i < 200000; ++i) {
            int64_t x0=P(rng),y0=P(rng),z0=P(rng),w0=w();
            int64_t x1=P(rng),y1=P(rng),z1=P(rng),w1=w();
            int64_t x2=P(rng),y2=P(rng),z2=P(rng),w2=w();
            int64_t x3=P(rng),y3=P(rng),z3=P(rng),w3=w();
            Vec3HE p0{double(x0),double(y0),double(z0),double(w0)};
            Vec3HE p1{double(x1),double(y1),double(z1),double(w1)};
            Vec3HE p2{double(x2),double(y2),double(z2),double(w2)};
            Vec3HE p3{double(x3),double(y3),double(z3),double(w3)};
            int got = orient3d(p0,p1,p2,p3);
            int exp = oracle_orient3d_h(x0,y0,z0,w0,x1,y1,z1,w1,x2,y2,z2,w2,x3,y3,z3,w3);
            if (got != exp) ++mism;
        }
        CHECK(mism == 0, "orient3d(homogeneous) matches oracle on 200k configs (+/- w)");
    }

    // ===== convention cross-check: w==1 equals Day-3 Cartesian ===============
    {
        std::uniform_int_distribution<int64_t> P(-(1<<18), (1<<18));
        int mism2 = 0, mism3 = 0;
        for (int i = 0; i < 100000; ++i) {
            int64_t ax=P(rng),ay=P(rng),bx=P(rng),by=P(rng),cx=P(rng),cy=P(rng);
            int h = orient2d(make_vec2he(double(ax),double(ay)),
                             make_vec2he(double(bx),double(by)),
                             make_vec2he(double(cx),double(cy)));
            int cart = sm::exact::orient2d(double(ax),double(ay),double(bx),double(by),
                                           double(cx),double(cy));
            if (h != cart) ++mism2;

            int64_t az=P(rng),bz=P(rng),cz=P(rng),dx=P(rng),dy=P(rng),dz=P(rng);
            int h3 = orient3d(make_vec3he(double(ax),double(ay),double(az)),
                              make_vec3he(double(bx),double(by),double(bz)),
                              make_vec3he(double(cx),double(cy),double(cz)),
                              make_vec3he(double(dx),double(dy),double(dz)));
            int cart3 = sm::exact::orient3d(double(ax),double(ay),double(az),
                                            double(bx),double(by),double(bz),
                                            double(cx),double(cy),double(cz),
                                            double(dx),double(dy),double(dz));
            if (h3 != cart3) ++mism3;
        }
        CHECK(mism2 == 0, "orient2d(w=1) == Cartesian orient2d on 100k configs");
        CHECK(mism3 == 0, "orient3d(w=1) == Cartesian orient3d on 100k configs");
    }

    // ===== real rational intersection point (no division) ====================
    // Segment A: (0,0)->(4,4). Segment B: (4,0)->(0,4). They cross at (2,2).
    // The exact crossing of two integer segments is rational; represent it
    // homogeneously. Here it happens to be (2,2) but we build it as x=8,y=8,w=4
    // (= 2,2) WITHOUT dividing, to mimic how intersection points are stored.
    {
        Vec2HE X(8.0, 8.0, 4.0);          // == (2,2), denominator 4, undivided
        // Triangle (0,0),(4,0),(4,4) is CCW; (2,2) lies on its hypotenuse edge
        // (0,0)-(4,4): orient2d of that edge with X must be exactly 0.
        int on_edge = orient2d(make_vec2he(0,0), make_vec2he(4,4), X);
        CHECK(on_edge == 0, "rational point (8,8,4)=(2,2) is exactly on line (0,0)-(4,4)");
        // It is strictly left of edge (0,0)->(4,0):
        int side = orient2d(make_vec2he(0,0), make_vec2he(4,0), X);
        CHECK(side > 0, "rational point (2,2) is left of edge (0,0)->(4,0)");
    }

    if (failures) {
        std::cerr << "\n" << failures << " homogeneous check(s) FAILED\n";
        return 1;
    }
    std::cout << "\nOK: homogeneous predicates exact (no division).\n";
    return 0;
}
