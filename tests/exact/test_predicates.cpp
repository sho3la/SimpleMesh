// ============================================================================
//  test_predicates.cpp - golden test for exact geometric predicates.
// ----------------------------------------------------------------------------
//  THE critical numeric checkpoint. Every predicate sign is compared to an
//  exact 128-bit integer exact-determinant oracle on integer coordinates (so the oracle is
//  itself exact). We test:
//    * generic random configs               -> predicate == oracle  (zero misses)
//    * exactly-degenerate families           -> predicate == 0
//        collinear (orient2d), coplanar (orient3d), cocircular (in_circle)
//    * geometric sanity (known orientations / inside-outside)
//    * teeth: on degenerate families with large coords naive double is WRONG.
//
//  Zero tolerance: any single mismatch fails the build.
// ============================================================================
#include "simplemesh/exact/Predicates.h"
#include "Int128.h"

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

// ---- exact integer oracles -------------------------------------------------
static int oracle_orient2d(i128 ax, i128 ay, i128 bx, i128 by, i128 cx, i128 cy) {
    return isign((bx - ax) * (cy - ay) - (by - ay) * (cx - ax));
}
static int oracle_orient3d(i128 ax, i128 ay, i128 az, i128 bx, i128 by, i128 bz,
                           i128 cx, i128 cy, i128 cz, i128 dx, i128 dy, i128 dz) {
    i128 bax = bx-ax, bay = by-ay, baz = bz-az;
    i128 cax = cx-ax, cay = cy-ay, caz = cz-az;
    i128 dax = dx-ax, day = dy-ay, daz = dz-az;
    i128 det = bax * (cay*daz - caz*day)
             - bay * (cax*daz - caz*dax)
             + baz * (cax*day - cay*dax);
    return isign(det);
}
static int oracle_in_circle(i128 ax, i128 ay, i128 bx, i128 by,
                            i128 cx, i128 cy, i128 dx, i128 dy) {
    i128 adx=ax-dx, ady=ay-dy, bdx=bx-dx, bdy=by-dy, cdx=cx-dx, cdy=cy-dy;
    i128 al = adx*adx + ady*ady, bl = bdx*bdx + bdy*bdy, cl = cdx*cdx + cdy*cdy;
    i128 det = al*(bdx*cdy - cdx*bdy) + bl*(cdx*ady - adx*cdy) + cl*(adx*bdy - bdx*ady);
    return isign(det);
}

int main() {
    std::mt19937_64 rng(20260628);

    // ===== geometric sanity ==================================================
    CHECK(orient2d(0,0, 1,0, 0,1) > 0, "orient2d CCW triangle positive");
    CHECK(orient2d(0,0, 0,1, 1,0) < 0, "orient2d CW triangle negative");
    CHECK(orient2d(0,0, 1,1, 2,2) == 0, "orient2d collinear zero");
    CHECK(orient3d(0,0,0, 1,0,0, 0,1,0, 0,0,1) != 0, "orient3d non-flat tetra");
    CHECK(orient3d(0,0,0, 1,0,0, 0,1,0, 2,2,0) == 0, "orient3d coplanar zero");
    CHECK(in_circle(0,0, 4,0, 0,4, 1,1) > 0, "in_circle interior point positive");
    CHECK(in_circle(0,0, 4,0, 0,4, 10,10) < 0, "in_circle exterior point negative");
    CHECK(in_circle(5,0, 0,5, -5,0, 0,-5) == 0, "in_circle cocircular zero");

    // ===== orient2d: random vs oracle ========================================
    {
        std::uniform_int_distribution<int64_t> P(-(1<<20), (1<<20));
        int mism = 0;
        for (int i = 0; i < 200000; ++i) {
            int64_t ax=P(rng),ay=P(rng),bx=P(rng),by=P(rng),cx=P(rng),cy=P(rng);
            int got = orient2d(double(ax),double(ay),double(bx),double(by),double(cx),double(cy));
            int exp = oracle_orient2d(ax,ay,bx,by,cx,cy);
            if (got != exp) ++mism;
        }
        CHECK(mism == 0, "orient2d matches oracle on 200k random configs");
    }
    // orient2d: collinear family must be EXACTLY 0. Note: on *exactly* collinear
    // integer points naive double is also 0 (the two 2x2-determinant terms are
    // roundings of the identical exact product P, and round(P)-round(P)==0), so
    // there are no "teeth" to show here - double is correct. The meaningful fact
    // is that the interval filter cannot certify 0, so EVERY one of these goes
    // through the exact Expansion fallback and still returns 0. orient2d's truly
    // catastrophic cases are *nearly* collinear (fractional) configs, exercised
    // by the random-vs-oracle test above whenever the filter abstains.
    {
        std::uniform_int_distribution<int64_t> P(-(1<<26), (1<<26));
        std::uniform_int_distribution<int64_t> K(2, 9);
        int pred_nonzero = 0, trials = 200000;
        for (int i = 0; i < trials; ++i) {
            int64_t ax=P(rng),ay=P(rng),dx=P(rng),dy=P(rng);
            int64_t k = K(rng);
            int64_t bx=ax+dx, by=ay+dy;          // b = a + d
            int64_t cx=ax+k*dx, cy=ay+k*dy;      // c = a + k*d  (collinear)
            if (orient2d(double(ax),double(ay),double(bx),double(by),double(cx),double(cy)) != 0)
                ++pred_nonzero;
        }
        CHECK(pred_nonzero == 0,
              "orient2d exactly 0 on all collinear configs (all via exact fallback)");
    }

    // ===== orient3d: random vs oracle ========================================
    {
        std::uniform_int_distribution<int64_t> P(-(1<<15), (1<<15));
        int mism = 0;
        for (int i = 0; i < 200000; ++i) {
            int64_t ax=P(rng),ay=P(rng),az=P(rng), bx=P(rng),by=P(rng),bz=P(rng);
            int64_t cx=P(rng),cy=P(rng),cz=P(rng), dx=P(rng),dy=P(rng),dz=P(rng);
            int got = orient3d(double(ax),double(ay),double(az), double(bx),double(by),double(bz),
                               double(cx),double(cy),double(cz), double(dx),double(dy),double(dz));
            int exp = oracle_orient3d(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz);
            if (got != exp) ++mism;
        }
        CHECK(mism == 0, "orient3d matches oracle on 200k random configs");
    }
    // orient3d: coplanar family (exact 0) + teeth
    {
        std::uniform_int_distribution<int64_t> P(-(1<<24), (1<<24));
        std::uniform_int_distribution<int64_t> K(-4, 4);
        int pred_nonzero = 0, naive_nonzero = 0, trials = 200000;
        for (int i = 0; i < trials; ++i) {
            int64_t ax=P(rng),ay=P(rng),az=P(rng);
            int64_t ux=P(rng),uy=P(rng),uz=P(rng), vx=P(rng),vy=P(rng),vz=P(rng);
            int64_t p=K(rng), q=K(rng);
            int64_t bx=ax+ux,by=ay+uy,bz=az+uz;
            int64_t cx=ax+vx,cy=ay+vy,cz=az+vz;
            int64_t dx=ax+p*ux+q*vx, dy=ay+p*uy+q*vy, dz=az+p*uz+q*vz;  // d in plane(a,u,v)
            if (orient3d(double(ax),double(ay),double(az), double(bx),double(by),double(bz),
                         double(cx),double(cy),double(cz), double(dx),double(dy),double(dz)) != 0)
                ++pred_nonzero;
            // naive double 3x3
            double Bx=double(bx)-double(ax),By=double(by)-double(ay),Bz=double(bz)-double(az);
            double Cx=double(cx)-double(ax),Cy=double(cy)-double(ay),Cz=double(cz)-double(az);
            double Dx=double(dx)-double(ax),Dy=double(dy)-double(ay),Dz=double(dz)-double(az);
            double nd = Bx*(Cy*Dz-Cz*Dy) - By*(Cx*Dz-Cz*Dx) + Bz*(Cx*Dy-Cy*Dx);
            if (nd != 0.0) ++naive_nonzero;
        }
        CHECK(pred_nonzero == 0, "orient3d exactly 0 on all coplanar configs");
        CHECK(naive_nonzero > 0,
              "naive double orient3d WRONG on " + std::to_string(naive_nonzero) +
              "/" + std::to_string(trials) + " coplanar configs (teeth)");
    }

    // ===== in_circle: random vs oracle =======================================
    {
        std::uniform_int_distribution<int64_t> P(-(1<<12), (1<<12));
        int mism = 0;
        for (int i = 0; i < 200000; ++i) {
            int64_t ax=P(rng),ay=P(rng),bx=P(rng),by=P(rng),cx=P(rng),cy=P(rng),dx=P(rng),dy=P(rng);
            int got = in_circle(double(ax),double(ay),double(bx),double(by),
                                double(cx),double(cy),double(dx),double(dy));
            int exp = oracle_in_circle(ax,ay,bx,by,cx,cy,dx,dy);
            if (got != exp) ++mism;
        }
        CHECK(mism == 0, "in_circle matches oracle on 200k random configs");
    }
    // in_circle: cocircular family (Gaussian-integer circle r^2 = 5^2*t^2) -> 0
    {
        // Points on circle of radius 5*s centered at origin: (3,4),(4,3),(5,0),(0,5)...
        const int base[][2] = {{5,0},{4,3},{3,4},{0,5},{-3,4},{-4,3},{-5,0},{0,-5},{3,-4},{4,-3}};
        int pred_nonzero = 0, trials = 0;
        std::uniform_int_distribution<int> S(1, 1000);
        std::uniform_int_distribution<int> I(0, 9);
        for (int t = 0; t < 50000; ++t) {
            int s = S(rng);
            int i0=I(rng),i1=I(rng),i2=I(rng),i3=I(rng);
            if (i0==i1||i0==i2||i0==i3||i1==i2||i1==i3||i2==i3) continue;
            ++trials;
            auto X=[&](int i){return double(base[i][0]*s);};
            auto Y=[&](int i){return double(base[i][1]*s);};
            if (in_circle(X(i0),Y(i0),X(i1),Y(i1),X(i2),Y(i2),X(i3),Y(i3)) != 0)
                ++pred_nonzero;
        }
        CHECK(pred_nonzero == 0,
              "in_circle exactly 0 on " + std::to_string(trials) + " cocircular configs");
    }

    if (failures) {
        std::cerr << "\n" << failures << " predicate check(s) FAILED\n";
        return 1;
    }
    std::cout << "\nOK: exact predicates match the integer oracle (zero misses).\n";
    return 0;
}
