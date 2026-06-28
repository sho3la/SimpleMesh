// ============================================================================
//  test_cdt.cpp - test for the CDT combinatorial core (point insertion).
// ----------------------------------------------------------------------------
//  The spec for an *unconstrained* Delaunay triangulation is self-validating, so
//  this needs no external golden file: after inserting points into a big
//  enclosing triangle, the result must be
//    (1) a structurally valid triangulation: every triangle consistently
//        oriented, adjacency symmetric, and the shared edge of adjacent
//        triangles uses the same two vertices;
//    (2) Delaunay: for every interior edge whose four involved vertices are all
//        real (not the 3 enclosing super-vertices), the vertex opposite across
//        the edge is NOT strictly inside the circumcircle of the triangle
//        (checked with the Day-3 EXACT in_circle).
//  Both are checked exactly. A grid input also exercises cocircular degeneracies.
// ============================================================================
#include "simplemesh/exact/CDT2d.h"
#include "simplemesh/exact/Predicates.h"

#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace sm::exact;
using index_t = CDT2d::index_t;

static int failures = 0;
#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } \
        else         { std::cout << "ok:   " << msg << "\n"; }             \
    } while (0)

// Structural validity + exact Delaunay property. Returns counts via refs.
static void audit(const CDT2d& T, int& bad_orient, int& bad_adj, int& non_delaunay) {
    bad_orient = bad_adj = non_delaunay = 0;
    Sign g = T.global_orientation();
    for (index_t t = 0; t < T.nT(); ++t) {
        index_t a = T.tri_v(t,0), b = T.tri_v(t,1), c = T.tri_v(t,2);
        // (1a) consistent orientation
        if (Sign(orient2d(T.px(a),T.py(a),T.px(b),T.py(b),T.px(c),T.py(c))) != g)
            ++bad_orient;
        // (1b) + (2) per edge
        for (index_t le = 0; le < 3; ++le) {
            index_t t2 = T.tri_adj(t, le);
            if (t2 == CDT_NO_INDEX) continue;
            // adjacency symmetry: t2 must point back to t
            index_t back = CDT_NO_INDEX;
            for (index_t e2 = 0; e2 < 3; ++e2) if (T.tri_adj(t2,e2)==t) back = e2;
            if (back == CDT_NO_INDEX) { ++bad_adj; continue; }
            if (t2 < t) continue;  // visit each edge once
            // shared edge vertices must match (t's edge le is opposite vertex le)
            index_t u = T.tri_v(t,(le+1)%3), w = T.tri_v(t,(le+2)%3);
            index_t u2 = T.tri_v(t2,(back+1)%3), w2 = T.tri_v(t2,(back+2)%3);
            if (!((u==w2 && w==u2))) { ++bad_adj; continue; }
            // Delaunay: opposite vertex of t2 not inside circumcircle(a,b,c).
            index_t d = T.tri_v(t2, back);
            index_t aa=T.tri_v(t,0), bb=T.tri_v(t,1), cc=T.tri_v(t,2);
            // skip edges touching the enclosing super-triangle (verts 0,1,2)
            if (aa<3||bb<3||cc<3||d<3) continue;
            int ic = in_circle(T.px(aa),T.py(aa),T.px(bb),T.py(bb),
                               T.px(cc),T.py(cc),T.px(d),T.py(d));
            // triangles are CCW (g==POSITIVE) -> inside == +1 is a violation
            if (ic * g > 0) ++non_delaunay;
        }
    }
}

static CDT2d make_enclosing() {
    CDT2d T;
    // CCW, large enough to contain all test points (|coord| <= ~1e5).
    T.create_enclosing_triangle(-1e7,-1e7,  1e7,-1e7,  0.0,1e7);
    return T;
}

int main() {
    // ===== random points =====================================================
    {
        std::mt19937_64 rng(20260628);
        std::uniform_real_distribution<double> U(-1e5, 1e5);
        CDT2d T = make_enclosing();
        const int N = 3000;
        for (int i = 0; i < N; ++i) T.insert_point(U(rng), U(rng));
        int bo, ba, nd; audit(T, bo, ba, nd);
        std::cout << "      random: " << T.nv() << " verts, " << T.nT()
                  << " triangles\n";
        CHECK(bo == 0, "random: all triangles consistently oriented");
        CHECK(ba == 0, "random: adjacency symmetric & edges consistent");
        CHECK(nd == 0, "random: Delaunay empty-circumcircle holds (exact)");
    }

    // ===== integer grid (heavy cocircular degeneracy) ========================
    {
        CDT2d T = make_enclosing();
        for (int x = 0; x < 25; ++x)
            for (int y = 0; y < 25; ++y)
                T.insert_point(double(x) * 100.0, double(y) * 100.0);
        int bo, ba, nd; audit(T, bo, ba, nd);
        std::cout << "      grid:   " << T.nv() << " verts, " << T.nT()
                  << " triangles\n";
        CHECK(bo == 0, "grid: all triangles consistently oriented");
        CHECK(ba == 0, "grid: adjacency symmetric & edges consistent");
        CHECK(nd == 0, "grid: Delaunay holds even with cocircular points (exact)");
    }

    // ===== duplicate points are folded =======================================
    {
        CDT2d T = make_enclosing();
        index_t a = T.insert_point(10.0, 20.0);
        index_t b = T.insert_point(10.0, 20.0);  // same point
        CHECK(a == b, "duplicate point returns the existing vertex");
        index_t base_nv = T.nv();
        T.insert_point(10.0, 20.0);
        CHECK(T.nv() == base_nv, "re-inserting duplicate does not grow nv");
    }

    if (failures) {
        std::cerr << "\n" << failures << " CDT check(s) FAILED\n";
        return 1;
    }
    std::cout << "\nOK: Delaunay point insertion is valid & exact.\n";
    return 0;
}
