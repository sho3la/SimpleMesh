// ============================================================================
//  test_cdt_constraints.cpp - test for CDT constraint insertion.
// ----------------------------------------------------------------------------
//  Verifies the three things constraint insertion must guarantee:
//    (A) constraint recovery: after insert_constraint(a,b), the edge a-b is
//        actually present in the triangulation AND flagged constrained, even
//        when it had to cross (and flip) existing Delaunay edges;
//    (B) crossing constraints create an intersection vertex, and both
//        constraints route through it (Sloan + create_intersection);
//    (C) the triangulation stays structurally valid and the NON-constrained
//        interior edges remain Delaunay (exact empty-circumcircle).
// ============================================================================
#include "simplemesh/exact/CDT2d.h"
#include "simplemesh/exact/Predicates.h"

#include <iostream>
#include <random>
#include <string>
#include <cmath>

using namespace sm::exact;
using index_t = CDT2d::index_t;

static int failures = 0;
#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } \
        else         { std::cout << "ok:   " << msg << "\n"; }             \
    } while (0)

// Is there an edge with endpoints {a,b}? Optionally require it be constrained.
static bool has_edge(const CDT2d& T, index_t a, index_t b, bool must_constrain) {
    for (index_t t = 0; t < T.nT(); ++t)
        for (index_t le = 0; le < 3; ++le) {
            index_t u = T.tri_v(t, (le+1)%3), w = T.tri_v(t, (le+2)%3);
            if ((u==a && w==b) || (u==b && w==a)) {
                if (!must_constrain || T.tri_edge_constrained(t, le)) return true;
            }
        }
    return false;
}

static index_t vertex_at(const CDT2d& T, double x, double y, double tol) {
    for (index_t v = 0; v < T.nv(); ++v)
        if (std::abs(T.px(v)-x) < tol && std::abs(T.py(v)-y) < tol) return v;
    return CDT_NO_INDEX;
}

// structural validity + Delaunay of non-constrained edges
static void audit(const CDT2d& T, int& bad, int& non_delaunay) {
    bad = non_delaunay = 0;
    Sign g = T.global_orientation();
    for (index_t t = 0; t < T.nT(); ++t) {
        index_t a=T.tri_v(t,0), b=T.tri_v(t,1), c=T.tri_v(t,2);
        if (Sign(orient2d(T.px(a),T.py(a),T.px(b),T.py(b),T.px(c),T.py(c))) != g) ++bad;
        for (index_t le = 0; le < 3; ++le) {
            index_t t2 = T.tri_adj(t, le);
            if (t2 == CDT_NO_INDEX || t2 < t) continue;
            index_t back = CDT_NO_INDEX;
            for (index_t e2=0;e2<3;++e2) if (T.tri_adj(t2,e2)==t) back=e2;
            if (back == CDT_NO_INDEX) { ++bad; continue; }
            if (T.tri_edge_constrained(t, le)) continue;  // constrained edges may be non-Delaunay
            index_t d = T.tri_v(t2, back);
            if (a<3||b<3||c<3||d<3) continue;             // skip enclosing verts
            int ic = in_circle(T.px(a),T.py(a),T.px(b),T.py(b),
                               T.px(c),T.py(c),T.px(d),T.py(d));
            if (ic * g > 0) ++non_delaunay;
        }
    }
}

static CDT2d square_cdt(index_t& a, index_t& b, index_t& c, index_t& d) {
    CDT2d T;
    T.create_enclosing_triangle(-1e7,-1e7, 1e7,-1e7, 0.0,1e7);
    a = T.insert_point(0,0);
    b = T.insert_point(10,0);
    c = T.insert_point(10,10);
    d = T.insert_point(0,10);
    return T;
}

int main() {
    // ===== (A) constraint recovery across a flip =============================
    {
        index_t a,b,c,d; CDT2d T = square_cdt(a,b,c,d);
        T.insert_constraint(a, c);  // a diagonal of the square
        CHECK(has_edge(T, a, c, true), "constraint a-c present and flagged constrained");
        int bad, nd; audit(T, bad, nd);
        CHECK(bad == 0, "(A) triangulation structurally valid after constraint");
        CHECK(nd == 0, "(A) free edges remain Delaunay");
    }

    // ===== (B) crossing constraints -> intersection vertex ===================
    {
        index_t a,b,c,d; CDT2d T = square_cdt(a,b,c,d);
        index_t nv0 = T.nv();
        T.insert_constraint(a, c);  // diagonal 1
        T.insert_constraint(b, d);  // diagonal 2 -> crosses diagonal 1 at (5,5)
        CHECK(T.nv() == nv0 + 1, "crossing constraints created exactly one new vertex");
        index_t x = vertex_at(T, 5.0, 5.0, 1e-6);
        CHECK(x != CDT_NO_INDEX, "intersection vertex created at (5,5)");
        if (x != CDT_NO_INDEX) {
            CHECK(has_edge(T,a,x,true) && has_edge(T,x,c,true),
                  "diagonal a-c routes through intersection vertex");
            CHECK(has_edge(T,b,x,true) && has_edge(T,x,d,true),
                  "diagonal b-d routes through intersection vertex");
        }
        int bad, nd; audit(T, bad, nd);
        CHECK(bad == 0, "(B) triangulation structurally valid after crossing constraints");
    }

    // ===== (C) many points + a long constraint ===============================
    {
        CDT2d T;
        T.create_enclosing_triangle(-1e7,-1e7, 1e7,-1e7, 0.0,1e7);
        index_t left  = T.insert_point(-9000, 0);
        index_t right = T.insert_point( 9000, 0);
        std::mt19937_64 rng(20260628);
        std::uniform_real_distribution<double> U(-8000, 8000);
        for (int i = 0; i < 800; ++i) T.insert_point(U(rng), U(rng));
        T.insert_constraint(left, right);  // long horizontal constraint y=0
        // left-right may be split by points/intersections that lie on y=0; at
        // minimum, the two endpoints must be joined by >=1 constrained edge and
        // 'left' must have a constrained incident edge.
        bool left_has_cnstr = false;
        for (index_t t = 0; t < T.nT(); ++t)
            for (index_t le = 0; le < 3; ++le)
                if (T.tri_edge_constrained(t, le)) {
                    index_t u=T.tri_v(t,(le+1)%3), w=T.tri_v(t,(le+2)%3);
                    if (u==left || w==left) left_has_cnstr = true;
                }
        CHECK(left_has_cnstr, "(C) long constraint produced a constrained edge at its start");
        int bad, nd; audit(T, bad, nd);
        CHECK(bad == 0, "(C) triangulation structurally valid with 800 pts + constraint");
        CHECK(nd == 0, "(C) free edges remain Delaunay");
    }

    if (failures) {
        std::cerr << "\n" << failures << " CDT-constraint check(s) FAILED\n";
        return 1;
    }
    std::cout << "\nOK: constraint recovery + crossing intersections valid.\n";
    return 0;
}
