// ============================================================================
//  test_mesh_in_triangle.cpp - exact-coordinate CDT + per-facet remesh.
// ----------------------------------------------------------------------------
//  The decisive property of the exact path: when two constraints cross, the
//  created intersection vertex lies EXACTLY on both constraint lines (exact
//  homogeneous orient2d == 0), with NO division. A double create_intersection
//  rounds the crossing and fails this. We verify:
//    (A) ExactCDT2d: crossing constraints whose true crossing is a non-dyadic
//        rational -> the exact vertex is exactly on both lines; the rounded
//        double point is NOT; the triangulation is a valid exact CDT.
//    (B) MeshInTriangle: remesh a tilted 3D facet with two crossing constraints
//        -> structurally valid retriangulation, intersection vertex created.
// ============================================================================
#include "simplemesh/exact/MeshInTriangle.h"
#include "simplemesh/exact/HomogeneousGeometry.h"

#include <iostream>
#include <string>

using namespace sm::exact;
using index_t = ExactCDT2d::index_t;

static int failures = 0;
#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } \
        else         { std::cout << "ok:   " << msg << "\n"; }             \
    } while (0)

// structural validity (+ exact Delaunay for ExactCDT2d when check_delaunay).
template <class CDT>
static void audit(const CDT& T, int& bad, int& nd, bool check_delaunay) {
    bad = nd = 0;
    Sign g = T.global_orientation();
    for (index_t t = 0; t < T.nT(); ++t) {
        index_t a=T.tri_v(t,0), b=T.tri_v(t,1), c=T.tri_v(t,2);
        if (Sign(orient2d(T.point(a),T.point(b),T.point(c))) != g) ++bad;
        for (index_t le = 0; le < 3; ++le) {
            index_t t2 = T.tri_adj(t, le);
            if (t2 == CDT_NO_INDEX || t2 < t) continue;
            index_t back = CDT_NO_INDEX;
            for (index_t e2=0;e2<3;++e2) if (T.tri_adj(t2,e2)==t) back=e2;
            if (back == CDT_NO_INDEX) { ++bad; continue; }
            if (!check_delaunay || T.tri_edge_constrained(t, le)) continue;
            index_t d = T.tri_v(t2, back);
            if (a<3||b<3||c<3||d<3) continue;
            int ic = incircle_2d(T.point(a),T.point(b),T.point(c),T.point(d),
                                 vec2h_length(T.point(a)), vec2h_length(T.point(b)),
                                 vec2h_length(T.point(c)), vec2h_length(T.point(d)));
            if (ic * g > 0) ++nd;
        }
    }
}

static Vec2HE P(double x, double y) { return Vec2HE(x, y, 1.0); }

int main() {
    std::cout << std::unitbuf;  // flush each write so a crash shows progress
    // ===== (A) ExactCDT2d: exact crossing vertex =============================
    {
        ExactCDT2d T;
        T.create_enclosing_triangle(P(-100,-100), P(100,-100), P(0,100));
        // Quad whose diagonals cross at (2/3, 14/9) - NOT a dyadic rational.
        index_t a = T.insert_point(P(0,0));
        index_t b = T.insert_point(P(3,0));
        index_t c = T.insert_point(P(3,7));
        index_t d = T.insert_point(P(0,2));
        index_t nv0 = T.nv();
        T.insert_constraint(a, c);     // line (0,0)-(3,7)
        T.insert_constraint(b, d);     // line (3,0)-(0,2) -> crosses a-c
        CHECK(T.nv() == nv0 + 1, "exactly one intersection vertex created");
        index_t x = nv0;               // the created vertex
        // EXACT: x lies on both constraint lines.
        CHECK(orient2d(T.point(a), T.point(c), T.point(x)) == 0,
              "intersection vertex is EXACTLY on line a-c");
        CHECK(orient2d(T.point(b), T.point(d), T.point(x)) == 0,
              "intersection vertex is EXACTLY on line b-d");
        // Teeth: the rounded double crossing point is NOT exactly on the lines.
        // Solve in double: t=2/9 along a->c gives (6/9,14/9).
        Vec2HE xd = P(6.0/9.0, 14.0/9.0);
        bool dbl_off = (orient2d(T.point(a), T.point(c), xd) != 0) ||
                       (orient2d(T.point(b), T.point(d), xd) != 0);
        CHECK(dbl_off, "the rounded double crossing point is OFF the lines (teeth)");
        // Valid exact CDT.
        int bad, nd; audit(T, bad, nd, true);
        CHECK(bad == 0, "(A) exact CDT structurally valid");
        CHECK(nd == 0, "(A) exact homogeneous-incircle Delaunay holds on free edges");
        // Constraint recovery through the intersection vertex.
        auto has = [&](index_t u, index_t w) {
            for (index_t t=0;t<T.nT();++t) for (index_t le=0;le<3;++le) {
                index_t p=T.tri_v(t,(le+1)%3), q=T.tri_v(t,(le+2)%3);
                if (((p==u&&q==w)||(p==w&&q==u)) && T.tri_edge_constrained(t,le)) return true;
            } return false; };
        CHECK(has(a,x)&&has(x,c)&&has(b,x)&&has(x,d),
              "(A) both constraints recovered through the intersection vertex");
    }

    // ===== (B) MeshInTriangle: remesh a tilted 3D facet ======================
    {
        // Facet in a tilted plane (forces a non-z normal_axis projection).
        double C0[3]={0,0,0}, C1[3]={6,0,3}, C2[3]={0,6,3};
        MeshInTriangle M;
        M.begin_facet(C0, C1, C2);
        // Four points strictly INSIDE the facet (x>=0, y>=0, x+y<6), on its
        // plane z = (x+y)/2. Constraints A-B (line x+y=4) and E-F (line y=x)
        // cross at (2,2), inside the facet -> exact intersection vertex.
        double A[3]={1,3,2}, B[3]={3,1,2}, E[3]={1,1,1}, F[3]={2.5,2.5,2.5};
        index_t a = M.add_vertex(A);
        index_t b = M.add_vertex(B);
        index_t e = M.add_vertex(E);
        index_t f = M.add_vertex(F);
        index_t nv0 = M.cdt().nv();
        M.add_constraint(a, b);    // crosses
        M.add_constraint(e, f);    // -> exact intersection inside the facet
        CHECK(M.cdt().nv() == nv0 + 1, "(B) facet remesh created one intersection vertex");
        int bad, nd; audit(M.cdt(), bad, nd, false);
        CHECK(bad == 0, "(B) facet retriangulation structurally valid (projected, exact)");
    }

    if (failures) {
        std::cerr << "\n" << failures << " mesh-in-triangle check(s) FAILED\n";
        return 1;
    }
    std::cout << "\nOK: exact per-facet remesh with division-free intersections.\n";
    return 0;
}
