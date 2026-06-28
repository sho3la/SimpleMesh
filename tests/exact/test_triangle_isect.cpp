// ============================================================================
//  test_triangle_isect.cpp - test for symbolic triangle/triangle isect.
// ----------------------------------------------------------------------------
//  Validates the symbolic region pairs against hand-derived golden sets, plus
//  the properties the arrangement engine relies on:
//    * a proper non-coplanar crossing yields the exact (T1-region, T2-region)
//      pairs naming the two endpoints of the intersection segment,
//    * COPLANAR overlap is detected (the case the old double Moller test in
//      MeshChecker.cpp:324 silently skips),
//    * a single shared vertex is a degenerate (dim-0) contact, not a crossing,
//    * disjoint triangles report nothing,
//    * swapping T1<->T2 yields the swapped region set (swap_T1_T2 symmetry).
// ============================================================================
#include "simplemesh/exact/TriangleIntersection.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

using namespace sm::exact;

static int failures = 0;
#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } \
        else         { std::cout << "ok:   " << msg << "\n"; }             \
    } while (0)

using V = double[3];

static std::vector<TriangleIsect> sorted(std::vector<TriangleIsect> v) {
    std::sort(v.begin(), v.end());
    return v;
}
static bool eq(std::vector<TriangleIsect> a, std::vector<TriangleIsect> b) {
    return sorted(std::move(a)) == sorted(std::move(b));
}
static std::string show(const std::vector<TriangleIsect>& v) {
    static const char* n[T_RGN_NB] = {"T1.P0","T1.P1","T1.P2","T2.P0","T2.P1","T2.P2",
        "T1.E0","T1.E1","T1.E2","T2.E0","T2.E1","T2.E2","T1.T","T2.T"};
    std::string s;
    for (auto& p : v) s += std::string("(")+n[p.first]+","+n[p.second]+") ";
    return s;
}

int main() {
    // ===== A: proper non-coplanar crossing (corpus two_tris_cross) ===========
    // T1 in z=0; T2 in plane y=0.5. Edge q0-q1 (=T2.E2) pierces T1 interior at
    // (0.5,0.5,0); vertex q2=(1.5,0.5,0) (=T2.P2) lies exactly on T1 edge p1-p2
    // (=T1.E0). EXACT golden: {(T1.T,T2.E2),(T1.E0,T2.P2)}.
    {
        V p0={0,0,0}, p1={2,0,0}, p2={0,2,0};
        V q0={0.5,0.5,-1}, q1={0.5,0.5,1}, q2={1.5,0.5,0};
        std::vector<TriangleIsect> r;
        bool hit = triangles_intersections(p0,p1,p2, q0,q1,q2, r);
        std::cout << "      A crossing: " << show(r) << "\n";
        CHECK(hit, "A: crossing reported as non-degenerate");
        CHECK(eq(r, {{T1_RGN_T,T2_RGN_E2}, {T1_RGN_E0,T2_RGN_P2}}),
              "A: exact symbolic region set matches golden");
        // bool-only entry point agrees
        CHECK(triangles_intersections(p0,p1,p2,q0,q1,q2) == hit,
              "A: predicate-only entry agrees with symbolic");
        // swap symmetry: result(T2,T1) == swap of result(T1,T2)
        std::vector<TriangleIsect> rs;
        triangles_intersections(q0,q1,q2, p0,p1,p2, rs);
        std::vector<TriangleIsect> rs_back;
        for (auto& pr : rs) {
            TriangleRegion a = swap_T1_T2(pr.first), b = swap_T1_T2(pr.second);
            if (!is_in_T1(a)) std::swap(a, b);
            rs_back.push_back({a, b});
        }
        CHECK(eq(rs_back, r), "A: swap_T1_T2(result(T2,T1)) == result(T1,T2)");
    }

    // ===== B: coplanar overlap (the case the double Moller test skips) =======
    {
        V p0={0,0,0}, p1={2,0,0}, p2={0,2,0};
        V q0={1.5,-0.5,0}, q1={1.5,1.5,0}, q2={-0.5,1.5,0};
        std::vector<TriangleIsect> r;
        bool hit = triangles_intersections(p0,p1,p2, q0,q1,q2, r);
        std::cout << "      B coplanar: " << show(r) << "\n";
        CHECK(hit, "B: coplanar overlap detected (exact succeeds where Moller skips)");
        CHECK(!r.empty(), "B: coplanar overlap produced intersection records");
    }

    // ===== C: single shared vertex -> degenerate dim-0 contact ===============
    {
        V p0={0,0,0}, p1={1,0,0}, p2={0,1,0};
        V q0={0,0,0}, q1={-1,0,0.5}, q2={0,-1,1};   // share (0,0,0)=p0=q0
        std::vector<TriangleIsect> r;
        bool hit = triangles_intersections(p0,p1,p2, q0,q1,q2, r);
        std::cout << "      C shared-vertex: " << show(r) << "\n";
        CHECK(!hit, "C: single shared vertex is NOT a non-degenerate intersection");
        CHECK(eq(r, {{T1_RGN_P0, T2_RGN_P0}}), "C: result is exactly the shared vertex pair");
    }

    // ===== D: disjoint, far apart ============================================
    {
        V p0={0,0,0}, p1={1,0,0}, p2={0,1,0};
        V q0={10,10,10}, q1={11,10,10}, q2={10,11,10};
        std::vector<TriangleIsect> r;
        bool hit = triangles_intersections(p0,p1,p2, q0,q1,q2, r);
        CHECK(!hit && r.empty(), "D: disjoint triangles report nothing");
    }

    // ===== E: shared full edge = ADJACENCY, not a transversal intersection ===
    // Two triangles sharing edge p0-p1, folded apart. The intersection IS the
    // shared edge, but it is detected as two COINCIDENT-VERTEX pairs (dim 0), so
    // has_non_degenerate is false: a shared edge is normal mesh
    // adjacency, NOT a self-intersection to remesh. (Same rule as a shared
    // vertex in case C.)
    {
        V p0={0,0,0}, p1={1,0,0}, p2={0,1,0};
        V q0={0,0,0}, q1={1,0,0}, q2={0,0,1};   // share edge (0,0,0)-(1,0,0)
        std::vector<TriangleIsect> r;
        bool hit = triangles_intersections(p0,p1,p2, q0,q1,q2, r);
        std::cout << "      E shared-edge: " << show(r) << "\n";
        CHECK(!hit, "E: a perfectly shared edge is adjacency, not a transversal isect");
        CHECK(eq(r, {{T1_RGN_P0,T2_RGN_P0}, {T1_RGN_P1,T2_RGN_P1}}),
              "E: shared edge reported as its two coincident-vertex pairs");
    }

    if (failures) {
        std::cerr << "\n" << failures << " triangle-isect check(s) FAILED\n";
        return 1;
    }
    std::cout << "\nOK: symbolic triangle intersection matches golden sets.\n";
    return 0;
}
