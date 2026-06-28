// ============================================================================
//  test_surface_intersection.cpp - arrangement assembly.
// ----------------------------------------------------------------------------
//  The acceptance property of the whole port: after the exact arrangement, no
//  two facets interpenetrate any more - every intersection curve has become a
//  shared mesh edge. We verify, on each corpus mesh:
//
//    (1) BEFORE: the input soup has the documented number of intersecting pairs
//        (exact Day-6 detector among non-adjacent triangles).
//    (2) resolve_intersections() runs and (for the crossing cases) creates new
//        intersection vertices and sub-triangles.
//    (3) AFTER: ZERO non-adjacent triangle pairs still have a non-degenerate
//        (dim>=1) intersection - i.e. all interpenetration is gone, only shared
//        edges/vertices remain (which are adjacency, skipped exactly as
//        MeshChecker.cpp:392 skips shared-vertex pairs).
//
//  The oracle is the EXACT detector, which - unlike the tolerant double
//  Moller test - also sees coplanar overlap, so coplanar_overlap is a real test
//  here, not a skipped case.
// ============================================================================
#include "simplemesh/exact/SurfaceIntersection.h"

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace sm::exact;

#ifndef SM_ISECT_DATA_DIR
#define SM_ISECT_DATA_DIR "."
#endif

static int failures = 0;
#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } \
        else         { std::cout << "ok:   " << msg << "\n"; }             \
    } while (0)

using Pts = std::vector<std::array<double,3>>;
using Tris = std::vector<std::array<std::size_t,3>>;

// Minimal OBJ soup reader (triangulates polygons by fan). We read a SOUP, not a
// half-edge mesh, because the input may already be two separate components and
// the output is non-manifold.
static bool load_obj_soup(const std::string& path, Pts& V, Tris& F) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() < 2) continue;
        std::istringstream ss(line);
        std::string tag; ss >> tag;
        if (tag == "v") {
            std::array<double,3> p; ss >> p[0] >> p[1] >> p[2];
            V.push_back(p);
        } else if (tag == "f") {
            std::vector<std::size_t> idx;
            std::string tok;
            while (ss >> tok) {
                std::size_t slash = tok.find('/');
                if (slash != std::string::npos) tok = tok.substr(0, slash);
                idx.push_back(std::size_t(std::stol(tok)) - 1);  // 1-based -> 0
            }
            for (std::size_t k = 2; k < idx.size(); ++k)
                F.push_back({ idx[0], idx[k-1], idx[k] });
        }
    }
    return true;
}

static bool share_vertex(const std::array<std::size_t,3>& a,
                         const std::array<std::size_t,3>& b) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) if (a[i] == b[j]) return true;
    return false;
}

// Count non-adjacent triangle pairs that have a non-degenerate intersection
// (exact Day-6 oracle). Adjacency (shared vertex/edge) is skipped.
static int count_interpenetrations(const Pts& V, const Tris& F) {
    int n = 0;
    for (std::size_t i = 0; i < F.size(); ++i)
        for (std::size_t j = i+1; j < F.size(); ++j) {
            if (share_vertex(F[i], F[j])) continue;
            if (triangles_intersections(
                    V[F[i][0]].data(), V[F[i][1]].data(), V[F[i][2]].data(),
                    V[F[j][0]].data(), V[F[j][1]].data(), V[F[j][2]].data()))
                ++n;
        }
    return n;
}

struct Case { const char* name; bool expect_crossing; };

int main() {
    std::cout << std::unitbuf;
    const std::string dir = SM_ISECT_DATA_DIR;

    // coplanar_overlap: the double detector reports 0, but the EXACT detector
    // sees the area overlap, so it IS a crossing case here.
    // tangent_touch: only a shared point -> exact detector treats as adjacency.
    Case cases[] = {
        { "two_tris_cross",     true  },
        { "two_quads_overlap",  true  },
        { "cube_x_cube",        true  },
        { "self_fold",          true  },
        { "coplanar_overlap",   true  },
        { "tangent_touch",      false },
    };

    for (const Case& c : cases) {
        std::string path = dir + "/" + c.name + ".obj";
        Pts V; Tris F;
        if (!load_obj_soup(path, V, F)) {
            std::cerr << "FAIL: cannot open " << path << "\n"; ++failures; continue;
        }
        std::cout << "\n=== " << c.name << " (" << V.size() << " v, "
                  << F.size() << " t) ===\n";

        int before = count_interpenetrations(V, F);
        std::cout << "    interpenetrating pairs before: " << before << "\n";
        if (c.expect_crossing)
            CHECK(before > 0, std::string(c.name) + ": input has interpenetration to resolve");
        else
            CHECK(before == 0, std::string(c.name) + ": input has no interpenetration (adjacency only)");

        ArrangementResult R = resolve_intersections(V, F);
        std::cout << "    intersecting pairs: " << R.n_intersection_pairs
                  << ", output: " << R.points.size() << " v, "
                  << R.triangles.size() << " t\n";

        int after = count_interpenetrations(R.points, R.triangles);
        CHECK(after == 0, std::string(c.name) + ": ZERO interpenetration after arrangement");

        if (c.expect_crossing) {
            CHECK(R.n_intersection_pairs > 0, std::string(c.name) + ": found intersecting facet pairs");
            CHECK(R.points.size() > V.size(),
                  std::string(c.name) + ": new intersection vertices were created");
            CHECK(R.triangles.size() > F.size(),
                  std::string(c.name) + ": facets were subdivided");
        }
    }

    if (failures) {
        std::cerr << "\n" << failures << " arrangement check(s) FAILED\n";
        return 1;
    }
    std::cout << "\nOK: exact arrangement removes all interpenetration on the corpus.\n";
    return 0;
}
