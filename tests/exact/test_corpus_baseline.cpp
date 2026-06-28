// ============================================================================
//  test_corpus_baseline.cpp - baseline oracle for the exact-arrangement
//  self-intersection port.
// ----------------------------------------------------------------------------
//  This does NOT test the (not-yet-ported) exact kernel. Its job is to:
//    1. prove every corpus mesh loads,
//    2. record the "before" self-intersection count for each one using
//       SimpleMesh's EXISTING detector (MeshChecker, src/MeshChecker.cpp), and
//    3. assert that count matches what the corpus was DESIGNED to contain.
//
//  The end-to-end arrangement test loads
//  the same files, run resolve_self_intersections(), and assert the detector
//  then reports ZERO intersections. So this file fixes the regression oracle.
//
//  Usage: test_corpus_baseline [data_dir]
//         data_dir defaults to the compiled-in SM_ISECT_DATA_DIR.
// ============================================================================
#include "simplemesh/Mesh.h"
#include "simplemesh/MeshChecker.h"

#include <iostream>
#include <string>
#include <vector>

#ifndef SM_ISECT_DATA_DIR
#define SM_ISECT_DATA_DIR "."
#endif

static int failures = 0;
#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } \
        else         { std::cout << "ok:   " << msg << "\n"; }             \
    } while (0)

// Detect tri-tri self-intersections with the existing checker.
static size_t count_self_intersections(const sm::Mesh& m) {
    sm::CheckOptions opt;
    opt.check_self_intersections = true;
    sm::MeshChecker checker(m);
    return checker.check(opt).self_intersections.size();
}

// How the EXISTING double-precision detector is expected to behave on each
// case. Some of these intentionally differ from the exact-arrangement semantics;
// differences are the gaps the ported exact kernel must close (tracked in
// the exact path must handle).
enum class Behavior {
    CROSS,            // genuine interior crossing; detector reports > 0
    COPLANAR_SKIPPED, // coplanar overlap; double detector SKIPS it (MeshChecker.cpp:324)
                      //   -> exact port MUST detect & resolve this.
    BOUNDARY_TOUCH,   // single shared point; tolerant detector reports it as a touch
                      //   -> exact port treats it symbolically (vertex region), no shell removal.
};

struct Case {
    const char* file;
    Behavior    behavior;
    const char* note;
};

int main(int argc, char** argv) {
    std::string dir = (argc > 1) ? argv[1] : SM_ISECT_DATA_DIR;
    if (!dir.empty() && dir.back() != '/' && dir.back() != '\\') dir += '/';

    // Expectations are BY CONSTRUCTION (see generate_corpus.py) and match the
    // CURRENT double-precision detector, not yet the exact kernel.
    const std::vector<Case> cases = {
        {"two_tris_cross.obj",    Behavior::CROSS,            "X-crossing triangles"},
        {"two_quads_overlap.obj", Behavior::CROSS,            "crossing quad sheets"},
        {"cube_x_cube.obj",       Behavior::CROSS,            "overlapping cubes (boolean input)"},
        {"self_fold.obj",         Behavior::CROSS,            "folded ribbon pierces itself"},
        {"coplanar_overlap.obj",  Behavior::COPLANAR_SKIPPED, "coplanar area overlap"},
        {"tangent_touch.obj",     Behavior::BOUNDARY_TOUCH,   "single-point boundary touch"},
    };

    for (const Case& c : cases) {
        sm::Mesh m;
        bool loaded = m.read_obj(dir + c.file);
        CHECK(loaded, std::string("load ") + c.file);
        if (!loaded) continue;

        size_t n = count_self_intersections(m);
        std::cout << "      " << c.file << ": " << m.n_faces()
                  << " faces, " << n << " intersecting pair(s)  (" << c.note << ")\n";

        switch (c.behavior) {
        case Behavior::CROSS:
            CHECK(n > 0, std::string(c.file) + " has interior crossings as designed");
            break;
        case Behavior::COPLANAR_SKIPPED:
            CHECK(n == 0, std::string(c.file) +
                  " coplanar pair skipped by double detector (exact port must fix)");
            break;
        case Behavior::BOUNDARY_TOUCH:
            CHECK(n > 0, std::string(c.file) +
                  " boundary touch reported by tolerant detector (exact port treats symbolically)");
            break;
        }
    }

    if (failures) {
        std::cerr << "\n" << failures << " baseline check(s) FAILED\n";
        return 1;
    }
    std::cout << "\nBaseline OK: corpus loads and detector matches design.\n";
    return 0;
}
