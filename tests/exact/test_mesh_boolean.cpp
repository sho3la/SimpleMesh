// ============================================================================
//  test_mesh_boolean.cpp - shell classification & boolean ops.
// ----------------------------------------------------------------------------
//  The end-to-end payoff of the port: feed two overlapping closed solids into
//  the arrangement + classifier and get back a CLEAN boolean
//  result. For the two unit cubes A=[0,1]^3 and B=[0.5,1.5]^3 we verify, for
//  union / intersection / difference, that the result is:
//    * watertight & manifold & consistently oriented (MeshChecker),
//    * free of self-intersections (the whole point),
//    * geometrically correct: exact enclosed VOLUME and bounding box.
//        union        : V = 1 + 1 - 0.125 = 1.875,  bbox [0,1.5]^3
//        intersection : V = 0.125,                  bbox [0.5,1]^3
//        difference   : V = 1 - 0.125 = 0.875,      bbox [0,1]^3  (A \ B)
//  These invariants pin the classifier exactly.
// ============================================================================
#include "simplemesh/exact/MeshBoolean.h"
#include "simplemesh/Mesh.h"
#include "simplemesh/MeshChecker.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace sm::exact;
using sm::Mesh;
using sm::Vec3;
using sm::VertexHandle;
using sm::MeshChecker;
using sm::CheckOptions;

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

static bool load_obj_soup(const std::string& path, Pts& V, Tris& F) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() < 2) continue;
        std::istringstream ss(line);
        std::string tag; ss >> tag;
        if (tag == "v") { std::array<double,3> p; ss >> p[0] >> p[1] >> p[2]; V.push_back(p); }
        else if (tag == "f") {
            std::vector<std::size_t> idx; std::string tok;
            while (ss >> tok) {
                auto s = tok.find('/'); if (s != std::string::npos) tok = tok.substr(0, s);
                idx.push_back(std::size_t(std::stol(tok)) - 1);
            }
            for (std::size_t k = 2; k < idx.size(); ++k) F.push_back({ idx[0], idx[k-1], idx[k] });
        }
    }
    return true;
}

// Signed enclosed volume (divergence theorem). Positive for outward normals.
static double enclosed_volume(const BooleanResult& R) {
    double vol = 0.0;
    for (const auto& t : R.triangles) {
        const auto& a = R.points[t[0]]; const auto& b = R.points[t[1]]; const auto& c = R.points[t[2]];
        double cx = b[1]*c[2]-b[2]*c[1], cy = b[2]*c[0]-b[0]*c[2], cz = b[0]*c[1]-b[1]*c[0];
        vol += (a[0]*cx + a[1]*cy + a[2]*cz);
    }
    return vol / 6.0;
}

static void bbox(const BooleanResult& R, double lo[3], double hi[3]) {
    for (int a = 0; a < 3; ++a) { lo[a] = 1e300; hi[a] = -1e300; }
    for (const auto& p : R.points)
        for (int a = 0; a < 3; ++a) { lo[a] = std::min(lo[a], p[a]); hi[a] = std::max(hi[a], p[a]); }
}

// Build a half-edge mesh from the (manifold) boolean result, so MeshChecker can
// certify watertight/manifold/no-self-intersection.
static bool build_mesh(const BooleanResult& R, Mesh& M) {
    std::vector<VertexHandle> vh;
    for (const auto& p : R.points) vh.push_back(M.add_vertex(Vec3(p[0], p[1], p[2])));
    bool all = true;
    for (const auto& t : R.triangles) {
        if (t[0]==t[1] || t[1]==t[2] || t[2]==t[0]) continue;
        auto f = M.add_face({ vh[t[0]], vh[t[1]], vh[t[2]] });
        if (!f.is_valid()) all = false;
    }
    return all;
}

// Append an axis-aligned cube [lo,hi]^3 (outward-oriented) to a soup.
static void add_cube(Pts& V, Tris& F, double lx, double ly, double lz,
                     double hx, double hy, double hz) {
    std::size_t b = V.size();
    V.push_back({lx,ly,lz}); V.push_back({hx,ly,lz});
    V.push_back({hx,hy,lz}); V.push_back({lx,hy,lz});
    V.push_back({lx,ly,hz}); V.push_back({hx,ly,hz});
    V.push_back({hx,hy,hz}); V.push_back({lx,hy,hz});
    // 6 faces, outward CCW, each split into 2 triangles.
    int q[6][4] = {
        {0,3,2,1}, // bottom z=lz (normal -z)
        {4,5,6,7}, // top    z=hz (+z)
        {0,1,5,4}, // y=ly (-y)
        {2,3,7,6}, // y=hy (+y)
        {1,2,6,5}, // x=hx (+x)
        {0,4,7,3}  // x=lx (-x)
    };
    for (auto& f : q) {
        F.push_back({ b+std::size_t(f[0]), b+std::size_t(f[1]), b+std::size_t(f[2]) });
        F.push_back({ b+std::size_t(f[0]), b+std::size_t(f[2]), b+std::size_t(f[3]) });
    }
}

struct Expect { BoolOp op; const char* name; double vol; double lo[3]; double hi[3]; };

static void run_case(const Pts& V, const Tris& F, const Expect& e) {
    std::cout << "\n--- " << e.name << " ---\n";
    BooleanResult R = mesh_boolean(V, F, e.op);
    std::cout << "    operands=" << R.n_operands << "  result: "
              << R.points.size() << " v, " << R.triangles.size() << " t\n";
    CHECK(R.triangles.size() > 0, std::string(e.name) + ": produced a non-empty result");

    Mesh M;
    bool built = build_mesh(R, M);
    CHECK(built, std::string(e.name) + ": result is a valid (manifold) mesh - all faces added");

    CheckOptions opt; opt.check_self_intersections = true;
    auto rep = MeshChecker(M).check(opt);
    CHECK(rep.is_manifold, std::string(e.name) + ": manifold");
    CHECK(rep.is_closed,   std::string(e.name) + ": watertight (closed)");
    CHECK(rep.is_oriented, std::string(e.name) + ": consistently oriented");
    CHECK(rep.self_intersections.empty(), std::string(e.name) + ": no self-intersections");

    double vol = enclosed_volume(R);
    std::cout << "    volume = " << vol << " (expected " << e.vol << ")\n";
    CHECK(std::fabs(vol - e.vol) < 1e-9, std::string(e.name) + ": correct enclosed volume");

    double lo[3], hi[3]; bbox(R, lo, hi);
    bool box_ok = true;
    for (int a = 0; a < 3; ++a)
        box_ok = box_ok && std::fabs(lo[a]-e.lo[a]) < 1e-12 && std::fabs(hi[a]-e.hi[a]) < 1e-12;
    std::cout << "    bbox = [" << lo[0]<<","<<lo[1]<<","<<lo[2] << "]..["
              << hi[0]<<","<<hi[1]<<","<<hi[2] << "]\n";
    CHECK(box_ok, std::string(e.name) + ": correct bounding box");
}

int main() {
    std::cout << std::unitbuf;
    const std::string path = std::string(SM_ISECT_DATA_DIR) + "/cube_x_cube.obj";
    Pts V; Tris F;
    if (!load_obj_soup(path, V, F)) {
        std::cerr << "FAIL: cannot open " << path << "\n"; return 1;
    }
    std::cout << "cube_x_cube: " << V.size() << " v, " << F.size() << " t\n";

    Expect cases[] = {
        { BoolOp::Union,        "union",        1.875, {0,0,0},          {1.5,1.5,1.5} },
        { BoolOp::Intersection, "intersection", 0.125, {0.5,0.5,0.5},    {1,1,1}       },
        { BoolOp::Difference,   "difference",   0.875, {0,0,0},          {1,1,1}       },
    };
    for (const auto& e : cases) run_case(V, F, e);

    // ---- Nested cubes: small fully inside big, DISJOINT surfaces. ----------
    // Two connected components -> exercises classify's ray-traced component
    // inclusion (ray-traced per object), the part a
    // naive per-facet probe gets wrong. big=A=[0,4]^3, small=B=[1,2]^3 (B subset A).
    std::cout << "\n========== nested cubes (small inside big) ==========\n";
    Pts V2; Tris F2;
    add_cube(V2, F2, 0,0,0, 4,4,4);   // operand 0 = big (A)
    add_cube(V2, F2, 1,1,1, 2,2,2);   // operand 1 = small (B)
    Expect nested[] = {
        { BoolOp::Union,        "nested-union",        64.0, {0,0,0}, {4,4,4} },
        { BoolOp::Intersection, "nested-intersection",  1.0, {1,1,1}, {2,2,2} },
        { BoolOp::Difference,   "nested-difference",   63.0, {0,0,0}, {4,4,4} }, // hollow
    };
    for (const auto& e : nested) run_case(V2, F2, e);

    if (failures) {
        std::cerr << "\n" << failures << " boolean check(s) FAILED\n";
        return 1;
    }
    std::cout << "\nOK: exact boolean ops (union/intersection/difference) are clean & correct.\n";
    return 0;
}
