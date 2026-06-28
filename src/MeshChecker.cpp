// ============================================================================
//  SimpleMesh - MeshChecker.cpp : implementation of the validity checker
// ============================================================================
#include "simplemesh/MeshChecker.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <unordered_map>

namespace sm {

// ============================================================================
//  Report helpers
// ============================================================================
bool MeshCheckReport::ok() const {
    return bad_halfedges.empty() && nonmanifold_vertices.empty() &&
           nonmanifold_edges.empty() && degenerate_faces.empty() &&
           duplicate_vertices.empty() && duplicate_faces.empty() &&
           self_intersections.empty();
}

static void append_offenders(std::ostringstream& os, const char* label,
                             size_t count, const std::vector<int>& ex) {
    os << "  " << label << ": " << count;
    if (count) {
        os << "  e.g. [";
        for (size_t i = 0; i < ex.size() && i < 8; ++i) os << (i ? "," : "") << ex[i];
        os << (ex.size() > 8 ? ",..." : "") << "]";
    }
    os << "\n";
}

std::string MeshCheckReport::summary() const {
    std::ostringstream os;
    os << (ok() ? "MeshCheck: PASS\n" : "MeshCheck: FAIL\n");
    os << "  components: " << n_components
       << "  euler(V-E+F): " << euler
       << "  genus: " << genus
       << "  boundary loops: " << n_boundary_loops << "\n";
    os << "  closed: "   << (is_closed   ? "yes" : "no")
       << "  manifold: " << (is_manifold ? "yes" : "no")
       << "  oriented: " << (is_oriented ? "yes" : "no") << "\n";

    append_offenders(os, "bad halfedges",        bad_halfedges.size(), bad_halfedges);
    append_offenders(os, "non-manifold vertices", nonmanifold_vertices.size(), nonmanifold_vertices);
    append_offenders(os, "non-manifold edges",    nonmanifold_edges.size(), nonmanifold_edges);
    append_offenders(os, "isolated vertices",     isolated_vertices.size(), isolated_vertices);
    append_offenders(os, "degenerate faces",      degenerate_faces.size(), degenerate_faces);
    os << "  boundary edges: " << boundary_edges.size() << "\n";
    os << "  duplicate vertices: " << duplicate_vertices.size() << "\n";
    os << "  duplicate faces: "    << duplicate_faces.size() << "\n";
    if (!self_intersections.empty())
        os << "  self-intersections: " << self_intersections.size() << "\n";
    return os.str();
}

// ============================================================================
//  A. Half-edge connectivity integrity
// ----------------------------------------------------------------------------
//  These are the structural invariants every later algorithm assumes. We test
//  each live half-edge independently so one broken link doesn't hide others.
// ============================================================================
static void check_connectivity(const Mesh& m, MeshCheckReport& r) {
    for (auto h : m.all_halfedges()) {
        bool bad = false;

        // opposite is an involution, encoded as h XOR 1
        if (m.opposite_halfedge(m.opposite_halfedge(h)) != h) bad = true;

        // next / prev are inverses
        if (m.prev_halfedge(m.next_halfedge(h)) != h) bad = true;
        if (m.next_halfedge(m.prev_halfedge(h)) != h) bad = true;

        // the chain is vertex-consistent: h ends where next(h) starts
        if (m.to_vertex(h) != m.from_vertex(m.next_halfedge(h))) bad = true;

        // face is constant around the loop next(h) follows
        if (m.face(h) != m.face(m.next_halfedge(h))) bad = true;

        // a boundary half-edge's loop must stay on the boundary
        if (m.is_boundary(h) && !m.is_boundary(m.next_halfedge(h))) bad = true;

        if (bad) r.bad_halfedges.push_back(h.idx());
    }
}

// ============================================================================
//  B. Manifoldness, boundary loops, components, Euler / genus
// ============================================================================

// A vertex is non-manifold ("bow-tie") if its incident faces form more than one
// fan. Two independent signals catch every case:
//   (1) MORE THAN ONE outgoing BOUNDARY half-edge -> the fans meet the boundary
//       in several separate sectors (the usual bow-tie built via add_face; a
//       single fan walk can't see this because the boundary stitches the sectors
//       into one loop).
//   (2) the single fan reachable by circulating next(opposite(.)) is smaller
//       than the vertex's degree -> a fully-interior pinch.
static void check_nonmanifold_vertices(const Mesh& m, MeshCheckReport& r,
                                       const std::vector<int>& outgoing_count,
                                       const std::vector<int>& boundary_out_count) {
    for (auto v : m.all_vertices()) {
        if (m.is_isolated(v)) continue;
        int total = outgoing_count[v.idx()];
        if (total <= 1) continue;

        if (boundary_out_count[v.idx()] > 1) {        // (1) multiple boundary sectors
            r.nonmanifold_vertices.push_back(v.idx());
            continue;
        }

        // (2) walk one fan and compare to the degree
        HalfedgeHandle start = m.halfedge(v);
        HalfedgeHandle h = start;
        int fan = 0;
        const int cap = total + 2;
        do {
            ++fan;
            h = m.next_halfedge(m.opposite_halfedge(h));
            if (fan > cap) break;                     // broken link guard
        } while (h != start);

        if (fan != total) r.nonmanifold_vertices.push_back(v.idx());
    }
}

// Union-find over vertices, joined by edges -> connected components.
namespace {
struct DSU {
    std::vector<int> p;
    explicit DSU(int n) : p(n) { for (int i = 0; i < n; ++i) p[i] = i; }
    int find(int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
    void join(int a, int b) { p[find(a)] = find(b); }
};
} // namespace

static void analyze_topology(const Mesh& m, MeshCheckReport& r) {
    // ---- per-vertex outgoing half-edge counts (one scan) -----------------
    std::vector<int> outgoing(m.n_vertices(), 0), boundary_out(m.n_vertices(), 0);
    for (auto h : m.all_halfedges()) {
        int fv = m.from_vertex(h).idx();
        outgoing[fv]++;
        if (m.is_boundary(h)) boundary_out[fv]++;
    }

    check_nonmanifold_vertices(m, r, outgoing, boundary_out);

    // ---- isolated vertices ------------------------------------------------
    for (auto v : m.all_vertices())
        if (m.is_isolated(v)) r.isolated_vertices.push_back(v.idx());

    // ---- non-manifold edges + boundary edges ------------------------------
    // In our pair-encoded kernel an edge can never be shared by >2 faces, but we
    // still verify it (a real check in array-of-faces libraries) and tally the
    // boundary. A "non-manifold edge" here = both half-edges boundary yet the
    // edge is not isolated, or the two faces coincide.
    int live_edges = 0;
    for (auto e : m.all_edges()) {
        ++live_edges;
        HalfedgeHandle h0 = m.halfedge(e, 0), h1 = m.halfedge(e, 1);
        if (m.is_boundary(e)) r.boundary_edges.push_back(e.idx());
        FaceHandle f0 = m.face(h0), f1 = m.face(h1);
        if (f0.is_valid() && f0 == f1)            // same face on both sides
            r.nonmanifold_edges.push_back(e.idx());
    }
    r.is_closed = r.boundary_edges.empty();

    // ---- boundary loops ---------------------------------------------------
    std::vector<char> seen(m.n_halfedges(), 0);
    int loops = 0;
    for (auto h : m.all_halfedges()) {
        if (!m.is_boundary(h) || seen[h.idx()]) continue;
        ++loops;
        HalfedgeHandle c = h;
        int cap = static_cast<int>(m.n_halfedges()) + 1;
        do { seen[c.idx()] = 1; c = m.next_halfedge(c); }
        while (c != h && cap-- > 0);
    }
    r.n_boundary_loops = loops;

    // ---- connected components (union vertices via edges) ------------------
    int nV = static_cast<int>(m.n_vertices());
    DSU dsu(nV);
    for (auto e : m.all_edges()) {
        HalfedgeHandle h = m.halfedge(e, 0);
        dsu.join(m.from_vertex(h).idx(), m.to_vertex(h).idx());
    }
    int comps = 0;
    for (auto v : m.all_vertices())
        if (dsu.find(v.idx()) == v.idx()) ++comps;     // count live roots
    r.n_components = comps;

    // ---- Euler characteristic & genus ------------------------------------
    int V = 0; for (auto v : m.all_vertices()) { (void)v; ++V; }
    int F = 0; for (auto f : m.all_faces())    { (void)f; ++F; }
    int E = live_edges;
    r.euler = V - E + F;
    // For a closed orientable manifold of C components: chi = 2C - 2g, so
    // g = (2C - chi)/2.  With boundary loops B: chi = 2C - 2g - B.
    r.genus = (2 * r.n_components - r.euler - r.n_boundary_loops) / 2;

    // ---- orientation: every edge's two half-edges run opposite ------------
    bool oriented = true;
    for (auto e : m.all_edges()) {
        HalfedgeHandle h0 = m.halfedge(e, 0), h1 = m.halfedge(e, 1);
        if (m.to_vertex(h0) != m.from_vertex(h1) ||
            m.to_vertex(h1) != m.from_vertex(h0)) { oriented = false; break; }
    }
    r.is_oriented = oriented;

    r.is_manifold = r.bad_halfedges.empty() && r.nonmanifold_vertices.empty() &&
                    r.nonmanifold_edges.empty();
}

// ============================================================================
//  C. Degeneracy & duplication
// ============================================================================
static void check_degenerate_faces(const Mesh& m, MeshCheckReport& r) {
    for (auto f : m.all_faces()) {
        auto fv = m.face_vertices(f);
        bool bad = fv.size() < 3;
        for (size_t i = 0; i < fv.size() && !bad; ++i)        // repeated vertex?
            for (size_t j = i + 1; j < fv.size(); ++j)
                if (fv[i] == fv[j]) { bad = true; break; }
        if (!bad && m.calc_face_area(f) < 1e-18) bad = true;  // ~zero area
        if (bad) r.degenerate_faces.push_back(f.idx());
    }
}

// FNV-1a hashers for the spatial grid / face-set lookups (O(1), neighbor-aware).
namespace {
struct Cell { int x, y, z; bool operator==(const Cell& o) const { return x == o.x && y == o.y && z == o.z; } };
struct CellHash { size_t operator()(const Cell& c) const {
    size_t h = 1469598103934665603ULL;
    for (int v : { c.x, c.y, c.z }) { h ^= (size_t)(unsigned)v; h *= 1099511628211ULL; }
    return h; } };
struct VecHash { size_t operator()(const std::vector<int>& v) const {
    size_t h = 1469598103934665603ULL;
    for (int x : v) { h ^= (size_t)(unsigned)x; h *= 1099511628211ULL; }
    return h; } };
} // namespace

// Neighbor-aware spatial hash: a vertex within tol of an already-seen one (in
// its cell OR a neighbour cell) is flagged duplicate. Catches cross-cell pairs.
static void check_duplicate_vertices(const Mesh& m, MeshCheckReport& r, double tol) {
    if (tol <= 0) return;
    const double inv = 1.0 / tol, tol2 = tol * tol;
    auto cell_of = [&](const Vec3& p) {
        return Cell{ (int)std::floor(p.x * inv), (int)std::floor(p.y * inv), (int)std::floor(p.z * inv) };
    };
    std::unordered_map<Cell, std::vector<int>, CellHash> grid;
    for (auto v : m.all_vertices()) {
        const Vec3& p = m.point(v);
        Cell c = cell_of(p);
        int match = -1;
        for (int dz = -1; dz <= 1 && match < 0; ++dz)
        for (int dy = -1; dy <= 1 && match < 0; ++dy)
        for (int dx = -1; dx <= 1 && match < 0; ++dx) {
            auto it = grid.find(Cell{ c.x + dx, c.y + dy, c.z + dz });
            if (it == grid.end()) continue;
            for (int q : it->second)
                if ((m.point(VertexHandle(q)) - p).sqrnorm() <= tol2) { match = q; break; }
        }
        if (match >= 0) r.duplicate_vertices.emplace_back(match, v.idx());
        else grid[c].push_back(v.idx());
    }
}

// Two faces are duplicates if their vertex SETS match (order/winding ignored).
static void check_duplicate_faces(const Mesh& m, MeshCheckReport& r) {
    std::unordered_map<std::vector<int>, int, VecHash> seen;   // sorted verts -> first face
    for (auto f : m.all_faces()) {
        std::vector<int> key;
        for (auto v : m.face_vertices(f)) key.push_back(v.idx());
        std::sort(key.begin(), key.end());
        auto it = seen.find(key);
        if (it != seen.end()) r.duplicate_faces.emplace_back(it->second, f.idx());
        else seen.emplace(std::move(key), f.idx());
    }
}

// ============================================================================
//  D. Triangle/triangle self-intersection (Moller's interval-overlap test)
// ============================================================================
namespace {

// signed distances of triangle V's vertices to plane of triangle U
inline void plane_dists(const Vec3 V[3], const Vec3& N, double d, double dv[3]) {
    for (int i = 0; i < 3; ++i) dv[i] = N.dot(V[i]) + d;
}

// Compute the [min,max] interval where triangle (with vertex param coords pv on
// the intersection line L, and signed plane distances d) crosses L.
inline bool tri_interval(const double pv[3], const double d[3], double out[2]) {
    // arrange so the lone vertex on one side is isolated
    int i0 = 0, i1 = 1, i2 = 2;
    // pick the vertex whose sign differs (product of the other two > 0)
    if (d[0] * d[1] > 0)      { i0 = 2; i1 = 0; i2 = 1; }
    else if (d[0] * d[2] > 0) { i0 = 1; i1 = 0; i2 = 2; }
    double t0 = pv[i1] + (pv[i0] - pv[i1]) * d[i1] / (d[i1] - d[i0]);
    double t1 = pv[i2] + (pv[i0] - pv[i2]) * d[i2] / (d[i2] - d[i0]);
    out[0] = std::min(t0, t1);
    out[1] = std::max(t0, t1);
    return true;
}

bool tri_tri_intersect(const Vec3 U[3], const Vec3 V[3]) {
    const double EPS = 1e-12;
    Vec3 N1 = (U[1] - U[0]).cross(U[2] - U[0]);
    double d1 = -N1.dot(U[0]);
    double dv[3]; plane_dists(V, N1, d1, dv);
    if (dv[0] >  EPS && dv[1] >  EPS && dv[2] >  EPS) return false;
    if (dv[0] < -EPS && dv[1] < -EPS && dv[2] < -EPS) return false;

    Vec3 N2 = (V[1] - V[0]).cross(V[2] - V[0]);
    double d2 = -N2.dot(V[0]);
    double du[3]; plane_dists(U, N2, d2, du);
    if (du[0] >  EPS && du[1] >  EPS && du[2] >  EPS) return false;
    if (du[0] < -EPS && du[1] < -EPS && du[2] < -EPS) return false;

    Vec3 D = N1.cross(N2);
    if (D.sqrnorm() < EPS) return false;          // coplanar: skip (teaching simplification)

    // project triangle vertices onto the intersection line direction D
    int axis = 0; double a = std::abs(D.x);
    if (std::abs(D.y) > a) { a = std::abs(D.y); axis = 1; }
    if (std::abs(D.z) > a) { axis = 2; }
    auto coord = [&](const Vec3& p) { return p[axis]; };

    double pu[3] = { coord(U[0]), coord(U[1]), coord(U[2]) };
    double pv[3] = { coord(V[0]), coord(V[1]), coord(V[2]) };
    double iu[2], iv[2];
    tri_interval(pu, du, iu);
    tri_interval(pv, dv, iv);
    return iu[0] <= iv[1] + EPS && iv[0] <= iu[1] + EPS;   // intervals overlap
}

} // namespace

// Sweep-and-prune broad phase. Naive all-pairs is O(F^2) in BOTH time and the
// size of the result vector, which on a large self-intersecting mesh thrashes
// memory. We instead sort triangles by their AABB min-x and only test pairs
// whose x-intervals overlap, then reject by full AABB before the exact test.
// Two hard caps (face count + hit count) guarantee bounded work and memory.
static void check_self_intersections(const Mesh& m, MeshCheckReport& r,
                                     const CheckOptions& opt) {
    struct Tri { int f; Vec3 p[3]; int v[3]; Vec3 lo, hi; };
    std::vector<Tri> tris;
    for (auto f : m.all_faces()) {
        auto fv = m.face_vertices(f);
        if (fv.size() != 3) continue;
        Tri t; t.f = f.idx();
        t.lo = Vec3( 1e300,  1e300,  1e300);
        t.hi = Vec3(-1e300, -1e300, -1e300);
        for (int k = 0; k < 3; ++k) {
            t.p[k] = m.point(fv[k]); t.v[k] = fv[k].idx();
            for (int a = 0; a < 3; ++a) {
                t.lo[a] = std::min(t.lo[a], t.p[k][a]);
                t.hi[a] = std::max(t.hi[a], t.p[k][a]);
            }
        }
        tris.push_back(t);
    }

    // safety guard: refuse the check on meshes too large to bound (caller can
    // raise the cap explicitly). 0 means "no limit".
    if (opt.self_intersection_max_faces &&
        tris.size() > opt.self_intersection_max_faces)
        return;

    std::sort(tris.begin(), tris.end(),
              [](const Tri& a, const Tri& b) { return a.lo.x < b.lo.x; });

    auto aabb_overlap = [](const Tri& a, const Tri& b) {
        return a.lo.x <= b.hi.x && b.lo.x <= a.hi.x &&
               a.lo.y <= b.hi.y && b.lo.y <= a.hi.y &&
               a.lo.z <= b.hi.z && b.lo.z <= a.hi.z;
    };
    auto shares_vertex = [](const Tri& a, const Tri& b) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) if (a.v[i] == b.v[j]) return true;
        return false;
    };

    const size_t cap = opt.self_intersection_max_hits;
    for (size_t i = 0; i < tris.size(); ++i) {
        for (size_t j = i + 1; j < tris.size(); ++j) {
            if (tris[j].lo.x > tris[i].hi.x) break;     // sweep: no later x can overlap
            if (!aabb_overlap(tris[i], tris[j])) continue;
            if (shares_vertex(tris[i], tris[j])) continue;
            if (tri_tri_intersect(tris[i].p, tris[j].p)) {
                r.self_intersections.emplace_back(tris[i].f, tris[j].f);
                if (cap && r.self_intersections.size() >= cap) return;   // bound memory
            }
        }
    }
}

// ============================================================================
//  Driver
// ============================================================================
MeshCheckReport MeshChecker::check(const CheckOptions& opt) const {
    MeshCheckReport r;
    check_connectivity(mesh_, r);
    analyze_topology(mesh_, r);
    check_degenerate_faces(mesh_, r);
    if (opt.check_duplicate_vertices) check_duplicate_vertices(mesh_, r, opt.vertex_merge_tol);
    if (opt.check_duplicate_faces)    check_duplicate_faces(mesh_, r);
    if (opt.check_self_intersections) check_self_intersections(mesh_, r, opt);
    return r;
}

bool MeshChecker::is_valid() const {
    MeshCheckReport r;
    check_connectivity(mesh_, r);
    analyze_topology(mesh_, r);
    check_degenerate_faces(mesh_, r);
    return r.bad_halfedges.empty() && r.nonmanifold_vertices.empty() &&
           r.nonmanifold_edges.empty() && r.degenerate_faces.empty();
}

} // namespace sm
