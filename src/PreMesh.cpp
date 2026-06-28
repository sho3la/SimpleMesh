// ============================================================================
//  SimpleMesh - PreMesh.cpp : the pre-processing repair layer
// ----------------------------------------------------------------------------
//  Cleans a raw face soup (weld, drop degenerate/duplicate, reorient) and
//  separates non-manifold edges and bow-tie vertices into manifold fans before
//  building the halfedge mesh. As a final safety net, any face the kernel still
//  rejects at emit() is re-emitted on duplicated vertices so nothing is dropped.
// ============================================================================
#include "simplemesh/PreMesh.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace {
// ---- hashing helpers (replace std::map/std::tuple hot paths) ---------------
// FNV-1a over a sequence of ints.
inline size_t fnv(std::initializer_list<int> vals) {
    size_t h = 1469598103934665603ULL;
    for (int v : vals) { h ^= (size_t)(unsigned)v; h *= 1099511628211ULL; }
    return h;
}
struct Cell { int x, y, z; bool operator==(const Cell& o) const { return x == o.x && y == o.y && z == o.z; } };
struct CellHash { size_t operator()(const Cell& c) const { return fnv({ c.x, c.y, c.z }); } };
struct VecHash  { size_t operator()(const std::vector<int>& v) const {
    size_t h = 1469598103934665603ULL;
    for (int x : v) { h ^= (size_t)(unsigned)x; h *= 1099511628211ULL; }
    return h; } };
} // namespace

namespace sm {

uint64_t PreMesh::ekey(int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
            static_cast<uint32_t>(b);
}

// ----------------------------------------------------------------------------
//  Construction
// ----------------------------------------------------------------------------
PreMesh PreMesh::from_soup(std::vector<Vec3> pos, std::vector<std::vector<int>> f) {
    PreMesh p;
    p.positions = std::move(pos);
    p.faces     = std::move(f);
    return p;
}

PreMesh PreMesh::from_mesh(const Mesh& m) {
    PreMesh p;
    std::vector<int> remap(m.n_vertices(), -1);
    for (auto v : m.all_vertices()) {
        remap[v.idx()] = static_cast<int>(p.positions.size());
        p.positions.push_back(m.point(v));
    }
    for (auto fh : m.all_faces()) {
        std::vector<int> face;
        for (auto v : m.face_vertices(fh)) face.push_back(remap[v.idx()]);
        p.faces.push_back(std::move(face));
    }
    return p;
}

// ----------------------------------------------------------------------------
//  Radial map (the radial cycle, flattened)
// ----------------------------------------------------------------------------
int PreMesh::build_radial() {
    radial_.clear();
    for (int f = 0; f < static_cast<int>(faces.size()); ++f) {
        const auto& fv = faces[f];
        const int n = static_cast<int>(fv.size());
        for (int s = 0; s < n; ++s)
            radial_[ekey(fv[s], fv[(s + 1) % n])].push_back({f, s});
    }
    radial_built_ = true;
    int nm = 0;
    for (auto& kv : radial_) if (kv.second.size() >= 3) ++nm;
    return nm;
}

int PreMesh::edge_valence(int a, int b) const {
    auto it = radial_.find(ekey(a, b));
    return it == radial_.end() ? 0 : static_cast<int>(it->second.size());
}

// ----------------------------------------------------------------------------
//  Weld colocated vertices - neighbor-aware spatial hash, O(V).
// ----------------------------------------------------------------------------
//  Cells are tol-sized; a representative is matched by scanning the vertex's own
//  cell AND its 26 neighbours, so two points closer than tol always merge even
//  when they straddle a cell boundary (the bug in a single-cell lookup). Uses an
//  unordered_map, not std::map, so this is genuinely O(V) not O(V log V).
// ----------------------------------------------------------------------------
int PreMesh::weld_(double tol) {
    if (tol <= 0) return 0;
    const double inv = 1.0 / tol, tol2 = tol * tol;
    auto cell_of = [&](const Vec3& p) {
        return Cell{ (int)std::floor(p.x * inv), (int)std::floor(p.y * inv), (int)std::floor(p.z * inv) };
    };
    std::unordered_map<Cell, std::vector<int>, CellHash> grid;
    std::vector<int> vmap(positions.size());
    int merged = 0;
    for (int i = 0; i < static_cast<int>(positions.size()); ++i) {
        const Vec3& p = positions[i];
        Cell c = cell_of(p);
        int found = -1;
        for (int dz = -1; dz <= 1 && found < 0; ++dz)
        for (int dy = -1; dy <= 1 && found < 0; ++dy)
        for (int dx = -1; dx <= 1 && found < 0; ++dx) {
            auto it = grid.find(Cell{ c.x + dx, c.y + dy, c.z + dz });
            if (it == grid.end()) continue;
            for (int r : it->second)
                if ((positions[r] - p).sqrnorm() <= tol2) { found = r; break; }
        }
        if (found >= 0) { vmap[i] = found; ++merged; }
        else { grid[c].push_back(i); vmap[i] = i; }
    }
    for (auto& f : faces) for (int& v : f) v = vmap[v];
    radial_built_ = false;
    return merged;
}

// ----------------------------------------------------------------------------
//  Drop degenerate faces (< 3 / repeated vertex / area below min_area).
//  min_area > 0 also removes thin slivers and caps that survive exact tests.
// ----------------------------------------------------------------------------
int PreMesh::drop_degenerate_(double min_area) {
    const double area_thresh = std::max(min_area, 1e-18);
    std::vector<std::vector<int>> keep;
    int removed = 0;
    for (auto& f : faces) {
        bool bad = f.size() < 3;
        for (size_t i = 0; i < f.size() && !bad; ++i)
            for (size_t j = i + 1; j < f.size(); ++j)
                if (f[i] == f[j]) { bad = true; break; }
        if (!bad) {                              // area test (triangulated fan)
            double area = 0;
            for (size_t i = 1; i + 1 < f.size(); ++i)
                area += 0.5 * (positions[f[i]] - positions[f[0]])
                                  .cross(positions[f[i + 1]] - positions[f[0]]).norm();
            if (area < area_thresh) bad = true;
        }
        if (bad) ++removed; else keep.push_back(std::move(f));
    }
    faces.swap(keep);
    radial_built_ = false;
    return removed;
}

// ----------------------------------------------------------------------------
//  Drop duplicate faces (same vertex set, any winding). unordered_map.
// ----------------------------------------------------------------------------
int PreMesh::drop_duplicate_() {
    std::unordered_map<std::vector<int>, char, VecHash> seen;
    seen.reserve(faces.size() * 2);
    std::vector<std::vector<int>> keep;
    int removed = 0;
    for (auto& f : faces) {
        std::vector<int> key = f;
        std::sort(key.begin(), key.end());
        if (seen.count(key)) { ++removed; continue; }
        seen.emplace(std::move(key), 1);
        keep.push_back(std::move(f));
    }
    faces.swap(keep);
    radial_built_ = false;
    return removed;
}

// ----------------------------------------------------------------------------
//  Coherent orientation via BFS over manifold (radial==2) adjacency.
// ----------------------------------------------------------------------------
int PreMesh::reorient_() {
    if (!radial_built_) build_radial();
    const int nf = static_cast<int>(faces.size());
    std::vector<char> visited(nf, 0), flip(nf, 0);
    int flips = 0;

    auto directed = [&](int f, int s, int& da, int& db) {
        const auto& fv = faces[f];
        int a = fv[s], b = fv[(s + 1) % fv.size()];
        if (flip[f]) std::swap(a, b);
        da = a; db = b;
    };

    for (int seed = 0; seed < nf; ++seed) {
        if (visited[seed]) continue;
        visited[seed] = 1;
        std::vector<int> stack{ seed };
        while (!stack.empty()) {
            int f = stack.back(); stack.pop_back();
            const auto& fv = faces[f];
            for (int s = 0; s < static_cast<int>(fv.size()); ++s) {
                auto it = radial_.find(ekey(fv[s], fv[(s + 1) % fv.size()]));
                if (it == radial_.end() || it->second.size() != 2) continue;  // boundary/non-manifold
                int da, db; directed(f, s, da, db);
                for (auto& c : it->second) {
                    if (c.face == f || visited[c.face]) continue;
                    int ga, gb; directed(c.face, c.slot, ga, gb);
                    bool consistent = (da == gb && db == ga);
                    visited[c.face] = 1;
                    if (!consistent) { flip[c.face] = 1; ++flips; }
                    stack.push_back(c.face);
                }
            }
        }
    }
    for (int f = 0; f < nf; ++f) if (flip[f]) std::reverse(faces[f].begin(), faces[f].end());
    radial_built_ = false;
    return flips;
}

// ----------------------------------------------------------------------------
//  Separate non-manifold edges and bow-tie vertices into manifold fans.
// ----------------------------------------------------------------------------
//  This separates the loop-fans meeting at a
//  vertex and re-joins ("splices") only those that share a consistent manifold
//  edge. We do it globally with union-find over face-CORNERS:
//
//    * a corner is a (face, slot); it references one vertex faces[face][slot].
//    * across every MANIFOLD edge (exactly two corners, oppositely wound) we
//      UNION the two corners that reference each shared endpoint.
//    * connected components of corners = the output vertices.
//
//  A bow-tie vertex's fans never get unioned -> they become separate vertices
//  (dissociation). The faces around a non-manifold edge (radial >= 3) aren't
//  unioned across it -> the edge separates into manifold sheets. After this pass
//  every edge has <= 2 faces and every vertex is a single fan, so emit can't fail.
// ----------------------------------------------------------------------------
namespace {
struct DSU {
    std::vector<int> p;
    explicit DSU(int n) : p(n) { for (int i = 0; i < n; ++i) p[i] = i; }
    int find(int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
    void join(int a, int b) { p[find(a)] = find(b); }
};
} // namespace

int PreMesh::split_fans_() {
    if (!radial_built_) build_radial();

    // corner ids: base[f] + slot
    std::vector<int> base(faces.size() + 1, 0);
    for (size_t f = 0; f < faces.size(); ++f)
        base[f + 1] = base[f] + static_cast<int>(faces[f].size());
    const int total = base.back();
    if (total == 0) return 0;
    DSU dsu(total);
    auto cid = [&](int f, int s) { return base[f] + s; };

    for (auto& kv : radial_) {
        auto& cs = kv.second;
        if (cs.size() != 2) continue;                    // boundary / non-manifold -> no union
        int f0 = cs[0].face, s0 = cs[0].slot, n0 = static_cast<int>(faces[f0].size());
        int f1 = cs[1].face, s1 = cs[1].slot, n1 = static_cast<int>(faces[f1].size());
        int a0 = faces[f0][s0], b0 = faces[f0][(s0 + 1) % n0];
        int a1 = faces[f1][s1], b1 = faces[f1][(s1 + 1) % n1];
        if (a0 == b1 && b0 == a1) {                      // consistently wound -> splice
            dsu.join(cid(f0, s0),            cid(f1, (s1 + 1) % n1));   // endpoint a0==b1
            dsu.join(cid(f0, (s0 + 1) % n0), cid(f1, s1));             // endpoint b0==a1
        }
    }

    // component -> new vertex; position copied from the corner's original vertex
    std::vector<int> comp(total, -1);
    std::vector<Vec3> np;
    for (size_t f = 0; f < faces.size(); ++f)
        for (int s = 0; s < static_cast<int>(faces[f].size()); ++s) {
            int r = dsu.find(cid(static_cast<int>(f), s));
            if (comp[r] < 0) { comp[r] = static_cast<int>(np.size()); np.push_back(positions[faces[f][s]]); }
        }
    // count distinct original vertices used (to report how many copies were added)
    std::vector<char> usedOrig(positions.size(), 0);
    for (auto& f : faces) for (int v : f) usedOrig[v] = 1;
    int used = 0; for (char c : usedOrig) used += c;

    for (size_t f = 0; f < faces.size(); ++f)
        for (int s = 0; s < static_cast<int>(faces[f].size()); ++s)
            faces[f][s] = comp[dsu.find(cid(static_cast<int>(f), s))];

    const int newV = static_cast<int>(np.size());
    positions.swap(np);
    radial_built_ = false;
    return newV - used;                                  // extra copies created
}

// ----------------------------------------------------------------------------
//  Strip tiny islands: drop faces in connected components smaller than
//  min_faces (scan/print debris). Returns the number of components removed.
// ----------------------------------------------------------------------------
int PreMesh::drop_small_components_(int min_faces) {
    if (min_faces <= 1 || faces.empty()) return 0;
    const int nf = static_cast<int>(faces.size());
    DSU dsu(nf);
    std::unordered_map<uint64_t, int> edge2face;
    edge2face.reserve(nf * 3);
    for (int f = 0; f < nf; ++f) {
        const auto& fv = faces[f]; int n = static_cast<int>(fv.size());
        for (int s = 0; s < n; ++s) {
            uint64_t k = ekey(fv[s], fv[(s + 1) % n]);
            auto it = edge2face.find(k);
            if (it == edge2face.end()) edge2face.emplace(k, f);
            else dsu.join(f, it->second);
        }
    }
    std::unordered_map<int, int> size;             // root -> face count
    for (int f = 0; f < nf; ++f) size[dsu.find(f)]++;
    std::vector<std::vector<int>> keep;
    std::unordered_map<int, int> removed_roots;
    for (int f = 0; f < nf; ++f) {
        int root = dsu.find(f);
        if (size[root] < min_faces) removed_roots[root]++;
        else keep.push_back(std::move(faces[f]));
    }
    faces.swap(keep);
    radial_built_ = false;
    return static_cast<int>(removed_roots.size());
}

// ----------------------------------------------------------------------------
//  Drop unused vertices
// ----------------------------------------------------------------------------
int PreMesh::drop_unused_() {
    std::vector<int> idx(positions.size(), -1);
    for (auto& f : faces) for (int v : f) idx[v] = 0;
    std::vector<Vec3> np;
    int removed = 0;
    for (int i = 0; i < static_cast<int>(positions.size()); ++i) {
        if (idx[i] == 0) { idx[i] = static_cast<int>(np.size()); np.push_back(positions[i]); }
        else ++removed;
    }
    for (auto& f : faces) for (int& v : f) v = idx[v];
    positions.swap(np);
    radial_built_ = false;
    return removed;
}

// ----------------------------------------------------------------------------
//  Driver
// ----------------------------------------------------------------------------
RepairReport PreMesh::repair(const RepairOptions& opt) {
    RepairReport r;
    if (opt.weld_vertices)        r.vertices_merged    = weld_(opt.vertex_merge_tol);
    if (opt.remove_degenerate || opt.min_face_area > 0)
        r.degenerate_removed = drop_degenerate_(opt.min_face_area);
    if (opt.remove_duplicate)     r.duplicate_removed  = drop_duplicate_();
    if (opt.min_component_faces > 1) r.components_removed = drop_small_components_(opt.min_component_faces);
    r.nm_edges = build_radial();
    if (opt.reorient)          r.faces_flipped      = reorient_();
    if (opt.split_nonmanifold) r.vertices_split     = split_fans_();
    if (opt.remove_unused)     r.vertices_removed    = drop_unused_();
    return r;
}

Mesh PreMesh::to_mesh(const RepairOptions& opt, RepairReport* report) const {
    Mesh out;
    std::vector<VertexHandle> vh(positions.size());
    for (size_t i = 0; i < positions.size(); ++i) vh[i] = out.add_vertex(positions[i]);

    int split = 0, failed = 0;
    for (const auto& f : faces) {
        std::vector<VertexHandle> fv;
        fv.reserve(f.size());
        for (int v : f) fv.push_back(vh[v]);
        if (out.add_face(fv).is_valid()) continue;

        // Rejected (non-manifold edge / bow-tie). Re-emit on its own vertices so
        // it survives as a separate sheet rather than being dropped.
        if (opt.split_nonmanifold) {
            std::vector<VertexHandle> dup;
            dup.reserve(f.size());
            for (int v : f) dup.push_back(out.add_vertex(positions[v]));
            if (out.add_face(dup).is_valid()) ++split; else ++failed;
        } else ++failed;
    }
    if (report) { report->faces_split = split; report->faces_failed = failed; }
    return out;
}

// ============================================================================
//  Soup loaders (parse file -> soup; no halfedge build, nothing dropped)
// ============================================================================
namespace {

template <class T> T read_le(std::istream& is) { T v; is.read(reinterpret_cast<char*>(&v), sizeof(T)); return v; }

bool load_stl_soup(const std::string& path, PreMesh& p) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    in.seekg(0, std::ios::end); std::streamsize size = in.tellg(); in.seekg(0);

    // binary STL is exactly 84 + 50*count bytes
    char header[80] = {0};
    in.read(header, 80);
    uint32_t count = read_le<uint32_t>(in);
    bool binary = (std::streamsize)(84 + 50ull * count) == size;

    auto add_tri = [&](const Vec3& a, const Vec3& b, const Vec3& c) {
        int base = static_cast<int>(p.positions.size());
        p.positions.push_back(a); p.positions.push_back(b); p.positions.push_back(c);
        p.faces.push_back({ base, base + 1, base + 2 });   // weld() will merge shared verts
    };

    if (binary) {
        for (uint32_t i = 0; i < count; ++i) {
            read_le<float>(in); read_le<float>(in); read_le<float>(in);   // normal
            Vec3 v[3];
            for (int k = 0; k < 3; ++k)
                v[k] = Vec3(read_le<float>(in), read_le<float>(in), read_le<float>(in));
            read_le<uint16_t>(in);                                         // attr
            add_tri(v[0], v[1], v[2]);
        }
        return true;
    }
    // ASCII
    in.clear(); in.seekg(0);
    std::string tok; std::vector<Vec3> v;
    while (in >> tok) {
        if (tok == "vertex") {
            double x, y, z; in >> x >> y >> z; v.push_back(Vec3(x, y, z));
            if (v.size() == 3) { add_tri(v[0], v[1], v[2]); v.clear(); }
        }
    }
    return !p.faces.empty();
}

bool load_obj_soup(const std::string& path, PreMesh& p) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream ls(line); std::string tag; ls >> tag;
        if (tag == "v") { double x, y, z; ls >> x >> y >> z; p.positions.push_back(Vec3(x, y, z)); }
        else if (tag == "f") {
            std::vector<int> fv; std::string t;
            while (ls >> t) {
                int vi = std::stoi(t.substr(0, t.find('/')));
                if (vi < 0) vi = static_cast<int>(p.positions.size()) + vi; else vi -= 1;
                fv.push_back(vi);
            }
            if (fv.size() >= 3) p.faces.push_back(std::move(fv));
        }
    }
    return true;
}

bool load_off_soup(const std::string& path, PreMesh& p) {
    std::ifstream in(path);
    if (!in) return false;
    std::string tok; in >> tok;
    if (tok.rfind("OFF", 0) != 0) return false;
    if (tok != "OFF") return false;                 // ascii OFF only
    int nv, nf, ne; in >> nv >> nf >> ne;
    for (int i = 0; i < nv; ++i) { double x, y, z; in >> x >> y >> z; p.positions.push_back(Vec3(x, y, z)); }
    for (int i = 0; i < nf; ++i) {
        int k; in >> k; std::vector<int> fv(k);
        for (int j = 0; j < k; ++j) in >> fv[j];
        if (k >= 3) p.faces.push_back(std::move(fv));
    }
    return true;
}

std::string ext_of(const std::string& path) {
    auto dot = path.find_last_of('.');
    std::string e = dot == std::string::npos ? "" : path.substr(dot + 1);
    for (auto& c : e) c = static_cast<char>(std::tolower(c));
    return e;
}

} // namespace

PreMesh load_soup(const std::string& path) {
    PreMesh p;
    const std::string e = ext_of(path);
    if      (e == "stl") load_stl_soup(path, p);
    else if (e == "obj") load_obj_soup(path, p);
    else if (e == "off") load_off_soup(path, p);
    // (ply soup loader is a follow-up; use read_ply then PreMesh::from_mesh.)
    return p;
}

Mesh load_and_repair(const std::string& path, const RepairOptions& opt) {
    PreMesh p = load_soup(path);
    RepairReport r = p.repair(opt);
    return p.to_mesh(opt, &r);
}

} // namespace sm
