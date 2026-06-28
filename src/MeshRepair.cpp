// ============================================================================
//  SimpleMesh - MeshRepair.cpp : rebuild-based mesh cleanup
// ============================================================================
#include "simplemesh/MeshRepair.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace sm {

std::string RepairReport::summary() const {
    std::ostringstream os;
    os << "MeshRepair: " << (changed() ? "changes applied\n" : "nothing to fix\n");
    os << "  vertices merged:   " << vertices_merged    << "\n";
    os << "  faces flipped:     " << faces_flipped      << "\n";
    os << "  degenerate faces:  " << degenerate_removed << "\n";
    os << "  duplicate faces:   " << duplicate_removed  << "\n";
    os << "  vertices removed:  " << vertices_removed   << "\n";
    if (vertices_split)     os << "  vertices split:    " << vertices_split     << "\n";
    if (components_removed) os << "  components removed:" << components_removed  << "\n";
    if (faces_split)        os << "  faces split:       " << faces_split        << "\n";
    if (faces_failed)       os << "  faces rejected:    " << faces_failed       << "\n";
    return os.str();
}

namespace {

// FNV-1a hashers (replace the std::map/std::tuple hot paths with O(1) lookups).
struct Cell { int x, y, z; bool operator==(const Cell& o) const { return x == o.x && y == o.y && z == o.z; } };
struct CellHash { size_t operator()(const Cell& c) const {
    size_t h = 1469598103934665603ULL;
    for (int v : { c.x, c.y, c.z }) { h ^= (size_t)(unsigned)v; h *= 1099511628211ULL; }
    return h; } };
struct VecHash { size_t operator()(const std::vector<int>& v) const {
    size_t h = 1469598103934665603ULL;
    for (int x : v) { h ^= (size_t)(unsigned)x; h *= 1099511628211ULL; }
    return h; } };

// ----- the cleaned soup we rebuild from ------------------------------------
struct Soup {
    std::vector<Vec3>             pos;     // vertex positions
    std::vector<std::vector<int>> faces;   // faces as vertex-index lists
};

// Extract live geometry + connectivity from the mesh as a flat soup.
Soup extract(const Mesh& m) {
    Soup s;
    std::vector<int> remap(m.n_vertices(), -1);   // old vertex idx -> soup idx
    for (auto v : m.all_vertices()) {
        remap[v.idx()] = static_cast<int>(s.pos.size());
        s.pos.push_back(m.point(v));
    }
    for (auto f : m.all_faces()) {
        std::vector<int> face;
        for (auto v : m.face_vertices(f)) face.push_back(remap[v.idx()]);
        s.faces.push_back(std::move(face));
    }
    return s;
}

// ----- 1. weld colocated vertices (neighbor-aware spatial hash) -------------
// tol-sized cells; a vertex matches a representative in its own cell OR any of
// the 26 neighbours, so points closer than tol always merge even across a cell
// boundary. unordered_map keeps it O(V). Fills vmap (soup idx -> representative).
int weld(const Soup& s, double tol, std::vector<int>& vmap) {
    vmap.assign(s.pos.size(), -1);
    if (tol <= 0) { for (size_t i = 0; i < s.pos.size(); ++i) vmap[i] = (int)i; return 0; }

    const double inv = 1.0 / tol, tol2 = tol * tol;
    auto cell_of = [&](const Vec3& p) {
        return Cell{ (int)std::floor(p.x * inv), (int)std::floor(p.y * inv), (int)std::floor(p.z * inv) };
    };
    std::unordered_map<Cell, std::vector<int>, CellHash> grid;
    int merged = 0;
    for (size_t i = 0; i < s.pos.size(); ++i) {
        const Vec3& p = s.pos[i];
        Cell c = cell_of(p);
        int found = -1;
        for (int dz = -1; dz <= 1 && found < 0; ++dz)
        for (int dy = -1; dy <= 1 && found < 0; ++dy)
        for (int dx = -1; dx <= 1 && found < 0; ++dx) {
            auto it = grid.find(Cell{ c.x + dx, c.y + dy, c.z + dz });
            if (it == grid.end()) continue;
            for (int r : it->second)
                if ((s.pos[r] - p).sqrnorm() <= tol2) { found = r; break; }
        }
        if (found >= 0) { vmap[i] = found; ++merged; }
        else { grid[c].push_back((int)i); vmap[i] = (int)i; }
    }
    return merged;
}

// ----- 2. reorient facets coherently ---------------------------------------
// BFS over face adjacency (faces sharing an undirected edge). Two faces are
// consistently oriented when they traverse the shared edge in OPPOSITE
// directions; if not, the neighbour is flipped. Seeds one face per component.
// Returns the number of faces flipped. (Edges shared by != 2 faces are skipped -
// non-manifold seams can't be coherently oriented.)
int reorient(std::vector<std::vector<int>>& faces) {
    // Canonical undirected-edge key. Note: std::minmax returns a pair of
    // *references* to its arguments, so returning it from a lambda would dangle;
    // build the pair from values instead.
    auto ekey = [](int a, int b) {
        return std::make_pair(std::min(a, b), std::max(a, b));
    };
    // undirected edge -> the (face, directed a->b) records touching it
    std::map<std::pair<int,int>, std::vector<std::tuple<int,int,int>>> edge2;
    for (int f = 0; f < (int)faces.size(); ++f) {
        const auto& fv = faces[f];
        for (size_t i = 0; i < fv.size(); ++i) {
            int a = fv[i], b = fv[(i + 1) % fv.size()];
            edge2[ekey(a, b)].emplace_back(f, a, b);
        }
    }

    std::vector<char> visited(faces.size(), 0), flip(faces.size(), 0);
    int flips = 0;
    for (int seed = 0; seed < (int)faces.size(); ++seed) {
        if (visited[seed]) continue;
        visited[seed] = 1;
        std::vector<int> stack{ seed };
        while (!stack.empty()) {
            int f = stack.back(); stack.pop_back();
            const auto& fv = faces[f];
            for (size_t i = 0; i < fv.size(); ++i) {
                int a = fv[i], b = fv[(i + 1) % fv.size()];
                // f currently traverses a->b (accounting for its own flip)
                int da = a, db = b;
                if (flip[f]) std::swap(da, db);
                auto& recs = edge2[ekey(a, b)];
                if (recs.size() != 2) continue;            // boundary or non-manifold
                for (auto& rec : recs) {
                    int g = std::get<0>(rec);
                    if (g == f) continue;
                    int ga = std::get<1>(rec), gb = std::get<2>(rec);
                    if (flip[g]) std::swap(ga, gb);
                    bool consistent = (da == gb && db == ga);  // opposite traversal
                    if (!visited[g]) {
                        visited[g] = 1;
                        if (!consistent) { flip[g] = 1; ++flips; }
                        stack.push_back(g);
                    }
                }
            }
        }
    }
    for (int f = 0; f < (int)faces.size(); ++f)
        if (flip[f]) std::reverse(faces[f].begin(), faces[f].end());
    return flips;
}

} // namespace

// ============================================================================
//  Driver
// ============================================================================
RepairReport repair_mesh(Mesh& m, const RepairOptions& opt) {
    RepairReport rep;
    Soup s = extract(m);

    // 1. weld -------------------------------------------------------------
    if (opt.weld_vertices) {
        std::vector<int> vmap;
        rep.vertices_merged = weld(s, opt.vertex_merge_tol, vmap);
        for (auto& f : s.faces)
            for (int& v : f) v = vmap[v];
    }

    // 3. drop degenerate facets (also catches faces that collapsed after weld)
    if (opt.remove_degenerate) {
        std::vector<std::vector<int>> keep;
        for (auto& f : s.faces) {
            bool bad = f.size() < 3;
            for (size_t i = 0; i < f.size() && !bad; ++i)
                for (size_t j = i + 1; j < f.size(); ++j)
                    if (f[i] == f[j]) { bad = true; break; }
            if (bad) ++rep.degenerate_removed; else keep.push_back(std::move(f));
        }
        s.faces.swap(keep);
    }

    // 4. drop duplicate facets (same vertex set, any winding) -------------
    if (opt.remove_duplicate) {
        std::unordered_map<std::vector<int>, char, VecHash> seen;
        seen.reserve(s.faces.size() * 2);
        std::vector<std::vector<int>> keep;
        for (auto& f : s.faces) {
            std::vector<int> key = f;
            std::sort(key.begin(), key.end());
            if (seen.count(key)) { ++rep.duplicate_removed; continue; }
            seen.emplace(std::move(key), 1);
            keep.push_back(std::move(f));
        }
        s.faces.swap(keep);
    }

    // 2. reorient (after dropping junk so adjacency is clean) -------------
    if (opt.reorient)
        rep.faces_flipped = reorient(s.faces);

    // 5. drop unused vertices --------------------------------------------
    std::vector<int> final_index(s.pos.size(), -1);
    if (opt.remove_unused) {
        for (auto& f : s.faces) for (int v : f) final_index[v] = 0;   // mark used
        int next = 0;
        for (size_t i = 0; i < s.pos.size(); ++i)
            if (final_index[i] == 0) final_index[i] = next++;
            else { final_index[i] = -1; ++rep.vertices_removed; }
    } else {
        for (size_t i = 0; i < s.pos.size(); ++i) final_index[i] = (int)i;
    }

    // ---- rebuild --------------------------------------------------------
    Mesh out;
    std::vector<VertexHandle> vh(s.pos.size());
    for (size_t i = 0; i < s.pos.size(); ++i)
        if (final_index[i] >= 0) vh[i] = out.add_vertex(s.pos[i]);
    for (auto& f : s.faces) {
        std::vector<VertexHandle> fv;
        fv.reserve(f.size());
        for (int v : f) fv.push_back(vh[v]);
        if (!out.add_face(fv).is_valid()) ++rep.faces_failed;
    }

    m = std::move(out);   // Mesh is move-only; replace in place
    return rep;
}

} // namespace sm
