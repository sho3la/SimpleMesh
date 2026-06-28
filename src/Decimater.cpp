// ============================================================================
//  SimpleMesh - Decimater.cpp : the modular decimation engine + stock modules
// ============================================================================
#include "simplemesh/Decimater.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>

namespace sm {

// ============================================================================
//  Small geometry helper: unit normal of a triangle from three points.
// ============================================================================
static Vec3 tri_normal(const Vec3& a, const Vec3& b, const Vec3& c) {
    return (b - a).cross(c - a).normalized();
}

// ============================================================================
//  ModQuadric (priority) - per-vertex Garland-Heckbert error quadrics.
// ============================================================================
void ModQuadric::initialize(Mesh& m) {
    Q_.assign(m.n_vertices(), Quadric());
    for (auto f : m.all_faces()) {
        Vec3 n = m.calc_face_normal(f);
        auto fv = m.face_vertices(f);
        if (fv.empty()) continue;
        const Vec3& p0 = m.point(fv[0]);
        double pd = -(n.x * p0.x + n.y * p0.y + n.z * p0.z);
        Quadric K = Quadric::from_plane(n.x, n.y, n.z, pd);
        for (auto v : fv) Q_[v.idx()] = Q_[v.idx()] + K;
    }
}

float ModQuadric::collapse_priority(const Mesh&, const CollapseInfo& ci) {
    // Error of merging both quadrics, evaluated at the surviving position p1.
    Quadric q = Q_[ci.v0.idx()] + Q_[ci.v1.idx()];
    double err = q.eval(ci.p1);
    if (err < 0) err = 0;                         // guard tiny negatives
    return static_cast<float>(err);
}

void ModQuadric::postprocess_collapse(Mesh&, const CollapseInfo& ci) {
    Q_[ci.v1.idx()] = Q_[ci.v1.idx()] + Q_[ci.v0.idx()];   // accumulate onto v1
}

// ============================================================================
//  ModEdgeLength (binary) - veto if any edge newly incident to v1 is too long.
// ============================================================================
float ModEdgeLength::collapse_priority(const Mesh& m, const CollapseInfo& ci) {
    // After the collapse, v1 inherits v0's neighbours: check those edge lengths.
    for (VertexHandle w : m.vertex_vertices(ci.v0)) {
        if (w == ci.v1) continue;
        if ((m.point(w) - ci.p1).norm() > max_length_) return ILLEGAL_COLLAPSE;
    }
    return LEGAL_COLLAPSE;
}

// ============================================================================
//  ModNormalFlipping (binary) - veto collapses that fold the surface.
// ----------------------------------------------------------------------------
//  For every triangle around v0 that SURVIVES the collapse (i.e. is not one of
//  the <=2 triangles on the collapsed edge), compare its normal now vs. after
//  v0 is moved onto p1. A large swing means a foldover.
// ============================================================================
float ModNormalFlipping::collapse_priority(const Mesh& m, const CollapseInfo& ci) {
    const double cos_max = std::cos(max_angle_);
    FaceHandle fl = m.face(ci.heh);
    FaceHandle fr = m.face(m.opposite_halfedge(ci.heh));

    for (HalfedgeHandle oh : m.vertex_outgoing_halfedges(ci.v0)) {
        FaceHandle f = m.face(oh);
        if (!f.is_valid() || f == fl || f == fr) continue;   // boundary / dies
        auto fv = m.face_vertices(f);
        if (fv.size() != 3) continue;

        // the triangle's three points, with v0 swapped for p1 in the "after"
        Vec3 a = m.point(fv[0]), b = m.point(fv[1]), c = m.point(fv[2]);
        Vec3 a2 = (fv[0] == ci.v0) ? ci.p1 : a;
        Vec3 b2 = (fv[1] == ci.v0) ? ci.p1 : b;
        Vec3 c2 = (fv[2] == ci.v0) ? ci.p1 : c;

        Vec3 n_old = tri_normal(a, b, c);
        Vec3 n_new = tri_normal(a2, b2, c2);
        if (n_old.dot(n_new) < cos_max) return ILLEGAL_COLLAPSE;
    }
    return LEGAL_COLLAPSE;
}

// ============================================================================
//  ModAspectRatio (binary) - veto collapses that create sliver triangles.
//  aspect = longest_edge / shortest_altitude = longest_edge^2 / (2*area).
// ============================================================================
float ModAspectRatio::collapse_priority(const Mesh& m, const CollapseInfo& ci) {
    FaceHandle fl = m.face(ci.heh);
    FaceHandle fr = m.face(m.opposite_halfedge(ci.heh));

    for (HalfedgeHandle oh : m.vertex_outgoing_halfedges(ci.v0)) {
        FaceHandle f = m.face(oh);
        if (!f.is_valid() || f == fl || f == fr) continue;
        auto fv = m.face_vertices(f);
        if (fv.size() != 3) continue;

        Vec3 p[3];
        for (int k = 0; k < 3; ++k) p[k] = (fv[k] == ci.v0) ? ci.p1 : m.point(fv[k]);

        double e0 = (p[1] - p[0]).norm();
        double e1 = (p[2] - p[1]).norm();
        double e2 = (p[0] - p[2]).norm();
        double longest = std::max(e0, std::max(e1, e2));
        double area = 0.5 * (p[1] - p[0]).cross(p[2] - p[0]).norm();
        if (area < 1e-20) return ILLEGAL_COLLAPSE;     // degenerate
        double aspect = longest * longest / (2.0 * area);
        if (aspect > max_aspect_) return ILLEGAL_COLLAPSE;
    }
    return LEGAL_COLLAPSE;
}

// ============================================================================
//  Decimator engine
// ============================================================================

CollapseInfo Decimator::make_info(HalfedgeHandle h) const {
    CollapseInfo ci;
    ci.heh = h;
    ci.v0 = mesh_.from_vertex(h);
    ci.v1 = mesh_.to_vertex(h);
    ci.p0 = mesh_.point(ci.v0);
    ci.p1 = mesh_.point(ci.v1);
    return ci;
}

bool Decimator::passes_binary(const CollapseInfo& ci) const {
    for (const auto& mod : modules_)
        if (mod->is_binary() &&
            mod->collapse_priority(mesh_, ci) == ILLEGAL_COLLAPSE)
            return false;
    return true;
}

bool Decimator::best_collapse(EdgeHandle e, float& cost_out, HalfedgeHandle& heh_out) {
    bool found = false;
    float best = std::numeric_limits<float>::max();
    for (int i = 0; i < 2; ++i) {
        HalfedgeHandle h = mesh_.halfedge(e, i);
        if (!mesh_.is_collapse_ok(h)) continue;
        CollapseInfo ci = make_info(h);
        if (!passes_binary(ci)) continue;
        float c = priority_->collapse_priority(mesh_, ci);
        if (c == ILLEGAL_COLLAPSE) continue;
        if (c < best) { best = c; heh_out = h; found = true; }
    }
    cost_out = best;
    return found;
}

namespace {
struct DEntry { float cost; int edge; int version; };
struct DGreater { bool operator()(const DEntry& a, const DEntry& b) const { return a.cost > b.cost; } };
} // namespace

size_t Decimator::decimate_to_faces(size_t target_faces) {
    // locate the single priority module
    priority_ = nullptr;
    for (auto& mod : modules_) if (!mod->is_binary()) { priority_ = mod.get(); break; }
    if (!priority_) return 0;                          // nothing to score by

    for (auto& mod : modules_) mod->initialize(mesh_);

    std::vector<int> version(mesh_.n_edges(), 0);
    std::priority_queue<DEntry, std::vector<DEntry>, DGreater> pq;

    for (auto e : mesh_.all_edges()) {
        float cost; HalfedgeHandle h;
        if (best_collapse(e, cost, h))
            pq.push({ cost, e.idx(), version[e.idx()] });
    }

    size_t live_faces = 0;
    for (auto f : mesh_.all_faces()) { (void)f; ++live_faces; }

    size_t n_collapses = 0;
    while (live_faces > target_faces && !pq.empty()) {
        DEntry top = pq.top(); pq.pop();
        EdgeHandle e(top.edge);
        if (mesh_.is_deleted(e) || version[e.idx()] != top.version) continue;   // stale

        // Re-derive the cheapest legal orientation now (geometry may have moved).
        float cost; HalfedgeHandle h;
        if (!best_collapse(e, cost, h)) continue;

        CollapseInfo ci = make_info(h);
        size_t removed = (mesh_.face(h).is_valid() ? 1 : 0) +
                         (mesh_.face(mesh_.opposite_halfedge(h)).is_valid() ? 1 : 0);

        mesh_.collapse(h);
        for (auto& mod : modules_) mod->postprocess_collapse(mesh_, ci);
        live_faces -= removed;
        ++n_collapses;

        // refresh every edge still incident to the surviving vertex v1
        for (HalfedgeHandle oh : mesh_.vertex_outgoing_halfedges(ci.v1)) {
            EdgeHandle ie = mesh_.edge(oh);
            ++version[ie.idx()];
            float c2; HalfedgeHandle h2;
            if (best_collapse(ie, c2, h2))
                pq.push({ c2, ie.idx(), version[ie.idx()] });
        }
    }

    mesh_.garbage_collection();
    return n_collapses;
}

} // namespace sm
