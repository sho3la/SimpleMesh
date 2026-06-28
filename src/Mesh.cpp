// ============================================================================
//  SimpleMesh - Mesh.cpp : implementation of the algorithmic mesh operations
// ============================================================================
#include "simplemesh/Mesh.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <limits>
#include <sstream>
#include <iostream>

namespace sm {

// Silent by default - see Mesh.h. Bulk-importing messy real-world meshes
// would otherwise print one line per rejected non-manifold face.
bool Mesh::warn_invalid_face = false;

// ============================================================================
//  Element factories
// ============================================================================

VertexHandle Mesh::new_vertex(const Vec3& p) {
    vertices_.push_back(Vertex());
    points_.push_back(p);
    v_deleted_.push_back(false);
    vprops_.push_back();                // keep vertex properties in lock-step
    return VertexHandle(static_cast<int>(vertices_.size()) - 1);
}

VertexHandle Mesh::add_vertex(const Vec3& p) { return new_vertex(p); }

// Create a fresh edge = two opposite halfedges, appended as a 2e / 2e+1 pair.
// Returns the halfedge that points from -> to.
HalfedgeHandle Mesh::new_edge(VertexHandle from, VertexHandle to) {
    const int base = static_cast<int>(halfedges_.size());
    halfedges_.push_back(Halfedge());   // index base   : from -> to
    halfedges_.push_back(Halfedge());   // index base+1 : to   -> from
    e_deleted_.push_back(false);        // one flag per edge (pair)
    eprops_.push_back();                // one edge-property slot
    hprops_.push_back(); hprops_.push_back();  // two halfedge-property slots
    HalfedgeHandle h0(base), h1(base + 1);
    set_to_vertex(h0, to);
    set_to_vertex(h1, from);
    return h0;
}

FaceHandle Mesh::new_face() {
    faces_.push_back(Face());
    f_deleted_.push_back(false);
    fprops_.push_back();                // keep face properties in lock-step
    return FaceHandle(static_cast<int>(faces_.size()) - 1);
}

// ============================================================================
//  Boundary tests that need to circulate
// ============================================================================

bool Mesh::is_boundary(VertexHandle v) const {
    HalfedgeHandle h = halfedge(v);
    // An isolated vertex (no halfedge) counts as boundary.
    if (!h.is_valid()) return true;
    // Thanks to adjust_outgoing_halfedge, a boundary vertex always stores a
    // boundary outgoing halfedge, so this single check is sufficient.
    return is_boundary(h);
}

bool Mesh::is_boundary(FaceHandle f) const {
    for (HalfedgeHandle h : face_halfedges(f))
        if (is_boundary(opposite_halfedge(h)))
            return true;
    return false;
}

// ============================================================================
//  find_halfedge : search the outgoing halfedges of `from` for one ending at
//  `to`. O(valence), which is fine for mesh construction.
// ============================================================================

HalfedgeHandle Mesh::find_halfedge(VertexHandle from, VertexHandle to) const {
    HalfedgeHandle start = halfedge(from);
    if (!start.is_valid()) return HalfedgeHandle();  // isolated vertex

    HalfedgeHandle h = start;
    do {
        if (to_vertex(h) == to) return h;
        // Rotate to the next outgoing halfedge around `from`:
        //   opposite(h) ends at `from`, so next(opposite(h)) starts at `from`.
        h = next_halfedge(opposite_halfedge(h));
    } while (h != start);

    return HalfedgeHandle();
}

// ============================================================================
//  adjust_outgoing_halfedge : enforce the "boundary vertices point to a
//  boundary halfedge" invariant.
// ============================================================================

void Mesh::adjust_outgoing_halfedge(VertexHandle v) {
    for (HalfedgeHandle h : vertex_outgoing_halfedges(v)) {
        if (is_boundary(h)) {
            set_halfedge(v, h);
            return;
        }
    }
}

// ============================================================================
//  add_face : the core topological operation.
// ----------------------------------------------------------------------------
//  The hard part is
//  that, when a new face is sewn into an existing fan of faces around a vertex,
//  the chain of boundary halfedges has to be re-routed ("relinked") so that the
//  boundary stays a single consistent loop. The four `case` blocks below handle
//  the four combinations of "is this side a brand-new edge or a reused one".
// ============================================================================

FaceHandle Mesh::add_face(const std::vector<VertexHandle>& vhs) {
    const size_t n = vhs.size();
    assert(n > 2 && "a face needs at least 3 vertices");
    if (n < 3) return FaceHandle();

    // Reject DEGENERATE input: a valid polygon never repeats a vertex. Without
    // this guard a face like (v,v,v) builds self-loop edges whose next-cycle
    // never closes, and the face walkers (face_halfedges) would then spin
    // forever, appending until they exhaust memory. Importers rely on this to
    // skip the bad faces in malformed "soup" files.
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j)
            if (vhs[i] == vhs[j]) return FaceHandle();

    // Per-side scratch data, indexed like the input vertices.
    std::vector<HalfedgeHandle> edge_he(n);
    std::vector<bool>           is_new(n);
    std::vector<bool>           needs_adjust(n, false);

    // Deferred (halfedge -> next) assignments. We must NOT change next-pointers
    // while still inspecting the old topology, so we queue them and apply last.
    std::vector<std::pair<HalfedgeHandle, HalfedgeHandle>> next_cache;
    next_cache.reserve(6 * n);

    // --- 1. topology sanity check + locate existing halfedges --------------
    for (size_t i = 0, ii = 1; i < n; ++i, ++ii, ii %= n) {
        if (!is_boundary(vhs[i])) {
            if (warn_invalid_face) std::cerr << "add_face: complex vertex (vertex already fully surrounded)\n";
            return FaceHandle();
        }
        edge_he[i] = find_halfedge(vhs[i], vhs[ii]);
        is_new[i]  = !edge_he[i].is_valid();
        if (!is_new[i] && !is_boundary(edge_he[i])) {
            if (warn_invalid_face) std::cerr << "add_face: complex edge (edge already has two faces)\n";
            return FaceHandle();
        }
    }

    // --- 2. relink existing patches so a free gap opens up ------------------
    for (size_t i = 0, ii = 1; i < n; ++i, ++ii, ii %= n) {
        if (!is_new[i] && !is_new[ii]) {
            HalfedgeHandle inner_prev = edge_he[i];
            HalfedgeHandle inner_next = edge_he[ii];
            if (next_halfedge(inner_prev) != inner_next) {
                // Find a boundary gap to splice the existing patch into. The cap
                // stops a runaway walk on pathological non-manifold input (common
                // in raw scan/print meshes) rather than spinning forever.
                HalfedgeHandle outer_prev    = opposite_halfedge(inner_next);
                HalfedgeHandle boundary_prev = outer_prev;
                size_t guard = halfedges_.size() + 1;
                do {
                    boundary_prev = opposite_halfedge(next_halfedge(boundary_prev));
                } while (!is_boundary(boundary_prev) && guard-- > 0);
                if (!is_boundary(boundary_prev)) return FaceHandle();   // no gap found
                HalfedgeHandle boundary_next = next_halfedge(boundary_prev);

                if (boundary_prev == inner_prev) {
                    if (warn_invalid_face) std::cerr << "add_face: patch re-linking failed (non-manifold)\n";
                    return FaceHandle();
                }

                HalfedgeHandle patch_start = next_halfedge(inner_prev);
                HalfedgeHandle patch_end   = prev_halfedge(inner_next);

                next_cache.emplace_back(boundary_prev, patch_start);
                next_cache.emplace_back(patch_end,   boundary_next);
                next_cache.emplace_back(inner_prev,  inner_next);
            }
        }
    }

    // --- 3. create the edges that don't exist yet ---------------------------
    for (size_t i = 0, ii = 1; i < n; ++i, ++ii, ii %= n)
        if (is_new[i])
            edge_he[i] = new_edge(vhs[i], vhs[ii]);

    // --- 4. create the face record ------------------------------------------
    FaceHandle fh = new_face();
    set_halfedge(fh, edge_he[n - 1]);

    // --- 5. wire up the halfedges -------------------------------------------
    for (size_t i = 0, ii = 1; i < n; ++i, ++ii, ii %= n) {
        VertexHandle   vh         = vhs[ii];
        HalfedgeHandle inner_prev = edge_he[i];
        HalfedgeHandle inner_next = edge_he[ii];

        // 2-bit code: bit0 = "this side is new", bit1 = "next side is new".
        int id = (is_new[i] ? 1 : 0) | (is_new[ii] ? 2 : 0);

        if (id) {
            HalfedgeHandle outer_prev = opposite_halfedge(inner_next);
            HalfedgeHandle outer_next = opposite_halfedge(inner_prev);

            switch (id) {
                case 1: {  // prev side new, next side old
                    HalfedgeHandle boundary_prev = prev_halfedge(inner_next);
                    next_cache.emplace_back(boundary_prev, outer_next);
                    set_halfedge(vh, outer_next);
                    break;
                }
                case 2: {  // prev side old, next side new
                    HalfedgeHandle boundary_next = next_halfedge(inner_prev);
                    next_cache.emplace_back(outer_prev, boundary_next);
                    set_halfedge(vh, boundary_next);
                    break;
                }
                case 3: {  // both sides new
                    if (!halfedge(vh).is_valid()) {
                        set_halfedge(vh, outer_next);
                        next_cache.emplace_back(outer_prev, outer_next);
                    } else {
                        HalfedgeHandle boundary_next = halfedge(vh);
                        HalfedgeHandle boundary_prev = prev_halfedge(boundary_next);
                        next_cache.emplace_back(boundary_prev, outer_next);
                        next_cache.emplace_back(outer_prev,   boundary_next);
                    }
                    break;
                }
            }
            // Always link the two inner halfedges of the new face together.
            next_cache.emplace_back(inner_prev, inner_next);
        } else {
            // Both sides reused: the new face fills an existing corner. We may
            // need to re-pick the vertex's outgoing halfedge afterwards.
            needs_adjust[ii] = (halfedge(vh) == inner_next);
        }

        set_face(edge_he[i], fh);
    }

    // --- 6. apply all deferred next-pointer assignments ---------------------
    for (auto& kv : next_cache)
        set_next(kv.first, kv.second);

    // --- 7. restore the boundary-outgoing-halfedge invariant ----------------
    for (size_t i = 0; i < n; ++i)
        if (needs_adjust[i])
            adjust_outgoing_halfedge(vhs[i]);

    return fh;
}

FaceHandle Mesh::add_triangle(VertexHandle a, VertexHandle b, VertexHandle c) {
    return add_face({a, b, c});
}

void Mesh::triangulate() {
    // Snapshot the geometry + each face's vertex list, then rebuild from
    // scratch as triangles. Rebuilding sidesteps fragile in-place surgery and
    // guarantees a consistent result (same pattern as loop_subdivide).
    Mesh out;
    for (size_t i = 0; i < points_.size(); ++i)
        if (!v_deleted_[i]) out.add_vertex(points_[i]);
    // NOTE: assumes a compact mesh (call garbage_collection() first if needed),
    // so vertex indices carry over unchanged.
    for (size_t f = 0; f < faces_.size(); ++f) {
        if (f_deleted_[f]) continue;
        auto vs = face_vertices(FaceHandle((int)f));
        for (size_t k = 1; k + 1 < vs.size(); ++k)         // fan (0, k, k+1)
            out.add_face({ vs[0], vs[k], vs[k + 1] });
    }
    *this = std::move(out);
}

// ============================================================================
//  Circulators
// ============================================================================

std::vector<HalfedgeHandle> Mesh::vertex_outgoing_halfedges(VertexHandle v) const {
    std::vector<HalfedgeHandle> out;
    HalfedgeHandle start = halfedge(v);
    if (!start.is_valid()) return out;
    // Safety cap (see face_halfedges): bail rather than loop forever / OOM if the
    // one-ring fails to close on corrupt topology.
    const size_t cap = halfedges_.size() + 1;
    HalfedgeHandle h = start;
    do {
        out.push_back(h);
        h = next_halfedge(opposite_halfedge(h));
    } while (h != start && out.size() <= cap);
    return out;
}

std::vector<VertexHandle> Mesh::vertex_vertices(VertexHandle v) const {
    std::vector<VertexHandle> out;
    for (HalfedgeHandle h : vertex_outgoing_halfedges(v))
        out.push_back(to_vertex(h));
    return out;
}

std::vector<FaceHandle> Mesh::vertex_faces(VertexHandle v) const {
    std::vector<FaceHandle> out;
    for (HalfedgeHandle h : vertex_outgoing_halfedges(v))
        if (!is_boundary(h))
            out.push_back(face(h));
    return out;
}

size_t Mesh::valence(VertexHandle v) const {
    return vertex_outgoing_halfedges(v).size();
}

std::vector<HalfedgeHandle> Mesh::face_halfedges(FaceHandle f) const {
    std::vector<HalfedgeHandle> out;
    HalfedgeHandle start = halfedge(f);
    if (!start.is_valid()) return out;
    // Safety cap: a well-formed face has at most n_halfedges sides. If the
    // next-cycle fails to close (corrupt topology), bail out rather than loop
    // forever and exhaust memory. Costs one comparison per step on valid meshes.
    const size_t cap = halfedges_.size() + 1;
    HalfedgeHandle h = start;
    do {
        out.push_back(h);
        h = next_halfedge(h);
    } while (h != start && out.size() <= cap);
    return out;
}

std::vector<VertexHandle> Mesh::face_vertices(FaceHandle f) const {
    std::vector<VertexHandle> out;
    for (HalfedgeHandle h : face_halfedges(f))
        out.push_back(to_vertex(h));
    return out;
}

std::vector<EdgeHandle> Mesh::vertex_edges(VertexHandle v) const {
    std::vector<EdgeHandle> out;
    for (HalfedgeHandle h : vertex_outgoing_halfedges(v))
        out.push_back(edge(h));
    return out;
}

std::vector<HalfedgeHandle> Mesh::vertex_incoming_halfedges(VertexHandle v) const {
    std::vector<HalfedgeHandle> out;
    for (HalfedgeHandle h : vertex_outgoing_halfedges(v))
        out.push_back(opposite_halfedge(h));
    return out;
}

std::vector<EdgeHandle> Mesh::face_edges(FaceHandle f) const {
    std::vector<EdgeHandle> out;
    for (HalfedgeHandle h : face_halfedges(f))
        out.push_back(edge(h));
    return out;
}

std::vector<FaceHandle> Mesh::face_faces(FaceHandle f) const {
    std::vector<FaceHandle> out;
    for (HalfedgeHandle h : face_halfedges(f)) {
        HalfedgeHandle o = opposite_halfedge(h);
        if (!is_boundary(o)) out.push_back(face(o));
    }
    return out;
}

// ============================================================================
//  Geometry
// ============================================================================

Vec3 Mesh::calc_face_normal(FaceHandle f) const {
    // Newell's method: robust for arbitrary (even non-planar) polygons and
    // degenerates gracefully to the triangle cross product.
    Vec3 n{0, 0, 0};
    auto verts = face_vertices(f);
    const size_t m = verts.size();
    for (size_t i = 0; i < m; ++i) {
        const Vec3& a = point(verts[i]);
        const Vec3& b = point(verts[(i + 1) % m]);
        n.x += (a.y - b.y) * (a.z + b.z);
        n.y += (a.z - b.z) * (a.x + b.x);
        n.z += (a.x - b.x) * (a.y + b.y);
    }
    return n.normalized();
}

Vec3 Mesh::calc_face_centroid(FaceHandle f) const {
    Vec3 c{0, 0, 0};
    auto verts = face_vertices(f);
    for (VertexHandle v : verts) c += point(v);
    return verts.empty() ? c : c / static_cast<double>(verts.size());
}

double Mesh::calc_edge_length(EdgeHandle e) const {
    HalfedgeHandle h = halfedge(e, 0);
    return (point(to_vertex(h)) - point(from_vertex(h))).norm();
}

double Mesh::calc_face_area(FaceHandle f) const {
    // The Newell vector's length is twice the (signed) polygon area.
    Vec3 n{0, 0, 0};
    auto verts = face_vertices(f);
    const size_t m = verts.size();
    for (size_t i = 0; i < m; ++i) {
        const Vec3& a = point(verts[i]);
        const Vec3& b = point(verts[(i + 1) % m]);
        n.x += (a.y - b.y) * (a.z + b.z);
        n.y += (a.z - b.z) * (a.x + b.x);
        n.z += (a.x - b.x) * (a.y + b.y);
    }
    return 0.5 * n.norm();
}

Vec3 Mesh::calc_edge_vector(HalfedgeHandle h) const {
    return point(to_vertex(h)) - point(from_vertex(h));
}

Vec3 Mesh::calc_edge_midpoint(EdgeHandle e) const {
    HalfedgeHandle h = halfedge(e, 0);
    return (point(to_vertex(h)) + point(from_vertex(h))) * 0.5;
}

Vec3 Mesh::calc_vertex_normal(VertexHandle v) const {
    Vec3 n{0, 0, 0};
    for (FaceHandle f : vertex_faces(v)) n += calc_face_normal(f);
    return n.normalized();
}

double Mesh::calc_dihedral_angle(EdgeHandle e) const {
    if (is_boundary(e)) return 0.0;
    HalfedgeHandle h = halfedge(e, 0);
    Vec3 n0 = calc_face_normal(face(h));
    Vec3 n1 = calc_face_normal(face(opposite_halfedge(h)));
    double d = std::max(-1.0, std::min(1.0, n0.dot(n1)));
    double angle = std::acos(d);
    // sign: convex vs concave, using the edge direction
    Vec3 ev = calc_edge_vector(h).normalized();
    return (n0.cross(n1).dot(ev) < 0.0) ? -angle : angle;
}

double Mesh::calc_sector_angle(HalfedgeHandle h) const {
    // Angle at the vertex `to_vertex(prev)` == from_vertex(h)'s corner inside
    // the face: between the incoming edge (-prev) and the outgoing edge (h).
    Vec3 v1 = calc_edge_vector(h).normalized();
    Vec3 v0 = (calc_edge_vector(prev_halfedge(h)) * -1.0).normalized();
    double d = std::max(-1.0, std::min(1.0, v0.dot(v1)));
    return std::acos(d);
}

double Mesh::surface_area() const {
    double total = 0.0;
    for (size_t f = 0; f < faces_.size(); ++f)
        if (!f_deleted_[f]) total += calc_face_area(FaceHandle((int)f));
    return total;
}

Vec3 Mesh::center_of_mass() const {
    Vec3 c{0, 0, 0};
    size_t n = 0;
    for (size_t i = 0; i < points_.size(); ++i)
        if (!v_deleted_[i]) { c += points_[i]; ++n; }
    return n ? c / (double)n : c;
}

std::pair<Vec3, Vec3> Mesh::bounding_box() const {
    const double inf = std::numeric_limits<double>::infinity();
    Vec3 lo{inf, inf, inf}, hi{-inf, -inf, -inf};
    for (size_t i = 0; i < points_.size(); ++i) {
        if (v_deleted_[i]) continue;
        const Vec3& p = points_[i];
        lo.x = std::min(lo.x, p.x); lo.y = std::min(lo.y, p.y); lo.z = std::min(lo.z, p.z);
        hi.x = std::max(hi.x, p.x); hi.y = std::max(hi.y, p.y); hi.z = std::max(hi.z, p.z);
    }
    return {lo, hi};
}

// ============================================================================
//  Deletion (lazy) - elements are flagged, not removed. See garbage_collection.
// ============================================================================

void Mesh::delete_face(FaceHandle fh, bool delete_isolated) {
    assert(fh.is_valid() && !is_deleted(fh));
    f_deleted_[fh.idx()] = true;

    // Collect this face's boundary edges (those whose opposite is already a
    // boundary halfedge) - deleting our face turns them fully unused.
    std::vector<EdgeHandle>   deleted_edges;
    std::vector<VertexHandle> vhandles;

    for (HalfedgeHandle hh : face_halfedges(fh)) {
        set_boundary(hh);  // detach: this halfedge now borders a hole
        if (is_boundary(opposite_halfedge(hh)))
            deleted_edges.push_back(edge(hh));
        vhandles.push_back(to_vertex(hh));
    }

    // Remove the now-isolated boundary edges and patch the boundary loop.
    for (EdgeHandle e : deleted_edges) {
        HalfedgeHandle h0 = halfedge(e, 0);
        VertexHandle   v0 = to_vertex(h0);
        HalfedgeHandle next0 = next_halfedge(h0), prev0 = prev_halfedge(h0);

        HalfedgeHandle h1 = halfedge(e, 1);
        VertexHandle   v1 = to_vertex(h1);
        HalfedgeHandle next1 = next_halfedge(h1), prev1 = prev_halfedge(h1);

        set_next(prev0, next1);
        set_next(prev1, next0);
        e_deleted_[e.idx()] = true;

        if (halfedge(v0) == h1) {                 // v0 pointed through this edge
            if (next0 == h1) {                    // v0 became isolated
                if (delete_isolated) v_deleted_[v0.idx()] = true;
                set_isolated(v0);
            } else set_halfedge(v0, next0);
        }
        if (halfedge(v1) == h0) {
            if (next1 == h0) {
                if (delete_isolated) v_deleted_[v1.idx()] = true;
                set_isolated(v1);
            } else set_halfedge(v1, next1);
        }
    }

    for (VertexHandle v : vhandles)
        if (!is_isolated(v)) adjust_outgoing_halfedge(v);
}

void Mesh::delete_vertex(VertexHandle vh, bool delete_isolated) {
    std::vector<FaceHandle> faces = vertex_faces(vh);  // snapshot before mutation
    for (FaceHandle f : faces)
        if (!is_deleted(f)) delete_face(f, delete_isolated);
    v_deleted_[vh.idx()] = true;
}

void Mesh::delete_edge(EdgeHandle e, bool delete_isolated) {
    FaceHandle f0 = face(halfedge(e, 0));
    FaceHandle f1 = face(halfedge(e, 1));
    if (f0.is_valid() && !is_deleted(f0)) delete_face(f0, delete_isolated);
    if (f1.is_valid() && !is_deleted(f1)) delete_face(f1, delete_isolated);
    e_deleted_[e.idx()] = true;
}

// ============================================================================
//  Garbage collection - compact the arrays, dropping deleted elements.
//  We rebuild fresh arrays and remap every internal handle. Edges are kept in
//  pairs so the opposite() == (idx XOR 1) invariant survives.
//  NOTE: every handle the caller holds is invalidated afterwards.
// ============================================================================

void Mesh::garbage_collection() {
    const int nV = static_cast<int>(vertices_.size());
    const int nE = static_cast<int>(n_edges());
    const int nF = static_cast<int>(faces_.size());

    // old index -> new index (or -1 if dropped)
    std::vector<int> vmap(nV, -1), emap(nE, -1), fmap(nF, -1);
    int newV = 0, newE = 0, newF = 0;
    for (int i = 0; i < nV; ++i) if (!v_deleted_[i]) vmap[i] = newV++;
    for (int i = 0; i < nE; ++i) if (!e_deleted_[i]) emap[i] = newE++;
    for (int i = 0; i < nF; ++i) if (!f_deleted_[i]) fmap[i] = newF++;

    auto map_v = [&](VertexHandle v) {
        return v.is_valid() ? VertexHandle(vmap[v.idx()]) : v; };
    auto map_f = [&](FaceHandle f) {
        return f.is_valid() ? FaceHandle(fmap[f.idx()]) : f; };
    auto map_h = [&](HalfedgeHandle h) {
        if (!h.is_valid()) return HalfedgeHandle();
        int ne = emap[h.idx() >> 1];
        return ne < 0 ? HalfedgeHandle() : HalfedgeHandle(ne * 2 + (h.idx() & 1));
    };

    std::vector<Vertex> nv(newV);
    std::vector<Vec3>   np(newV);
    for (int i = 0; i < nV; ++i) if (vmap[i] >= 0) {
        nv[vmap[i]].outgoing_halfedge = map_h(vertices_[i].outgoing_halfedge);
        np[vmap[i]] = points_[i];
    }

    std::vector<Halfedge> nh(newE * 2);
    for (int e = 0; e < nE; ++e) if (emap[e] >= 0) {
        for (int s = 0; s < 2; ++s) {
            const Halfedge& H = halfedges_[2 * e + s];
            Halfedge& D = nh[2 * emap[e] + s];
            D.to_vertex = map_v(H.to_vertex);
            D.face      = map_f(H.face);
            D.next      = map_h(H.next);
            D.prev      = map_h(H.prev);
        }
    }

    std::vector<Face> nf(newF);
    for (int i = 0; i < nF; ++i) if (fmap[i] >= 0)
        nf[fmap[i]].halfedge = map_h(faces_[i].halfedge);

    vertices_ = std::move(nv);
    points_   = std::move(np);
    halfedges_ = std::move(nh);
    faces_    = std::move(nf);
    v_deleted_.assign(newV, false);
    e_deleted_.assign(newE, false);
    f_deleted_.assign(newF, false);

    // Compact custom properties with the same maps. Halfedges follow their
    // edge: old halfedge 2e+s -> 2*emap[e]+s (or dropped if the edge is gone).
    std::vector<int> hmap_old_to_new(nE * 2, -1);
    for (int e = 0; e < nE; ++e) if (emap[e] >= 0) {
        hmap_old_to_new[2 * e]     = 2 * emap[e];
        hmap_old_to_new[2 * e + 1] = 2 * emap[e] + 1;
    }
    vprops_.compact(vmap, newV);
    eprops_.compact(emap, newE);
    fprops_.compact(fmap, newF);
    hprops_.compact(hmap_old_to_new, newE * 2);
}

// ============================================================================
//  Edge flip (triangle meshes)
// ============================================================================

bool Mesh::is_flip_ok(EdgeHandle e) const {
    if (is_boundary(e)) return false;             // boundary edges can't flip

    HalfedgeHandle h = halfedge(e, 0), o = halfedge(e, 1);
    // The two "tip" vertices that would become the new edge's endpoints.
    VertexHandle a = to_vertex(next_halfedge(h));
    VertexHandle b = to_vertex(next_halfedge(o));
    if (a == b) return false;                     // degenerate

    // If a-b already exists, flipping would duplicate it -> non-manifold.
    for (VertexHandle nb : vertex_vertices(a))
        if (nb == b) return false;
    return true;
}

void Mesh::flip(EdgeHandle e) {
    assert(is_flip_ok(e));
    HalfedgeHandle a0 = halfedge(e, 0), b0 = halfedge(e, 1);
    HalfedgeHandle a1 = next_halfedge(a0), a2 = next_halfedge(a1);
    HalfedgeHandle b1 = next_halfedge(b0), b2 = next_halfedge(b1);

    VertexHandle va0 = to_vertex(a0), va1 = to_vertex(a1);
    VertexHandle vb0 = to_vertex(b0), vb1 = to_vertex(b1);
    FaceHandle   fa = face(a0), fb = face(b0);

    set_to_vertex(a0, va1);
    set_to_vertex(b0, vb1);

    set_next(a0, a2); set_next(a2, b1); set_next(b1, a0);
    set_next(b0, b2); set_next(b2, a1); set_next(a1, b0);

    set_face(a1, fb); set_face(b1, fa);
    set_halfedge(fa, a0); set_halfedge(fb, b0);

    if (halfedge(va0) == b0) set_halfedge(va0, a1);
    if (halfedge(vb0) == a0) set_halfedge(vb0, b1);
}

// ============================================================================
//  Edge split (triangle meshes): insert v on edge e, subdividing neighbours.
// ============================================================================

void Mesh::split(EdgeHandle e, VertexHandle vh) {
    HalfedgeHandle h0 = halfedge(e, 0), o0 = halfedge(e, 1);
    VertexHandle   v2 = to_vertex(o0);

    HalfedgeHandle e1 = new_edge(vh, v2);
    HalfedgeHandle t1 = opposite_halfedge(e1);

    FaceHandle f0 = face(h0), f3 = face(o0);
    set_halfedge(vh, h0);
    set_to_vertex(o0, vh);

    if (!is_boundary(h0)) {
        HalfedgeHandle h1 = next_halfedge(h0), h2 = next_halfedge(h1);
        VertexHandle   v1 = to_vertex(h1);
        HalfedgeHandle e0 = new_edge(vh, v1);
        HalfedgeHandle t0 = opposite_halfedge(e0);
        FaceHandle     f1 = new_face();
        set_halfedge(f0, h0); set_halfedge(f1, h2);
        set_face(h1, f0); set_face(t0, f0); set_face(h0, f0);
        set_face(h2, f1); set_face(t1, f1); set_face(e0, f1);
        set_next(h0, h1); set_next(h1, t0); set_next(t0, h0);
        set_next(e0, h2); set_next(h2, t1); set_next(t1, e0);
    } else {
        set_next(prev_halfedge(h0), t1);
        set_next(t1, h0);
    }

    if (!is_boundary(o0)) {
        HalfedgeHandle o1 = next_halfedge(o0), o2 = next_halfedge(o1);
        VertexHandle   v3 = to_vertex(o1);
        HalfedgeHandle e2 = new_edge(vh, v3);
        HalfedgeHandle t2 = opposite_halfedge(e2);
        FaceHandle     f2 = new_face();
        set_halfedge(f2, o1); set_halfedge(f3, o0);
        set_face(o1, f2); set_face(t2, f2); set_face(e1, f2);
        set_face(o2, f3); set_face(o0, f3); set_face(e2, f3);
        set_next(e1, o1); set_next(o1, t2); set_next(t2, e1);
        set_next(o0, e2); set_next(e2, o2); set_next(o2, o0);
    } else {
        set_next(e1, next_halfedge(o0));
        set_next(o0, e1);
        set_halfedge(vh, e1);
    }

    if (halfedge(v2) == h0) set_halfedge(v2, t1);
}

// ============================================================================
//  Face split: fan the polygon `f` around an interior vertex `v`.
// ============================================================================

void Mesh::split(FaceHandle fh, VertexHandle vh) {
    HalfedgeHandle hend = halfedge(fh);
    HalfedgeHandle hh   = next_halfedge(hend);

    HalfedgeHandle hold = new_edge(to_vertex(hend), vh);
    set_next(hend, hold);
    set_face(hold, fh);
    hold = opposite_halfedge(hold);

    while (hh != hend) {
        HalfedgeHandle hnext = next_halfedge(hh);
        FaceHandle fnew = new_face();
        set_halfedge(fnew, hh);
        HalfedgeHandle hnew = new_edge(to_vertex(hh), vh);
        set_next(hnew, hold);
        set_next(hold, hh);
        set_next(hh, hnew);
        set_face(hnew, fnew);
        set_face(hold, fnew);
        set_face(hh, fnew);
        hold = opposite_halfedge(hnew);
        hh = hnext;
    }

    set_next(hold, hend);
    set_next(next_halfedge(hend), hold);
    set_face(hold, fh);
    set_halfedge(vh, hold);
}

// ============================================================================
//  Edge collapse: merge from_vertex(h) into to_vertex(h).
// ============================================================================

bool Mesh::is_collapse_ok(HalfedgeHandle v0v1) {
    if (is_deleted(edge(v0v1))) return false;

    HalfedgeHandle v1v0 = opposite_halfedge(v0v1);
    VertexHandle   v0 = to_vertex(v1v0);   // the vertex being removed
    VertexHandle   v1 = to_vertex(v0v1);   // the vertex kept
    if (is_deleted(v0) || is_deleted(v1)) return false;

    bool v0v1_tri = !is_boundary(v0v1) && face_halfedges(face(v0v1)).size() == 3;
    bool v1v0_tri = !is_boundary(v1v0) && face_halfedges(face(v1v0)).size() == 3;

    VertexHandle v_01_n = to_vertex(next_halfedge(v0v1));
    VertexHandle v_10_n = to_vertex(next_halfedge(v1v0));

    // The two side edges of each collapsing triangle must not both be boundary.
    VertexHandle vl, vr;
    if (v0v1_tri) {
        HalfedgeHandle h1 = next_halfedge(v0v1), h2 = next_halfedge(h1);
        vl = to_vertex(h1);
        if (is_boundary(opposite_halfedge(h1)) && is_boundary(opposite_halfedge(h2)))
            return false;
    }
    if (v1v0_tri) {
        HalfedgeHandle h1 = next_halfedge(v1v0), h2 = next_halfedge(h1);
        vr = to_vertex(h1);
        if (is_boundary(opposite_halfedge(h1)) && is_boundary(opposite_halfedge(h2)))
            return false;
    }
    if (vl.is_valid() && vl == vr) return false;

    // An interior edge between two boundary vertices can't be collapsed.
    if (is_boundary(v0) && is_boundary(v1) && !is_boundary(v0v1) && !is_boundary(v1v0))
        return false;

    // The 1-rings of v0 and v1 may only share the two tip vertices, otherwise
    // the collapse would fold the surface onto itself.
    std::vector<VertexHandle> ring1 = vertex_vertices(v1);
    for (VertexHandle vv : vertex_vertices(v0)) {
        bool shared = false;
        for (VertexHandle w : ring1) if (w == vv) { shared = true; break; }
        if (shared &&
            !(vv == v_01_n && v0v1_tri) &&
            !(vv == v_10_n && v1v0_tri))
            return false;
    }
    return true;
}

void Mesh::collapse(HalfedgeHandle hh) {
    HalfedgeHandle h0 = hh;
    HalfedgeHandle h1 = next_halfedge(h0);
    HalfedgeHandle o0 = opposite_halfedge(h0);
    HalfedgeHandle o1 = next_halfedge(o0);

    collapse_edge(h0);

    // If either incident face was a triangle it has degenerated into a 2-edge
    // "loop" that must be collapsed away too.
    if (next_halfedge(next_halfedge(h1)) == h1) collapse_loop(next_halfedge(h1));
    if (next_halfedge(next_halfedge(o1)) == o1) collapse_loop(o1);
}

void Mesh::collapse_edge(HalfedgeHandle h) {
    HalfedgeHandle hn = next_halfedge(h), hp = prev_halfedge(h);
    HalfedgeHandle o  = opposite_halfedge(h);
    HalfedgeHandle on = next_halfedge(o), op = prev_halfedge(o);
    FaceHandle     fh = face(h), fo = face(o);
    VertexHandle   vh = to_vertex(h);   // kept
    VertexHandle   vo = to_vertex(o);   // removed

    // Re-target every halfedge that pointed at vo so it points at vh.
    for (HalfedgeHandle outg : vertex_outgoing_halfedges(vo))
        set_to_vertex(opposite_halfedge(outg), vh);

    set_next(hp, hn);
    set_next(op, on);

    if (fh.is_valid()) set_halfedge(fh, hn);
    if (fo.is_valid()) set_halfedge(fo, on);

    if (halfedge(vh) == o) set_halfedge(vh, hn);
    adjust_outgoing_halfedge(vh);
    set_isolated(vo);

    e_deleted_[edge(h).idx()] = true;
    v_deleted_[vo.idx()] = true;
}

void Mesh::collapse_loop(HalfedgeHandle h0) {
    HalfedgeHandle h1 = next_halfedge(h0);
    HalfedgeHandle o0 = opposite_halfedge(h0);
    HalfedgeHandle o1 = opposite_halfedge(h1);
    VertexHandle   v0 = to_vertex(h0), v1 = to_vertex(h1);
    FaceHandle     fh = face(h0), fo = face(o0);

    set_next(h1, next_halfedge(o0));
    set_next(prev_halfedge(o0), h1);

    set_face(h1, fo);

    set_halfedge(v0, h1); adjust_outgoing_halfedge(v0);
    set_halfedge(v1, o1); adjust_outgoing_halfedge(v1);

    if (fo.is_valid() && halfedge(fo) == o0) set_halfedge(fo, h1);

    if (fh.is_valid()) {
        set_halfedge(fh, HalfedgeHandle());
        f_deleted_[fh.idx()] = true;
    }
    e_deleted_[edge(h0).idx()] = true;
}

// ============================================================================
//  Vertex split: the inverse of collapse (triangle meshes).
// ----------------------------------------------------------------------------
//  Built from two
//  primitives that are themselves the inverses of the collapse helpers:
//
//    * insert_loop(h) - splits the open region left of `h` off into its own
//      triangle by stitching in a parallel halfedge (inverse of collapse_loop).
//    * insert_edge(v0, h0, h1) - re-introduces the edge v0->v1 between the two
//      wings, re-homing v1's halfedges that now belong to v0 (inverse of
//      collapse_edge).
// ============================================================================

HalfedgeHandle Mesh::insert_loop(HalfedgeHandle h0) {
    HalfedgeHandle o0 = opposite_halfedge(h0);

    VertexHandle v0 = to_vertex(o0);   // from-vertex of h0
    VertexHandle v1 = to_vertex(h0);   // to-vertex   of h0

    HalfedgeHandle h1 = new_edge(v1, v0);
    HalfedgeHandle o1 = opposite_halfedge(h1);

    FaceHandle f0 = face(h0);
    FaceHandle f1 = new_face();

    // halfedge -> halfedge: splice o1 into h0's old loop, and close h0/h1 into
    // a new 2-gon that new_face() will own.
    set_next(prev_halfedge(h0), o1);
    set_next(o1, next_halfedge(h0));
    set_next(h1, h0);
    set_next(h0, h1);

    // halfedge -> face
    set_face(o1, f0);
    set_face(h0, f1);
    set_face(h1, f1);

    // face -> halfedge
    set_halfedge(f1, h0);
    if (f0.is_valid()) set_halfedge(f0, o1);

    adjust_outgoing_halfedge(v0);
    adjust_outgoing_halfedge(v1);
    return h1;
}

HalfedgeHandle Mesh::insert_edge(VertexHandle v0, HalfedgeHandle h0, HalfedgeHandle h1) {
    assert(h0.is_valid() && h1.is_valid());
    VertexHandle v1 = to_vertex(h0);
    assert(v1 == to_vertex(h1));

    HalfedgeHandle v0v1 = new_edge(v0, v1);
    HalfedgeHandle v1v0 = opposite_halfedge(v0v1);

    // vertex -> halfedge
    set_halfedge(v0, v0v1);
    set_halfedge(v1, v1v0);

    // halfedge -> halfedge: weave the new edge between the two wing halfedges.
    set_next(v0v1, next_halfedge(h0));
    set_next(h0, v0v1);
    set_next(v1v0, next_halfedge(h1));
    set_next(h1, v1v0);

    // every halfedge now reachable around v0 must point at v0 (was v1 before).
    for (HalfedgeHandle oh : vertex_outgoing_halfedges(v0))
        set_to_vertex(opposite_halfedge(oh), v0);

    // halfedge -> face
    set_face(v0v1, face(h0));
    set_face(v1v0, face(h1));
    if (face(v0v1).is_valid()) set_halfedge(face(v0v1), v0v1);
    if (face(v1v0).is_valid()) set_halfedge(face(v1v0), v1v0);

    adjust_outgoing_halfedge(v0);
    adjust_outgoing_halfedge(v1);
    return v0v1;
}

HalfedgeHandle Mesh::vertex_split(VertexHandle v0, VertexHandle v1,
                                  VertexHandle vl, VertexHandle vr) {
    HalfedgeHandle v1vl, vlv1, vrv1;

    // build loop from halfedge v1->vl
    if (vl.is_valid()) {
        v1vl = find_halfedge(v1, vl);
        assert(v1vl.is_valid());
        vlv1 = insert_loop(v1vl);
    }
    // build loop from halfedge vr->v1
    if (vr.is_valid()) {
        vrv1 = find_halfedge(vr, v1);
        assert(vrv1.is_valid());
        insert_loop(vrv1);
    }
    // boundary cases: the missing wing is the boundary halfedge into v1
    if (!vl.is_valid()) vlv1 = prev_halfedge(halfedge(v1));
    if (!vr.is_valid()) vrv1 = prev_halfedge(halfedge(v1));

    // split vertex v1 into the edge v0-v1
    return insert_edge(v0, vlv1, vrv1);
}

// ============================================================================
//  Minimal OBJ I/O
// ============================================================================

bool Mesh::read_obj(const std::string& path) {
    std::ifstream in(path);
    if (!in) { std::cerr << "read_obj: cannot open " << path << "\n"; return false; }

    std::vector<VertexHandle> vmap;  // OBJ index (1-based) -> handle
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream ls(line);
        std::string tag;
        ls >> tag;
        if (tag == "v") {
            double x, y, z; ls >> x >> y >> z;
            vmap.push_back(add_vertex({x, y, z}));
        } else if (tag == "f") {
            std::vector<VertexHandle> fv;
            std::string tok;
            while (ls >> tok) {
                // handle "v", "v/vt", "v/vt/vn" - we only need the vertex index.
                int vi = std::stoi(tok.substr(0, tok.find('/')));
                fv.push_back(vmap[vi - 1]);
            }
            if (fv.size() >= 3) add_face(fv);
        }
    }
    return true;
}

bool Mesh::write_obj(const std::string& path) const {
    std::ofstream out(path);
    if (!out) { std::cerr << "write_obj: cannot open " << path << "\n"; return false; }

    out << "# written by SimpleMesh\n";
    // Skip deleted vertices and build old-index -> 1-based-OBJ-index map.
    std::vector<int> obj_index(points_.size(), 0);
    int next = 1;
    for (size_t i = 0; i < points_.size(); ++i) {
        if (v_deleted_[i]) continue;
        out << "v " << points_[i].x << ' ' << points_[i].y << ' ' << points_[i].z << '\n';
        obj_index[i] = next++;
    }
    for (size_t fi = 0; fi < faces_.size(); ++fi) {
        if (f_deleted_[fi]) continue;
        out << 'f';
        for (VertexHandle v : face_vertices(FaceHandle(static_cast<int>(fi))))
            out << ' ' << obj_index[v.idx()];  // remapped, OBJ is 1-based
        out << '\n';
    }
    return true;
}

} // namespace sm
