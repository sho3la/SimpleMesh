// ============================================================================
//  SimpleMesh - Algorithms.cpp : Loop subdivision + QEM decimation
// ============================================================================
#include "simplemesh/Algorithms.h"

#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <vector>

namespace sm {

// ============================================================================
//  Loop subdivision
// ----------------------------------------------------------------------------
//  We REBUILD a fresh mesh rather than surgically splitting in place - far
//  clearer and impossible to leave in an inconsistent state. Each old triangle
//  becomes four:
//
//            c                     c
//           / \                   / \
//          /   \      ===>     mca---mbc
//         /     \             / \   / \
//        a-------b           a---mab---b
//
//  New positions use the standard Loop weights:
//    * old vertex (interior): (1 - n*beta)*v + beta*sum(neighbours)
//    * old vertex (boundary): 3/4*v + 1/8*(two boundary neighbours)
//    * edge vertex (interior): 3/8*(a+b) + 1/8*(c+d)   (c,d = opposite apexes)
//    * edge vertex (boundary): 1/2*(a+b)
// ============================================================================

static Mesh loop_subdivide_once(const Mesh& in) {
    const double PI = 3.14159265358979323846;
    Mesh out;

    const size_t nV = in.n_vertices();
    const size_t nE = in.n_edges();

    // --- 1. repositioned old vertices (indices 0..nV-1 carry over) ---------
    for (size_t vi = 0; vi < nV; ++vi) {
        VertexHandle v(static_cast<int>(vi));

        if (in.is_boundary(v)) {
            // boundary neighbours = endpoints reached along boundary edges
            std::vector<VertexHandle> bn;
            for (auto h : in.voh_range(v))
                if (in.is_boundary(in.edge(h))) bn.push_back(in.to_vertex(h));

            if (bn.size() == 2)
                out.add_vertex(in.point(v) * 0.75 +
                               (in.point(bn[0]) + in.point(bn[1])) * 0.125);
            else
                out.add_vertex(in.point(v));   // corner / non-manifold: keep
        } else {
            auto nb = in.vertex_vertices(v);
            const int n = static_cast<int>(nb.size());
            const double t = 3.0 / 8.0 + 0.25 * std::cos(2.0 * PI / n);
            const double beta = (1.0 / n) * (5.0 / 8.0 - t * t);
            Vec3 sum{0, 0, 0};
            for (auto w : nb) sum += in.point(w);
            out.add_vertex(in.point(v) * (1.0 - n * beta) + sum * beta);
        }
    }

    // --- 2. one new vertex per edge ----------------------------------------
    std::vector<VertexHandle> emid(nE);
    for (size_t ei = 0; ei < nE; ++ei) {
        EdgeHandle e(static_cast<int>(ei));
        HalfedgeHandle h = in.halfedge(e, 0);
        VertexHandle a = in.from_vertex(h), b = in.to_vertex(h);

        if (in.is_boundary(e)) {
            emid[ei] = out.add_vertex((in.point(a) + in.point(b)) * 0.5);
        } else {
            VertexHandle c = in.to_vertex(in.next_halfedge(h));
            VertexHandle d = in.to_vertex(in.next_halfedge(in.opposite_halfedge(h)));
            emid[ei] = out.add_vertex((in.point(a) + in.point(b)) * (3.0 / 8.0) +
                                      (in.point(c) + in.point(d)) * (1.0 / 8.0));
        }
    }

    // --- 3. four faces per old triangle ------------------------------------
    for (size_t fi = 0; fi < in.n_faces(); ++fi) {
        FaceHandle f(static_cast<int>(fi));
        auto he = in.face_halfedges(f);
        if (he.size() != 3) continue;          // Loop is defined for triangles

        VertexHandle a(in.from_vertex(he[0]).idx());   // he[0]: a->b
        VertexHandle b(in.from_vertex(he[1]).idx());   // he[1]: b->c
        VertexHandle c(in.from_vertex(he[2]).idx());   // he[2]: c->a
        VertexHandle mab = emid[in.edge(he[0]).idx()];
        VertexHandle mbc = emid[in.edge(he[1]).idx()];
        VertexHandle mca = emid[in.edge(he[2]).idx()];

        out.add_face({a,   mab, mca});
        out.add_face({mab, b,   mbc});
        out.add_face({mca, mbc, c});
        out.add_face({mab, mbc, mca});
    }

    return out;
}

Mesh loop_subdivide(const Mesh& in, int iterations) {
    Mesh m = loop_subdivide_once(in);
    for (int k = 1; k < iterations; ++k)
        m = loop_subdivide_once(m);
    return m;
}

// ============================================================================
//  Catmull-Clark subdivision (primal, quad-producing)
// ----------------------------------------------------------------------------
//  Like Loop, we REBUILD rather than edit in place. Three kinds of new point:
//
//    * face point   F = centroid of the face's vertices.
//    * edge point   E = average of the edge's two endpoints and the two
//                       adjacent face points (interior); midpoint (boundary).
//    * vertex point V = (F_avg + 2*R_avg + (n-3)*P) / n   for an interior
//                       vertex of valence n, where F_avg / R_avg are the
//                       averages of the incident face points / edge midpoints;
//                       boundary vertices follow the cubic B-spline rule
//                       3/4 P + 1/8 (nbr1 + nbr2).
//
//  Then each old face of k sides becomes k quads, one per corner:
//        [ V(corner), E(outgoing edge), F, E(incoming edge) ].
// ============================================================================

static Mesh catmull_clark_once(const Mesh& in) {
    Mesh out;
    const size_t nV = in.n_vertices();
    const size_t nF = in.n_faces();
    const size_t nE = in.n_edges();

    // --- face points -------------------------------------------------------
    std::vector<Vec3> fpoint(nF);
    for (size_t fi = 0; fi < nF; ++fi)
        fpoint[fi] = in.calc_face_centroid(FaceHandle(static_cast<int>(fi)));

    // --- 1. updated old vertices occupy indices 0..nV-1 --------------------
    for (size_t vi = 0; vi < nV; ++vi) {
        VertexHandle v(static_cast<int>(vi));
        const Vec3 P = in.point(v);

        if (in.is_boundary(v)) {
            std::vector<VertexHandle> bn;     // boundary neighbours
            for (auto h : in.voh_range(v))
                if (in.is_boundary(in.edge(h))) bn.push_back(in.to_vertex(h));
            if (bn.size() == 2)
                out.add_vertex(P * 0.75 + (in.point(bn[0]) + in.point(bn[1])) * 0.125);
            else
                out.add_vertex(P);            // corner / non-manifold: keep
        } else {
            Vec3 Favg{0, 0, 0}, Ravg{0, 0, 0};
            int nf = 0, ne = 0;
            for (auto f : in.vertex_faces(v))  { Favg += fpoint[f.idx()]; ++nf; }
            for (auto e : in.vertex_edges(v))  { Ravg += in.calc_edge_midpoint(e); ++ne; }
            const double n = static_cast<double>(ne);   // valence
            Favg = Favg / (double)nf;
            Ravg = Ravg / (double)ne;
            out.add_vertex((Favg + Ravg * 2.0 + P * (n - 3.0)) / n);
        }
    }

    // --- 2. face points appended (index fp_idx[f]) -------------------------
    std::vector<VertexHandle> fp_idx(nF);
    for (size_t fi = 0; fi < nF; ++fi)
        fp_idx[fi] = out.add_vertex(fpoint[fi]);

    // --- 3. edge points appended (index ep_idx[e]) -------------------------
    std::vector<VertexHandle> ep_idx(nE);
    for (size_t ei = 0; ei < nE; ++ei) {
        EdgeHandle e(static_cast<int>(ei));
        HalfedgeHandle h = in.halfedge(e, 0);
        Vec3 a = in.point(in.from_vertex(h)), b = in.point(in.to_vertex(h));
        if (in.is_boundary(e)) {
            ep_idx[ei] = out.add_vertex((a + b) * 0.5);
        } else {
            Vec3 f0 = fpoint[in.face(h).idx()];
            Vec3 f1 = fpoint[in.face(in.opposite_halfedge(h)).idx()];
            ep_idx[ei] = out.add_vertex((a + b + f0 + f1) * 0.25);
        }
    }

    // --- 4. one quad per corner of every old face --------------------------
    for (size_t fi = 0; fi < nF; ++fi) {
        FaceHandle f(static_cast<int>(fi));
        auto he = in.face_halfedges(f);            // h_i : corner_i -> corner_{i+1}
        const size_t k = he.size();
        for (size_t i = 0; i < k; ++i) {
            HalfedgeHandle h    = he[i];
            HalfedgeHandle hprev = he[(i + k - 1) % k];   // incoming edge at corner
            VertexHandle corner(in.from_vertex(h).idx());
            out.add_face({ corner,
                           ep_idx[in.edge(h).idx()],
                           fp_idx[fi],
                           ep_idx[in.edge(hprev).idx()] });
        }
    }
    return out;
}

Mesh catmull_clark(const Mesh& in, int iterations) {
    Mesh m = catmull_clark_once(in);
    for (int k = 1; k < iterations; ++k)
        m = catmull_clark_once(m);
    return m;
}

// ============================================================================
//  Sqrt(3) subdivision (Kobbelt 2000)
// ----------------------------------------------------------------------------
//  Triangle scheme with three moves per step:
//    1. insert a centroid vertex in every face and split it 1 -> 3,
//    2. flip every ORIGINAL interior edge,
//    3. relax old interior vertices: P' = (1 - a_n) P + a_n * mean(neighbours),
//       a_n = (4 - 2 cos(2*pi/n)) / 9.
//
//  We build the 1 -> 3 split as a fresh mesh and then reuse the tested edge
//  flip() for step 2 - flipping each original edge a-b turns it into the edge
//  between the two new centroid vertices. Boundary vertices are kept fixed and
//  boundary edges are not flipped (the proper boundary scheme alternates rules
//  every other step; we keep it simple for teaching).
// ============================================================================

static Mesh sqrt3_subdivide_once(const Mesh& in) {
    const double PI = 3.14159265358979323846;
    Mesh out;
    const size_t nV = in.n_vertices();
    const size_t nF = in.n_faces();

    // --- 1. relaxed old vertices (indices 0..nV-1 carry over) --------------
    for (size_t vi = 0; vi < nV; ++vi) {
        VertexHandle v(static_cast<int>(vi));
        if (in.is_boundary(v)) {
            out.add_vertex(in.point(v));          // boundary kept fixed
        } else {
            auto nb = in.vertex_vertices(v);
            const int n = static_cast<int>(nb.size());
            const double a = (4.0 - 2.0 * std::cos(2.0 * PI / n)) / 9.0;
            Vec3 sum{0, 0, 0};
            for (auto w : nb) sum += in.point(w);
            out.add_vertex(in.point(v) * (1.0 - a) + (sum / (double)n) * a);
        }
    }

    // --- 2. one centroid vertex per face, then 1 -> 3 split ----------------
    std::vector<VertexHandle> cp(nF);
    for (size_t fi = 0; fi < nF; ++fi)
        cp[fi] = out.add_vertex(in.calc_face_centroid(FaceHandle(static_cast<int>(fi))));

    for (size_t fi = 0; fi < nF; ++fi) {
        FaceHandle f(static_cast<int>(fi));
        auto he = in.face_halfedges(f);
        if (he.size() != 3) continue;             // sqrt3 is defined for triangles
        VertexHandle a(in.from_vertex(he[0]).idx());
        VertexHandle b(in.from_vertex(he[1]).idx());
        VertexHandle c(in.from_vertex(he[2]).idx());
        out.add_face({ a, b, cp[fi] });
        out.add_face({ b, c, cp[fi] });
        out.add_face({ c, a, cp[fi] });
    }

    // --- 3. flip every original interior edge ------------------------------
    for (size_t ei = 0; ei < in.n_edges(); ++ei) {
        EdgeHandle e(static_cast<int>(ei));
        if (in.is_boundary(e)) continue;          // leave the boundary alone
        HalfedgeHandle h = in.halfedge(e, 0);
        VertexHandle a(in.from_vertex(h).idx()), b(in.to_vertex(h).idx());
        HalfedgeHandle nh = out.find_halfedge(a, b);   // same a-b edge in `out`
        if (!nh.is_valid()) continue;
        EdgeHandle ne = out.edge(nh);
        if (out.is_flip_ok(ne)) out.flip(ne);
    }
    return out;
}

Mesh sqrt3_subdivide(const Mesh& in, int iterations) {
    Mesh m = sqrt3_subdivide_once(in);
    for (int k = 1; k < iterations; ++k)
        m = sqrt3_subdivide_once(m);
    return m;
}

// ============================================================================
//  QEM decimation (Garland & Heckbert)
// ----------------------------------------------------------------------------
//  1. Give every vertex a quadric = sum of the plane-quadrics of its faces.
//  2. For every edge, the contraction error is (Q0+Q1) evaluated at the best of
//     three candidate targets {endpoint0, endpoint1, midpoint}. (The full method
//     solves a 3x3 system for the optimal point; the candidate set is a robust,
//     allocation-light approximation that needs no matrix inverse.)
//  3. Repeatedly collapse the cheapest valid edge, merging quadrics, until the
//     face budget is met. A priority queue with per-edge version stamps discards
//     entries made stale by earlier collapses.
// ============================================================================

namespace {

struct QEntry {
    double cost;
    Vec3   target;
    int    edge;
    int    version;
};
struct QGreater { bool operator()(const QEntry& a, const QEntry& b) const { return a.cost > b.cost; } };

} // namespace

void quadric_decimate(Mesh& m, size_t target_faces) {
    // --- per-vertex quadrics ----------------------------------------------
    std::vector<Quadric> Q(m.n_vertices());
    for (auto f : m.all_faces()) {
        Vec3 n = m.calc_face_normal(f);
        auto fv = m.face_vertices(f);
        if (fv.empty()) continue;
        double pd = -(n.x * m.point(fv[0]).x + n.y * m.point(fv[0]).y + n.z * m.point(fv[0]).z);
        Quadric K = Quadric::from_plane(n.x, n.y, n.z, pd);
        for (auto v : fv) Q[v.idx()] = Q[v.idx()] + K;
    }

    // cost+target for an edge given the current geometry & quadrics
    auto evaluate = [&](EdgeHandle e) -> QEntry {
        HalfedgeHandle h = m.halfedge(e, 0);
        VertexHandle a = m.from_vertex(h), b = m.to_vertex(h);
        Quadric q = Q[a.idx()] + Q[b.idx()];
        Vec3 cand[3] = { m.point(a), m.point(b),
                         (m.point(a) + m.point(b)) * 0.5 };
        QEntry best{ std::numeric_limits<double>::max(), cand[0], e.idx(), 0 };
        for (const Vec3& p : cand) {
            double err = q.eval(p);
            if (err < best.cost) { best.cost = err; best.target = p; }
        }
        return best;
    };

    std::vector<int> version(m.n_edges(), 0);
    std::priority_queue<QEntry, std::vector<QEntry>, QGreater> pq;
    for (auto e : m.all_edges()) pq.push(evaluate(e));

    size_t live_faces = 0;
    for (auto f : m.all_faces()) { (void)f; ++live_faces; }

    // --- greedy collapse loop ---------------------------------------------
    while (live_faces > target_faces && !pq.empty()) {
        QEntry top = pq.top(); pq.pop();
        EdgeHandle e(top.edge);

        if (m.is_deleted(e) || version[e.idx()] != top.version) continue;  // stale

        // Pick an orientation we are allowed to collapse. collapse(h) removes
        // from_vertex(h) and keeps to_vertex(h).
        HalfedgeHandle h0 = m.halfedge(e, 0), h1 = m.halfedge(e, 1);
        HalfedgeHandle h;
        if      (m.is_collapse_ok(h0)) h = h0;
        else if (m.is_collapse_ok(h1)) h = h1;
        else continue;

        // how many faces this collapse will remove (triangles on each side)
        size_t removed = (m.face(h).is_valid() ? 1 : 0) +
                         (m.face(m.opposite_halfedge(h)).is_valid() ? 1 : 0);

        VertexHandle kept    = m.to_vertex(h);
        VertexHandle removed_v = m.from_vertex(h);
        Quadric merged = Q[kept.idx()] + Q[removed_v.idx()];

        m.collapse(h);
        m.set_point(kept, top.target);
        Q[kept.idx()] = merged;
        live_faces -= removed;

        // refresh every edge still incident to the kept vertex
        for (auto oh : m.voh_range(kept)) {
            EdgeHandle ie = m.edge(oh);
            ++version[ie.idx()];
            QEntry ne = evaluate(ie);
            ne.version = version[ie.idx()];
            pq.push(ne);
        }
    }

    m.garbage_collection();
}

// ============================================================================
//  Laplacian (umbrella) smoothing
// ----------------------------------------------------------------------------
//  Each interior vertex drifts toward the centroid of its neighbours. We
//  double-buffer (compute all new positions from the OLD ones, then commit) so
//  the result is independent of vertex ordering - a Jacobi iteration.
// ============================================================================

void laplacian_smooth(Mesh& m, int iterations, double lambda) {
    std::vector<Vec3> next(m.n_vertices());
    for (int it = 0; it < iterations; ++it) {
        for (auto v : m.all_vertices()) {
            if (m.is_boundary(v)) { next[v.idx()] = m.point(v); continue; }
            Vec3 avg{0, 0, 0};
            int n = 0;
            for (auto w : m.vv_range(v)) { avg += m.point(w); ++n; }
            next[v.idx()] = n ? m.point(v) + (avg / (double)n - m.point(v)) * lambda
                              : m.point(v);
        }
        for (auto v : m.all_vertices()) m.set_point(v, next[v.idx()]);
    }
}

// ============================================================================
//  Cotangent-weighted Laplacian smoothing
// ----------------------------------------------------------------------------
//  The uniform Laplacian above weights every neighbour equally, so it distorts
//  irregular triangulations. The cotangent weights w_ij = cot(alpha) + cot(beta)
//  - alpha,beta the angles OPPOSITE edge ij in its two triangles - make the
//  operator agree with the smooth Laplace-Beltrami operator, so smoothing no
//  longer depends on how the surface happens to be meshed.
// ============================================================================

// cot of the angle at the apex of h's triangle (the vertex that next(h) points
// to), i.e. the angle subtended opposite the edge of h. 0 on a boundary side.
static double cot_opposite_angle(const Mesh& m, HalfedgeHandle h) {
    if (m.is_boundary(h)) return 0.0;
    VertexHandle apex = m.to_vertex(m.next_halfedge(h));
    Vec3 a = m.point(m.from_vertex(h)) - m.point(apex);   // apex -> from
    Vec3 b = m.point(m.to_vertex(h))   - m.point(apex);   // apex -> to
    double denom = a.cross(b).norm();
    if (denom < 1e-20) return 0.0;
    return a.dot(b) / denom;                               // cos/sin = cot
}

void cotan_smooth(Mesh& m, int iterations, double lambda) {
    std::vector<Vec3> next(m.n_vertices());
    for (int it = 0; it < iterations; ++it) {
        for (auto v : m.all_vertices()) {
            if (m.is_boundary(v)) { next[v.idx()] = m.point(v); continue; }
            Vec3 acc{0, 0, 0};
            double wsum = 0.0;
            for (HalfedgeHandle h : m.vertex_outgoing_halfedges(v)) {
                // edge v->w has the two opposite angles in h's and opp(h)'s faces
                double w = cot_opposite_angle(m, h) +
                           cot_opposite_angle(m, m.opposite_halfedge(h));
                acc  += m.point(m.to_vertex(h)) * w;
                wsum += w;
            }
            next[v.idx()] = (std::abs(wsum) > 1e-20)
                          ? m.point(v) + (acc / wsum - m.point(v)) * lambda
                          : m.point(v);
        }
        for (auto v : m.all_vertices()) m.set_point(v, next[v.idx()]);
    }
}

// ============================================================================
//  Hole filling
// ----------------------------------------------------------------------------
//  A "hole" is a closed loop of boundary half-edges, which chain via
//  next_halfedge. We walk the loop, collect its vertices, and fan-triangulate
//  it. Boundary half-edges already wind with the hole on their left - exactly
//  how a face winds around its interior - so the fill faces reuse them directly,
//  in the SAME order (no reversal): add_face turns each boundary half-edge into
//  a real face half-edge.
// ============================================================================

bool fill_hole(Mesh& m, HalfedgeHandle boundary_h, int max_edges) {
    if (!m.is_boundary(boundary_h)) return false;

    std::vector<VertexHandle> loop;          // in boundary-traversal order
    HalfedgeHandle h = boundary_h;
    do {
        loop.push_back(m.from_vertex(h));
        h = m.next_halfedge(h);
    } while (h != boundary_h && loop.size() < m.n_halfedges() + 1);
    if (loop.size() < 3) return false;
    if (max_edges > 0 && static_cast<int>(loop.size()) > max_edges)
        return false;                        // loop too big - not a hole we fill

    bool any = false;
    for (size_t i = 1; i + 1 < loop.size(); ++i)
        if (m.add_face({ loop[0], loop[i], loop[i + 1] }).is_valid()) any = true;
    return any;
}

int fill_holes(Mesh& m, int max_edges) {
    int filled = 0;
    // Snapshot the half-edge count: filling adds faces but never removes the
    // boundary half-edges we still need to test (their face becomes valid, so
    // is_boundary turns false and we won't re-enter a filled loop).
    size_t nH = m.n_halfedges();
    for (size_t i = 0; i < nH; ++i) {
        HalfedgeHandle h(static_cast<int>(i));
        if (m.is_deleted(h) || !m.is_boundary(h)) continue;
        if (fill_hole(m, h, max_edges)) ++filled;
    }
    return filled;
}

// ============================================================================
//  Remove small connected components (stray islands / scan-print debris).
//  Union faces sharing an edge, count faces per component, delete the small
//  ones. In place; runs garbage_collection. Returns #components removed.
// ============================================================================
int remove_small_components(Mesh& m, size_t min_faces) {
    if (min_faces <= 1) return 0;
    const int nf = static_cast<int>(m.n_faces());
    std::vector<int> parent(nf);
    for (int i = 0; i < nf; ++i) parent[i] = i;
    std::function<int(int)> find = [&](int x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };

    for (int f = 0; f < nf; ++f) {
        FaceHandle fh(f);
        if (m.is_deleted(fh)) continue;
        for (FaceHandle g : m.face_faces(fh))         // edge-adjacent faces
            if (!m.is_deleted(g)) parent[find(f)] = find(g.idx());
    }
    std::vector<int> count(nf, 0);
    for (int f = 0; f < nf; ++f)
        if (!m.is_deleted(FaceHandle(f))) count[find(f)]++;

    int removed = 0;
    std::vector<char> root_removed(nf, 0);
    for (int f = 0; f < nf; ++f) {
        FaceHandle fh(f);
        if (m.is_deleted(fh)) continue;
        int r = find(f);
        if (count[r] < static_cast<int>(min_faces)) {
            if (!root_removed[r]) { root_removed[r] = 1; ++removed; }
            m.delete_face(fh, true);
        }
    }
    if (removed) m.garbage_collection();
    return removed;
}

// ============================================================================
//  Midpoint & Butterfly subdivision (1 -> 4, interpolating)
// ----------------------------------------------------------------------------
//  Both keep the OLD vertices exactly (interpolating) and add one vertex per
//  edge, then emit four sub-triangles per face - the same connectivity as Loop.
//  Only the edge-vertex rule differs, so they share one builder.
// ============================================================================

static Mesh build_1to4(const Mesh& in, const std::vector<Vec3>& edge_pos) {
    Mesh out;
    const size_t nV = in.n_vertices();
    const size_t nE = in.n_edges();

    for (size_t vi = 0; vi < nV; ++vi)                    // old vertices kept
        out.add_vertex(in.point(VertexHandle(static_cast<int>(vi))));

    std::vector<VertexHandle> emid(nE);
    for (size_t ei = 0; ei < nE; ++ei)
        emid[ei] = out.add_vertex(edge_pos[ei]);

    for (size_t fi = 0; fi < in.n_faces(); ++fi) {
        FaceHandle f(static_cast<int>(fi));
        auto he = in.face_halfedges(f);
        if (he.size() != 3) continue;
        VertexHandle a(in.from_vertex(he[0]).idx());
        VertexHandle b(in.from_vertex(he[1]).idx());
        VertexHandle c(in.from_vertex(he[2]).idx());
        VertexHandle mab = emid[in.edge(he[0]).idx()];
        VertexHandle mbc = emid[in.edge(he[1]).idx()];
        VertexHandle mca = emid[in.edge(he[2]).idx()];
        out.add_face({a,   mab, mca});
        out.add_face({mab, b,   mbc});
        out.add_face({mca, mbc, c});
        out.add_face({mab, mbc, mca});
    }
    return out;
}

static Mesh midpoint_once(const Mesh& in) {
    std::vector<Vec3> ep(in.n_edges());
    for (size_t ei = 0; ei < in.n_edges(); ++ei) {
        EdgeHandle e(static_cast<int>(ei));
        HalfedgeHandle h = in.halfedge(e, 0);
        ep[ei] = (in.point(in.from_vertex(h)) + in.point(in.to_vertex(h))) * 0.5;
    }
    return build_1to4(in, ep);
}

Mesh midpoint_subdivide(const Mesh& in, int iterations) {
    Mesh m = midpoint_once(in);
    for (int k = 1; k < iterations; ++k) m = midpoint_once(m);
    return m;
}

// Butterfly stencil for one interior edge h: a->b.
//   new = 1/2 (a+b) + 1/8 (c+d) - 1/16 (e1+e2+e3+e4)
// where c,d are the two opposite apexes and e1..e4 are the outer "wing" apexes.
// Falls back to the midpoint if any wing is missing (boundary / open fan).
static Vec3 butterfly_edge_point(const Mesh& m, HalfedgeHandle h) {
    HalfedgeHandle o = m.opposite_halfedge(h);
    Vec3 a = m.point(m.from_vertex(h));
    Vec3 b = m.point(m.to_vertex(h));
    if (m.is_boundary(h) || m.is_boundary(o))
        return (a + b) * 0.5;                          // boundary -> midpoint

    HalfedgeHandle t1 = m.next_halfedge(h);            // b->c
    HalfedgeHandle t2 = m.next_halfedge(t1);           // c->a
    HalfedgeHandle s1 = m.next_halfedge(o);            // a->d
    HalfedgeHandle s2 = m.next_halfedge(s1);           // d->b
    VertexHandle c = m.to_vertex(t1), d = m.to_vertex(s1);

    HalfedgeHandle w1 = m.opposite_halfedge(t1);       // wing across (b,c)
    HalfedgeHandle w2 = m.opposite_halfedge(t2);       // wing across (c,a)
    HalfedgeHandle w3 = m.opposite_halfedge(s1);       // wing across (a,d)
    HalfedgeHandle w4 = m.opposite_halfedge(s2);       // wing across (d,b)
    if (m.is_boundary(w1) || m.is_boundary(w2) ||
        m.is_boundary(w3) || m.is_boundary(w4))        // irregular -> drop wings
        return (a + b) * 0.5 + (m.point(c) + m.point(d)) * 0.125;

    Vec3 e1 = m.point(m.to_vertex(m.next_halfedge(w1)));
    Vec3 e2 = m.point(m.to_vertex(m.next_halfedge(w2)));
    Vec3 e3 = m.point(m.to_vertex(m.next_halfedge(w3)));
    Vec3 e4 = m.point(m.to_vertex(m.next_halfedge(w4)));

    return (a + b) * 0.5 + (m.point(c) + m.point(d)) * 0.125
         - (e1 + e2 + e3 + e4) * (1.0 / 16.0);
}

static Mesh butterfly_once(const Mesh& in) {
    std::vector<Vec3> ep(in.n_edges());
    for (size_t ei = 0; ei < in.n_edges(); ++ei)
        ep[ei] = butterfly_edge_point(in, in.halfedge(EdgeHandle(static_cast<int>(ei)), 0));
    return build_1to4(in, ep);
}

Mesh butterfly_subdivide(const Mesh& in, int iterations) {
    Mesh m = butterfly_once(in);
    for (int k = 1; k < iterations; ++k) m = butterfly_once(m);
    return m;
}

// ============================================================================
//  Longest-edge (Rivara) refinement, in place
// ----------------------------------------------------------------------------
//  Repeatedly bisect every edge longer than the bound. split(edge, v) cuts BOTH
//  incident triangles, so the result is always conformal (no hanging nodes). We
//  rescan after each pass because splitting creates new edges; each split halves
//  an edge, so the loop terminates once all edges are within the bound.
// ============================================================================

void longest_edge_subdivide(Mesh& m, double max_edge_length) {
    const double thr2 = max_edge_length * max_edge_length;
    bool changed = true;
    while (changed) {
        changed = false;
        size_t nE = m.n_edges();                 // only scan edges that exist now
        for (size_t ei = 0; ei < nE; ++ei) {
            EdgeHandle e(static_cast<int>(ei));
            if (m.is_deleted(e)) continue;
            HalfedgeHandle h = m.halfedge(e, 0);
            Vec3 d = m.point(m.to_vertex(h)) - m.point(m.from_vertex(h));
            if (d.sqrnorm() <= thr2) continue;
            VertexHandle mid = m.add_vertex(m.calc_edge_midpoint(e));
            m.split(e, mid);
            changed = true;
        }
    }
}

} // namespace sm
