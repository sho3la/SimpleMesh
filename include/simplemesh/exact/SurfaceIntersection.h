// ============================================================================
//  simplemesh/exact/SurfaceIntersection.h - surface arrangement assembly
// ----------------------------------------------------------------------------
//  Resolves all triangle/triangle intersections of a triangle SOUP into an exact
//  surface arrangement, where every intersection curve has become a shared mesh
//  edge. The pipeline:
//
//    1. Broad phase  - AABB sweep-and-prune candidate pairs (the same approach
//       as the self-intersection detector in MeshChecker.cpp).
//    2. Symbolic intersection (TriangleTriangleIntersection) for each candidate
//       pair -> a set of (TriangleRegion,TriangleRegion) records that name the
//       intersection points combinatorially.
//    3. Per-facet exact remesh: a FacetMesher inserts those records as
//       constraints into an exact-coordinate constrained Delaunay triangulation;
//       constraint crossings become new vertices whose 3D coordinates are
//       computed EXACTLY (plane-line / three-planes) with NO DIVISION
//       (homogeneous coordinates).
//    4. Global merge: every output vertex (original corner or intersection
//       point) is looked up in a std::map keyed by its EXACT 3D coordinate, so
//       coincident points - including the triple points where three facets meet,
//       computed independently on each facet - collapse to one index by exact
//       equality (NOT tolerance). This is what turns interpenetration into
//       shared edges.
//
//  The result is a triangle soup that may be NON-manifold along intersection
//  curves (that is the correct arrangement; shell selection lives in
//  MeshBoolean.h). Because it can be non-manifold it is returned as a soup, not a
//  half-edge mesh.
//
//  Compile strict-FP (-ffp-contract=off / /fp:strict).
// ============================================================================
#pragma once

#include "CDT2d.h"
#include "HomogeneousGeometry.h"
#include "Predicates.h"
#include "TriangleIntersection.h"

#include <array>
#include <map>
#include <vector>
#include <algorithm>
#include <cstddef>
#include <functional>
#include <numeric>

namespace sm {
namespace exact {

// ----------------------------------------------------------------------------
//  Exact 3D arithmetic helpers (a plain ExpansionNt x/y/z triple).
// ----------------------------------------------------------------------------
struct EVec3 { ExpansionNt x, y, z; };

inline EVec3 ev_from(const double* p) { return EVec3{ ExpansionNt(p[0]), ExpansionNt(p[1]), ExpansionNt(p[2]) }; }
// difference b - a as an exact vector
inline EVec3 ev_sub(const double* a, const double* b) {
    return EVec3{ ExpansionNt(b[0]) - ExpansionNt(a[0]),
                  ExpansionNt(b[1]) - ExpansionNt(a[1]),
                  ExpansionNt(b[2]) - ExpansionNt(a[2]) };
}
inline EVec3 ev_cross(const EVec3& a, const EVec3& b) {
    return EVec3{ a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x };
}
inline ExpansionNt ev_dot(const EVec3& a, const EVec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

// Coordinate accessor on a homogeneous 3D point.
inline const ExpansionNt& he_coord(const Vec3HE& p, int a) {
    return (a == 0) ? p.x : (a == 1) ? p.y : p.z;
}
inline Vec3HE he_from_double(const double* p) {
    return Vec3HE(ExpansionNt(p[0]), ExpansionNt(p[1]), ExpansionNt(p[2]), ExpansionNt(1.0));
}

// Homogeneous point at rational parameter t = tnum/tden along p1->p2, with p1,p2
// given as Cartesian doubles. No division (homogeneous result).
//   P = ((tden - tnum)*p1 + tnum*p2) / tden
inline Vec3HE mix3(const ExpansionNt& tnum, const ExpansionNt& tden,
                   const double* p1, const double* p2) {
    ExpansionNt sn = tden - tnum;
    return Vec3HE(sn*ExpansionNt(p1[0]) + tnum*ExpansionNt(p2[0]),
                  sn*ExpansionNt(p1[1]) + tnum*ExpansionNt(p2[1]),
                  sn*ExpansionNt(p1[2]) + tnum*ExpansionNt(p2[2]),
                  tden);
}

// Exact plane(p1,p2,p3) /\ line(q1,q2), Moller-Trumbore form.
// Precondition: the intersection exists (d != 0).
inline Vec3HE plane_line_intersection(const double* p1, const double* p2, const double* p3,
                                      const double* q1, const double* q2) {
    EVec3 D  = ev_sub(q1, q2);   // q2 - q1
    EVec3 E1 = ev_sub(p1, p2);   // p2 - p1
    EVec3 E2 = ev_sub(p1, p3);   // p3 - p1
    EVec3 AO = ev_sub(p1, q1);   // q1 - p1
    EVec3 N  = ev_cross(E1, E2);
    ExpansionNt d = -ev_dot(D, N);
    ExpansionNt tnum = ev_dot(AO, N);
    return mix3(tnum, d, q1, q2);
}

// Exact intersection of three supporting planes. Returns false iff the planes
// do not meet in a single point (then the caller falls back to a 2D edge-edge
// intersection).
inline bool three_planes_intersection(Vec3HE& result,
    const double* p1, const double* p2, const double* p3,
    const double* q1, const double* q2, const double* q3,
    const double* r1, const double* r2, const double* r3) {
    EVec3 N1 = ev_cross(ev_sub(p1, p2), ev_sub(p1, p3));
    EVec3 N2 = ev_cross(ev_sub(q1, q2), ev_sub(q1, q3));
    EVec3 N3 = ev_cross(ev_sub(r1, r2), ev_sub(r1, r3));
    EVec3 B{ ev_dot(N1, ev_from(p1)), ev_dot(N2, ev_from(q1)), ev_dot(N3, ev_from(r1)) };

    ExpansionNt w = det3x3<ExpansionNt>(N1.x, N1.y, N1.z,
                                        N2.x, N2.y, N2.z,
                                        N3.x, N3.y, N3.z);
    if (w.sign() == 0) return false;

    ExpansionNt x = det3x3<ExpansionNt>(B.x, N1.y, N1.z,
                                        B.y, N2.y, N2.z,
                                        B.z, N3.y, N3.z);
    ExpansionNt y = det3x3<ExpansionNt>(N1.x, B.x, N1.z,
                                        N2.x, B.y, N2.z,
                                        N3.x, B.z, N3.z);
    ExpansionNt z = det3x3<ExpansionNt>(N1.x, N1.y, B.x,
                                        N2.x, N2.y, B.y,
                                        N3.x, N3.y, B.z);
    result = Vec3HE(x, y, z, w);
    return true;
}

// ----------------------------------------------------------------------------
//  Exact lexicographic comparison of homogeneous 3D points: orders by
//  (x/w, y/w, z/w) using cross-multiplication.
// ----------------------------------------------------------------------------
// sign of n1/d1 - n2/d2.
inline int ratio_compare(const ExpansionNt& n1, const ExpansionNt& d1,
                         const ExpansionNt& n2, const ExpansionNt& d2) {
    ExpansionNt num = n1*d2 - n2*d1;
    return num.sign() * d1.sign() * d2.sign();
}
struct ExactPointCompare {
    bool operator()(const Vec3HE& a, const Vec3HE& b) const {
        int s = ratio_compare(a.x, a.w, b.x, b.w); if (s) return s < 0;
        s = ratio_compare(a.y, a.w, b.y, b.w);     if (s) return s < 0;
        s = ratio_compare(a.z, a.w, b.z, b.w);     return s < 0;
    }
};

// ----------------------------------------------------------------------------
//  Symbolic intersection record: an intersection of facet f1 with facet f2,
//  either a single point (A==B) or a segment with endpoints A,B.
// ----------------------------------------------------------------------------
struct IsectInfo {
    std::size_t f1, f2;
    TriangleRegion A_rgn_f1, A_rgn_f2;
    TriangleRegion B_rgn_f1, B_rgn_f2;
    bool is_point() const { return A_rgn_f1 == B_rgn_f1 && A_rgn_f2 == B_rgn_f2; }
    void flip() {
        std::swap(f1, f2);
        A_rgn_f1 = swap_T1_T2(A_rgn_f1); A_rgn_f2 = swap_T1_T2(A_rgn_f2);
        B_rgn_f1 = swap_T1_T2(B_rgn_f1); B_rgn_f2 = swap_T1_T2(B_rgn_f2);
        std::swap(A_rgn_f1, A_rgn_f2);
        std::swap(B_rgn_f1, B_rgn_f2);
    }
};

// The lowest-dimensional region of a single triangle that contains both R1 and
// R2 (both regions on the same triangle).
inline TriangleRegion regions_convex_hull(TriangleRegion R1, TriangleRegion R2) {
    if (R1 == R2) return R1;
    TriangleRegion R = is_in_T1(R1) ? T1_RGN_T : T2_RGN_T;
    if (region_dim(R1) == 1 && region_dim(R2) == 0) {
        TriangleRegion v1, v2; get_edge_vertices(R1, v1, v2);
        if (R2 == v1 || R2 == v2) R = R1;
    } else if (region_dim(R2) == 1 && region_dim(R1) == 0) {
        TriangleRegion v1, v2; get_edge_vertices(R2, v1, v2);
        if (R1 == v1 || R1 == v2) R = R2;
    } else if (region_dim(R1) == 0 && region_dim(R2) == 0) {
        for (TriangleRegion E : { T1_RGN_E0, T1_RGN_E1, T1_RGN_E2,
                                  T2_RGN_E0, T2_RGN_E1, T2_RGN_E2 }) {
            TriangleRegion v1, v2; get_edge_vertices(E, v1, v2);
            if ((R1 == v1 && R2 == v2) || (R1 == v2 && R2 == v1)) { R = E; break; }
        }
    }
    return R;
}

// ----------------------------------------------------------------------------
//  FacetMesher - remesh one facet f1 with the constraints coming from its
//  intersections with other facets. Exact-coordinate CDT; every vertex carries
//  its EXACT 3D point so the global merge is exact.
// ----------------------------------------------------------------------------
class FacetMesher : public CDTBase2d {
public:
    using Soup_pts = std::vector<std::array<double,3>>;
    using Soup_tri = std::vector<std::array<std::size_t,3>>;

    FacetMesher(const Soup_pts& V, const Soup_tri& F) : V_(V), F_(F) {
        exact_incircle_ = true;
        exact_intersections_ = true;
    }

    void clear() override {
        CDTBase2d::clear();
        pt3_.clear(); l_.clear(); edges_.clear();
    }

    void begin_facet(std::size_t f) {
        f1_ = f;
        const double* P0 = fv(f, 0);
        const double* P1 = fv(f, 1);
        const double* P2 = fv(f, 2);
        nax_ = triangle_normal_axis(P0, P1, P2);
        u_ = (nax_ + 1) % 3;
        v_ = (nax_ + 2) % 3;
        push_vertex(he_from_double(P0));
        push_vertex(he_from_double(P1));
        push_vertex(he_from_double(P2));
        // Keep the macro triangle CCW in the projection so orient_012_>0
        // swap the projection axes if the macro triangle comes out CW.
        if (exact::orient2d(proj(0), proj(1), proj(2)) < 0) {
            std::swap(u_, v_);
            recompute_lengths();
        }
        CDTBase2d::create_enclosing_triangle(0, 1, 2);
        // The 3 facet boundary edges become constraints (ids 0,1,2). They carry
        // no f2 (NO_INDEX): a constraint never crosses them in a facet interior.
        edges_.push_back(Edge{1, 2, NO_F, T2_RGN_T});
        edges_.push_back(Edge{2, 0, NO_F, T2_RGN_T});
        edges_.push_back(Edge{0, 1, NO_F, T2_RGN_T});
    }

    // One intersection point named combinatorially.
    index_t add_vertex(std::size_t f2, TriangleRegion R1, TriangleRegion R2) {
        if (region_dim(R1) == 0) return index_t(R1);   // a corner of f1 (0,1,2)
        push_vertex(compute_geometry(f2, R1, R2));
        index_t v = CDTBase2d::insert(pt3_.size() - 1);
        if (pt3_.size() > nv()) pop_vertex();           // was a duplicate
        return v;
    }

    // One intersection segment.
    void add_edge(std::size_t f2,
                  TriangleRegion AR1, TriangleRegion AR2,
                  TriangleRegion BR1, TriangleRegion BR2) {
        index_t v1 = add_vertex(f2, AR1, AR2);
        index_t v2 = add_vertex(f2, BR1, BR2);
        // If both extremities are on the same edge of f1, the edge is generated
        // when that boundary edge is remeshed; don't add it as a constraint.
        if (region_dim(regions_convex_hull(AR1, BR1)) == 1) return;
        edges_.push_back(Edge{v1, v2, f2, regions_convex_hull(AR2, BR2)});
    }

    // Insert all constraints (ids match edges_ order) into the triangulation.
    void commit() {
        for (const Edge& E : edges_) CDTBase2d::insert_constraint(E.v1, E.v2);
    }

    // ---- read-back ----------------------------------------------------------
    using CDTBase2d::nT;
    using CDTBase2d::nv;
    index_t tri_v(index_t t, index_t lv) const { return Tv(t, lv); }
    const Vec3HE& point(index_t v) const { return pt3_[v]; }

protected:
    Sign orient2d(index_t i, index_t j, index_t k) const override {
        return Sign(exact::orient2d(proj(i), proj(j), proj(k)));
    }
    Sign incircle(index_t i, index_t j, index_t k, index_t l) const override {
        return Sign(exact::incircle_2d(proj(i), proj(j), proj(k), proj(l),
                                       l_[i], l_[j], l_[k], l_[l]));
    }
    // A constraint crossing inside the facet: a SECONDARY intersection vertex.
    // It is the meeting of the three facet planes f1, f2(E1), f3(E2) -> exact 3D
    // with no division.
    index_t create_intersection(index_t E1, index_t /*i*/, index_t /*j*/,
                                index_t E2, index_t /*k*/, index_t /*l*/) override {
        Vec3HE I;
        get_edge_edge_intersection(E1, E2, I);
        push_vertex(I);
        index_t v = nv_;
        v2T_.push_back(CDT_NO_INDEX);
        ++nv_;
        return v;
    }

private:
    static constexpr std::size_t NO_F = std::size_t(-1);
    struct Edge { index_t v1, v2; std::size_t f2; TriangleRegion R2; };

    const double* fv(std::size_t f, int lv) const { return V_[F_[f][lv]].data(); }

    Vec2HE proj(index_t v) const {
        return Vec2HE(he_coord(pt3_[v], u_), he_coord(pt3_[v], v_), pt3_[v].w);
    }
    double proj_length(const Vec3HE& p) const {
        double x = he_coord(p, u_).estimate();
        double y = he_coord(p, v_).estimate();
        double w = p.w.estimate();
        return (x*x + y*y) / (w*w);
    }
    void push_vertex(const Vec3HE& p) { pt3_.push_back(p); l_.push_back(proj_length(p)); }
    void pop_vertex() { pt3_.pop_back(); l_.pop_back(); }
    void recompute_lengths() { for (std::size_t i = 0; i < pt3_.size(); ++i) l_[i] = proj_length(pt3_[i]); }

    // The exact 3D location of a primary intersection vertex (on facet f1_, from
    // its intersection with facet f2, named by regions R1 on f1 and R2 on f2).
    Vec3HE compute_geometry(std::size_t f2, TriangleRegion R1, TriangleRegion R2) const {
        // case 1: f1 vertex (handled by add_vertex short-circuit, kept for safety)
        if (region_dim(R1) == 0) {
            int lv = int(R1) - int(T1_RGN_P0);
            return he_from_double(fv(f1_, lv));
        }
        // case 2: f2 vertex
        if (region_dim(R2) == 0) {
            int lv = int(R2) - int(T2_RGN_P0);
            return he_from_double(fv(f2, lv));
        }
        // case 3: f1 (face or edge) /\ f2 edge, in 3D
        if ((region_dim(R1) == 2 || region_dim(R1) == 1) && region_dim(R2) == 1) {
            const double* p1 = fv(f1_, 0);
            const double* p2 = fv(f1_, 1);
            const double* p3 = fv(f1_, 2);
            int e = int(R2) - int(T2_RGN_E0);
            const double* q1 = fv(f2, (e+1)%3);
            const double* q2 = fv(f2, (e+2)%3);
            bool seg_seg_2D = (region_dim(R1) == 1 &&
                               orient3d(p1, p2, p3, q1) == 0 &&
                               orient3d(p1, p2, p3, q2) == 0);
            if (!seg_seg_2D) return plane_line_intersection(p1, p2, p3, q1, q2);
        }
        // case 4: f1 edge /\ f2 face, in 3D
        if (region_dim(R1) == 1 && region_dim(R2) == 2) {
            const double* p1 = fv(f2, 0);
            const double* p2 = fv(f2, 1);
            const double* p3 = fv(f2, 2);
            int e = int(R1) - int(T1_RGN_E0);
            const double* q1 = fv(f1_, (e+1)%3);
            const double* q2 = fv(f1_, (e+2)%3);
            return plane_line_intersection(p1, p2, p3, q1, q2);
        }
        // case 5: f1 edge /\ f2 edge, in 2D (coplanar) - intersect in the (u,v)
        // projection, mix along the 3D edge of f1.
        {
            int e1 = int(R1) - int(T1_RGN_E0);
            int e2 = int(R2) - int(T2_RGN_E0);
            const double* P1 = fv(f1_, (e1+1)%3);
            const double* P2 = fv(f1_, (e1+2)%3);
            const double* Q1 = fv(f2, (e2+1)%3);
            const double* Q2 = fv(f2, (e2+2)%3);
            ExpansionNt tnum, tden;
            seg_seg_2D_param(P1, P2, Q1, Q2, tnum, tden);
            return mix3(tnum, tden, P1, P2);
        }
    }

    // Exact location of a secondary (constraint x constraint) vertex: the meeting
    // of three facet planes, with a 2D fallback when they do not meet in a point.
    void get_edge_edge_intersection(index_t e1, index_t e2, Vec3HE& I) const {
        std::size_t f1 = f1_;
        std::size_t f2 = edges_[e1].f2;
        std::size_t f3 = edges_[e2].f2;
        if (f2 != NO_F && f3 != NO_F &&
            three_planes_intersection(I,
                fv(f1,0), fv(f1,1), fv(f1,2),
                fv(f2,0), fv(f2,1), fv(f2,2),
                fv(f3,0), fv(f3,1), fv(f3,2))) {
            return;
        }
        get_edge_edge_intersection_2D(e1, e2, I);
    }

    // 2D fallback for a coplanar secondary crossing.
    void get_edge_edge_intersection_2D(index_t e1, index_t e2, Vec3HE& I) const {
        const Edge& E1 = edges_[e1];
        const Edge& E2 = edges_[e2];
        if (region_dim(E1.R2) == 1 && region_dim(E2.R2) == 1) {
            int le1 = int(E1.R2) - int(T2_RGN_E0);
            int le2 = int(E2.R2) - int(T2_RGN_E0);
            const double* p1 = fv(E1.f2, (le1+1)%3);
            const double* p2 = fv(E1.f2, (le1+2)%3);
            const double* q1 = fv(E2.f2, (le2+1)%3);
            const double* q2 = fv(E2.f2, (le2+2)%3);
            // C1 = p2-p1, C2 = q1-q2, B = q1-p1; t = det(B,C2)/det(C1,C2)
            ExpansionNt c1x = ExpansionNt(p2[u_])-ExpansionNt(p1[u_]);
            ExpansionNt c1y = ExpansionNt(p2[v_])-ExpansionNt(p1[v_]);
            ExpansionNt c2x = ExpansionNt(q1[u_])-ExpansionNt(q2[u_]);
            ExpansionNt c2y = ExpansionNt(q1[v_])-ExpansionNt(q2[v_]);
            ExpansionNt bx  = ExpansionNt(q1[u_])-ExpansionNt(p1[u_]);
            ExpansionNt by  = ExpansionNt(q1[v_])-ExpansionNt(p1[v_]);
            ExpansionNt tden = det2x2(c1x, c1y, c2x, c2y);
            ExpansionNt tnum = det2x2(bx, by, c2x, c2y);
            I = mix3(tnum, tden, p1, p2);
        } else {
            std::size_t fa = E1.f2; TriangleRegion Ra = E1.R2;
            std::size_t fb = E2.f2; TriangleRegion Rb = E2.R2;
            if (region_dim(Ra) == 1) { std::swap(fa, fb); std::swap(Ra, Rb); }
            int e = int(Rb) - int(T2_RGN_E0);
            I = plane_line_intersection(fv(fa,0), fv(fa,1), fv(fa,2),
                                        fv(fb,(e+1)%3), fv(fb,(e+2)%3));
        }
    }

    // 2D segment-segment crossing parameter t = tnum/tden along P1->P2, in the
    // facet (u,v) projection.
    void seg_seg_2D_param(const double* P1, const double* P2,
                          const double* Q1, const double* Q2,
                          ExpansionNt& tnum, ExpansionNt& tden) const {
        ExpansionNt d1x = ExpansionNt(P2[u_])-ExpansionNt(P1[u_]);
        ExpansionNt d1y = ExpansionNt(P2[v_])-ExpansionNt(P1[v_]);
        ExpansionNt d2x = ExpansionNt(Q2[u_])-ExpansionNt(Q1[u_]);
        ExpansionNt d2y = ExpansionNt(Q2[v_])-ExpansionNt(Q1[v_]);
        ExpansionNt aox = ExpansionNt(Q1[u_])-ExpansionNt(P1[u_]);
        ExpansionNt aoy = ExpansionNt(Q1[v_])-ExpansionNt(P1[v_]);
        tden = det2x2(d1x, d1y, d2x, d2y);
        tnum = det2x2(aox, aoy, d2x, d2y);
    }

    const Soup_pts& V_;
    const Soup_tri& F_;
    std::size_t f1_ = NO_F;
    int nax_ = 2, u_ = 0, v_ = 1;
    std::vector<Vec3HE> pt3_;     // exact 3D coordinate of every CDT vertex
    std::vector<double> l_;       // projected squared length (for incircle)
    std::vector<Edge>   edges_;   // constraints, indexed by constraint id
};

// ----------------------------------------------------------------------------
//  The arrangement driver.
// ----------------------------------------------------------------------------
struct ArrangementResult {
    std::vector<std::array<double,3>>      points;       // rounded output verts
    std::vector<std::array<std::size_t,3>> triangles;    // output soup
    std::vector<Vec3HE>                    exact_points;  // parallel to points
    std::vector<char>                      is_original;   // parallel to points:
                                                          //   1 = original input vertex
    std::vector<std::size_t>               triangle_source;  // parallel to triangles:
                                                             //   originating input facet
    std::size_t                            n_input_vertices = 0;
    std::size_t                            n_intersection_pairs = 0;
};

namespace detail {

// Exact-vertex deduplication map: find or create a vertex by exact coordinate.
class ExactVertexSet {
public:
    explicit ExactVertexSet(ArrangementResult& R) : R_(R) {}
    // 'original' marks an actual input-mesh vertex (as opposed to a created
    // intersection vertex); used later by shell classification.
    std::size_t find_or_create(const Vec3HE& p, bool original) {
        auto it = map_.find(p);
        if (it != map_.end()) {
            if (original) R_.is_original[it->second] = 1;
            return it->second;
        }
        std::size_t idx = R_.points.size();
        double w = p.w.estimate();
        R_.points.push_back({ p.x.estimate()/w, p.y.estimate()/w, p.z.estimate()/w });
        R_.exact_points.push_back(p);
        R_.is_original.push_back(original ? 1 : 0);
        map_.emplace(p, idx);
        return idx;
    }
private:
    ArrangementResult& R_;
    std::map<Vec3HE, std::size_t, ExactPointCompare> map_;
};

// Connected components of the input soup (= the boolean operands). Returns the
// operand id of each facet and sets nb_operands.
inline std::vector<int> operand_of_facet(
    const std::vector<std::array<double,3>>& V,
    const std::vector<std::array<std::size_t,3>>& F, int& nb_operands) {
    std::vector<std::size_t> parent(V.size());
    std::iota(parent.begin(), parent.end(), std::size_t(0));
    std::function<std::size_t(std::size_t)> find = [&](std::size_t x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    auto uni = [&](std::size_t a, std::size_t b) { parent[find(a)] = find(b); };
    for (const auto& f : F) { uni(f[0], f[1]); uni(f[1], f[2]); }
    std::map<std::size_t,int> root_id;
    std::vector<int> op(F.size());
    nb_operands = 0;
    for (std::size_t i = 0; i < F.size(); ++i) {
        std::size_t r = find(F[i][0]);
        auto it = root_id.find(r);
        if (it == root_id.end()) it = root_id.emplace(r, nb_operands++).first;
        op[i] = it->second;
    }
    return op;
}

// AABB sweep-and-prune candidate pairs (SimpleMesh MeshChecker.cpp:309 idea).
inline std::vector<std::pair<std::size_t,std::size_t>> broad_phase(
    const std::vector<std::array<double,3>>& V,
    const std::vector<std::array<std::size_t,3>>& F) {
    struct Box { std::size_t f; std::array<double,3> lo, hi; };
    std::vector<Box> b(F.size());
    for (std::size_t f = 0; f < F.size(); ++f) {
        b[f].f = f;
        b[f].lo = { 1e300, 1e300, 1e300 };
        b[f].hi = { -1e300, -1e300, -1e300 };
        for (int k = 0; k < 3; ++k)
            for (int a = 0; a < 3; ++a) {
                double c = V[F[f][k]][a];
                b[f].lo[a] = std::min(b[f].lo[a], c);
                b[f].hi[a] = std::max(b[f].hi[a], c);
            }
    }
    std::sort(b.begin(), b.end(), [](const Box& a, const Box& c){ return a.lo[0] < c.lo[0]; });
    std::vector<std::pair<std::size_t,std::size_t>> pairs;
    for (std::size_t i = 0; i < b.size(); ++i)
        for (std::size_t j = i+1; j < b.size(); ++j) {
            if (b[j].lo[0] > b[i].hi[0]) break;
            bool ov = true;
            for (int a = 0; a < 3 && ov; ++a)
                ov = b[i].lo[a] <= b[j].hi[a] && b[j].lo[a] <= b[i].hi[a];
            if (!ov) continue;
            std::size_t f1 = b[i].f, f2 = b[j].f;
            pairs.emplace_back(std::min(f1,f2), std::max(f1,f2));
        }
    return pairs;
}

} // namespace detail

inline ArrangementResult resolve_intersections(
    const std::vector<std::array<double,3>>& V,
    const std::vector<std::array<std::size_t,3>>& F) {

    ArrangementResult R;
    R.n_input_vertices = V.size();

    // ---- phase 1+2: broad phase, then exact symbolic intersection ----------
    auto pairs = detail::broad_phase(V, F);
    std::vector<IsectInfo> intersections;
    std::vector<char> facet_has_isect(F.size(), 0);

    for (auto pr : pairs) {
        std::size_t f1 = pr.first, f2 = pr.second;
        const double* p0 = V[F[f1][0]].data();
        const double* p1 = V[F[f1][1]].data();
        const double* p2 = V[F[f1][2]].data();
        const double* q0 = V[F[f2][0]].data();
        const double* q1 = V[F[f2][1]].data();
        const double* q2 = V[F[f2][2]].data();
        std::vector<TriangleIsect> I;
        bool nd = triangles_intersections(p0, p1, p2, q0, q1, q2, I);
        if (!nd || I.empty()) continue;   // adjacency / no proper intersection
        ++R.n_intersection_pairs;
        facet_has_isect[f1] = facet_has_isect[f2] = 1;

        auto push_both = [&](IsectInfo II) {
            intersections.push_back(II);
            II.flip();
            intersections.push_back(II);
        };

        if (I.size() > 2) {
            // Coplanar overlap: edges are valid pairs whose convex hull stays
            // 1D on f1 or on f2 (combinatorial convex hull).
            for (std::size_t i1 = 0; i1 < I.size(); ++i1)
                for (std::size_t i2 = 0; i2 < i1; ++i2) {
                    IsectInfo II{ f1, f2, I[i1].first, I[i1].second,
                                          I[i2].first, I[i2].second };
                    TriangleRegion AB1 = regions_convex_hull(II.A_rgn_f1, II.B_rgn_f1);
                    TriangleRegion AB2 = regions_convex_hull(II.A_rgn_f2, II.B_rgn_f2);
                    if (region_dim(AB1) == 1 || region_dim(AB2) == 1) push_both(II);
                }
        } else {
            TriangleRegion A1 = I[0].first,  A2 = I[0].second;
            TriangleRegion B1 = A1,          B2 = A2;
            if (I.size() == 2) { B1 = I[1].first; B2 = I[1].second; }
            push_both(IsectInfo{ f1, f2, A1, A2, B1, B2 });
        }
    }

    // ---- group intersections by f1 -----------------------------------------
    std::sort(intersections.begin(), intersections.end(),
              [](const IsectInfo& a, const IsectInfo& b){ return a.f1 < b.f1; });

    // ---- output: original (non-intersected) facets pass through ------------
    detail::ExactVertexSet vset(R);
    auto emit_original = [&](std::size_t f) {
        std::array<std::size_t,3> t;
        for (int k = 0; k < 3; ++k)
            t[k] = vset.find_or_create(he_from_double(V[F[f][k]].data()), true);
        R.triangles.push_back(t);
        R.triangle_source.push_back(f);
    };
    for (std::size_t f = 0; f < F.size(); ++f)
        if (!facet_has_isect[f]) emit_original(f);

    // ---- phase 3: remesh each intersected facet ----------------------------
    FacetMesher MIT(V, F);
    std::size_t i = 0;
    while (i < intersections.size()) {
        std::size_t f1 = intersections[i].f1;
        std::size_t j = i;
        MIT.clear();
        MIT.begin_facet(f1);
        while (j < intersections.size() && intersections[j].f1 == f1) {
            const IsectInfo& II = intersections[j];
            if (II.is_point())
                MIT.add_vertex(II.f2, II.A_rgn_f1, II.A_rgn_f2);
            else
                MIT.add_edge(II.f2, II.A_rgn_f1, II.A_rgn_f2,
                                    II.B_rgn_f1, II.B_rgn_f2);
            ++j;
        }
        MIT.commit();
        // Map this facet's sub-triangles into the global exact-merged soup.
        // Vertices 0,1,2 are the facet's original corners (is_original).
        std::vector<std::size_t> g(MIT.nv());
        for (FacetMesher::index_t v = 0; v < MIT.nv(); ++v)
            g[v] = vset.find_or_create(MIT.point(v), v < 3);
        // Source-facet normal: every sub-triangle must inherit the SAME winding
        // as the input facet (the shell-classification model relies on it). The CDT
        // emits triangles in its own projected orientation, so flip any whose 3D
        // normal opposes the source facet's.
        const double* S0 = V[F[f1][0]].data();
        const double* S1 = V[F[f1][1]].data();
        const double* S2 = V[F[f1][2]].data();
        double ns[3] = {
            (S1[1]-S0[1])*(S2[2]-S0[2]) - (S1[2]-S0[2])*(S2[1]-S0[1]),
            (S1[2]-S0[2])*(S2[0]-S0[0]) - (S1[0]-S0[0])*(S2[2]-S0[2]),
            (S1[0]-S0[0])*(S2[1]-S0[1]) - (S1[1]-S0[1])*(S2[0]-S0[0]) };
        for (FacetMesher::index_t t = 0; t < MIT.nT(); ++t) {
            std::array<std::size_t,3> tri = { g[MIT.tri_v(t,0)],
                                              g[MIT.tri_v(t,1)],
                                              g[MIT.tri_v(t,2)] };
            if (tri[0] == tri[1] || tri[1] == tri[2] || tri[2] == tri[0]) continue;
            const auto& a = R.points[tri[0]];
            const auto& b = R.points[tri[1]];
            const auto& c = R.points[tri[2]];
            double nt[3] = {
                (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]),
                (b[2]-a[2])*(c[0]-a[0]) - (b[0]-a[0])*(c[2]-a[2]),
                (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]) };
            if (ns[0]*nt[0] + ns[1]*nt[1] + ns[2]*nt[2] < 0.0)
                std::swap(tri[1], tri[2]);
            R.triangles.push_back(tri);
            R.triangle_source.push_back(f1);
        }
        i = j;
    }

    return R;
}

} // namespace exact
} // namespace sm
