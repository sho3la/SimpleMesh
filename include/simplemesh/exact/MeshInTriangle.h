// ============================================================================
//  simplemesh/exact/MeshInTriangle.h - exact-coordinate CDT + per-facet remesh
// ----------------------------------------------------------------------------
//  Two pieces used by the per-facet exact remesh:
//
//   * ExactCDT2d : CDTBase2d - a constrained Delaunay triangulation whose
//     vertices carry EXACT homogeneous 2D coordinates (Vec2HE). It overrides the
//     three geometry hooks with the exact homogeneous predicates and builds
//     intersection vertices with NO DIVISION (exact 2D line-line cross + mix).
//     exact_incircle_ and exact_intersections_ are both true.
//
//   * MeshInTriangle - remeshes ONE 3D triangle facet: it projects the facet to
//     exact 2D (dropping the normal axis), uses the facet's 3 corners as the CDT
//     macro-triangle, inserts the intersection vertices/constraints that fall on
//     that facet, and reads back the exact retriangulation.
//
//  Compile strict-FP (-ffp-contract=off / /fp:strict).
// ============================================================================
#pragma once

#include "CDT2d.h"
#include "HomogeneousGeometry.h"
#include "Predicates.h"

#include <vector>
#include <array>

namespace sm {
namespace exact {

// ----------------------------------------------------------------------------
//  ExactCDT2d - CDT on exact homogeneous coordinates.
// ----------------------------------------------------------------------------
class ExactCDT2d : public CDTBase2d {
public:
    ExactCDT2d() {
        exact_incircle_ = true;
        exact_intersections_ = true;
    }

    void clear() override {
        CDTBase2d::clear();
        point_.clear(); length_.clear(); cnstr_.clear();
    }

    const Vec2HE& point(index_t v) const { return point_[v]; }

    // The 3 macro-triangle corners (must be CCW so orient_012_ == POSITIVE).
    void create_enclosing_triangle(const Vec2HE& p0, const Vec2HE& p1, const Vec2HE& p2) {
        add_point(p0); add_point(p1); add_point(p2);
        CDTBase2d::create_enclosing_triangle(0, 1, 2);
    }

    index_t insert_point(const Vec2HE& p, index_t hint = CDT_NO_INDEX) {
        index_t v = nv();
        if (v >= point_.size()) add_point(p);
        index_t w = CDTBase2d::insert(v, hint);
        if (w != v) { point_.pop_back(); length_.pop_back(); }  // duplicate
        return w;
    }

    // Constrain edge (v1,v2); remember the *original* endpoints so intersection
    // points keep simple coordinates.
    void insert_constraint(index_t v1, index_t v2) {
        cnstr_.push_back({v1, v2});
        CDTBase2d::insert_constraint(v1, v2);
    }

    using CDTBase2d::nT;
    using CDTBase2d::nv;
    index_t tri_v(index_t t, index_t lv) const { return Tv(t, lv); }
    index_t tri_adj(index_t t, index_t le) const { return Tadj(t, le); }
    bool tri_edge_constrained(index_t t, index_t le) const { return Tedge_is_constrained(t, le); }
    Sign global_orientation() const { return orient_012_; }

protected:
    Sign orient2d(index_t i, index_t j, index_t k) const override {
        return Sign(exact::orient2d(point_[i], point_[j], point_[k]));
    }
    Sign incircle(index_t i, index_t j, index_t k, index_t l) const override {
        return Sign(exact::incircle_2d(point_[i], point_[j], point_[k], point_[l],
                                       length_[i], length_[j], length_[k], length_[l]));
    }
    // Exact intersection of constraints E1 and E2, no division.
    index_t create_intersection(index_t E1, index_t /*i*/, index_t /*j*/,
                                index_t E2, index_t /*k*/, index_t /*l*/) override {
        index_t i = cnstr_[E1][0], j = cnstr_[E1][1];
        index_t k = cnstr_[E2][0], l = cnstr_[E2][1];
        Vec2HE U = point_[j] - point_[i];
        Vec2HE V = point_[l] - point_[k];
        Vec2HE D = point_[k] - point_[i];
        ExpansionNt tnum = det2x2(D.x, D.y, V.x, V.y) * U.w;
        ExpansionNt tden = det2x2(U.x, U.y, V.x, V.y) * D.w;
        Vec2HE P = mix(tnum, tden, point_[i], point_[j]);
        index_t v = nv_;
        add_point(P);
        v2T_.push_back(CDT_NO_INDEX);
        ++nv_;
        return v;
    }

private:
    void add_point(const Vec2HE& p) {
        point_.push_back(p);
        length_.push_back(vec2h_length(p));
    }
    std::vector<Vec2HE> point_;
    std::vector<double>  length_;
    std::vector<std::array<index_t,2>> cnstr_;
};

// ----------------------------------------------------------------------------
//  MeshInTriangle - remesh one 3D facet exactly.
// ----------------------------------------------------------------------------
class MeshInTriangle {
public:
    using index_t = CDTBase2d::index_t;

    // Set the facet's three 3D corners; projects to exact 2D and seeds the CDT.
    void begin_facet(const double* P0, const double* P1, const double* P2) {
        nax_ = triangle_normal_axis(P0, P1, P2);
        u_ = (nax_ + 1) % 3;
        v_ = (nax_ + 2) % 3;
        Vec2HE p0 = project(P0), p1 = project(P1), p2 = project(P2);
        // Macro-triangle must be CCW for orient_012_ == POSITIVE; if the chosen
        // projection flipped it, swap the (u,v) axes.
        if (exact::orient2d(p0, p1, p2) < 0) {
            std::swap(u_, v_);
            p0 = project(P0); p1 = project(P1); p2 = project(P2);
        }
        cdt_.create_enclosing_triangle(p0, p1, p2);
    }

    // Add an intersection vertex given by its exact 3D coordinates (must lie on
    // the facet). Returns its CDT vertex index.
    index_t add_vertex(const double* P) { return cdt_.insert_point(project(P)); }

    // Constrain the edge between two existing CDT vertices.
    void add_constraint(index_t a, index_t b) { cdt_.insert_constraint(a, b); }

    const ExactCDT2d& cdt() const { return cdt_; }
    index_t corner(int i) const { return index_t(i); }  // corners are verts 0,1,2

private:
    Vec2HE project(const double* P) const { return Vec2HE(P[u_], P[v_], 1.0); }
    ExactCDT2d cdt_;
    int nax_ = 2, u_ = 0, v_ = 1;
};

} // namespace exact
} // namespace sm
