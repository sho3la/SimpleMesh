// ============================================================================
//  simplemesh/exact/TriangleIntersection.h - symbolic triangle/triangle isect
// ----------------------------------------------------------------------------
//  Symbolic triangle/triangle intersection. Given two
//  triangles T1=(p0,p1,p2), T2=(q0,q1,q2) it returns the SYMBOLIC description of
//  their intersection: a set of pairs (R1,R2) of TriangleRegion, where each pair
//  names one intersection vertex by the region of T1 and the region of T2 it
//  belongs to (a vertex P*, an edge E*, or the interior T). Working symbolically
//  (rather than with coordinates) is what keeps the later intersection-point
//  coordinates low-degree and exact.
//
//  Every geometric decision is an EXACT predicate: orient3d, orient2d
//  (projected), dot3d, plus aligned_3d / triangle_normal_axis.
//
//  Compile strict-FP (-ffp-contract=off / /fp:strict).
// ============================================================================
#pragma once

#include "Predicates.h"

#include <vector>
#include <utility>
#include <algorithm>
#include <cstdint>

namespace sm {
namespace exact {

// TriangleRegion: the 14 regions of the (T1,T2) pair.
enum TriangleRegion {
    T1_RGN_P0 = 0, T1_RGN_P1 = 1, T1_RGN_P2 = 2,
    T2_RGN_P0 = 3, T2_RGN_P1 = 4, T2_RGN_P2 = 5,
    T1_RGN_E0 = 6, T1_RGN_E1 = 7, T1_RGN_E2 = 8,
    T2_RGN_E0 = 9, T2_RGN_E1 = 10, T2_RGN_E2 = 11,
    T1_RGN_T  = 12, T2_RGN_T  = 13,
    T_RGN_NB  = 14
};

using TriangleIsect = std::pair<TriangleRegion, TriangleRegion>;

inline bool is_in_T1(TriangleRegion R) {
    return R==T1_RGN_P0 || R==T1_RGN_P1 || R==T1_RGN_P2 ||
           R==T1_RGN_E0 || R==T1_RGN_E1 || R==T1_RGN_E2 || R==T1_RGN_T;
}

inline int region_dim(TriangleRegion r) {
    static const int d[T_RGN_NB] = {0,0,0,0,0,0,1,1,1,1,1,1,2,2};
    return d[int(r)];
}

inline void get_triangle_vertices(TriangleRegion T, TriangleRegion& p0,
                                  TriangleRegion& p1, TriangleRegion& p2) {
    if (T == T1_RGN_T) { p0=T1_RGN_P0; p1=T1_RGN_P1; p2=T1_RGN_P2; }
    else               { p0=T2_RGN_P0; p1=T2_RGN_P1; p2=T2_RGN_P2; }
}
inline void get_triangle_edges(TriangleRegion T, TriangleRegion& e0,
                               TriangleRegion& e1, TriangleRegion& e2) {
    if (T == T1_RGN_T) { e0=T1_RGN_E0; e1=T1_RGN_E1; e2=T1_RGN_E2; }
    else               { e0=T2_RGN_E0; e1=T2_RGN_E1; e2=T2_RGN_E2; }
}
inline void get_edge_vertices(TriangleRegion E, TriangleRegion& q0, TriangleRegion& q1) {
    switch (E) {
    case T1_RGN_E0: q0=T1_RGN_P1; q1=T1_RGN_P2; break;
    case T1_RGN_E1: q0=T1_RGN_P2; q1=T1_RGN_P0; break;
    case T1_RGN_E2: q0=T1_RGN_P0; q1=T1_RGN_P1; break;
    case T2_RGN_E0: q0=T2_RGN_P1; q1=T2_RGN_P2; break;
    case T2_RGN_E1: q0=T2_RGN_P2; q1=T2_RGN_P0; break;
    case T2_RGN_E2: q0=T2_RGN_P0; q1=T2_RGN_P1; break;
    default: q0=q1=T_RGN_NB; break;
    }
}

inline TriangleRegion swap_T1_T2(TriangleRegion R) {
    if (R == T_RGN_NB) return T_RGN_NB;
    if (is_in_T1(R)) return TriangleRegion(int(R) + ((region_dim(R)==2) ? 1 : 3));
    return TriangleRegion(int(R) - ((region_dim(R)==2) ? 1 : 3));
}

// ----------------------------------------------------------------------------
//  TriangleTriangleIntersection - the symbolic intersection worker.
// ----------------------------------------------------------------------------
class TriangleTriangleIntersection {
public:
    static constexpr int CACHE_UNINITIALIZED = -2;

    TriangleTriangleIntersection(const double* p0, const double* p1, const double* p2,
                                 const double* q0, const double* q1, const double* q2,
                                 std::vector<TriangleIsect>* result)
        : result_(result) {
        set(0,p0); set(1,p1); set(2,p2); set(3,q0); set(4,q1); set(5,q2);
        for (int i = 0; i < 64; ++i) o3d_cache_[i] = CACHE_UNINITIALIZED;
        has_non_degenerate_intersection_ = false;
    }

    void compute() {
        if (result_) result_->clear();
        // degenerate triangle test
        if (triangle_dim(T1_RGN_P0,T1_RGN_P1,T1_RGN_P2) != 2 ||
            triangle_dim(T2_RGN_P0,T2_RGN_P1,T2_RGN_P2) != 2) return;

        // If T1 strictly on one side of T2's plane -> no intersection.
        {
            TriangleRegion p1,p2,p3,q1,q2,q3;
            get_triangle_vertices(T1_RGN_T, p1,p2,p3);
            get_triangle_vertices(T2_RGN_T, q1,q2,q3);
            int o1 = orient3d(q1,q2,q3,p1);
            int o2 = orient3d(q1,q2,q3,p2);
            int o3 = orient3d(q1,q2,q3,p3);
            if (o1*o2 == 1 && o2*o3 == 1 && o3*o1 == 1) return;
        }

        intersect_edge_triangle(T1_RGN_E0, T2_RGN_T); if (finished()) return;
        intersect_edge_triangle(T1_RGN_E1, T2_RGN_T); if (finished()) return;
        intersect_edge_triangle(T1_RGN_E2, T2_RGN_T); if (finished()) return;
        intersect_edge_triangle(T2_RGN_E0, T1_RGN_T); if (finished()) return;
        intersect_edge_triangle(T2_RGN_E1, T1_RGN_T); if (finished()) return;
        intersect_edge_triangle(T2_RGN_E2, T1_RGN_T); if (finished()) return;

        if (result_) {
            std::sort(result_->begin(), result_->end());
            result_->erase(std::unique(result_->begin(), result_->end()), result_->end());
        }
    }

    bool has_non_degenerate_intersection() const { return has_non_degenerate_intersection_; }

protected:
    bool finished() const { return result_ == nullptr && has_non_degenerate_intersection_; }

    void intersect_edge_triangle(TriangleRegion E, TriangleRegion T) {
        TriangleRegion R1 = E, R2 = T;
        TriangleRegion p1,p2,p3, e1,e2,e3, q1,q2;
        get_triangle_vertices(T,p1,p2,p3);
        get_triangle_edges(T,e1,e2,e3);
        get_edge_vertices(E,q1,q2);

        int o1 = orient3d(p1,p2,p3,q1);
        int o2 = orient3d(p1,p2,p3,q2);
        if (o1*o2 == 1) return;  // both extremities same side -> no isect

        if (o1 == 0 && o2 == 0) {
            // coplanar segment / triangle
            int nax = normal_axis(p1,p2,p3);
            int a1=orient2d(q1,p1,p2,nax), a2=orient2d(q1,p2,p3,nax), a3=orient2d(q1,p3,p1,nax);
            int b1=orient2d(q2,p1,p2,nax), b2=orient2d(q2,p2,p3,nax), b3=orient2d(q2,p3,p1,nax);
            if (a1*a2 > 0 && a2*a3 > 0 && a3*a1 > 0) { add_intersection(q1,T); if(finished())return; }
            if (b1*b2 > 0 && b2*b3 > 0 && b3*b1 > 0) { add_intersection(q2,T); if(finished())return; }
            intersect_edge_edge_2d(E,e1,nax); if (finished()) return;
            intersect_edge_edge_2d(E,e2,nax); if (finished()) return;
            intersect_edge_edge_2d(E,e3,nax); if (finished()) return;
        } else {
            if (o1 == 0) R1 = q1; else if (o2 == 0) R1 = q2;
            int oo1 = orient3d(p2,p3,q1,q2);
            int oo2 = orient3d(p3,p1,q1,q2);
            if (oo1*oo2 == -1) return;
            int oo3 = orient3d(p1,p2,q1,q2);
            int nb_zeros = (oo1==0) + (oo2==0) + (oo3==0);
            if (nb_zeros == 1) {
                R2 = (oo1==0) ? e1 : (oo2==0) ? e2 : e3;
            } else if (nb_zeros == 2) {
                R2 = (oo1!=0) ? p1 : (oo2!=0) ? p2 : p3;
            }
            bool outside = (oo1*oo2 == -1) || (oo2*oo3 == -1) || (oo3*oo1 == -1);
            if (!outside) add_intersection(R1,R2);
        }
    }

    void intersect_edge_edge_2d(TriangleRegion E1, TriangleRegion E2, int nax) {
        TriangleRegion R1 = E1, R2 = E2;
        TriangleRegion p1,p2,q1,q2;
        get_edge_vertices(E1,p1,p2);
        get_edge_vertices(E2,q1,q2);
        int a1 = orient2d(q1,q2,p1,nax), a2 = orient2d(q1,q2,p2,nax);
        if (a1 == 0 && a2 == 0) {
            intersect_edge_edge_1d(E1,E2);
        } else {
            if (a1 == 0) R1 = p1; else if (a2 == 0) R1 = p2;
            int b1 = orient2d(p1,p2,q1,nax), b2 = orient2d(p1,p2,q2,nax);
            if (b1 == 0) R2 = q1; else if (b2 == 0) R2 = q2;
            if (a1*a2 != 1 && b1*b2 != 1) add_intersection(R1,R2);
        }
    }

    void intersect_edge_edge_1d(TriangleRegion E1, TriangleRegion E2) {
        TriangleRegion p1,p2,q1,q2;
        get_edge_vertices(E1,p1,p2);
        get_edge_vertices(E2,q1,q2);
        int d1 = dot3d(p1,q1,q2), d2 = dot3d(p2,q1,q2);
        int d3 = dot3d(q1,p1,p2), d4 = dot3d(q2,p1,p2);
        if (d1==0 && d3==0 && points_are_identical(p1,q1)) add_intersection(p1,q1);
        if (d2==0 && d3==0 && points_are_identical(p2,q1)) add_intersection(p2,q1);
        if (d1==0 && d4==0 && points_are_identical(p1,q2)) add_intersection(p1,q2);
        if (d2==0 && d4==0 && points_are_identical(p2,q2)) add_intersection(p2,q2);
        if (d1 < 0) add_intersection(p1,E2);
        if (d2 < 0) add_intersection(p2,E2);
        if (d3 < 0) add_intersection(E1,q1);
        if (d4 < 0) add_intersection(E1,q2);
    }

    void add_intersection(TriangleRegion R1, TriangleRegion R2) {
        if (region_dim(R1) >= 1 || region_dim(R2) >= 1)
            has_non_degenerate_intersection_ = true;
        if (is_in_T1(R1)) {
            if (result_) result_->push_back(std::make_pair(R1,R2));
        } else {
            if (result_) result_->push_back(std::make_pair(R2,R1));
        }
    }

    // --- exact geometry on the six stored points (regions index p_[]) --------
    int orient3d(TriangleRegion i, TriangleRegion j, TriangleRegion k, TriangleRegion l) const {
        index_t idx = (1u<<int(i)) | (1u<<int(j)) | (1u<<int(k)) | (1u<<int(l));
        bool flip = odd_order(int(i),int(j),int(k),int(l));
        if (o3d_cache_[idx] == CACHE_UNINITIALIZED) {
            int o = exact::orient3d(p_[i],p_[j],p_[k],p_[l]);
            o3d_cache_[idx] = std::int8_t(flip ? -o : o);
        }
        return flip ? -int(o3d_cache_[idx]) : int(o3d_cache_[idx]);
    }

    static bool odd_order(int i, int j, int k, int l) {
        int tab[4] = {i,j,k,l}; bool r = false;
        for (int I = 0; I < 3; ++I)
            for (int J = 0; J < 3-I; ++J)
                if (tab[J] > tab[J+1]) { std::swap(tab[J],tab[J+1]); r = !r; }
        return r;
    }

    int orient2d(TriangleRegion i, TriangleRegion j, TriangleRegion k, int nax) const {
        double pi[2], pj[2], pk[2];
        for (int c = 0; c < 2; ++c) {
            int cc = (nax + 1 + c) % 3;
            pi[c] = p_[i][cc]; pj[c] = p_[j][cc]; pk[c] = p_[k][cc];
        }
        return exact::orient2d(pi, pj, pk);
    }

    int dot3d(TriangleRegion i, TriangleRegion j, TriangleRegion k) const {
        return exact::dot_3d(p_[i], p_[j], p_[k]);
    }

    int normal_axis(TriangleRegion v1, TriangleRegion v2, TriangleRegion v3) const {
        return exact::triangle_normal_axis(p_[v1], p_[v2], p_[v3]);
    }

    bool points_are_identical(TriangleRegion i, TriangleRegion j) const {
        return p_[i][0]==p_[j][0] && p_[i][1]==p_[j][1] && p_[i][2]==p_[j][2];
    }

    int triangle_dim(TriangleRegion i, TriangleRegion j, TriangleRegion k) {
        if (!exact::aligned_3d(p_[i], p_[j], p_[k])) return 2;
        if (points_are_identical(i,j) && points_are_identical(j,k)) return 0;
        return 1;
    }

private:
    using index_t = std::size_t;
    void set(int i, const double* p) { p_[i][0]=p[0]; p_[i][1]=p[1]; p_[i][2]=p[2]; }
    double p_[6][3];
    std::vector<TriangleIsect>* result_;
    bool has_non_degenerate_intersection_;
    mutable std::int8_t o3d_cache_[64];
};

// ----------------------------------------------------------------------------
//  Public entry points.
// ----------------------------------------------------------------------------
inline bool triangles_intersections(
    const double* p0, const double* p1, const double* p2,
    const double* q0, const double* q1, const double* q2,
    std::vector<TriangleIsect>& result
) {
    result.clear();
    TriangleTriangleIntersection I(p0,p1,p2,q0,q1,q2,&result);
    I.compute();
    return I.has_non_degenerate_intersection();
}

inline bool triangles_intersections(
    const double* p0, const double* p1, const double* p2,
    const double* q0, const double* q1, const double* q2
) {
    TriangleTriangleIntersection I(p0,p1,p2,q0,q1,q2,nullptr);
    I.compute();
    return I.has_non_degenerate_intersection();
}

} // namespace exact
} // namespace sm
