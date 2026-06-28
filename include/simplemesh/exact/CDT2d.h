// ============================================================================
//  simplemesh/exact/CDT2d.h - constrained Delaunay triangulation (combinatorics)
// ----------------------------------------------------------------------------
//  CDTBase2d is the abstract COMBINATORIAL constrained-Delaunay triangulation
//  engine; the concrete double-precision CDT2d builds on the exact predicates.
//  The geometry is reached only through the pure-virtual hooks orient2d(i,j,k) /
//  incircle(i,j,k,l) / create_intersection(...); that seam is what lets the
//  exact-coordinate per-triangle remesher plug in the exact homogeneous
//  predicates instead.
//
//  Covers: triangle storage + navigation, create_enclosing_triangle/quad,
//  Delaunay point insertion (locate-by-walking, split, edge flips), and
//  constraint insertion (insert_constraint / find_intersected_edges /
//  constrain_edges) with intersection-vertex creation.
// ============================================================================
#pragma once

#include "Predicates.h"

#include <vector>
#include <cstdint>
#include <random>
#include <functional>

namespace sm {
namespace exact {

// Orientation/sign result.
enum Sign { NEGATIVE = -1, ZERO = 0, POSITIVE = 1 };

static constexpr std::size_t CDT_NO_INDEX = std::size_t(-1);

// ----------------------------------------------------------------------------
//  CDTBase2d - abstract combinatorial CDT
// ----------------------------------------------------------------------------
class CDTBase2d {
public:
    using index_t = std::size_t;

    CDTBase2d() : nv_(0), ncnstr_(0), delaunay_(true), exact_incircle_(true),
                  exact_intersections_(false), orient_012_(POSITIVE),
                  rng_(0x9E3779B9u) {}
    virtual ~CDTBase2d() {}

    index_t nT() const { return T_.size() / 3; }
    index_t nv() const { return nv_; }

    virtual void clear() {
        nv_ = 0; ncnstr_ = 0;
        T_.clear(); Tadj_.clear(); v2T_.clear(); Tflags_.clear();
        Tecnstr_first_.clear(); ecnstr_val_.clear(); ecnstr_next_.clear();
        Tnext_.clear(); Tprev_.clear();
    }

protected:
    // ---- pure-virtual geometry hooks (the exact-predicate seam) -------------
    virtual Sign orient2d(index_t i, index_t j, index_t k) const = 0;
    virtual Sign incircle(index_t i, index_t j, index_t k, index_t l) const = 0;
    // Create the vertex where constraint E1=(i,j) crosses constraint E2=(k,l),
    // append it (coords supplied by the subclass), and return its index.
    virtual index_t create_intersection(index_t E1, index_t i, index_t j,
                                        index_t E2, index_t k, index_t l) = 0;

    virtual void begin_insert_transaction() {}
    virtual void commit_insert_transaction() {}
    virtual void rollback_insert_transaction() {}

    // ---- triangle storage accessors -----------------------------------------
    index_t Tv(index_t t, index_t lv) const { return T_[3*t+lv]; }
    index_t Tadj(index_t t, index_t le) const { return Tadj_[3*t+le]; }
    index_t vT(index_t v) const { return v2T_[v]; }
    static index_t find_3(const index_t* p, index_t x) {
        return (p[0]==x) ? 0 : (p[1]==x) ? 1 : (p[2]==x) ? 2 : CDT_NO_INDEX;
    }
    index_t Tv_find(index_t t, index_t v) const { return find_3(&T_[3*t], v); }
    index_t Tadj_find(index_t t1, index_t t2) const { return find_3(&Tadj_[3*t1], t2); }

    index_t Tedge_cnstr_first(index_t t, index_t le) const { return Tecnstr_first_[3*t+le]; }
    void Tset_edge_cnstr_first(index_t t, index_t le, index_t ec) { Tecnstr_first_[3*t+le] = ec; }
    bool Tedge_is_constrained(index_t t, index_t le) const {
        return Tedge_cnstr_first(t, le) != CDT_NO_INDEX;
    }
    index_t ncnstr() const { return ncnstr_; }
    index_t edge_cnstr(index_t ecit) const { return ecnstr_val_[ecit]; }
    index_t edge_cnstr_next(index_t ecit) const { return ecnstr_next_[ecit]; }

    // Append a constraint id to a triangle edge's constraint list.
    void Tadd_edge_cnstr(index_t t, index_t le, index_t cnstr_id) {
        for (index_t ec = Tedge_cnstr_first(t, le); ec != CDT_NO_INDEX;
             ec = edge_cnstr_next(ec))
            if (edge_cnstr(ec) == cnstr_id) return;
        ecnstr_val_.push_back(cnstr_id);
        ecnstr_next_.push_back(Tedge_cnstr_first(t, le));
        Tset_edge_cnstr_first(t, le, ecnstr_val_.size() - 1);
    }
    void Tadd_edge_cnstr_with_neighbor(index_t t, index_t le, index_t cnstr_id) {
        Tadd_edge_cnstr(t, le, cnstr_id);
        index_t t2 = Tadj(t, le);
        if (t2 != CDT_NO_INDEX) {
            index_t le2 = Tadj_find(t2, t);
            Tset_edge_cnstr_first(t2, le2, Tedge_cnstr_first(t, le));
        }
    }

    // Visit every triangle around v, calling doit(t, lv); stops when doit
    // returns true. Handles both interior and border fans.
    void for_each_T_around_v(index_t v, std::function<bool(index_t, index_t)> doit) {
        index_t t = vT(v);
        index_t lv = CDT_NO_INDEX;
        do {
            lv = Tv_find(t, v);
            if (doit(t, lv)) return;
            t = Tadj(t, (lv + 1) % 3);
        } while (t != vT(v) && t != CDT_NO_INDEX);
        if (t != CDT_NO_INDEX) return;  // interior vertex, done
        t = vT(v);
        lv = Tv_find(t, v);
        t = Tadj(t, (lv + 2) % 3);
        while (t != CDT_NO_INDEX) {
            lv = Tv_find(t, v);
            if (doit(t, lv)) return;
            t = Tadj(t, (lv + 2) % 3);
        }
    }

    index_t Tnew() {
        index_t t = nT();
        index_t nc = (t+1)*3;
        T_.resize(nc, CDT_NO_INDEX);
        Tadj_.resize(nc, CDT_NO_INDEX);
        Tecnstr_first_.resize(nc, CDT_NO_INDEX);
        Tflags_.resize(t+1, 0);
        Tnext_.resize(t+1, CDT_NO_INDEX);
        Tprev_.resize(t+1, CDT_NO_INDEX);
        return t;
    }

    void Tset(index_t t, index_t v1, index_t v2, index_t v3,
              index_t adj1, index_t adj2, index_t adj3,
              index_t e1 = CDT_NO_INDEX, index_t e2 = CDT_NO_INDEX,
              index_t e3 = CDT_NO_INDEX) {
        T_[3*t]=v1; T_[3*t+1]=v2; T_[3*t+2]=v3;
        Tadj_[3*t]=adj1; Tadj_[3*t+1]=adj2; Tadj_[3*t+2]=adj3;
        Tecnstr_first_[3*t]=e1; Tecnstr_first_[3*t+1]=e2; Tecnstr_first_[3*t+2]=e3;
        v2T_[v1]=t; v2T_[v2]=t; v2T_[v3]=t;
    }

    void Trot(index_t t, index_t lv) {
        if (lv != 0) {
            index_t i=3*t+lv, j=3*t+((lv+1)%3), k=3*t+((lv+2)%3);
            Tset(t, T_[i],T_[j],T_[k], Tadj_[i],Tadj_[j],Tadj_[k],
                 Tecnstr_first_[i],Tecnstr_first_[j],Tecnstr_first_[k]);
        }
    }

    void Tadj_set(index_t t, index_t le, index_t adj) { Tadj_[3*t+le] = adj; }

    // After re-linking t1->neighbor across le1, re-link neighbor->t1.
    void Tadj_back_connect(index_t t1, index_t le1, index_t prev_t2_adj_e2) {
        index_t t2 = Tadj(t1, le1);
        if (t2 == CDT_NO_INDEX) return;
        index_t le2 = Tadj_find(t2, prev_t2_adj_e2);
        Tadj_set(t2, le2, t1);
        Tset_edge_cnstr_first(t1, le1, Tedge_cnstr_first(t2, le2));
    }

    // ---- triangle flags / DList plumbing ------------------------------------
    enum { DLIST_S_ID = 0, DLIST_Q_ID = 1, DLIST_N_ID = 2, DLIST_NB = 3 };
    void Tset_flag(index_t t, index_t f) { Tflags_[t] |= std::uint8_t(1u << f); }
    void Treset_flag(index_t t, index_t f) { Tflags_[t] &= std::uint8_t(~(1u << f)); }
    bool Tflag_is_set(index_t t, index_t f) const { return (Tflags_[t] & (1u << f)) != 0; }
    bool Tis_in_list(index_t t) const {
        return (Tflags_[t] & std::uint8_t((1 << DLIST_NB) - 1)) != 0;
    }

    // Doubly-connected triangle list backed by Tnext_/Tprev_/Tflags_.
    struct DList {
        CDTBase2d& cdt_; index_t list_id_; index_t back_, front_;
        DList(CDTBase2d& cdt, index_t id) : cdt_(cdt), list_id_(id),
            back_(CDT_NO_INDEX), front_(CDT_NO_INDEX) {}
        explicit DList(CDTBase2d& cdt) : cdt_(cdt), list_id_(CDT_NO_INDEX),
            back_(CDT_NO_INDEX), front_(CDT_NO_INDEX) {}
        ~DList() { if (initialized()) clear(); }
        void initialize(index_t id) { list_id_ = id; }
        bool initialized() const { return list_id_ != CDT_NO_INDEX; }
        bool empty() const { return back_ == CDT_NO_INDEX; }
        bool contains(index_t t) const { return cdt_.Tflag_is_set(t, list_id_); }
        index_t front() const { return front_; }
        index_t back() const { return back_; }
        index_t next(index_t t) const { return cdt_.Tnext_[t]; }
        index_t prev(index_t t) const { return cdt_.Tprev_[t]; }
        void clear() {
            for (index_t t = front_; t != CDT_NO_INDEX; t = cdt_.Tnext_[t])
                cdt_.Treset_flag(t, list_id_);
            back_ = front_ = CDT_NO_INDEX;
        }
        void push_back(index_t t) {
            cdt_.Tset_flag(t, list_id_);
            if (empty()) {
                front_ = back_ = t;
                cdt_.Tnext_[t] = cdt_.Tprev_[t] = CDT_NO_INDEX;
            } else {
                cdt_.Tnext_[back_] = t; cdt_.Tprev_[t] = back_;
                cdt_.Tnext_[t] = CDT_NO_INDEX; back_ = t;
            }
        }
        index_t pop_back() {
            index_t t = back_;
            cdt_.Treset_flag(t, list_id_);
            back_ = cdt_.Tprev_[t];
            if (back_ == CDT_NO_INDEX) front_ = CDT_NO_INDEX;
            else cdt_.Tnext_[back_] = CDT_NO_INDEX;
            return t;
        }
        void push_front(index_t t) {
            cdt_.Tset_flag(t, list_id_);
            if (empty()) {
                back_ = front_ = t;
                cdt_.Tnext_[t] = cdt_.Tprev_[t] = CDT_NO_INDEX;
            } else {
                cdt_.Tprev_[t] = CDT_NO_INDEX;
                cdt_.Tprev_[front_] = t;
                cdt_.Tnext_[t] = front_;
                front_ = t;
            }
        }
    };

    // ---- enclosing simplex --------------------------------------------------
    void create_enclosing_triangle(index_t v0, index_t v1, index_t v2) {
        nv_ = 3; v2T_.resize(3, CDT_NO_INDEX);
        index_t t0 = Tnew();
        Tset(t0, v0, v1, v2, CDT_NO_INDEX, CDT_NO_INDEX, CDT_NO_INDEX);
        orient_012_ = orient2d(0, 1, 2);
    }

    void create_enclosing_quad(index_t v0, index_t v1, index_t v2, index_t v3) {
        nv_ = 4; v2T_.resize(4, CDT_NO_INDEX);
        index_t t0 = Tnew(), t1 = Tnew();
        Tset(t0, v0, v1, v3, t1, CDT_NO_INDEX, CDT_NO_INDEX);
        Tset(t1, v3, v1, v2, CDT_NO_INDEX, CDT_NO_INDEX, t0);
        orient_012_ = orient2d(0, 1, 2);
        if (Sign(incircle(v0, v1, v2, v3) * orient_012_) == POSITIVE) swap_edge(t0);
    }

    // ---- Delaunay point insertion -------------------------------------------
    index_t insert(index_t v, index_t hint = CDT_NO_INDEX) {
        bool keep_duplicates = false;
        if (v == nv()) { v2T_.push_back(CDT_NO_INDEX); ++nv_; }
        else { keep_duplicates = true; }

        begin_insert_transaction();
        Sign o[3];
        index_t t = locate(v, hint, o);
        int nb_z = (o[0]==ZERO) + (o[1]==ZERO) + (o[2]==ZERO);

        if (nb_z == 2) {  // duplicated vertex
            v = (o[0] != ZERO) ? Tv(t,0) : (o[1] != ZERO) ? Tv(t,1) : Tv(t,2);
            if (!keep_duplicates) { v2T_.pop_back(); --nv_; }
            rollback_insert_transaction();
            return v;
        }
        commit_insert_transaction();

        DList S(*this);
        if (delaunay_) S.initialize(DLIST_S_ID);

        if (nb_z == 1) {
            index_t le = (o[0]==ZERO) ? 0 : (o[1]==ZERO) ? 1 : 2;
            insert_vertex_in_edge(v, t, le, S);
        } else {
            insert_vertex_in_triangle(v, t, S);
        }

        if (delaunay_) Delaunayize_vertex_neighbors(v, S);
        return v;
    }

    index_t locate(index_t v, index_t hint, Sign* o) const {
        index_t t_pred = nT() + 1;
        index_t t = (hint == CDT_NO_INDEX) ? (index_t(rng_()) % nT()) : hint;
    still_walking:
        {
            index_t tv[3] = { Tv(t,0), Tv(t,1), Tv(t,2) };
            index_t e0 = index_t(rng_()) % 3;
            for (index_t de = 0; de < 3; ++de) {
                index_t le = (e0 + de) % 3;
                index_t t_next = Tadj(t, le);
                if (t_next == t_pred) { o[le] = POSITIVE; continue; }
                index_t v_bkp = tv[le];
                tv[le] = v;
                o[le] = Sign(orient_012_ * orient2d(tv[0], tv[1], tv[2]));
                if (o[le] != NEGATIVE) { tv[le] = v_bkp; continue; }
                // o[le] == NEGATIVE: we must cross edge le. If there is no
                // triangle there, the query point is OUTSIDE the triangulation
                // (a caller bug - all inserted points must lie within it). Don't
                // walk into nothing (would corrupt memory); abandon this edge.
                if (t_next == CDT_NO_INDEX) { tv[le] = v_bkp; continue; }
                t_pred = t; t = t_next; goto still_walking;
            }
        }
        return t;
    }

    void insert_vertex_in_triangle(index_t v, index_t t, DList& S) {
        index_t t1 = t;
        index_t v1=Tv(t1,0), v2=Tv(t1,1), v3=Tv(t1,2);
        index_t adj1=Tadj(t1,0), adj2=Tadj(t1,1), adj3=Tadj(t1,2);
        index_t t2 = Tnew(), t3 = Tnew();
        Tset(t1, v, v2, v3, adj1, t2, t3);
        Tset(t2, v, v3, v1, adj2, t3, t1);
        Tset(t3, v, v1, v2, adj3, t1, t2);
        Tadj_back_connect(t1, 0, t1);
        Tadj_back_connect(t2, 0, t1);
        Tadj_back_connect(t3, 0, t1);
        if (S.initialized()) { S.push_back(t1); S.push_back(t2); S.push_back(t3); }
    }

    void insert_vertex_in_edge(index_t v, index_t t, index_t le1) {
        DList S(*this);  // uninitialized -> no flip stack
        insert_vertex_in_edge(v, t, le1, S);
    }

    void insert_vertex_in_edge(index_t v, index_t t, index_t le1, DList& S) {
        index_t cnstr_first = Tedge_cnstr_first(t, le1);
        index_t t1 = t;
        index_t t2 = Tadj(t1, le1);
        index_t v1=Tv(t1,le1), v2=Tv(t1,(le1+1)%3), v3=Tv(t1,(le1+2)%3);
        index_t t1_adj2=Tadj(t1,(le1+1)%3), t1_adj3=Tadj(t1,(le1+2)%3);
        if (t2 != CDT_NO_INDEX) {
            index_t le2 = Tadj_find(t2, t1);
            index_t v4 = Tv(t2, le2);
            index_t t2_adj2=Tadj(t2,(le2+1)%3), t2_adj3=Tadj(t2,(le2+2)%3);
            index_t t3 = Tnew(), t4 = Tnew();
            Tset(t1, v, v1, v2, t1_adj3, t2, t4);
            Tset(t2, v, v2, v4, t2_adj2, t3, t1);
            Tset(t3, v, v4, v3, t2_adj3, t4, t2);
            Tset(t4, v, v3, v1, t1_adj2, t1, t3);
            Tadj_back_connect(t1, 0, t1);
            Tadj_back_connect(t2, 0, t2);
            Tadj_back_connect(t3, 0, t2);
            Tadj_back_connect(t4, 0, t1);
            Tset_edge_cnstr_first(t1,1,cnstr_first); Tset_edge_cnstr_first(t2,2,cnstr_first);
            Tset_edge_cnstr_first(t3,1,cnstr_first); Tset_edge_cnstr_first(t4,2,cnstr_first);
            if (S.initialized()) { S.push_back(t1); S.push_back(t2); S.push_back(t3); S.push_back(t4); }
        } else {
            t2 = Tnew();
            Tset(t1, v, v1, v2, t1_adj3, CDT_NO_INDEX, t2);
            Tset(t2, v, v3, v1, t1_adj2, t1, CDT_NO_INDEX);
            Tadj_back_connect(t1, 0, t1);
            Tadj_back_connect(t2, 0, t1);
            Tset_edge_cnstr_first(t1,1,cnstr_first); Tset_edge_cnstr_first(t2,2,cnstr_first);
            if (S.initialized()) { S.push_back(t1); S.push_back(t2); }
        }
    }

    void swap_edge(index_t t1, bool swap_t1_t2 = false) {
        index_t v1=Tv(t1,0), v2=Tv(t1,1), v3=Tv(t1,2);
        index_t t1_adj2=Tadj(t1,1), t1_adj3=Tadj(t1,2);
        index_t t2 = Tadj(t1,0);
        index_t le2 = Tadj_find(t2, t1);
        index_t v4 = Tv(t2, le2);
        index_t t2_adj2=Tadj(t2,(le2+1)%3), t2_adj3=Tadj(t2,(le2+2)%3);
        if (swap_t1_t2) {
            Tset(t2, v1, v4, v3, t2_adj3, t1_adj2, t1);
            Tset(t1, v1, v2, v4, t2_adj2, t2, t1_adj3);
            Tadj_back_connect(t2, 0, t2);
            Tadj_back_connect(t2, 1, t1);
            Tadj_back_connect(t1, 0, t2);
            Tadj_back_connect(t1, 2, t1);
        } else {
            Tset(t1, v1, v4, v3, t2_adj3, t1_adj2, t2);
            Tset(t2, v1, v2, v4, t2_adj2, t1, t1_adj3);
            Tadj_back_connect(t1, 0, t2);
            Tadj_back_connect(t1, 1, t1);
            Tadj_back_connect(t2, 0, t2);
            Tadj_back_connect(t2, 2, t1);
        }
    }

    bool is_convex_quad(index_t t) const {
        index_t v1=Tv(t,0), v2=Tv(t,1), v3=Tv(t,2);
        index_t t2 = Tadj(t,0);
        index_t le2 = Tadj_find(t2, t);
        index_t v4 = Tv(t2, le2);
        return orient2d(v1,v4,v3) == orient_012_ &&
               orient2d(v4,v1,v2) == orient_012_;
    }

    void Delaunayize_vertex_neighbors(index_t v, DList& S) {
        while (!S.empty()) {
            index_t t1 = S.pop_back();
            if (Tedge_is_constrained(t1, 0)) continue;
            index_t t2 = Tadj(t1, 0);
            if (t2 == CDT_NO_INDEX) continue;
            if (!exact_incircle_ && !is_convex_quad(t1)) continue;
            index_t v1=Tv(t2,0), v2=Tv(t2,1), v3=Tv(t2,2);
            if (Sign(incircle(v1, v2, v3, v) * orient_012_) == POSITIVE) {
                swap_edge(t1);
                S.push_back(t1);
                S.push_back(t2);
            }
        }
    }

    // ---- constraint insertion -----------------------------------------------
    struct ConstraintWalker {
        ConstraintWalker(index_t i_in, index_t j_in)
            : i(i_in), j(j_in), t_prev(CDT_NO_INDEX), v_prev(CDT_NO_INDEX),
              t(CDT_NO_INDEX), v(i_in) {}
        index_t i, j, t_prev, v_prev, t, v;
    };

public:
    // Constrain the triangulation to contain edge (i,j), creating intersection
    // vertices where it crosses existing constraints.
    void insert_constraint(index_t i, index_t j) {
        ++ncnstr_;
        index_t first_v_isect = nv_;
        DList Q(*this, DLIST_Q_ID);
        DList N(*this);
        if (delaunay_) N.initialize(DLIST_N_ID);
        while (i != j) {
            index_t k = find_intersected_edges(i, j, Q);
            if (delaunay_ && exact_intersections_ && k >= first_v_isect) {
                insert(k);
                Q.clear();
                Delaunayize_vertex_neighbors(k);
                index_t new_k = find_intersected_edges(i, j, Q);
                (void)new_k;
            }
            constrain_edges(i, k, Q, N);
            if (delaunay_) Delaunayize_new_edges(N);
            i = k;
        }
        if (delaunay_ && !exact_intersections_) {
            for (index_t v = first_v_isect; v < nv(); ++v)
                Delaunayize_vertex_neighbors(v);
        }
    }

protected:
    index_t find_intersected_edges(index_t i, index_t j, DList& Q) {
        ConstraintWalker W(i, j);
        while (W.v == i || W.v == CDT_NO_INDEX) {
            if (W.v != CDT_NO_INDEX) walk_constraint_v(W);
            else                     walk_constraint_t(W, Q);
        }
        return W.v;
    }

    void walk_constraint_v(ConstraintWalker& W) {
        index_t t_next = CDT_NO_INDEX, v_next = CDT_NO_INDEX;
        for_each_T_around_v(W.v, [&](index_t t_around_v, index_t le) -> bool {
            if (t_around_v == W.t_prev) return false;
            index_t v1 = Tv(t_around_v, (le+1)%3);
            index_t v2 = Tv(t_around_v, (le+2)%3);
            if (v1 == W.j || v2 == W.j) {
                v_next = W.j;
                index_t le_c = (v1 == W.j) ? (le+2)%3 : (le+1)%3;
                Tadd_edge_cnstr_with_neighbor(t_around_v, le_c, ncnstr_-1);
                return true;
            }
            Sign o1 = orient2d(W.i, W.j, v1);
            Sign o2 = orient2d(W.i, W.j, v2);
            Sign o3 = orient2d(v1, v2, W.j);
            Sign o4 = orient_012_;
            if (o1*o2 < 0 && o3*o4 < 0) {
                Trot(t_around_v, le);
                t_next = t_around_v;
                return true;
            } else {
                if (o1 == ZERO && o3*o4 < 0 && v1 != W.v_prev) {
                    v_next = v1;
                    Tadd_edge_cnstr_with_neighbor(t_around_v, (le+2)%3, ncnstr_-1);
                    return true;
                } else if (o2 == ZERO && o3*o4 < 0 && v2 != W.v_prev) {
                    v_next = v2;
                    Tadd_edge_cnstr_with_neighbor(t_around_v, (le+1)%3, ncnstr_-1);
                    return true;
                }
            }
            return false;
        });
        W.t_prev = W.t; W.v_prev = W.v; W.t = t_next; W.v = v_next;
    }

    void walk_constraint_t(ConstraintWalker& W, DList& Q) {
        index_t v_next = CDT_NO_INDEX, t_next = CDT_NO_INDEX;
        if (Tv(W.t,0)==W.j || Tv(W.t,1)==W.j || Tv(W.t,2)==W.j) {
            v_next = W.j;
        } else {
            for (index_t le = 0; le < 3; ++le) {
                if (Tadj(W.t, le) == W.t_prev) continue;
                index_t v1 = Tv(W.t, (le+1)%3);
                index_t v2 = Tv(W.t, (le+2)%3);
                Sign o1 = orient2d(W.i, W.j, v1);
                Sign o2 = orient2d(W.i, W.j, v2);
                if (o1*o2 < 0) {
                    Trot(W.t, le);
                    if (Tedge_is_constrained(W.t, 0)) {
                        v_next = create_intersection(
                            ncnstr()-1, W.i, W.j,
                            edge_cnstr(Tedge_cnstr_first(W.t,0)), v1, v2);
                        insert_vertex_in_edge(v_next, W.t, 0);
                        if (W.v_prev != CDT_NO_INDEX)
                            Tadd_edge_cnstr_with_neighbor(W.t, 2, ncnstr_-1);
                    } else {
                        Q.push_back(W.t);
                        t_next = Tadj(W.t, 0);
                    }
                    break;
                } else {
                    if (o1 == ZERO) { v_next = v1; break; }
                    else if (o2 == ZERO) { v_next = v2; break; }
                }
            }
        }
        W.t_prev = W.t; W.v_prev = W.v; W.t = t_next; W.v = v_next;
    }

    void constrain_edges(index_t i, index_t j, DList& Q, DList& N) {
        auto new_edge = [&](index_t t, index_t le) {
            Trot(t, le);
            if ((Tv(t,1)==i && Tv(t,2)==j) || (Tv(t,1)==j && Tv(t,2)==i))
                Tadd_edge_cnstr_with_neighbor(t, 0, ncnstr_-1);
            else if (N.initialized())
                N.push_back(t);
        };
        auto isect_edge = [&](index_t t, index_t le) { Trot(t, le); Q.push_front(t); };

        while (!Q.empty()) {
            index_t t1 = Q.pop_back();
            if (!is_convex_quad(t1)) {
                Q.push_front(t1);
            } else {
                index_t t2 = Tadj(t1, 0);
                bool no_isect = !Q.contains(t2);
                index_t v0 = Tv(t1, 0);
                bool t2v0_t1v2 = (Q.contains(t2) && Tv(t2,0) == Tv(t1,2));
                if (no_isect) {
                    swap_edge(t1);
                    new_edge(t1, 2);
                } else {
                    Sign o = Sign(orient2d(i, j, v0) * orient_012_);
                    if (t2v0_t1v2) {
                        swap_edge(t1, false);
                        if (o >= 0) new_edge(t1, 2);
                        else        isect_edge(t1, 2);
                    } else {
                        swap_edge(t1, true);
                        if (o > 0)  isect_edge(t1, 1);
                        else        new_edge(t1, 1);
                    }
                }
            }
        }
    }

    void Delaunayize_new_edges(DList& N) {
        bool swap_occured = true;
        while (swap_occured) {
            swap_occured = false;
            for (index_t t1 = N.front(); t1 != CDT_NO_INDEX; t1 = N.next(t1)) {
                if (Tedge_is_constrained(t1, 0)) continue;
                index_t v1 = Tv(t1,1), v2 = Tv(t1,2), v0 = Tv(t1,0);
                index_t t2 = Tadj(t1, 0);
                if (t2 == CDT_NO_INDEX) continue;
                if (!exact_incircle_ && !is_convex_quad(t1)) continue;
                index_t e2 = Tadj_find(t2, t1);
                index_t v3 = Tv(t2, e2);
                if (Sign(incircle(v0, v1, v2, v3) * orient_012_) == POSITIVE) {
                    if (Tv(t2,0) == Tv(t1,1)) { swap_edge(t1, true);  Trot(t1, 1); }
                    else                      { swap_edge(t1, false); Trot(t1, 2); }
                    swap_occured = true;
                }
            }
        }
        N.clear();
    }

    // Re-Delaunayize the fan around an interior vertex v.
    void Delaunayize_vertex_neighbors(index_t v) {
        DList S(*this, DLIST_S_ID);
        index_t t0 = vT(v), t = t0;
        do {
            index_t lv = Tv_find(t, v);
            Trot(t, lv);
            S.push_back(t);
            t = Tadj(t, 1);
        } while (t != t0 && t != CDT_NO_INDEX);
        Delaunayize_vertex_neighbors(v, S);
    }

    // ---- storage ------------------------------------------------------------
    index_t nv_, ncnstr_;
    bool delaunay_, exact_incircle_, exact_intersections_;
    Sign orient_012_;
    std::vector<index_t> T_, Tadj_, v2T_;
    std::vector<std::uint8_t> Tflags_;
    std::vector<index_t> Tecnstr_first_, ecnstr_val_, ecnstr_next_;
    std::vector<index_t> Tnext_, Tprev_;
    mutable std::mt19937 rng_;  // replaces Numeric::random_int32 (deterministic)
};

// ----------------------------------------------------------------------------
//  CDT2d - concrete double-precision CDT using the Day-3 exact predicates.
// ----------------------------------------------------------------------------
class CDT2d : public CDTBase2d {
public:
    // Establish the big enclosing triangle (must contain every point inserted
    // later, and be given CCW so orient_012_ == POSITIVE).
    void create_enclosing_triangle(double x0, double y0, double x1, double y1,
                                   double x2, double y2) {
        px_ = {x0, x1, x2}; py_ = {y0, y1, y2};
        CDTBase2d::create_enclosing_triangle(0, 1, 2);
    }

    // Insert a point; returns its vertex index (an existing index if duplicate).
    index_t insert_point(double x, double y, index_t hint = CDT_NO_INDEX) {
        index_t v = nv();
        if (v >= px_.size()) { px_.resize(v + 1); py_.resize(v + 1); }
        px_[v] = x; py_[v] = y;
        index_t w = CDTBase2d::insert(v, hint);
        if (w != v) { px_.resize(v); py_.resize(v); }  // duplicate: drop coords
        return w;
    }

    double px(index_t v) const { return px_[v]; }
    double py(index_t v) const { return py_[v]; }
    // Expose read-only combinatorics for tests / consumers.
    using CDTBase2d::nT;
    using CDTBase2d::nv;
    index_t tri_v(index_t t, index_t lv) const { return Tv(t, lv); }
    index_t tri_adj(index_t t, index_t le) const { return Tadj(t, le); }
    bool tri_edge_constrained(index_t t, index_t le) const { return Tedge_is_constrained(t, le); }
    Sign global_orientation() const { return orient_012_; }

    void clear() override { CDTBase2d::clear(); px_.clear(); py_.clear(); }

protected:
    Sign orient2d(index_t i, index_t j, index_t k) const override {
        return Sign(exact::orient2d(px_[i], py_[i], px_[j], py_[j], px_[k], py_[k]));
    }
    Sign incircle(index_t i, index_t j, index_t k, index_t l) const override {
        return Sign(exact::in_circle(px_[i], py_[i], px_[j], py_[j],
                                     px_[k], py_[k], px_[l], py_[l]));
    }
    // Intersection of constraint (i,j) with (k,l), computed in double precision
    // Not bit-exact - the EXACT version
    // lives in the Day-7 homogeneous subclass; this is the plain double CDT.
    index_t create_intersection(index_t /*E1*/, index_t i, index_t j,
                                index_t /*E2*/, index_t k, index_t l) override {
        double Ux = px_[j]-px_[i], Uy = py_[j]-py_[i];
        double Vx = px_[l]-px_[k], Vy = py_[l]-py_[k];
        double Dx = px_[k]-px_[i], Dy = py_[k]-py_[i];
        double delta = Ux*Vy - Uy*Vx;
        double t = (Dx*Vy - Dy*Vx) / delta;
        index_t v = nv_;
        px_.push_back(px_[i] + t*Ux);
        py_.push_back(py_[i] + t*Uy);
        v2T_.push_back(CDT_NO_INDEX);
        ++nv_;
        return v;
    }

private:
    std::vector<double> px_, py_;
};

} // namespace exact
} // namespace sm
