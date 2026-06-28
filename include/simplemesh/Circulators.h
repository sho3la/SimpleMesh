// ============================================================================
//  SimpleMesh - Circulators.h : lazy iterators & ranges
// ----------------------------------------------------------------------------
//  The vector-returning circulators (e.g. vertex_vertices) are simple and
//  Pythonic but allocate a heap buffer on every call. This
//  header adds *lazy* iterators that walk the connectivity on
//  demand with zero allocation, and wraps them in range objects so they work in
//  C++11 range-based for loops:
//
//      for (auto h : mesh.voh_range(v))   // outgoing halfedges of v
//          ...
//      for (auto v : mesh.all_vertices()) // every (non-deleted) vertex
//          ...
//
//  This file is included at the BOTTOM of Mesh.h, so the Mesh class is already
//  complete here and the iterators can call its public navigation methods. It
//  also #includes Mesh.h so it can be used stand-alone; the include guards make
//  the mutual inclusion safe.
// ============================================================================
#pragma once

#include "Mesh.h"

namespace sm {

// ============================================================================
//  1. Cyclic circulators (around a vertex or around a face)
// ----------------------------------------------------------------------------
//  A circulator walks a *cycle* of halfedges, so begin() and end() refer to the
//  same start halfedge - we can't tell them apart by handle alone. The classic
//  fix is an `active_` flag: begin() starts active, and once
//  ++ brings us back to the start we go inactive, becoming equal to end().
//
//  The cyclic logic lives once in this CRTP base; each concrete circulator only
//  says how to ROTATE to the next halfedge and how to turn the current halfedge
//  into the VALUE it yields.
// ============================================================================
template <class Derived, class Value>
class CircBase {
public:
    CircBase(const Mesh* m, HalfedgeHandle start, bool active)
        : mesh_(m), start_(start), cur_(start), active_(active) {}

    bool operator!=(const CircBase& o) const { return active_ != o.active_ || cur_ != o.cur_; }
    bool operator==(const CircBase& o) const { return !(*this != o); }

    Derived& operator++() {
        cur_ = self().rotate(cur_);
        if (cur_ == start_) active_ = false;   // completed the loop
        return self();
    }
    Value operator*() const { return self().value(cur_); }

    // current raw halfedge (handy for debugging / mixed traversal)
    HalfedgeHandle halfedge() const { return cur_; }

protected:
    Derived&       self()       { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }

    const Mesh*    mesh_;
    HalfedgeHandle start_, cur_;
    bool           active_;
};

/// Outgoing halfedges around a vertex.
class VOHCirculator : public CircBase<VOHCirculator, HalfedgeHandle> {
public:
    using CircBase::CircBase;
    // rotate: opposite(h) ends at our vertex, so next(opposite(h)) starts there.
    HalfedgeHandle rotate(HalfedgeHandle h) const { return mesh_->next_halfedge(mesh_->opposite_halfedge(h)); }
    HalfedgeHandle value(HalfedgeHandle h) const { return h; }
};

/// Neighbour vertices (1-ring) of a vertex.
class VVCirculator : public CircBase<VVCirculator, VertexHandle> {
public:
    using CircBase::CircBase;
    HalfedgeHandle rotate(HalfedgeHandle h) const { return mesh_->next_halfedge(mesh_->opposite_halfedge(h)); }
    VertexHandle   value(HalfedgeHandle h) const { return mesh_->to_vertex(h); }
};

/// Halfedges around a face.
class FHCirculator : public CircBase<FHCirculator, HalfedgeHandle> {
public:
    using CircBase::CircBase;
    HalfedgeHandle rotate(HalfedgeHandle h) const { return mesh_->next_halfedge(h); }
    HalfedgeHandle value(HalfedgeHandle h) const { return h; }
};

/// Vertices around a face (in order).
class FVCirculator : public CircBase<FVCirculator, VertexHandle> {
public:
    using CircBase::CircBase;
    HalfedgeHandle rotate(HalfedgeHandle h) const { return mesh_->next_halfedge(h); }
    VertexHandle   value(HalfedgeHandle h) const { return mesh_->to_vertex(h); }
};

/// A range adapter so a circulator works in `for (x : range)`.
template <class Circ>
class CircRange {
public:
    CircRange(const Mesh* m, HalfedgeHandle start) : mesh_(m), start_(start) {}
    Circ begin() const { return Circ(mesh_, start_, start_.is_valid()); }
    Circ end()   const { return Circ(mesh_, start_, false); }
private:
    const Mesh*    mesh_;
    HalfedgeHandle start_;
};

// ============================================================================
//  2. Linear element iterators (every vertex / edge / halfedge / face)
// ----------------------------------------------------------------------------
//  These walk the element arrays by index and SKIP elements flagged deleted,
//  so you can iterate a mesh that hasn't been garbage-collected
//  yet and only see live elements. `is_deleted` is overloaded per handle type,
//  so the right check is picked automatically from the Handle template arg.
// ============================================================================
template <class Handle>
class ElementIter {
public:
    ElementIter(const Mesh* m, int idx, int n) : mesh_(m), idx_(idx), n_(n) { skip_deleted(); }

    bool operator!=(const ElementIter& o) const { return idx_ != o.idx_; }
    bool operator==(const ElementIter& o) const { return idx_ == o.idx_; }  // py::make_iterator needs ==
    ElementIter& operator++() { ++idx_; skip_deleted(); return *this; }
    Handle operator*() const { return Handle(idx_); }

private:
    void skip_deleted() { while (idx_ < n_ && mesh_->is_deleted(Handle(idx_))) ++idx_; }
    const Mesh* mesh_;
    int idx_, n_;
};

template <class Handle>
class ElementRange {
public:
    ElementRange(const Mesh* m, int n) : mesh_(m), n_(n) {}
    ElementIter<Handle> begin() const { return ElementIter<Handle>(mesh_, 0,  n_); }
    ElementIter<Handle> end()   const { return ElementIter<Handle>(mesh_, n_, n_); }
private:
    const Mesh* mesh_;
    int n_;
};

// ============================================================================
//  3. The Mesh member functions that hand out these ranges (defined inline
//     here, declared in Mesh.h). `this` is a complete Mesh at this point.
// ============================================================================
inline CircRange<VOHCirculator> Mesh::voh_range(VertexHandle v) const { return CircRange<VOHCirculator>(this, halfedge(v)); }
inline CircRange<VVCirculator>  Mesh::vv_range(VertexHandle v)  const { return CircRange<VVCirculator>(this, halfedge(v)); }
inline CircRange<FHCirculator>  Mesh::fh_range(FaceHandle f)    const { return CircRange<FHCirculator>(this, halfedge(f)); }
inline CircRange<FVCirculator>  Mesh::fv_range(FaceHandle f)    const { return CircRange<FVCirculator>(this, halfedge(f)); }

inline ElementRange<VertexHandle>   Mesh::all_vertices()  const { return ElementRange<VertexHandle>(this,   static_cast<int>(n_vertices())); }
inline ElementRange<EdgeHandle>     Mesh::all_edges()     const { return ElementRange<EdgeHandle>(this,     static_cast<int>(n_edges())); }
inline ElementRange<HalfedgeHandle> Mesh::all_halfedges() const { return ElementRange<HalfedgeHandle>(this, static_cast<int>(n_halfedges())); }
inline ElementRange<FaceHandle>     Mesh::all_faces()     const { return ElementRange<FaceHandle>(this,     static_cast<int>(n_faces())); }

} // namespace sm
