// ============================================================================
//  SimpleMesh - SmartHandles.h : ergonomic "smart" handles
// ----------------------------------------------------------------------------
//  A plain handle (Handles.h) is just an index: to navigate you must keep
//  passing it back to the mesh -
//
//      auto h2 = mesh.opposite_halfedge(mesh.next_halfedge(h));
//      auto v  = mesh.to_vertex(h2);
//
//  A *smart* handle bundles the index together with a back-pointer to the mesh,
//  so the same walk reads as a fluent chain:
//
//      auto v = h.next().opp().to();
//
//  That is the ONLY thing smart handles add - pure ergonomics over the existing
//  navigation primitives. No new data-structure concept, no new storage.
//
//  Design notes
//  ------------
//    * Each SmartXHandle PUBLICLY inherits its plain XHandle, so a smart handle
//      is-a plain handle and can be passed anywhere a plain one is expected
//      (it simply slices off the mesh pointer).
//    * Smart handles BORROW the mesh (const Mesh*). They are valid only while
//      that mesh is alive and unmodified - just like an iterator.
//    * The range accessors (.vertices(), .faces(), ...) return std::vector of
//      smart handles. We return plain vectors (rather than lazy ranges) to stay
//      readable and trivially bindable to Python, consistent with the rest of
//      the library's vector-returning circulators.
//
//  Included at the BOTTOM of Mesh.h (after Circulators.h) so the Mesh class is
//  complete and these inline methods can call its public API.
// ============================================================================
#pragma once

#include "Mesh.h"
#include <vector>

namespace sm {

struct SmartVertexHandle;
struct SmartHalfedgeHandle;
struct SmartEdgeHandle;
struct SmartFaceHandle;

// ----------------------------------------------------------------------------
//  Common base: stores the borrowed mesh pointer. Status / boundary / deleted
//  predicates are shared by all four smart handle kinds, so they live here as
//  templates over the concrete derived type (CRTP).
// ----------------------------------------------------------------------------
class SmartBaseHandle {
public:
    explicit SmartBaseHandle(const Mesh* m = nullptr) : mesh_(m) {}
    /// The mesh this handle navigates (nullptr if default-constructed).
    const Mesh* mesh() const { return mesh_; }
protected:
    const Mesh* mesh_;
};

// make_smart: wrap a plain handle + mesh into the matching smart handle.
inline SmartVertexHandle   make_smart(VertexHandle   h, const Mesh* m);
inline SmartHalfedgeHandle make_smart(HalfedgeHandle h, const Mesh* m);
inline SmartEdgeHandle     make_smart(EdgeHandle     h, const Mesh* m);
inline SmartFaceHandle     make_smart(FaceHandle     h, const Mesh* m);

// ----------------------------------------------------------------------------
//  Smart vertex
// ----------------------------------------------------------------------------
struct SmartVertexHandle : public SmartBaseHandle, public VertexHandle {
    explicit SmartVertexHandle(int idx = -1, const Mesh* m = nullptr)
        : SmartBaseHandle(m), VertexHandle(idx) {}

    // --- single-step navigation (return smart handles, so you can chain) ---
    SmartHalfedgeHandle out()      const;   ///< an outgoing halfedge
    SmartHalfedgeHandle halfedge() const;   ///< alias for out()
    SmartHalfedgeHandle in()       const;   ///< an incoming halfedge (== out().opp())

    // --- neighbourhoods (smart-handle vectors) -----------------------------
    std::vector<SmartVertexHandle>   vertices()            const;  ///< 1-ring neighbour vertices
    std::vector<SmartHalfedgeHandle> outgoing_halfedges()  const;
    std::vector<SmartHalfedgeHandle> incoming_halfedges()  const;
    std::vector<SmartEdgeHandle>     edges()               const;
    std::vector<SmartFaceHandle>     faces()               const;

    // --- queries -----------------------------------------------------------
    const Vec3& point()     const { return mesh_->point(*this); }
    Vec3        normal()    const { return mesh_->calc_vertex_normal(*this); }
    size_t      valence()   const { return mesh_->valence(*this); }
    bool is_boundary() const { return mesh_->is_boundary(static_cast<const VertexHandle&>(*this)); }
    bool is_isolated() const { return mesh_->is_isolated(*this); }
    bool deleted()     const { return mesh_->is_deleted(static_cast<const VertexHandle&>(*this)); }
    bool selected() const { return mesh_->is_selected(static_cast<const VertexHandle&>(*this)); }
    bool tagged()   const { return mesh_->is_tagged(static_cast<const VertexHandle&>(*this)); }
    bool locked()   const { return mesh_->is_locked(static_cast<const VertexHandle&>(*this)); }
    bool feature()  const { return mesh_->is_feature(static_cast<const VertexHandle&>(*this)); }
    bool hidden()   const { return mesh_->is_hidden(static_cast<const VertexHandle&>(*this)); }
};

// ----------------------------------------------------------------------------
//  Smart halfedge - the workhorse for fluent walks.
// ----------------------------------------------------------------------------
struct SmartHalfedgeHandle : public SmartBaseHandle, public HalfedgeHandle {
    explicit SmartHalfedgeHandle(int idx = -1, const Mesh* m = nullptr)
        : SmartBaseHandle(m), HalfedgeHandle(idx) {}

    SmartHalfedgeHandle next() const;
    SmartHalfedgeHandle prev() const;
    SmartHalfedgeHandle opp()  const;   ///< opposite halfedge
    SmartVertexHandle   to()   const;   ///< vertex this halfedge points to
    SmartVertexHandle   from() const;   ///< vertex this halfedge starts at
    SmartEdgeHandle     edge() const;
    SmartFaceHandle     face() const;

    Vec3 vector() const { return mesh_->calc_edge_vector(*this); }

    bool is_boundary() const { return mesh_->is_boundary(static_cast<const HalfedgeHandle&>(*this)); }
    bool deleted()     const { return mesh_->is_deleted(static_cast<const HalfedgeHandle&>(*this)); }
    bool selected() const { return mesh_->is_selected(static_cast<const HalfedgeHandle&>(*this)); }
    bool tagged()   const { return mesh_->is_tagged(static_cast<const HalfedgeHandle&>(*this)); }
    bool locked()   const { return mesh_->is_locked(static_cast<const HalfedgeHandle&>(*this)); }
    bool feature()  const { return mesh_->is_feature(static_cast<const HalfedgeHandle&>(*this)); }
    bool hidden()   const { return mesh_->is_hidden(static_cast<const HalfedgeHandle&>(*this)); }
};

// ----------------------------------------------------------------------------
//  Smart edge
// ----------------------------------------------------------------------------
struct SmartEdgeHandle : public SmartBaseHandle, public EdgeHandle {
    explicit SmartEdgeHandle(int idx = -1, const Mesh* m = nullptr)
        : SmartBaseHandle(m), EdgeHandle(idx) {}

    SmartHalfedgeHandle halfedge(int i) const;  ///< the i-th (0/1) halfedge
    SmartHalfedgeHandle h0() const { return halfedge(0); }
    SmartHalfedgeHandle h1() const { return halfedge(1); }
    SmartVertexHandle   v0() const;             ///< from-vertex of h0
    SmartVertexHandle   v1() const;             ///< to-vertex of h0
    std::vector<SmartFaceHandle> faces() const; ///< the (1 or 2) incident faces

    double length()   const { return mesh_->calc_edge_length(*this); }
    Vec3   midpoint() const { return mesh_->calc_edge_midpoint(*this); }
    double dihedral_angle() const { return mesh_->calc_dihedral_angle(*this); }

    bool is_boundary() const { return mesh_->is_boundary(static_cast<const EdgeHandle&>(*this)); }
    bool deleted()     const { return mesh_->is_deleted(static_cast<const EdgeHandle&>(*this)); }
    bool selected() const { return mesh_->is_selected(static_cast<const EdgeHandle&>(*this)); }
    bool tagged()   const { return mesh_->is_tagged(static_cast<const EdgeHandle&>(*this)); }
    bool locked()   const { return mesh_->is_locked(static_cast<const EdgeHandle&>(*this)); }
    bool feature()  const { return mesh_->is_feature(static_cast<const EdgeHandle&>(*this)); }
    bool hidden()   const { return mesh_->is_hidden(static_cast<const EdgeHandle&>(*this)); }
};

// ----------------------------------------------------------------------------
//  Smart face
// ----------------------------------------------------------------------------
struct SmartFaceHandle : public SmartBaseHandle, public FaceHandle {
    explicit SmartFaceHandle(int idx = -1, const Mesh* m = nullptr)
        : SmartBaseHandle(m), FaceHandle(idx) {}

    SmartHalfedgeHandle halfedge() const;       ///< one halfedge of the face

    std::vector<SmartVertexHandle>   vertices()  const;
    std::vector<SmartHalfedgeHandle> halfedges() const;
    std::vector<SmartEdgeHandle>     edges()     const;
    std::vector<SmartFaceHandle>     faces()     const;  ///< edge-adjacent faces

    Vec3   normal()   const { return mesh_->calc_face_normal(*this); }
    Vec3   centroid() const { return mesh_->calc_face_centroid(*this); }
    double area()     const { return mesh_->calc_face_area(*this); }
    size_t valence()  const { return mesh_->face_vertices(*this).size(); }

    bool is_boundary() const { return mesh_->is_boundary(static_cast<const FaceHandle&>(*this)); }
    bool deleted()     const { return mesh_->is_deleted(static_cast<const FaceHandle&>(*this)); }
    bool selected() const { return mesh_->is_selected(static_cast<const FaceHandle&>(*this)); }
    bool tagged()   const { return mesh_->is_tagged(static_cast<const FaceHandle&>(*this)); }
    bool locked()   const { return mesh_->is_locked(static_cast<const FaceHandle&>(*this)); }
    bool feature()  const { return mesh_->is_feature(static_cast<const FaceHandle&>(*this)); }
    bool hidden()   const { return mesh_->is_hidden(static_cast<const FaceHandle&>(*this)); }
};

// ============================================================================
//  make_smart implementations (all four structs are complete now)
// ============================================================================
inline SmartVertexHandle   make_smart(VertexHandle   h, const Mesh* m) { return SmartVertexHandle  (h.idx(), m); }
inline SmartHalfedgeHandle make_smart(HalfedgeHandle h, const Mesh* m) { return SmartHalfedgeHandle(h.idx(), m); }
inline SmartEdgeHandle     make_smart(EdgeHandle     h, const Mesh* m) { return SmartEdgeHandle    (h.idx(), m); }
inline SmartFaceHandle     make_smart(FaceHandle     h, const Mesh* m) { return SmartFaceHandle    (h.idx(), m); }

// small helper: convert a vector<PlainHandle> to vector<SmartHandle>
template <class Smart, class Plain>
inline std::vector<Smart> smart_vector(const std::vector<Plain>& in, const Mesh* m) {
    std::vector<Smart> out;
    out.reserve(in.size());
    for (const Plain& h : in) out.push_back(make_smart(h, m));
    return out;
}

// ============================================================================
//  Inline navigation bodies. Each is a one-liner over the existing Mesh API,
//  re-wrapped as a smart handle so the result keeps chaining.
// ============================================================================

// --- SmartVertexHandle ---
inline SmartHalfedgeHandle SmartVertexHandle::out()      const { return make_smart(mesh_->halfedge(static_cast<const VertexHandle&>(*this)), mesh_); }
inline SmartHalfedgeHandle SmartVertexHandle::halfedge() const { return out(); }
inline SmartHalfedgeHandle SmartVertexHandle::in()       const { return out().opp(); }

inline std::vector<SmartVertexHandle>   SmartVertexHandle::vertices()           const { return smart_vector<SmartVertexHandle>(mesh_->vertex_vertices(*this), mesh_); }
inline std::vector<SmartHalfedgeHandle> SmartVertexHandle::outgoing_halfedges() const { return smart_vector<SmartHalfedgeHandle>(mesh_->vertex_outgoing_halfedges(*this), mesh_); }
inline std::vector<SmartHalfedgeHandle> SmartVertexHandle::incoming_halfedges() const { return smart_vector<SmartHalfedgeHandle>(mesh_->vertex_incoming_halfedges(*this), mesh_); }
inline std::vector<SmartEdgeHandle>     SmartVertexHandle::edges()              const { return smart_vector<SmartEdgeHandle>(mesh_->vertex_edges(*this), mesh_); }
inline std::vector<SmartFaceHandle>     SmartVertexHandle::faces()              const { return smart_vector<SmartFaceHandle>(mesh_->vertex_faces(*this), mesh_); }

// --- SmartHalfedgeHandle ---
inline SmartHalfedgeHandle SmartHalfedgeHandle::next() const { return make_smart(mesh_->next_halfedge(*this), mesh_); }
inline SmartHalfedgeHandle SmartHalfedgeHandle::prev() const { return make_smart(mesh_->prev_halfedge(*this), mesh_); }
inline SmartHalfedgeHandle SmartHalfedgeHandle::opp()  const { return make_smart(mesh_->opposite_halfedge(*this), mesh_); }
inline SmartVertexHandle   SmartHalfedgeHandle::to()   const { return make_smart(mesh_->to_vertex(*this), mesh_); }
inline SmartVertexHandle   SmartHalfedgeHandle::from() const { return make_smart(mesh_->from_vertex(*this), mesh_); }
inline SmartEdgeHandle     SmartHalfedgeHandle::edge() const { return make_smart(mesh_->edge(*this), mesh_); }
inline SmartFaceHandle     SmartHalfedgeHandle::face() const { return make_smart(mesh_->face(*this), mesh_); }

// --- SmartEdgeHandle ---
inline SmartHalfedgeHandle SmartEdgeHandle::halfedge(int i) const { return make_smart(mesh_->halfedge(static_cast<const EdgeHandle&>(*this), i), mesh_); }
inline SmartVertexHandle   SmartEdgeHandle::v0()           const { return h0().from(); }
inline SmartVertexHandle   SmartEdgeHandle::v1()           const { return h0().to(); }
inline std::vector<SmartFaceHandle> SmartEdgeHandle::faces() const {
    std::vector<SmartFaceHandle> out;
    for (int i = 0; i < 2; ++i) {
        SmartFaceHandle f = halfedge(i).face();
        if (f.is_valid()) out.push_back(f);
    }
    return out;
}

// --- SmartFaceHandle ---
inline SmartHalfedgeHandle SmartFaceHandle::halfedge() const { return make_smart(mesh_->halfedge(static_cast<const FaceHandle&>(*this)), mesh_); }
inline std::vector<SmartVertexHandle>   SmartFaceHandle::vertices()  const { return smart_vector<SmartVertexHandle>(mesh_->face_vertices(*this), mesh_); }
inline std::vector<SmartHalfedgeHandle> SmartFaceHandle::halfedges() const { return smart_vector<SmartHalfedgeHandle>(mesh_->face_halfedges(*this), mesh_); }
inline std::vector<SmartEdgeHandle>     SmartFaceHandle::edges()     const { return smart_vector<SmartEdgeHandle>(mesh_->face_edges(*this), mesh_); }
inline std::vector<SmartFaceHandle>     SmartFaceHandle::faces()     const { return smart_vector<SmartFaceHandle>(mesh_->face_faces(*this), mesh_); }

} // namespace sm
