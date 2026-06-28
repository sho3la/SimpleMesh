// ============================================================================
//  SimpleMesh - Mesh.h : the halfedge mesh class
// ----------------------------------------------------------------------------
//  This is the public C++ API of the library. It is a simplified, single-file
//  implementation of the core ideas of a half-edge polygon connectivity:
//
//      * elements stored in flat arrays, referenced by typed handles
//      * undirected edges = pairs of opposite halfedges (opposite = h XOR 1)
//      * faces of arbitrary polygon size (triangles, quads, ...)
//      * topologically-correct add_face that stitches halfedges together
//      * navigation primitives + convenience "circulator" collectors
//
//  Trivial one-line accessors are defined inline here. The interesting,
//  algorithmic methods (add_face, find_halfedge, OBJ I/O) live in Mesh.cpp so
//  the header stays readable.
// ============================================================================
#pragma once

#include "Vec3.h"
#include "Handles.h"
#include "Items.h"
#include "Properties.h"

#include <vector>
#include <string>
#include <utility>
#include <cassert>

namespace sm {

// Forward declarations for the lazy iterators defined in Circulators.h
// (included at the bottom of this file). They are only needed as return types
// of declarations here, so incomplete types are fine.
template <class Handle> class ElementRange;
template <class Circ>   class CircRange;
class VOHCirculator;
class VVCirculator;
class FHCirculator;
class FVCirculator;

// Smart handles - defined in SmartHandles.h, included at bottom.
struct SmartVertexHandle;
struct SmartHalfedgeHandle;
struct SmartEdgeHandle;
struct SmartFaceHandle;

class Mesh {
public:
    // A Mesh owns property arrays through unique_ptr, so it is a move-only
    // resource type: cheap to move, never silently copied. Algorithms like
    // loop_subdivide return a Mesh by value (moved), and pybind11 relies on the
    // move constructor existing for return-by-value.
    Mesh() = default;
    Mesh(Mesh&&) = default;
    Mesh& operator=(Mesh&&) = default;
    Mesh(const Mesh&) = delete;
    Mesh& operator=(const Mesh&) = delete;

    // ======================================================================
    //  Construction of geometry & topology
    // ======================================================================

    /// Add a new isolated vertex at position `p`. Returns its handle.
    VertexHandle add_vertex(const Vec3& p);

    /// Add a polygonal face spanning the given vertices (given in CCW order).
    /// Returns an invalid handle if the face would violate manifoldness.
    /// This is the central, non-trivial operation - see Mesh.cpp for the full,
    /// commented algorithm.
    FaceHandle add_face(const std::vector<VertexHandle>& vertices);

    /// Convenience overload for the common triangle case.
    FaceHandle add_triangle(VertexHandle a, VertexHandle b, VertexHandle c);

    /// When true, add_face prints to std::cerr explaining why it rejected a face
    /// (complex vertex/edge). Default false: bulk importers of messy meshes would
    /// otherwise emit one line per non-manifold face (gigabytes on a bad file).
    static bool warn_invalid_face;

    // ======================================================================
    //  Counts
    // ======================================================================
    size_t n_vertices()  const { return vertices_.size(); }
    size_t n_halfedges() const { return halfedges_.size(); }
    size_t n_edges()     const { return halfedges_.size() / 2; }
    size_t n_faces()     const { return faces_.size(); }

    // ======================================================================
    //  Geometry access
    // ======================================================================
    const Vec3& point(VertexHandle v) const { return points_[v.idx()]; }
    void set_point(VertexHandle v, const Vec3& p) { points_[v.idx()] = p; }

    /// Direct access to the contiguous position buffer. Because std::vector<Vec3>
    /// is laid out as 3*N packed doubles, this is exactly what numpy needs for a
    /// zero-copy (N,3) view - see the buffer-protocol bindings.
    const std::vector<Vec3>& points() const { return points_; }
    std::vector<Vec3>&       points()       { return points_; }

    // ======================================================================
    //  Navigation primitives - the vocabulary of the halfedge structure.
    //  Master these eight functions and you can traverse anything.
    // ======================================================================

    /// An (arbitrary) outgoing halfedge of a vertex.
    HalfedgeHandle halfedge(VertexHandle v) const { return vertices_[v.idx()].outgoing_halfedge; }
    /// An (arbitrary) halfedge belonging to a face.
    HalfedgeHandle halfedge(FaceHandle f) const { return faces_[f.idx()].halfedge; }
    /// The i-th (0 or 1) halfedge of an edge.
    HalfedgeHandle halfedge(EdgeHandle e, int i) const { return HalfedgeHandle(e.idx() * 2 + i); }

    /// Vertex a halfedge points to / comes from.
    VertexHandle to_vertex(HalfedgeHandle h)   const { return halfedges_[h.idx()].to_vertex; }
    VertexHandle from_vertex(HalfedgeHandle h) const { return to_vertex(opposite_halfedge(h)); }

    /// Walk around the incident face.
    HalfedgeHandle next_halfedge(HalfedgeHandle h) const { return halfedges_[h.idx()].next; }
    HalfedgeHandle prev_halfedge(HalfedgeHandle h) const { return halfedges_[h.idx()].prev; }

    /// The opposite halfedge. Encoded purely by index arithmetic: the two
    /// halfedges of edge e are 2e and 2e+1, so flipping the lowest bit swaps
    /// them. This is why opposite() needs no stored pointer.
    HalfedgeHandle opposite_halfedge(HalfedgeHandle h) const { return HalfedgeHandle(h.idx() ^ 1); }

    /// Incident face of a halfedge (invalid for boundary halfedges).
    FaceHandle face(HalfedgeHandle h) const { return halfedges_[h.idx()].face; }

    /// The edge a halfedge belongs to, and back.
    EdgeHandle edge(HalfedgeHandle h) const { return EdgeHandle(h.idx() >> 1); }

    // ======================================================================
    //  Boundary tests
    // ======================================================================
    bool is_boundary(HalfedgeHandle h) const { return !face(h).is_valid(); }
    bool is_boundary(EdgeHandle e)     const { return is_boundary(halfedge(e, 0)) || is_boundary(halfedge(e, 1)); }
    bool is_boundary(VertexHandle v)   const;   // see Mesh.cpp
    bool is_boundary(FaceHandle f)     const;   // see Mesh.cpp

    // ======================================================================
    //  Circulators - convenience collectors.
    //  These could be lazy iterator objects; here we return plain
    //  vectors. That is less efficient but far easier to read, to bind to
    //  Python, and to learn from. Profiling-sensitive code can always drop down
    //  to the navigation primitives above.
    // ======================================================================

    /// All vertices adjacent to `v` (its 1-ring).
    std::vector<VertexHandle> vertex_vertices(VertexHandle v) const;
    /// All outgoing halfedges of `v`.
    std::vector<HalfedgeHandle> vertex_outgoing_halfedges(VertexHandle v) const;
    /// All faces incident to `v`.
    std::vector<FaceHandle> vertex_faces(VertexHandle v) const;
    /// Number of edges incident to `v`.
    size_t valence(VertexHandle v) const;

    /// All edges incident to `v`.
    std::vector<EdgeHandle> vertex_edges(VertexHandle v) const;
    /// All incoming halfedges of `v` (i.e. those pointing at `v`).
    std::vector<HalfedgeHandle> vertex_incoming_halfedges(VertexHandle v) const;

    /// The vertices of a face, in order.
    std::vector<VertexHandle> face_vertices(FaceHandle f) const;
    /// The halfedges of a face, in order.
    std::vector<HalfedgeHandle> face_halfedges(FaceHandle f) const;
    /// The edges of a face, in order.
    std::vector<EdgeHandle> face_edges(FaceHandle f) const;
    /// The faces edge-adjacent to `f` (skips boundary edges).
    std::vector<FaceHandle> face_faces(FaceHandle f) const;

    /// Find the halfedge pointing from `from` to `to`, or an invalid handle.
    HalfedgeHandle find_halfedge(VertexHandle from, VertexHandle to) const;

    // ======================================================================
    //  Lazy ranges - zero-allocation alternatives to the
    //  vector-returning circulators above, usable in range-based for loops.
    //  Defined inline in Circulators.h (included at the bottom of this file).
    // ======================================================================
    CircRange<VOHCirculator> voh_range(VertexHandle v) const;  ///< outgoing halfedges of v
    CircRange<VVCirculator>  vv_range(VertexHandle v)  const;  ///< neighbour vertices of v
    CircRange<FHCirculator>  fh_range(FaceHandle f)    const;  ///< halfedges of face f
    CircRange<FVCirculator>  fv_range(FaceHandle f)    const;  ///< vertices of face f

    ElementRange<VertexHandle>   all_vertices()  const;  ///< every live vertex
    ElementRange<EdgeHandle>     all_edges()     const;  ///< every live edge
    ElementRange<HalfedgeHandle> all_halfedges() const;  ///< every live halfedge
    ElementRange<FaceHandle>     all_faces()     const;  ///< every live face

    // ======================================================================
    //  Smart handles - wrap a plain handle + this mesh so
    //  navigation can be chained fluently: `smart(h).next().opp().to()`.
    //  Defined inline in SmartHandles.h (included at the bottom of this file).
    // ======================================================================
    SmartVertexHandle   smart(VertexHandle v)   const;
    SmartHalfedgeHandle smart(HalfedgeHandle h) const;
    SmartEdgeHandle     smart(EdgeHandle e)     const;
    SmartFaceHandle     smart(FaceHandle f)     const;

    // ======================================================================
    //  Geometry queries
    // ======================================================================
    /// Unit face normal (Newell's method; robust for polygons).
    Vec3 calc_face_normal(FaceHandle f) const;
    /// Face area (half the magnitude of the Newell normal vector).
    double calc_face_area(FaceHandle f) const;
    /// Centroid (average of the face's vertices).
    Vec3 calc_face_centroid(FaceHandle f) const;
    /// Euclidean length of an edge.
    double calc_edge_length(EdgeHandle e) const;
    /// The vector along a halfedge: point(to) - point(from).
    Vec3 calc_edge_vector(HalfedgeHandle h) const;
    /// Midpoint of an edge.
    Vec3 calc_edge_midpoint(EdgeHandle e) const;
    /// Unit vertex normal = normalized sum of incident unit face normals.
    Vec3 calc_vertex_normal(VertexHandle v) const;
    /// Signed dihedral angle (radians) across an edge, between its two incident
    /// face normals. 0 for a flat edge; returns 0 on boundary edges.
    double calc_dihedral_angle(EdgeHandle e) const;
    /// Interior angle (radians) at the corner that halfedge `h` points to,
    /// inside h's face (between -h and next(h)).
    double calc_sector_angle(HalfedgeHandle h) const;

    // --- whole-mesh measures ---------------------------------------------
    /// Sum of all face areas.
    double surface_area() const;
    /// Average of all (live) vertex positions.
    Vec3 center_of_mass() const;
    /// Axis-aligned bounding box as (min_corner, max_corner).
    std::pair<Vec3, Vec3> bounding_box() const;

    // ======================================================================
    //  Status / deletion
    // ----------------------------------------------------------------------
    //  Deletion is "lazy": elements are flagged deleted but stay in the arrays
    //  so existing handles keep their indices. Call garbage_collection() to
    //  actually compact the arrays (which DOES invalidate handles).
    // ======================================================================
    bool is_deleted(VertexHandle v)   const { return v_deleted_[v.idx()]; }
    bool is_deleted(EdgeHandle e)     const { return e_deleted_[e.idx()]; }
    bool is_deleted(HalfedgeHandle h) const { return e_deleted_[edge(h).idx()]; }
    bool is_deleted(FaceHandle f)     const { return f_deleted_[f.idx()]; }

    /// A vertex with no incident halfedge (e.g. after its faces were deleted).
    bool is_isolated(VertexHandle v) const { return !halfedge(v).is_valid(); }

    /// Flag a face deleted, detaching it from the mesh. If `delete_isolated`
    /// is true, vertices that become isolated are flagged deleted too.
    void delete_face(FaceHandle f, bool delete_isolated = true);
    /// Flag a vertex and all its incident faces deleted.
    void delete_vertex(VertexHandle v, bool delete_isolated = true);
    /// Delete the faces on both sides of an edge (and the edge if it had none).
    void delete_edge(EdgeHandle e, bool delete_isolated = true);

    /// Physically remove all deleted elements and compact the arrays.
    /// WARNING: invalidates every handle the caller is holding.
    void garbage_collection();

    // ======================================================================
    //  Topological editing operators - triangle-oriented.
    // ======================================================================

    /// Can edge `e` be flipped without creating a non-manifold/degenerate mesh?
    bool is_flip_ok(EdgeHandle e) const;
    /// Flip the shared diagonal of the two triangles adjacent to `e`.
    void flip(EdgeHandle e);

    /// Split edge `e` by inserting existing vertex `v` at its midpoint,
    /// subdividing the (up to two) incident triangles. Set v's position first.
    void split(EdgeHandle e, VertexHandle v);
    /// Split face `f` into a triangle fan around interior vertex `v`.
    void split(FaceHandle f, VertexHandle v);

    /// Is collapsing halfedge `h` (merging from-vertex into to-vertex) legal?
    bool is_collapse_ok(HalfedgeHandle h);
    /// Collapse halfedge `h`: its from-vertex is merged into its to-vertex and
    /// the incident triangles degenerate away. Marks elements deleted (lazy).
    void collapse(HalfedgeHandle h);

    /// Vertex split - the exact inverse of `collapse` (triangle meshes), the key
    /// operator for progressive meshes / view-dependent LOD. Splits vertex `v1`
    /// into the new edge `v0`-`v1`, re-creating the two wing triangles toward
    /// `vl` (left) and `vr` (right). `v0` must be a freshly added, isolated
    /// vertex (use add_vertex to place it). Either of `vl` / `vr` may be invalid
    /// for a boundary split. Returns the new halfedge pointing v0 -> v1.
    ///
    /// Round-trip: `collapse(vertex_split(v0, v1, vl, vr))` restores the mesh.
    HalfedgeHandle vertex_split(VertexHandle v0, VertexHandle v1,
                                VertexHandle vl, VertexHandle vr);

    /// Fan-triangulate every non-triangle face in place (rebuilds the mesh).
    /// Positions are preserved; the result is a pure triangle mesh.
    void triangulate();

    // ======================================================================
    //  Custom properties
    // ----------------------------------------------------------------------
    //  Attach arbitrary named, typed data arrays to any element kind. The
    //  arrays stay in lock-step with the elements automatically: adding a
    //  vertex grows every vertex property, and garbage_collection compacts
    //  them too. See Properties.h for the type-erasure mechanism.
    //
    //  These are templates, so they live in the header. Usage:
    //      auto q = mesh.add_vertex_property<double>("quality", 0.0);
    //      mesh.property(q, v) = 1.5;
    // ======================================================================
    template <class T>
    VertexPropHandle<T> add_vertex_property(const std::string& name, const T& def = T())
    { return VertexPropHandle<T>{ vprops_.add<T>(name, def) }; }
    template <class T>
    HalfedgePropHandle<T> add_halfedge_property(const std::string& name, const T& def = T())
    { return HalfedgePropHandle<T>{ hprops_.add<T>(name, def) }; }
    template <class T>
    EdgePropHandle<T> add_edge_property(const std::string& name, const T& def = T())
    { return EdgePropHandle<T>{ eprops_.add<T>(name, def) }; }
    template <class T>
    FacePropHandle<T> add_face_property(const std::string& name, const T& def = T())
    { return FacePropHandle<T>{ fprops_.add<T>(name, def) }; }

    // Element/property access. Returns a reference so you can read and write.
    template <class T> T& property(VertexPropHandle<T> p, VertexHandle v)
    { return (*vprops_.get<T>(p.id))[v.idx()]; }
    template <class T> const T& property(VertexPropHandle<T> p, VertexHandle v) const
    { return (*vprops_.get<T>(p.id))[v.idx()]; }

    template <class T> T& property(HalfedgePropHandle<T> p, HalfedgeHandle h)
    { return (*hprops_.get<T>(p.id))[h.idx()]; }
    template <class T> const T& property(HalfedgePropHandle<T> p, HalfedgeHandle h) const
    { return (*hprops_.get<T>(p.id))[h.idx()]; }

    template <class T> T& property(EdgePropHandle<T> p, EdgeHandle e)
    { return (*eprops_.get<T>(p.id))[e.idx()]; }
    template <class T> const T& property(EdgePropHandle<T> p, EdgeHandle e) const
    { return (*eprops_.get<T>(p.id))[e.idx()]; }

    template <class T> T& property(FacePropHandle<T> p, FaceHandle f)
    { return (*fprops_.get<T>(p.id))[f.idx()]; }
    template <class T> const T& property(FacePropHandle<T> p, FaceHandle f) const
    { return (*fprops_.get<T>(p.id))[f.idx()]; }

    void remove_vertex_property_by_name(const std::string& n)   { int i = vprops_.find(n); if (i >= 0) vprops_.remove(i); }
    void remove_halfedge_property_by_name(const std::string& n) { int i = hprops_.find(n); if (i >= 0) hprops_.remove(i); }
    void remove_edge_property_by_name(const std::string& n)     { int i = eprops_.find(n); if (i >= 0) eprops_.remove(i); }
    void remove_face_property_by_name(const std::string& n)     { int i = fprops_.find(n); if (i >= 0) fprops_.remove(i); }
    size_t n_vertex_properties()   const { return vprops_.n_properties(); }
    size_t n_halfedge_properties() const { return hprops_.n_properties(); }
    size_t n_edge_properties()     const { return eprops_.n_properties(); }
    size_t n_face_properties()     const { return fprops_.n_properties(); }

    // ======================================================================
    //  Status bits
    // ----------------------------------------------------------------------
    //  These status flags are stored as an ordinary built-in `int` property per
    //  element kind, created lazily on first use.
    //  The payoff of riding on the property system: status arrays grow with new
    //  elements and get COMPACTED by garbage_collection automatically - no extra
    //  code. (`deleted` stays separate because it drives iteration and GC.)
    // ======================================================================
    enum StatusBits : int {
        SELECTED = 1, TAGGED = 2, LOCKED = 4, FEATURE = 8, HIDDEN = 16
    };

    // raw status word per element (0 if the status property doesn't exist yet)
    int  get_status(VertexHandle h)   const { return v_status_id_ < 0 ? 0 : (*vprops_.get<int>(v_status_id_))[h.idx()]; }
    int  get_status(HalfedgeHandle h) const { return h_status_id_ < 0 ? 0 : (*hprops_.get<int>(h_status_id_))[h.idx()]; }
    int  get_status(EdgeHandle h)     const { return e_status_id_ < 0 ? 0 : (*eprops_.get<int>(e_status_id_))[h.idx()]; }
    int  get_status(FaceHandle h)     const { return f_status_id_ < 0 ? 0 : (*fprops_.get<int>(f_status_id_))[h.idx()]; }

    void set_status(VertexHandle h, int s)   { ensure_v_status(); (*vprops_.get<int>(v_status_id_))[h.idx()] = s; }
    void set_status(HalfedgeHandle h, int s) { ensure_h_status(); (*hprops_.get<int>(h_status_id_))[h.idx()] = s; }
    void set_status(EdgeHandle h, int s)     { ensure_e_status(); (*eprops_.get<int>(e_status_id_))[h.idx()] = s; }
    void set_status(FaceHandle h, int s)     { ensure_f_status(); (*fprops_.get<int>(f_status_id_))[h.idx()] = s; }

    // generic bit test / change, usable with any element handle
    template <class H> bool test_status(H h, int bit) const { return (get_status(h) & bit) != 0; }
    template <class H> void change_status(H h, int bit, bool on) {
        int s = get_status(h);
        set_status(h, on ? (s | bit) : (s & ~bit));
    }

    // named convenience accessors (templated over the element handle type)
    template <class H> bool is_selected(H h) const { return test_status(h, SELECTED); }
    template <class H> bool is_tagged(H h)   const { return test_status(h, TAGGED); }
    template <class H> bool is_locked(H h)   const { return test_status(h, LOCKED); }
    template <class H> bool is_feature(H h)  const { return test_status(h, FEATURE); }
    template <class H> bool is_hidden(H h)   const { return test_status(h, HIDDEN); }
    template <class H> void set_selected(H h, bool v = true) { change_status(h, SELECTED, v); }
    template <class H> void set_tagged(H h, bool v = true)   { change_status(h, TAGGED, v); }
    template <class H> void set_locked(H h, bool v = true)   { change_status(h, LOCKED, v); }
    template <class H> void set_feature(H h, bool v = true)  { change_status(h, FEATURE, v); }
    template <class H> void set_hidden(H h, bool v = true)   { change_status(h, HIDDEN, v); }

    // ======================================================================
    //  File I/O. OBJ, PLY and STL.
    //  All writers skip deleted elements and remap indices, so they are safe to
    //  call before garbage_collection. Polygon faces are fan-triangulated when
    //  the format (STL) requires triangles.
    // ======================================================================
    bool read_obj(const std::string& path);
    bool write_obj(const std::string& path) const;

    bool read_ply(const std::string& path);
    bool write_ply(const std::string& path, bool binary = false) const;

    bool read_stl(const std::string& path);
    bool write_stl(const std::string& path, bool binary = true) const;

    bool read_off(const std::string& path);
    bool write_off(const std::string& path) const;

private:
    // ----- low-level element factories (used by add_face) -----------------
    VertexHandle   new_vertex(const Vec3& p);
    HalfedgeHandle new_edge(VertexHandle from, VertexHandle to);  // creates 2 halfedges, returns from->to
    FaceHandle     new_face();

    // ----- low-level setters ---------------------------------------------
    void set_next(HalfedgeHandle h, HalfedgeHandle nx) { halfedges_[h.idx()].next = nx; halfedges_[nx.idx()].prev = h; }
    void set_face(HalfedgeHandle h, FaceHandle f)      { halfedges_[h.idx()].face = f; }
    void set_to_vertex(HalfedgeHandle h, VertexHandle v) { halfedges_[h.idx()].to_vertex = v; }
    void set_halfedge(VertexHandle v, HalfedgeHandle h)  { vertices_[v.idx()].outgoing_halfedge = h; }
    void set_halfedge(FaceHandle f, HalfedgeHandle h)    { faces_[f.idx()].halfedge = h; }

    /// Mark a halfedge as a boundary halfedge (no incident face).
    void set_boundary(HalfedgeHandle h) { halfedges_[h.idx()].face = FaceHandle(); }
    /// Mark a vertex isolated (no outgoing halfedge).
    void set_isolated(VertexHandle v)   { vertices_[v.idx()].outgoing_halfedge = HalfedgeHandle(); }

    // collapse() helpers (see Mesh.cpp)
    void collapse_edge(HalfedgeHandle h);
    void collapse_loop(HalfedgeHandle h);

    // vertex_split() helpers - the inverses of collapse_loop / collapse_edge.
    // insert_loop turns the (open) face left of `h` into a triangle by adding a
    // matching halfedge; insert_edge re-introduces the split edge at vertex v0.
    HalfedgeHandle insert_loop(HalfedgeHandle h);
    HalfedgeHandle insert_edge(VertexHandle v0, HalfedgeHandle h0, HalfedgeHandle h1);

    /// Ensure the vertex's stored outgoing halfedge is a boundary one if any
    /// boundary halfedge exists. This invariant lets find_halfedge / add_face
    /// reason about boundaries correctly.
    void adjust_outgoing_halfedge(VertexHandle v);

    // ----- lazily-created status properties --------------------------------
    int v_status_id_ = -1, h_status_id_ = -1, e_status_id_ = -1, f_status_id_ = -1;
    void ensure_v_status() { if (v_status_id_ < 0) v_status_id_ = vprops_.add<int>("v:status", 0); }
    void ensure_h_status() { if (h_status_id_ < 0) h_status_id_ = hprops_.add<int>("h:status", 0); }
    void ensure_e_status() { if (e_status_id_ < 0) e_status_id_ = eprops_.add<int>("e:status", 0); }
    void ensure_f_status() { if (f_status_id_ < 0) f_status_id_ = fprops_.add<int>("f:status", 0); }

    // ----- storage --------------------------------------------------------
    std::vector<Vertex>   vertices_;
    std::vector<Vec3>     points_;      // parallel to vertices_
    std::vector<Halfedge> halfedges_;   // always an even count (edges = pairs)
    std::vector<Face>     faces_;

    // ----- per-element "deleted" status -----------------------------------
    // Parallel to the arrays above. Halfedges share their edge's flag.
    std::vector<bool>     v_deleted_;   // parallel to vertices_
    std::vector<bool>     e_deleted_;   // parallel to edges  (= halfedges_/2)
    std::vector<bool>     f_deleted_;   // parallel to faces_

    // ----- custom property containers -------------------------------------
    PropertyContainer     vprops_;      // one slot per vertex
    PropertyContainer     hprops_;      // one slot per halfedge
    PropertyContainer     eprops_;      // one slot per edge
    PropertyContainer     fprops_;      // one slot per face
};

} // namespace sm

// The lazy iterators need the fully-defined Mesh, so they are
// included here, after the class. See the header for the mutual-include note.
#include "Circulators.h"

// Smart handles likewise need the complete Mesh, and their
// factory methods are defined inline against the now-complete smart types.
#include "SmartHandles.h"

namespace sm {
inline SmartVertexHandle   Mesh::smart(VertexHandle v)   const { return make_smart(v, this); }
inline SmartHalfedgeHandle Mesh::smart(HalfedgeHandle h) const { return make_smart(h, this); }
inline SmartEdgeHandle     Mesh::smart(EdgeHandle e)     const { return make_smart(e, this); }
inline SmartFaceHandle     Mesh::smart(FaceHandle f)     const { return make_smart(f, this); }
} // namespace sm
