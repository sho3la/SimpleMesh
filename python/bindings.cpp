// ============================================================================
//  SimpleMesh - bindings.cpp : the pybind11 Python wrapper
// ----------------------------------------------------------------------------
//  This file is the entire "plugin" layer. It turns the C++ `sm::Mesh` class
//  into a Python module called `simplemesh`. The key concepts demonstrated:
//
//    * PYBIND11_MODULE   : defines the module entry point the Python import
//                          machinery looks for (PyInit_simplemesh).
//    * py::class_<T>     : exposes a C++ class, its constructor, methods and
//                          properties to Python.
//    * lambdas as glue   : where the C++ signature isn't directly Pythonic
//                          (e.g. Vec3 <-> tuple), we wrap it in a small lambda.
//    * automatic STL     : <pybind11/stl.h> auto-converts std::vector <-> list.
//
//  The handle types are exposed as opaque little classes carrying `.idx()` and
//  `.is_valid()`, mirroring the C++ API so the Python feels familiar.
// ============================================================================
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>   // py::array_t and the buffer protocol

#include "simplemesh/Mesh.h"
#include "simplemesh/Algorithms.h"
#include "simplemesh/Decimater.h"
#include "simplemesh/MeshChecker.h"
#include "simplemesh/MeshRepair.h"
#include "simplemesh/PreMesh.h"

#include <array>
#include <stdexcept>

namespace py = pybind11;
using namespace sm;

// The zero-copy numpy view below relies on std::vector<Vec3> being a tightly
// packed block of doubles (x0,y0,z0,x1,...). Guarantee that at compile time.
static_assert(sizeof(Vec3) == 3 * sizeof(double),
              "Vec3 must be exactly 3 packed doubles for the numpy buffer view");

// numpy arrays we accept as input: contiguous, and silently cast (e.g. int
// lists, float32) to the required dtype.
using ArrayD = py::array_t<double, py::array::c_style | py::array::forcecast>;
using ArrayI = py::array_t<int,    py::array::c_style | py::array::forcecast>;

PYBIND11_MODULE(simplemesh, m) {
    m.doc() = "SimpleMesh - a teaching halfedge mesh library (C++ core, Python API)";

    // ----- Vec3 -----------------------------------------------------------
    py::class_<Vec3>(m, "Vec3")
        .def(py::init<>())
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def_readwrite("x", &Vec3::x)
        .def_readwrite("y", &Vec3::y)
        .def_readwrite("z", &Vec3::z)
        .def("dot",   &Vec3::dot)
        .def("cross", &Vec3::cross)
        .def("norm",  &Vec3::norm)
        .def("normalized", &Vec3::normalized)
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def("__repr__", [](const Vec3& v) {
            return "Vec3(" + std::to_string(v.x) + ", " + std::to_string(v.y) +
                   ", " + std::to_string(v.z) + ")";
        });

    // ----- Handles --------------------------------------------------------
    // A tiny macro to bind each handle type identically.
#define BIND_HANDLE(CppType, PyName)                                          \
    py::class_<CppType>(m, PyName)                                            \
        .def(py::init<>())                                                    \
        .def(py::init<int>(), py::arg("idx"))                                 \
        .def("idx", &CppType::idx)                                            \
        .def("is_valid", &CppType::is_valid)                                  \
        .def("__eq__", [](const CppType& a, const CppType& b){ return a == b; }) \
        .def("__hash__", [](const CppType& h){ return h.idx(); })            \
        .def("__repr__", [](const CppType& h){                                \
            return std::string(PyName) + "(" + std::to_string(h.idx()) + ")"; })

    BIND_HANDLE(VertexHandle,   "VertexHandle");
    BIND_HANDLE(HalfedgeHandle, "HalfedgeHandle");
    BIND_HANDLE(EdgeHandle,     "EdgeHandle");
    BIND_HANDLE(FaceHandle,     "FaceHandle");
#undef BIND_HANDLE

    // ----- Smart handles --------------------------------------------------
    // A smart handle is a plain handle + a borrowed mesh pointer, so navigation
    // chains fluently: `mesh.smart_halfedge(h).next().opp().to()`. They inherit
    // the plain handle, so `.idx()` / `.is_valid()` still work, and they convert
    // implicitly to plain handles when passed back into Mesh methods.
    //
    // Lifetime: a smart handle BORROWS the mesh. The factory methods on Mesh use
    // keep_alive so the first smart handle keeps the mesh alive; keep the mesh in
    // a Python variable while you chain off it.
    py::class_<SmartVertexHandle>(m, "SmartVertexHandle")
        .def("idx", &SmartVertexHandle::idx)
        .def("is_valid", &SmartVertexHandle::is_valid)
        .def("out", &SmartVertexHandle::out)
        .def("halfedge", &SmartVertexHandle::halfedge)
        .def("in_", &SmartVertexHandle::in)
        .def("vertices", &SmartVertexHandle::vertices)
        .def("outgoing_halfedges", &SmartVertexHandle::outgoing_halfedges)
        .def("incoming_halfedges", &SmartVertexHandle::incoming_halfedges)
        .def("edges", &SmartVertexHandle::edges)
        .def("faces", &SmartVertexHandle::faces)
        .def("point", &SmartVertexHandle::point, py::return_value_policy::reference_internal)
        .def("normal", &SmartVertexHandle::normal)
        .def("valence", &SmartVertexHandle::valence)
        .def("is_boundary", &SmartVertexHandle::is_boundary)
        .def("is_isolated", &SmartVertexHandle::is_isolated)
        .def("deleted", &SmartVertexHandle::deleted)
        .def("selected", &SmartVertexHandle::selected)
        .def("tagged", &SmartVertexHandle::tagged)
        .def("locked", &SmartVertexHandle::locked)
        .def("feature", &SmartVertexHandle::feature)
        .def("hidden", &SmartVertexHandle::hidden)
        .def("__repr__", [](const SmartVertexHandle& h){ return "SmartVertexHandle(" + std::to_string(h.idx()) + ")"; });

    py::class_<SmartHalfedgeHandle>(m, "SmartHalfedgeHandle")
        .def("idx", &SmartHalfedgeHandle::idx)
        .def("is_valid", &SmartHalfedgeHandle::is_valid)
        .def("next", &SmartHalfedgeHandle::next)
        .def("prev", &SmartHalfedgeHandle::prev)
        .def("opp", &SmartHalfedgeHandle::opp)
        .def("to", &SmartHalfedgeHandle::to)
        .def("from_", &SmartHalfedgeHandle::from)
        .def("edge", &SmartHalfedgeHandle::edge)
        .def("face", &SmartHalfedgeHandle::face)
        .def("vector", &SmartHalfedgeHandle::vector)
        .def("is_boundary", &SmartHalfedgeHandle::is_boundary)
        .def("deleted", &SmartHalfedgeHandle::deleted)
        .def("selected", &SmartHalfedgeHandle::selected)
        .def("tagged", &SmartHalfedgeHandle::tagged)
        .def("locked", &SmartHalfedgeHandle::locked)
        .def("feature", &SmartHalfedgeHandle::feature)
        .def("hidden", &SmartHalfedgeHandle::hidden)
        .def("__repr__", [](const SmartHalfedgeHandle& h){ return "SmartHalfedgeHandle(" + std::to_string(h.idx()) + ")"; });

    py::class_<SmartEdgeHandle>(m, "SmartEdgeHandle")
        .def("idx", &SmartEdgeHandle::idx)
        .def("is_valid", &SmartEdgeHandle::is_valid)
        .def("halfedge", &SmartEdgeHandle::halfedge, py::arg("i"))
        .def("h0", &SmartEdgeHandle::h0)
        .def("h1", &SmartEdgeHandle::h1)
        .def("v0", &SmartEdgeHandle::v0)
        .def("v1", &SmartEdgeHandle::v1)
        .def("faces", &SmartEdgeHandle::faces)
        .def("length", &SmartEdgeHandle::length)
        .def("midpoint", &SmartEdgeHandle::midpoint)
        .def("dihedral_angle", &SmartEdgeHandle::dihedral_angle)
        .def("is_boundary", &SmartEdgeHandle::is_boundary)
        .def("deleted", &SmartEdgeHandle::deleted)
        .def("selected", &SmartEdgeHandle::selected)
        .def("tagged", &SmartEdgeHandle::tagged)
        .def("locked", &SmartEdgeHandle::locked)
        .def("feature", &SmartEdgeHandle::feature)
        .def("hidden", &SmartEdgeHandle::hidden)
        .def("__repr__", [](const SmartEdgeHandle& h){ return "SmartEdgeHandle(" + std::to_string(h.idx()) + ")"; });

    py::class_<SmartFaceHandle>(m, "SmartFaceHandle")
        .def("idx", &SmartFaceHandle::idx)
        .def("is_valid", &SmartFaceHandle::is_valid)
        .def("halfedge", &SmartFaceHandle::halfedge)
        .def("vertices", &SmartFaceHandle::vertices)
        .def("halfedges", &SmartFaceHandle::halfedges)
        .def("edges", &SmartFaceHandle::edges)
        .def("faces", &SmartFaceHandle::faces)
        .def("normal", &SmartFaceHandle::normal)
        .def("centroid", &SmartFaceHandle::centroid)
        .def("area", &SmartFaceHandle::area)
        .def("valence", &SmartFaceHandle::valence)
        .def("is_boundary", &SmartFaceHandle::is_boundary)
        .def("deleted", &SmartFaceHandle::deleted)
        .def("selected", &SmartFaceHandle::selected)
        .def("tagged", &SmartFaceHandle::tagged)
        .def("locked", &SmartFaceHandle::locked)
        .def("feature", &SmartFaceHandle::feature)
        .def("hidden", &SmartFaceHandle::hidden)
        .def("__repr__", [](const SmartFaceHandle& h){ return "SmartFaceHandle(" + std::to_string(h.idx()) + ")"; });

    // ----- Property handles -----------------------------------------------
    // Templates can't be bound generically, so we register the concrete
    // instantiations. For full C++ parity we cover every element kind
    // (vertex/halfedge/edge/face) x the three common value types
    // (float -> double, int, vec3 -> Vec3). Each handle is an opaque id holder.
#define REG_PROP_HANDLE(ELEMCAP, T, PYNAME) \
    py::class_<ELEMCAP##PropHandle<T>>(m, PYNAME).def("is_valid", &ELEMCAP##PropHandle<T>::is_valid)
    REG_PROP_HANDLE(Vertex,   double, "VertexFloatProperty");
    REG_PROP_HANDLE(Vertex,   int,    "VertexIntProperty");
    REG_PROP_HANDLE(Vertex,   Vec3,   "VertexVec3Property");
    REG_PROP_HANDLE(Halfedge, double, "HalfedgeFloatProperty");
    REG_PROP_HANDLE(Halfedge, int,    "HalfedgeIntProperty");
    REG_PROP_HANDLE(Halfedge, Vec3,   "HalfedgeVec3Property");
    REG_PROP_HANDLE(Edge,     double, "EdgeFloatProperty");
    REG_PROP_HANDLE(Edge,     int,    "EdgeIntProperty");
    REG_PROP_HANDLE(Edge,     Vec3,   "EdgeVec3Property");
    REG_PROP_HANDLE(Face,     double, "FaceFloatProperty");
    REG_PROP_HANDLE(Face,     int,    "FaceIntProperty");
    REG_PROP_HANDLE(Face,     Vec3,   "FaceVec3Property");
#undef REG_PROP_HANDLE

    // ----- Mesh -----------------------------------------------------------
    auto mesh_cls = py::class_<Mesh>(m, "Mesh")
        .def(py::init<>())

        // construction
        .def("add_vertex", &Mesh::add_vertex, py::arg("point"),
             "Add an isolated vertex and return its VertexHandle.")
        .def("add_face", &Mesh::add_face, py::arg("vertices"),
             "Add a polygon face from a list of VertexHandles (CCW order).")
        .def("add_triangle", &Mesh::add_triangle,
             py::arg("a"), py::arg("b"), py::arg("c"))

        // counts (exposed as read-only properties: mesh.n_vertices)
        .def_property_readonly("n_vertices",  &Mesh::n_vertices)
        .def_property_readonly("n_halfedges", &Mesh::n_halfedges)
        .def_property_readonly("n_edges",     &Mesh::n_edges)
        .def_property_readonly("n_faces",     &Mesh::n_faces)

        // geometry access
        .def("point", &Mesh::point, py::arg("vertex"),
             py::return_value_policy::reference_internal)
        .def("set_point", &Mesh::set_point, py::arg("vertex"), py::arg("point"))

        // navigation primitives
        .def("halfedge", py::overload_cast<VertexHandle>(&Mesh::halfedge, py::const_))
        .def("face_halfedge", py::overload_cast<FaceHandle>(&Mesh::halfedge, py::const_))
        .def("edge_halfedge", py::overload_cast<EdgeHandle, int>(&Mesh::halfedge, py::const_),
             py::arg("edge"), py::arg("i"), "The i-th (0 or 1) halfedge of an edge.")
        .def("to_vertex",   &Mesh::to_vertex)
        .def("from_vertex", &Mesh::from_vertex)
        .def("next_halfedge",     &Mesh::next_halfedge)
        .def("prev_halfedge",     &Mesh::prev_halfedge)
        .def("opposite_halfedge", &Mesh::opposite_halfedge)
        .def("face",  &Mesh::face)
        .def("edge",  &Mesh::edge)

        // boundary tests (overloaded -> disambiguate explicitly)
        .def("is_boundary_halfedge", py::overload_cast<HalfedgeHandle>(&Mesh::is_boundary, py::const_))
        .def("is_boundary_edge",     py::overload_cast<EdgeHandle>(&Mesh::is_boundary, py::const_))
        .def("is_boundary_vertex",   py::overload_cast<VertexHandle>(&Mesh::is_boundary, py::const_))
        .def("is_boundary_face",     py::overload_cast<FaceHandle>(&Mesh::is_boundary, py::const_))

        // circulators
        .def("vertex_vertices", &Mesh::vertex_vertices)
        .def("vertex_outgoing_halfedges", &Mesh::vertex_outgoing_halfedges)
        .def("vertex_incoming_halfedges", &Mesh::vertex_incoming_halfedges)
        .def("vertex_faces", &Mesh::vertex_faces)
        .def("vertex_edges", &Mesh::vertex_edges)
        .def("valence", &Mesh::valence)
        .def("face_vertices", &Mesh::face_vertices)
        .def("face_halfedges", &Mesh::face_halfedges)
        .def("face_edges", &Mesh::face_edges)
        .def("face_faces", &Mesh::face_faces)
        .def("find_halfedge", &Mesh::find_halfedge)

        // smart-handle factories. keep_alive<0,1> ties the
        // returned smart handle's lifetime to this mesh so its borrowed pointer
        // can't dangle while you chain navigation off it.
        .def("smart_vertex",   [](const Mesh& me, VertexHandle v)   { return me.smart(v); }, py::arg("vertex"),   py::keep_alive<0, 1>())
        .def("smart_halfedge", [](const Mesh& me, HalfedgeHandle h) { return me.smart(h); }, py::arg("halfedge"), py::keep_alive<0, 1>())
        .def("smart_edge",     [](const Mesh& me, EdgeHandle e)     { return me.smart(e); }, py::arg("edge"),     py::keep_alive<0, 1>())
        .def("smart_face",     [](const Mesh& me, FaceHandle f)     { return me.smart(f); }, py::arg("face"),     py::keep_alive<0, 1>())

        // geometry queries
        .def("calc_face_normal",   &Mesh::calc_face_normal)
        .def("calc_face_area",     &Mesh::calc_face_area)
        .def("calc_face_centroid", &Mesh::calc_face_centroid)
        .def("calc_edge_length",   &Mesh::calc_edge_length)
        .def("calc_edge_vector",   &Mesh::calc_edge_vector)
        .def("calc_edge_midpoint", &Mesh::calc_edge_midpoint)
        .def("calc_vertex_normal", &Mesh::calc_vertex_normal)
        .def("calc_dihedral_angle",&Mesh::calc_dihedral_angle)
        .def("calc_sector_angle",  &Mesh::calc_sector_angle)
        .def("surface_area",       &Mesh::surface_area)
        .def("center_of_mass",     &Mesh::center_of_mass)
        .def("bounding_box",       &Mesh::bounding_box)
        .def("triangulate",        &Mesh::triangulate)

        // status / deletion
        .def("is_deleted_vertex",   py::overload_cast<VertexHandle>(&Mesh::is_deleted, py::const_))
        .def("is_deleted_halfedge", py::overload_cast<HalfedgeHandle>(&Mesh::is_deleted, py::const_))
        .def("is_deleted_edge",     py::overload_cast<EdgeHandle>(&Mesh::is_deleted, py::const_))
        .def("is_deleted_face",     py::overload_cast<FaceHandle>(&Mesh::is_deleted, py::const_))
        .def("is_isolated", &Mesh::is_isolated)
        .def("delete_face",   &Mesh::delete_face,   py::arg("face"),   py::arg("delete_isolated") = true)
        .def("delete_vertex", &Mesh::delete_vertex, py::arg("vertex"), py::arg("delete_isolated") = true)
        .def("delete_edge",   &Mesh::delete_edge,   py::arg("edge"),   py::arg("delete_isolated") = true)
        .def("garbage_collection", &Mesh::garbage_collection,
             "Compact the arrays, removing deleted elements. Invalidates all handles.")

        // topological editing
        .def("is_flip_ok", &Mesh::is_flip_ok, py::arg("edge"))
        .def("flip",       &Mesh::flip,       py::arg("edge"))
        .def("split_edge", py::overload_cast<EdgeHandle, VertexHandle>(&Mesh::split),
             py::arg("edge"), py::arg("vertex"))
        .def("split_face", py::overload_cast<FaceHandle, VertexHandle>(&Mesh::split),
             py::arg("face"), py::arg("vertex"))
        .def("is_collapse_ok", &Mesh::is_collapse_ok, py::arg("halfedge"))
        .def("collapse",       &Mesh::collapse,       py::arg("halfedge"))
        .def("vertex_split",   &Mesh::vertex_split,
             py::arg("v0"), py::arg("v1"), py::arg("vl"), py::arg("vr"),
             "Inverse of collapse: split v1 into the new edge v0-v1 (v0 must be a "
             "fresh isolated vertex), recreating the wing triangles toward vl/vr "
             "(either may be an invalid handle for a boundary split).")

        // custom properties: the per-(element,type) add/get/set trios are
        // generated below the class via the BIND_PROP_METHODS macro. Here we
        // only expose the property-management helpers.
        .def("remove_vertex_property",   &Mesh::remove_vertex_property_by_name,   py::arg("name"))
        .def("remove_halfedge_property", &Mesh::remove_halfedge_property_by_name, py::arg("name"))
        .def("remove_edge_property",     &Mesh::remove_edge_property_by_name,     py::arg("name"))
        .def("remove_face_property",     &Mesh::remove_face_property_by_name,     py::arg("name"))
        .def_property_readonly("n_vertex_properties",   &Mesh::n_vertex_properties)
        .def_property_readonly("n_halfedge_properties", &Mesh::n_halfedge_properties)
        .def_property_readonly("n_edge_properties",     &Mesh::n_edge_properties)
        .def_property_readonly("n_face_properties",     &Mesh::n_face_properties)

        // ==================================================================
        //  numpy bulk / buffer-protocol bindings
        // ==================================================================

        // ZERO-COPY: a writable (N,3) view that aliases the mesh's own memory.
        // Mutating the returned array mutates the mesh in place. We take the
        // Python object `self` and pass it as the array's `base`, tying the
        // buffer's lifetime to the mesh so the memory can't be freed while
        // numpy still references it. (Using `self` avoids copying the mesh.)
        .def("points_view", [](py::object self) {
            Mesh& me = self.cast<Mesh&>();
            auto& pts = me.points();
            return py::array_t<double>(
                { (py::ssize_t)pts.size(), (py::ssize_t)3 },                       // shape
                { (py::ssize_t)sizeof(Vec3), (py::ssize_t)sizeof(double) },        // strides
                pts.empty() ? nullptr : &pts[0].x,                                 // data
                self);                                                             // base (no copy)
        }, "Zero-copy writable (N,3) view of the vertex positions.")

        // COPY-OUT / COPY-IN of positions.
        .def("get_vertices", [](const Mesh& me) {
            const auto& pts = me.points();
            py::array_t<double> arr({ (py::ssize_t)pts.size(), (py::ssize_t)3 });
            auto w = arr.mutable_unchecked<2>();
            for (py::ssize_t i = 0; i < (py::ssize_t)pts.size(); ++i) {
                w(i, 0) = pts[i].x; w(i, 1) = pts[i].y; w(i, 2) = pts[i].z;
            }
            return arr;
        }, "Return positions as a fresh (N,3) float64 array.")
        .def("set_vertices", [](Mesh& me, ArrayD arr) {
            if (arr.ndim() != 2 || arr.shape(1) != 3)
                throw std::runtime_error("set_vertices: expected an (N,3) array");
            if ((size_t)arr.shape(0) != me.n_vertices())
                throw std::runtime_error("set_vertices: row count must equal n_vertices");
            auto r = arr.unchecked<2>();
            for (py::ssize_t i = 0; i < arr.shape(0); ++i)
                me.set_point(VertexHandle((int)i), Vec3(r(i, 0), r(i, 1), r(i, 2)));
        }, py::arg("array"))
        .def("add_vertices", [](Mesh& me, ArrayD arr) {
            if (arr.ndim() != 2 || arr.shape(1) != 3)
                throw std::runtime_error("add_vertices: expected an (N,3) array");
            auto r = arr.unchecked<2>();
            for (py::ssize_t i = 0; i < arr.shape(0); ++i)
                me.add_vertex(Vec3(r(i, 0), r(i, 1), r(i, 2)));
            return (int)arr.shape(0);
        }, py::arg("array"), "Append vertices from an (N,3) array; returns the count added.")

        // Faces as an (M,3) int array (triangle meshes). Run garbage_collection
        // first if you have deleted elements, so indices stay dense.
        .def("get_triangles", [](const Mesh& me) {
            std::vector<std::array<int, 3>> tris;
            for (size_t f = 0; f < me.n_faces(); ++f) {
                FaceHandle fh((int)f);
                if (me.is_deleted(fh)) continue;
                auto vs = me.face_vertices(fh);
                if (vs.size() != 3)
                    throw std::runtime_error("get_triangles: mesh has a non-triangle face");
                tris.push_back({ vs[0].idx(), vs[1].idx(), vs[2].idx() });
            }
            py::array_t<int> arr({ (py::ssize_t)tris.size(), (py::ssize_t)3 });
            auto w = arr.mutable_unchecked<2>();
            for (py::ssize_t i = 0; i < (py::ssize_t)tris.size(); ++i)
                for (int k = 0; k < 3; ++k) w(i, k) = tris[i][k];
            return arr;
        }, "Return triangle connectivity as an (M,3) int32 array.")
        .def("add_triangles", [](Mesh& me, ArrayI arr) {
            if (arr.ndim() != 2 || arr.shape(1) != 3)
                throw std::runtime_error("add_triangles: expected an (M,3) array");
            auto r = arr.unchecked<2>();
            int added = 0;
            for (py::ssize_t i = 0; i < arr.shape(0); ++i)
                if (me.add_triangle(VertexHandle(r(i, 0)), VertexHandle(r(i, 1)),
                                    VertexHandle(r(i, 2))).is_valid()) ++added;
            return added;
        }, py::arg("array"), "Add triangles from an (M,3) index array; returns count added.")

        // Bulk transfer of whole property arrays.
        .def("vertex_float_to_numpy", [](Mesh& me, VertexPropHandle<double> p) {
            py::array_t<double> arr((py::ssize_t)me.n_vertices());
            auto w = arr.mutable_unchecked<1>();
            for (py::ssize_t i = 0; i < (py::ssize_t)me.n_vertices(); ++i)
                w(i) = me.property(p, VertexHandle((int)i));
            return arr;
        })
        .def("vertex_float_from_numpy", [](Mesh& me, VertexPropHandle<double> p, ArrayD arr) {
            if (arr.ndim() != 1 || (size_t)arr.shape(0) != me.n_vertices())
                throw std::runtime_error("expected a length-n_vertices 1D array");
            auto r = arr.unchecked<1>();
            for (py::ssize_t i = 0; i < arr.shape(0); ++i)
                me.property(p, VertexHandle((int)i)) = r(i);
        })
        .def("vertex_vec3_to_numpy", [](Mesh& me, VertexPropHandle<Vec3> p) {
            py::array_t<double> arr({ (py::ssize_t)me.n_vertices(), (py::ssize_t)3 });
            auto w = arr.mutable_unchecked<2>();
            for (py::ssize_t i = 0; i < (py::ssize_t)me.n_vertices(); ++i) {
                const Vec3& v = me.property(p, VertexHandle((int)i));
                w(i, 0) = v.x; w(i, 1) = v.y; w(i, 2) = v.z;
            }
            return arr;
        })
        .def("vertex_vec3_from_numpy", [](Mesh& me, VertexPropHandle<Vec3> p, ArrayD arr) {
            if (arr.ndim() != 2 || arr.shape(1) != 3 || (size_t)arr.shape(0) != me.n_vertices())
                throw std::runtime_error("expected an (n_vertices,3) array");
            auto r = arr.unchecked<2>();
            for (py::ssize_t i = 0; i < arr.shape(0); ++i)
                me.property(p, VertexHandle((int)i)) = Vec3(r(i, 0), r(i, 1), r(i, 2));
        })

        // ==================================================================
        //  Lazy iterators. py::make_iterator turns our C++
        //  begin/end pair into a Python iterator, so you can write
        //  `for v in mesh.vertices(): ...`. keep_alive<0,1> ties the
        //  iterator's lifetime to the mesh (arg 1 = self) so it can't dangle.
        // ==================================================================
        .def("vertices", [](const Mesh& me) {
            auto r = me.all_vertices();
            return py::make_iterator(r.begin(), r.end());
        }, py::keep_alive<0, 1>(), "Iterate live vertices lazily.")
        .def("faces", [](const Mesh& me) {
            auto r = me.all_faces();
            return py::make_iterator(r.begin(), r.end());
        }, py::keep_alive<0, 1>(), "Iterate live faces lazily.")
        .def("edges", [](const Mesh& me) {
            auto r = me.all_edges();
            return py::make_iterator(r.begin(), r.end());
        }, py::keep_alive<0, 1>(), "Iterate live edges lazily.")
        .def("halfedges", [](const Mesh& me) {
            auto r = me.all_halfedges();
            return py::make_iterator(r.begin(), r.end());
        }, py::keep_alive<0, 1>(), "Iterate live halfedges lazily.")
        .def("vv", [](const Mesh& me, VertexHandle v) {
            auto r = me.vv_range(v);
            return py::make_iterator(r.begin(), r.end());
        }, py::arg("vertex"), py::keep_alive<0, 1>(), "Lazily iterate the 1-ring neighbours of a vertex.")
        .def("voh", [](const Mesh& me, VertexHandle v) {
            auto r = me.voh_range(v);
            return py::make_iterator(r.begin(), r.end());
        }, py::arg("vertex"), py::keep_alive<0, 1>(), "Lazily iterate the outgoing halfedges of a vertex.")
        .def("fv", [](const Mesh& me, FaceHandle f) {
            auto r = me.fv_range(f);
            return py::make_iterator(r.begin(), r.end());
        }, py::arg("face"), py::keep_alive<0, 1>(), "Lazily iterate the vertices of a face.")
        .def("fh", [](const Mesh& me, FaceHandle f) {
            auto r = me.fh_range(f);
            return py::make_iterator(r.begin(), r.end());
        }, py::arg("face"), py::keep_alive<0, 1>(), "Lazily iterate the halfedges of a face.")

        // I/O (OBJ + PLY/STL)
        .def("read_obj",  &Mesh::read_obj,  py::arg("path"))
        .def("write_obj", &Mesh::write_obj, py::arg("path"))
        .def("read_ply",  &Mesh::read_ply,  py::arg("path"))
        .def("write_ply", &Mesh::write_ply, py::arg("path"), py::arg("binary") = false,
             "Write PLY; binary=True for binary_little_endian, else ASCII.")
        .def("read_stl",  &Mesh::read_stl,  py::arg("path"),
             "Read STL (auto-detects ASCII vs binary); welds coincident vertices.")
        .def("write_stl", &Mesh::write_stl, py::arg("path"), py::arg("binary") = true,
             "Write STL; binary=True (default) or ASCII. Polygons are fan-triangulated.")
        .def("read_off",  &Mesh::read_off,  py::arg("path"))
        .def("write_off", &Mesh::write_off, py::arg("path"))

        // bulk normals as numpy (N,3) / (F,3), computed on the fly
        .def("vertex_normals_to_numpy", [](const Mesh& me) {
            py::array_t<double> arr({ (py::ssize_t)me.n_vertices(), (py::ssize_t)3 });
            auto w = arr.mutable_unchecked<2>();
            for (py::ssize_t i = 0; i < (py::ssize_t)me.n_vertices(); ++i) {
                Vec3 n = me.calc_vertex_normal(VertexHandle((int)i));
                w(i, 0) = n.x; w(i, 1) = n.y; w(i, 2) = n.z;
            }
            return arr;
        }, "Per-vertex unit normals as an (N,3) array.")
        .def("face_normals_to_numpy", [](const Mesh& me) {
            py::array_t<double> arr({ (py::ssize_t)me.n_faces(), (py::ssize_t)3 });
            auto w = arr.mutable_unchecked<2>();
            for (py::ssize_t i = 0; i < (py::ssize_t)me.n_faces(); ++i) {
                Vec3 n = me.calc_face_normal(FaceHandle((int)i));
                w(i, 0) = n.x; w(i, 1) = n.y; w(i, 2) = n.z;
            }
            return arr;
        }, "Per-face unit normals as an (F,3) array.")

        .def("__repr__", [](const Mesh& mesh) {
            return "<simplemesh.Mesh V=" + std::to_string(mesh.n_vertices()) +
                   " E=" + std::to_string(mesh.n_edges()) +
                   " F=" + std::to_string(mesh.n_faces()) + ">";
        });

    // ----- per-(element,type) property accessors --------------------------
    // Generate add_<elem>_<type>_property / get_<elem>_<type> / set_<elem>_<type>
    // for every element kind and value type. Python can't assign through a C++
    // reference, hence explicit get/set rather than the reference-returning
    // C++ property(). This mirrors the templated C++ API one-to-one.
#define BIND_PROP_METHODS(elem, ELEMCAP, ELEMH, T, TYNAME)                                       \
    mesh_cls.def("add_" #elem "_" TYNAME "_property",                                            \
        [](Mesh& me, const std::string& n, T d) { return me.add_##elem##_property<T>(n, d); },   \
        py::arg("name"), py::arg("default") = T());                                              \
    mesh_cls.def("get_" #elem "_" TYNAME,                                                        \
        [](Mesh& me, ELEMCAP##PropHandle<T> p, ELEMH h) { return me.property(p, h); });          \
    mesh_cls.def("set_" #elem "_" TYNAME,                                                        \
        [](Mesh& me, ELEMCAP##PropHandle<T> p, ELEMH h, T v) { me.property(p, h) = v; })

    BIND_PROP_METHODS(vertex,   Vertex,   VertexHandle,   double, "float");
    BIND_PROP_METHODS(vertex,   Vertex,   VertexHandle,   int,    "int");
    BIND_PROP_METHODS(vertex,   Vertex,   VertexHandle,   Vec3,   "vec3");
    BIND_PROP_METHODS(halfedge, Halfedge, HalfedgeHandle, double, "float");
    BIND_PROP_METHODS(halfedge, Halfedge, HalfedgeHandle, int,    "int");
    BIND_PROP_METHODS(halfedge, Halfedge, HalfedgeHandle, Vec3,   "vec3");
    BIND_PROP_METHODS(edge,     Edge,     EdgeHandle,     double, "float");
    BIND_PROP_METHODS(edge,     Edge,     EdgeHandle,     int,    "int");
    BIND_PROP_METHODS(edge,     Edge,     EdgeHandle,     Vec3,   "vec3");
    BIND_PROP_METHODS(face,     Face,     FaceHandle,     double, "float");
    BIND_PROP_METHODS(face,     Face,     FaceHandle,     int,    "int");
    BIND_PROP_METHODS(face,     Face,     FaceHandle,     Vec3,   "vec3");
#undef BIND_PROP_METHODS

    // ----- status bits ----------------------------------------------------
    // bit constants on the module: simplemesh.SELECTED, etc.
    m.attr("SELECTED") = (int)Mesh::SELECTED;
    m.attr("TAGGED")   = (int)Mesh::TAGGED;
    m.attr("LOCKED")   = (int)Mesh::LOCKED;
    m.attr("FEATURE")  = (int)Mesh::FEATURE;
    m.attr("HIDDEN")   = (int)Mesh::HIDDEN;

    // raw status word per element
#define BIND_STATUS_RAW(elem, ELEMH)                                                       \
    mesh_cls.def(#elem "_status", [](const Mesh& me, ELEMH h) { return me.get_status(h); }); \
    mesh_cls.def("set_" #elem "_status", [](Mesh& me, ELEMH h, int s) { me.set_status(h, s); })
    BIND_STATUS_RAW(vertex,   VertexHandle);
    BIND_STATUS_RAW(halfedge, HalfedgeHandle);
    BIND_STATUS_RAW(edge,     EdgeHandle);
    BIND_STATUS_RAW(face,     FaceHandle);
#undef BIND_STATUS_RAW

    // named flag accessors: is_<elem>_<flag> / set_<elem>_<flag>
#define BIND_STATUS(elem, ELEMH, flag)                                                       \
    mesh_cls.def("is_" #elem "_" #flag, [](const Mesh& me, ELEMH h) { return me.is_##flag(h); }); \
    mesh_cls.def("set_" #elem "_" #flag, [](Mesh& me, ELEMH h, bool v) { me.set_##flag(h, v); },  \
                 py::arg("handle"), py::arg("value") = true)
#define BIND_STATUS_ELEM(elem, ELEMH)        \
    BIND_STATUS(elem, ELEMH, selected);      \
    BIND_STATUS(elem, ELEMH, tagged);        \
    BIND_STATUS(elem, ELEMH, locked);        \
    BIND_STATUS(elem, ELEMH, feature);       \
    BIND_STATUS(elem, ELEMH, hidden)
    BIND_STATUS_ELEM(vertex,   VertexHandle);
    BIND_STATUS_ELEM(halfedge, HalfedgeHandle);
    BIND_STATUS_ELEM(edge,     EdgeHandle);
    BIND_STATUS_ELEM(face,     FaceHandle);
#undef BIND_STATUS_ELEM
#undef BIND_STATUS

    // ----- algorithms, as free functions ----------------------------------
    m.def("loop_subdivide", &loop_subdivide, py::arg("mesh"), py::arg("iterations") = 1,
          "Return a new mesh that is `mesh` Loop-subdivided `iterations` times.");
    m.def("quadric_decimate", &quadric_decimate, py::arg("mesh"), py::arg("target_faces"),
          "Simplify `mesh` in place to ~target_faces via QEM edge collapses "
          "(runs garbage_collection, invalidating handles).");
    m.def("laplacian_smooth", &laplacian_smooth, py::arg("mesh"),
          py::arg("iterations") = 1, py::arg("lambda_") = 0.5,
          "Uniform-Laplacian smoothing in place; boundary vertices stay fixed.");
    m.def("catmull_clark", &catmull_clark, py::arg("mesh"), py::arg("iterations") = 1,
          "Return a new mesh that is `mesh` Catmull-Clark subdivided `iterations` "
          "times (every face becomes a fan of quads).");
    m.def("sqrt3_subdivide", &sqrt3_subdivide, py::arg("mesh"), py::arg("iterations") = 1,
          "Return a new mesh that is `mesh` sqrt(3)-subdivided `iterations` times "
          "(centroid split + edge flips; triangle meshes).");
    m.def("midpoint_subdivide", &midpoint_subdivide, py::arg("mesh"), py::arg("iterations") = 1,
          "Interpolating 1->4 subdivision; new vertices at edge midpoints, old "
          "vertices unchanged.");
    m.def("butterfly_subdivide", &butterfly_subdivide, py::arg("mesh"), py::arg("iterations") = 1,
          "Interpolating, smooth 1->4 subdivision using the 8-point butterfly "
          "stencil (midpoint fallback on boundaries / irregular fans).");
    m.def("longest_edge_subdivide", &longest_edge_subdivide, py::arg("mesh"), py::arg("max_edge_length"),
          "In-place Rivara refinement: bisect every edge longer than "
          "max_edge_length until none exceed it (conformal).");
    m.def("cotan_smooth", &cotan_smooth, py::arg("mesh"),
          py::arg("iterations") = 1, py::arg("lambda_") = 0.5,
          "Cotangent-weighted Laplacian smoothing in place; boundary fixed.");
    m.def("fill_holes", &fill_holes, py::arg("mesh"), py::arg("max_edges") = 0,
          "Fan-triangulate boundary loops with at most max_edges edges (0 = no "
          "limit); returns the number of holes filled. Use max_edges so a large "
          "open boundary is not wrongly capped.");
    m.def("fill_hole", &fill_hole, py::arg("mesh"), py::arg("boundary_halfedge"),
          py::arg("max_edges") = 0,
          "Fill the single hole that the given boundary half-edge bounds.");
    m.def("remove_small_components", &remove_small_components,
          py::arg("mesh"), py::arg("min_faces"),
          "Delete connected components with fewer than min_faces faces (debris/"
          "islands). Runs garbage_collection; returns #components removed.");

    // ----- modular decimation (Decimater + modules) -----------------------
    // The pluggable framework: register a priority module (ModQuadric) and any
    // number of binary veto modules, then decimate.
    py::class_<CollapseInfo>(m, "CollapseInfo")
        .def_readonly("heh", &CollapseInfo::heh)
        .def_readonly("v0", &CollapseInfo::v0)
        .def_readonly("v1", &CollapseInfo::v1)
        .def_readonly("p0", &CollapseInfo::p0)
        .def_readonly("p1", &CollapseInfo::p1);

    py::class_<Decimator>(m, "Decimator")
        .def(py::init<Mesh&>(), py::arg("mesh"), py::keep_alive<1, 2>())
        .def("add_quadric_module",
             [](Decimator& d) { d.add_module<ModQuadric>(); },
             "Add the Garland-Heckbert quadric-error priority module.")
        .def("add_edge_length_module",
             [](Decimator& d, double max_len) { d.add_module<ModEdgeLength>(max_len); },
             py::arg("max_length"),
             "Veto collapses that would create an edge longer than max_length.")
        .def("add_normal_flipping_module",
             [](Decimator& d, double max_angle) { d.add_module<ModNormalFlipping>(max_angle); },
             py::arg("max_angle_rad"),
             "Veto collapses that flip a face normal by more than max_angle radians.")
        .def("add_aspect_ratio_module",
             [](Decimator& d, double max_aspect) { d.add_module<ModAspectRatio>(max_aspect); },
             py::arg("max_aspect"),
             "Veto collapses that create a triangle with aspect ratio > max_aspect.")
        .def("decimate_to_faces", &Decimator::decimate_to_faces, py::arg("target_faces"),
             "Greedily collapse to <= target_faces (runs garbage_collection, "
             "invalidating handles). Returns the number of collapses performed.");

    // ----- MeshChecker ----------------------------------------------------
    py::class_<CheckOptions>(m, "CheckOptions")
        .def(py::init<>())
        .def_readwrite("check_duplicate_vertices", &CheckOptions::check_duplicate_vertices)
        .def_readwrite("vertex_merge_tol",         &CheckOptions::vertex_merge_tol)
        .def_readwrite("check_duplicate_faces",    &CheckOptions::check_duplicate_faces)
        .def_readwrite("check_self_intersections", &CheckOptions::check_self_intersections);

    py::class_<MeshCheckReport>(m, "MeshCheckReport")
        .def_readonly("bad_halfedges",        &MeshCheckReport::bad_halfedges)
        .def_readonly("nonmanifold_vertices", &MeshCheckReport::nonmanifold_vertices)
        .def_readonly("nonmanifold_edges",    &MeshCheckReport::nonmanifold_edges)
        .def_readonly("isolated_vertices",    &MeshCheckReport::isolated_vertices)
        .def_readonly("boundary_edges",       &MeshCheckReport::boundary_edges)
        .def_readonly("degenerate_faces",     &MeshCheckReport::degenerate_faces)
        .def_readonly("duplicate_vertices",   &MeshCheckReport::duplicate_vertices)
        .def_readonly("duplicate_faces",      &MeshCheckReport::duplicate_faces)
        .def_readonly("self_intersections",   &MeshCheckReport::self_intersections)
        .def_readonly("n_components",         &MeshCheckReport::n_components)
        .def_readonly("n_boundary_loops",     &MeshCheckReport::n_boundary_loops)
        .def_readonly("euler",                &MeshCheckReport::euler)
        .def_readonly("genus",                &MeshCheckReport::genus)
        .def_readonly("is_closed",            &MeshCheckReport::is_closed)
        .def_readonly("is_manifold",          &MeshCheckReport::is_manifold)
        .def_readonly("is_oriented",          &MeshCheckReport::is_oriented)
        .def("ok",      &MeshCheckReport::ok)
        .def("summary", &MeshCheckReport::summary)
        .def("__repr__", [](const MeshCheckReport& r) { return r.summary(); });

    py::class_<MeshChecker>(m, "MeshChecker")
        .def(py::init<const Mesh&>(), py::arg("mesh"), py::keep_alive<1, 2>())
        .def("check", &MeshChecker::check, py::arg("options") = CheckOptions(),
             "Run the full validity battery and return a MeshCheckReport.")
        .def("is_valid", &MeshChecker::is_valid,
             "Quick pass/fail: connectivity + manifoldness + no degenerate faces.");

    // ----- MeshRepair -----------------------------------------------------
    py::class_<RepairOptions>(m, "RepairOptions")
        .def(py::init<>())
        .def_readwrite("weld_vertices",     &RepairOptions::weld_vertices)
        .def_readwrite("vertex_merge_tol",  &RepairOptions::vertex_merge_tol)
        .def_readwrite("reorient",          &RepairOptions::reorient)
        .def_readwrite("remove_degenerate", &RepairOptions::remove_degenerate)
        .def_readwrite("remove_duplicate",  &RepairOptions::remove_duplicate)
        .def_readwrite("remove_unused",     &RepairOptions::remove_unused)
        .def_readwrite("split_nonmanifold", &RepairOptions::split_nonmanifold)
        .def_readwrite("min_face_area",     &RepairOptions::min_face_area)
        .def_readwrite("min_component_faces", &RepairOptions::min_component_faces);

    py::class_<RepairReport>(m, "RepairReport")
        .def_readonly("vertices_merged",   &RepairReport::vertices_merged)
        .def_readonly("faces_flipped",     &RepairReport::faces_flipped)
        .def_readonly("degenerate_removed",&RepairReport::degenerate_removed)
        .def_readonly("duplicate_removed", &RepairReport::duplicate_removed)
        .def_readonly("vertices_removed",  &RepairReport::vertices_removed)
        .def_readonly("faces_failed",      &RepairReport::faces_failed)
        .def_readonly("nm_edges",          &RepairReport::nm_edges)
        .def_readonly("faces_split",       &RepairReport::faces_split)
        .def_readonly("vertices_split",    &RepairReport::vertices_split)
        .def_readonly("components_removed",&RepairReport::components_removed)
        .def("changed",  &RepairReport::changed)
        .def("summary",  &RepairReport::summary)
        .def("__repr__", [](const RepairReport& r) { return r.summary(); });

    m.def("repair_mesh", &repair_mesh, py::arg("mesh"), py::arg("options") = RepairOptions(),
          "Rebuild-based cleanup (weld colocated vertices, reorient, drop "
          "degenerate/duplicate faces and unused vertices). Invalidates handles; "
          "returns a RepairReport.");

    // ----- PreMesh repair layer -------------------------------------------
    py::class_<PreMesh>(m, "PreMesh")
        .def_static("from_mesh", &PreMesh::from_mesh, py::arg("mesh"))
        .def_property_readonly("n_vertices", &PreMesh::n_vertices)
        .def_property_readonly("n_faces",    &PreMesh::n_faces)
        .def("build_radial", &PreMesh::build_radial,
             "Build the radial map; returns the number of non-manifold edges.")
        .def("edge_valence", &PreMesh::edge_valence, py::arg("a"), py::arg("b"),
             "Number of faces on undirected edge (a,b) after build_radial().")
        .def("repair", &PreMesh::repair, py::arg("options") = RepairOptions(),
             "Run the repair passes in place; returns a RepairReport.")
        .def("to_mesh", [](const PreMesh& p, const RepairOptions& opt) {
                 RepairReport r;
                 Mesh m = p.to_mesh(opt, &r);
                 return py::make_tuple(std::move(m), r);
             }, py::arg("options") = RepairOptions(),
             "Build the halfedge Mesh. Returns (mesh, RepairReport) where the "
             "report carries faces_split / faces_failed.");

    m.def("load_soup", &load_soup, py::arg("path"),
          "Parse a file (STL/OBJ/OFF) straight to a PreMesh soup (no halfedge "
          "build, so no faces are dropped).");
    m.def("load_and_repair", &load_and_repair, py::arg("path"), py::arg("options") = RepairOptions(),
          "Load a file, repair it through the PreMesh layer, and return a clean "
          "halfedge Mesh (load_soup -> repair -> to_mesh).");
}
