// ============================================================================
//  SimpleMesh - Items.h : the internal records stored in the mesh arrays
// ----------------------------------------------------------------------------
//  These are the actual structs the mesh keeps in std::vectors. The user never
//  touches them directly - they only ever see handles (see Handles.h). Keeping
//  them in their own header makes the storage layout of the halfedge data
//  structure explicit and easy to study.
//
//  The halfedge data structure represents each undirected edge as TWO directed
//  "halfedges" pointing in opposite directions. Almost all connectivity is
//  encoded on the halfedges:
//
//        v_from  --------- next ---------->  v_to
//               <--------- (opposite) ------
// ============================================================================
#pragma once

#include "Handles.h"

namespace sm {

/// A vertex stores a reference to ONE outgoing halfedge. From that single link
/// we can reach every neighbour by "circulating" around the vertex.
/// (The 3D position lives in a parallel array `points_`, not here, so geometry
///  and connectivity stay cleanly separated.)
struct Vertex {
    // Default-initialized to invalid. The explicit `= HalfedgeHandle()` is
    // needed because the handle's default constructor is `explicit`, so plain
    // aggregate `{}` initialization of members would be ill-formed under MSVC.
    HalfedgeHandle outgoing_halfedge = HalfedgeHandle();
};

/// A directed halfedge - the heart of the data structure.
struct Halfedge {
    VertexHandle   to_vertex = VertexHandle();   // the vertex this halfedge points TO
    FaceHandle     face      = FaceHandle();      // incident face, or invalid if boundary halfedge
    HalfedgeHandle next      = HalfedgeHandle();  // next halfedge around the incident face (CCW)
    HalfedgeHandle prev      = HalfedgeHandle();  // previous halfedge around the incident face
    // NOTE: the OPPOSITE halfedge is not stored. We use the classic trick of
    // pairing halfedges 2e and 2e+1 for edge e, so opposite(h) == (h XOR 1).
};

/// A face stores a reference to ONE of its boundary halfedges. Walking `next`
/// from there visits the whole face.
struct Face {
    HalfedgeHandle halfedge = HalfedgeHandle();  // one halfedge belonging to this face
};

} // namespace sm
