// ============================================================================
//  SimpleMesh - MeshRepair.h : fix the defects that MeshChecker reports
// ----------------------------------------------------------------------------
//  MeshChecker *finds* problems; MeshRepair *fixes* them. It is the cleanup pass
//  every importer runs on raw "triangle soup" before handing it to real
//  algorithms.
//
//  Strategy: REBUILD. We extract the raw (positions + index-faces) soup, clean
//  it on those flat arrays - where merging, dropping and reorienting are trivial
//  and can never leave the half-edge structure half-stitched - then rebuild the
//  mesh through the trusted add_vertex / add_face path. This mirrors how
//  loop_subdivide and triangulate already work.
//
//  The repair steps, in order (each optional):
//    1. weld colocated vertices
//    2. reorient facets coherently
//    3. drop degenerate facets          (repeated vertex / < 3 verts)
//    4. drop duplicate facets
//    5. drop unused (isolated) vertices
// ============================================================================
#pragma once

#include "Mesh.h"

#include <string>

namespace sm {

/// Which repairs to run, and the colocation tolerance.
struct RepairOptions {
    bool   weld_vertices      = true;
    double vertex_merge_tol   = 1e-10;
    bool   reorient           = true;   ///< make face winding globally coherent
    bool   remove_degenerate  = true;
    bool   remove_duplicate   = true;
    bool   remove_unused      = true;   ///< drop vertices referenced by no face
    /// PreMesh emit only: when a face can't be added to the manifold halfedge
    /// kernel (non-manifold edge / bow-tie), duplicate its vertices and emit it
    /// as a separate sheet instead of dropping it. Preserves all geometry.
    bool   split_nonmanifold  = true;
    /// Drop faces whose area is below this (0 = keep all). Removes slivers/caps
    /// that survive the exact-degeneracy test. PreMesh path only.
    double min_face_area      = 0.0;
    /// Drop connected components with fewer than this many faces (0 = keep all).
    /// Strips scan/print debris and stray islands. PreMesh path only.
    int    min_component_faces = 0;
};

/// What the repair actually did (all zero = the mesh was already clean).
struct RepairReport {
    int vertices_merged   = 0;   ///< colocated vertices folded away
    int faces_flipped     = 0;   ///< faces reoriented for coherence
    int degenerate_removed = 0;
    int duplicate_removed  = 0;
    int vertices_removed  = 0;   ///< unused vertices dropped
    int faces_failed      = 0;   ///< faces add_face still rejected (e.g. non-manifold)
    int nm_edges          = 0;   ///< non-manifold edges detected (radial >= 3)
    int faces_split       = 0;   ///< faces emitted as separate sheets (emit fallback)
    int vertices_split    = 0;   ///< new vertex copies from fan separation (PreMesh)
    int components_removed = 0;  ///< tiny islands dropped (min_component_faces)

    bool changed() const {
        return vertices_merged || faces_flipped || degenerate_removed ||
               duplicate_removed || vertices_removed || faces_failed ||
               faces_split || vertices_split || components_removed;
    }
    std::string summary() const;
};

/// Repair `m` in place, returning a report of what changed. All handles into `m`
/// are invalidated (the mesh is rebuilt).
RepairReport repair_mesh(Mesh& m, const RepairOptions& opt = RepairOptions());

} // namespace sm
