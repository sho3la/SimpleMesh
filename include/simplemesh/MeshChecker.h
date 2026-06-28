// ============================================================================
//  SimpleMesh - MeshChecker.h : an advanced mesh validity / sanity checker
// ----------------------------------------------------------------------------
//  A half-edge mesh has many invariants that, if violated, make every later
//  algorithm misbehave (often silently). This checker verifies them and reports
//  exactly which elements are wrong - the role a validity / repair pass plays in
//  production mesh libraries.
//
//  What it checks:
//
//   A. Half-edge CONNECTIVITY integrity (the structural invariants):
//        opposite involution, next/prev inverse, face constant around a loop,
//        to_vertex(h) == from_vertex(next(h)), boundary loops stay on boundary.
//      (the "pair up all the half-edges" idea, generalized.)
//
//   B. TOPOLOGY / manifoldness:
//        non-manifold vertices (more than one fan = "bow-tie"),
//        non-manifold edges (used by != 2 faces),
//        isolated vertices, boundary edges, connected components,
//        Euler characteristic & estimated genus, closed/watertight.
//
//   C. DEGENERACY & DUPLICATION (common repair targets):
//        faces with a repeated vertex or < 3 vertices or ~zero area,
//        colocated (duplicate) vertices within a tolerance,
//        duplicate faces (same vertex set).
//
//   D. GEOMETRY (optional, O(F^2)):
//        triangle/triangle self-intersections (Moller's test).
//
//  Nothing here mutates the mesh - it is a pure inspector.
// ============================================================================
#pragma once

#include "Mesh.h"

#include <string>
#include <utility>
#include <vector>

namespace sm {

/// Tunables for the more expensive / fuzzy checks.
struct CheckOptions {
    bool   check_duplicate_vertices  = true;
    double vertex_merge_tol          = 1e-10;  ///< colocated if closer than this
    bool   check_duplicate_faces     = true;
    bool   check_self_intersections  = false;  ///< opt-in (see caps below)
    /// Self-intersection uses a sweep-and-prune broad phase, but is still worst-
    /// case quadratic. These caps keep it bounded in time AND memory:
    ///   * skip the whole check above this face count (0 = no limit),
    ///   * stop after collecting this many intersecting pairs (0 = no limit).
    size_t self_intersection_max_faces = 500000;
    size_t self_intersection_max_hits  = 100000;
};

/// A structured report. Empty vectors / true flags mean "all good".
struct MeshCheckReport {
    // A. connectivity
    std::vector<int> bad_halfedges;          ///< halfedge indices failing an invariant
    // B. topology
    std::vector<int> nonmanifold_vertices;
    std::vector<int> nonmanifold_edges;
    std::vector<int> isolated_vertices;
    std::vector<int> boundary_edges;
    int  n_components = 0;
    int  n_boundary_loops = 0;
    int  euler = 0;                          ///< V - E + F (live elements)
    int  genus = 0;                          ///< (2*C - euler - B)/2, if manifold+closed-ish
    bool is_closed   = false;                ///< no boundary edges
    bool is_manifold = false;                ///< connectivity + no non-manifold elems
    bool is_oriented = false;                ///< every edge consistently oriented
    // C. degeneracy / duplication
    std::vector<int> degenerate_faces;
    std::vector<std::pair<int,int>> duplicate_vertices;  ///< colocated vertex pairs
    std::vector<std::pair<int,int>> duplicate_faces;     ///< same-vertex-set face pairs
    // D. geometry
    std::vector<std::pair<int,int>> self_intersections;  ///< intersecting face pairs

    /// True iff the mesh passed every check that was run.
    bool ok() const;
    /// Human-readable multi-line summary (counts + the first few offenders).
    std::string summary() const;
};

/// The checker. Construct around a mesh, then call check().
class MeshChecker {
public:
    explicit MeshChecker(const Mesh& m) : mesh_(m) {}

    /// Run the full battery and return the report.
    MeshCheckReport check(const CheckOptions& opt = CheckOptions()) const;

    /// Quick pass/fail: connectivity + manifoldness + no degenerate faces.
    /// (Skips the fuzzy duplicate / geometric checks.)
    bool is_valid() const;

private:
    const Mesh& mesh_;
};

} // namespace sm
