// ============================================================================
//  SimpleMesh - PreMesh.h : the pre-processing repair layer
// ----------------------------------------------------------------------------
//  A staging area BETWEEN raw input files and the manifold halfedge `Mesh`. Its
//  job: take arbitrary "triangle soup", repair it into halfedge-legal topology,
//  then emit the halfedge mesh - so the kernel never has to drop faces it can't
//  represent.
//
//  PreMesh is an indexed face set plus an on-demand RADIAL map (edge -> all the
//  face-corners on it), a static form of a radial cycle. The
//  radial count classifies every edge:
//      1  -> boundary,  2 -> manifold,  >=3 -> non-manifold.
//  Repairs (weld / drop degenerate+duplicate / reorient / split non-manifold)
//  run here, where non-manifold topology is representable; `to_mesh()` then
//  builds the halfedge mesh.
// ============================================================================
#pragma once

#include "Mesh.h"
#include "MeshRepair.h"   // RepairOptions, RepairReport

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace sm {

class PreMesh {
public:
    // ---- the soup ----------------------------------------------------------
    std::vector<Vec3>             positions;   ///< vertex coordinates
    std::vector<std::vector<int>> faces;       ///< polygons as vertex indices (order = winding)

    // ---- construction ------------------------------------------------------
    static PreMesh from_soup(std::vector<Vec3> pos, std::vector<std::vector<int>> f);
    static PreMesh from_mesh(const Mesh& m);   ///< extract a live mesh back to soup

    // ---- repair + emit -----------------------------------------------------
    /// Run the repair passes in place; returns what changed.
    RepairReport repair(const RepairOptions& opt = RepairOptions());
    /// Build the halfedge mesh. With opt.split_nonmanifold, faces the kernel
    /// rejects are re-emitted on duplicated vertices instead of dropped.
    Mesh to_mesh(const RepairOptions& opt, RepairReport* report = nullptr) const;
    Mesh to_mesh() const { return to_mesh(RepairOptions()); }

    size_t n_vertices() const { return positions.size(); }
    size_t n_faces()    const { return faces.size(); }

    // ---- radial structure (built on demand) --------------------------------
    /// A face-corner: the directed edge faces[face][slot] -> next-vertex.
    struct Corner { int face; int slot; };
    /// (Re)build the radial map (edge -> corners) and vertex -> corners. Returns
    /// the number of non-manifold edges (radial >= 3).
    int build_radial();
    /// Number of faces on the undirected edge (a,b) after build_radial().
    int edge_valence(int a, int b) const;

private:
    static uint64_t ekey(int a, int b);   // packed undirected-edge key (min,max)
    std::unordered_map<uint64_t, std::vector<Corner>> radial_;  // flat radial cycle
    bool radial_built_ = false;

    int  weld_(double tol);                 // merge colocated verts (neighbor-aware spatial hash)
    int  drop_degenerate_(double min_area); // drop zero-area / sliver faces
    int  drop_duplicate_();                 // drop repeated faces
    int  reorient_();                       // consistent winding (uses radial_)
    int  split_fans_();                     // split non-manifold edges & bow-ties into fans
    int  drop_small_components_(int min_faces);  // strip tiny islands
    int  drop_unused_();                    // drop unreferenced verts
};

// ----------------------------------------------------------------------------
//  Loaders + one-call convenience. load_soup parses a file straight to
//  soup (no halfedge build, so nothing is dropped). Supported: STL / OBJ / OFF.
// ----------------------------------------------------------------------------
PreMesh load_soup(const std::string& path);
Mesh    load_and_repair(const std::string& path, const RepairOptions& opt = RepairOptions());

} // namespace sm
