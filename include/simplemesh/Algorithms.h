// ============================================================================
//  SimpleMesh - Algorithms.h : mesh-processing algorithms
// ----------------------------------------------------------------------------
//  Built entirely on top of the core mesh primitives, these show that the data
//  structure is good for more than bookkeeping:
//
//    * loop_subdivide    - Loop subdivision surface refinement (1 -> 4 triangles
//                          per face, with the classic smoothing weights). Shows
//                          off construction + navigation.
//    * quadric_decimate  - Garland-Heckbert QEM mesh simplification by greedy
//                          edge collapse. Shows off the editing operators
//                          (collapse / is_collapse_ok), custom data per element,
//                          and garbage_collection.
//
//  Both are free functions in namespace sm so they stay decoupled from the Mesh
//  class - exactly how you would layer algorithms on a real library.
// ============================================================================
#pragma once

#include "Mesh.h"

namespace sm {

// ----------------------------------------------------------------------------
//  A 4x4 symmetric "error quadric" (Garland & Heckbert 1997). It encodes the
//  sum of squared distances from a point to a set of planes. Stored as the 10
//  unique entries of the symmetric matrix.
// ----------------------------------------------------------------------------
struct Quadric {
    double a = 0, b = 0, c = 0, d = 0,   //  [ a b c d ]
                  e = 0, f = 0, g = 0,   //  [ b e f g ]
                         h = 0, i = 0,   //  [ c f h i ]
                                j = 0;   //  [ d g i j ]

    Quadric operator+(const Quadric& o) const {
        return { a + o.a, b + o.b, c + o.c, d + o.d, e + o.e,
                 f + o.f, g + o.g, h + o.h, i + o.i, j + o.j };
    }

    /// Build the quadric of a single plane (nx,ny,nz,pd) with unit normal:
    /// K = p p^T where p = (nx,ny,nz,pd).
    static Quadric from_plane(double nx, double ny, double nz, double pd) {
        return { nx*nx, nx*ny, nx*nz, nx*pd,
                        ny*ny, ny*nz, ny*pd,
                               nz*nz, nz*pd,
                                      pd*pd };
    }

    /// Evaluate v^T Q v for v = (p, 1): the squared-distance error at point p.
    double eval(const Vec3& p) const {
        return a*p.x*p.x + 2*b*p.x*p.y + 2*c*p.x*p.z + 2*d*p.x
             +              e*p.y*p.y   + 2*f*p.y*p.z + 2*g*p.y
             +                            h*p.z*p.z   + 2*i*p.z
             +                                          j;
    }
};

/// Loop-subdivide a triangle mesh `in`, `iterations` times, returning the new
/// (finer, smoother) mesh. The input should be a compact triangle mesh
/// (call garbage_collection first if it has deleted elements).
Mesh loop_subdivide(const Mesh& in, int iterations = 1);

/// Simplify a triangle mesh in place down to (about) `target_faces` faces using
/// greedy quadric-error edge collapses. Runs garbage_collection at the end, so
/// all handles into `m` are invalidated.
void quadric_decimate(Mesh& m, size_t target_faces);

/// Catmull-Clark subdivision of an arbitrary polygon mesh, `iterations` times.
/// Every face becomes a fan of quads (so after one step the mesh is all-quads),
/// using the classic face-point / edge-point / vertex-point masks. A primal
/// quad scheme - a nice contrast to Loop (which is primal on triangles).
/// Boundaries follow the cubic B-spline rule (so they stay smooth curves).
Mesh catmull_clark(const Mesh& in, int iterations = 1);

/// Sqrt(3) subdivision (Kobbelt) of a triangle mesh, `iterations` times. Each
/// face gets a new centroid vertex (1 -> 3 triangles), every original interior
/// edge is flipped, and old interior vertices are relaxed with Kobbelt's mask.
/// Slower density growth than Loop (factor 3 per step, not 4) and the new
/// vertices interpolate face centroids. Boundary vertices are kept fixed.
Mesh sqrt3_subdivide(const Mesh& in, int iterations = 1);

/// Uniform-Laplacian (umbrella) smoothing: move each interior vertex a fraction
/// `lambda` of the way toward the average of its 1-ring neighbours, `iterations`
/// times. Boundary vertices are held fixed. Topology is unchanged.
void laplacian_smooth(Mesh& m, int iterations = 1, double lambda = 0.5);

/// Cotangent-weighted Laplacian smoothing - the discretization that respects
/// the surface metric (used for fairing / curvature flow). Each interior vertex
/// moves toward the cotan-weighted average of its neighbours; weights come from
/// the angles opposite each edge. Boundary vertices stay fixed. A Jacobi
/// (double-buffered) iteration, `iterations` times, step size `lambda`.
void cotan_smooth(Mesh& m, int iterations = 1, double lambda = 0.5);

// ----------------------------------------------------------------------------
//  Hole filling - close boundary loops with triangles. Returns how many holes
//  (boundary loops) were filled. `fill_hole` fills the single loop that the
//  given boundary half-edge belongs to (fan triangulation).
//
//  `max_edges` caps the boundary length of a loop that will be filled (0 = no
//  limit). Use it so only small holes are closed and the legitimate open
//  boundary of a surface is NOT "capped".
// ----------------------------------------------------------------------------
int  fill_holes(Mesh& m, int max_edges = 0);
bool fill_hole(Mesh& m, HalfedgeHandle boundary_h, int max_edges = 0);

/// Remove connected components with fewer than `min_faces` faces (stray islands
/// / scan or print debris). Operates in place (runs garbage_collection, so all
/// handles are invalidated). Returns the number of components removed.
int  remove_small_components(Mesh& m, size_t min_faces);

// ----------------------------------------------------------------------------
//  More subdivision schemes (joining loop / catmull_clark / sqrt3).
// ----------------------------------------------------------------------------

/// Midpoint (a.k.a. polyhedral) subdivision: 1 -> 4 triangles, new vertices at
/// edge midpoints, old vertices UNCHANGED. The simplest *interpolating* scheme -
/// a useful contrast to Loop (same connectivity, no smoothing).
Mesh midpoint_subdivide(const Mesh& in, int iterations = 1);

/// Butterfly subdivision: 1 -> 4 triangles, *interpolating* (old vertices kept)
/// but smooth, using the 8-point butterfly stencil on interior edges (falling
/// back to the midpoint on boundaries / irregular fans). The interpolating
/// counterpart to Loop.
Mesh butterfly_subdivide(const Mesh& in, int iterations = 1);

/// Longest-edge (Rivara) refinement, in place: repeatedly bisect every edge
/// longer than `max_edge_length` (splitting both incident triangles, so the
/// mesh stays conformal) until no edge exceeds the bound. Adaptive density -
/// the contrast to the uniform schemes. Topology grows; handles stay valid
/// (no garbage_collection is performed).
void longest_edge_subdivide(Mesh& m, double max_edge_length);

} // namespace sm
