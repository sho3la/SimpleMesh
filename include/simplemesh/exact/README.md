# `simplemesh/exact/` — exact-arithmetic kernel

This subtree is the **isolated**, exact-arithmetic foundation used to resolve
triangle/triangle self-intersections into an exact surface arrangement and to run
boolean operations on it. It is kept separate from SimpleMesh's deliberately
simple, `double`-precision public core so that the teaching API stays clean.

The layers build **bottom-up**, each with its own self-checking test:

| Header                   | Role                                                        |
|--------------------------|-------------------------------------------------------------|
| `Expansion.h`            | arbitrary-precision sums of doubles (exact arithmetic)      |
| `Interval.h`             | round-to-nearest interval arithmetic (predicate filter)     |
| `Predicates.h`           | exact orient2d / orient3d / in_circle (filtered)            |
| `HomogeneousGeometry.h`  | homogeneous points + exact projective predicates            |
| `CDT2d.h`                | constrained Delaunay triangulation (combinatorics)          |
| `TriangleIntersection.h` | symbolic triangle/triangle intersection                     |
| `MeshInTriangle.h`       | exact-coordinate CDT + per-facet remesh                     |
| `SurfaceIntersection.h`  | the arrangement assembly driver (soup in / soup out)        |
| `MeshBoolean.h`          | Weiler model + radial sort + boolean classification         |

**Invariant:** no tolerance constants in this subtree. Exact means exact.
