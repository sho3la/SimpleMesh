# SimpleMesh

**A compact half-edge mesh library for loading, repairing, editing, and processing
triangle & polygon meshes — with a clean C++ core and a one-to-one Python API.**

SimpleMesh takes messy real-world meshes (STL/OBJ/PLY/OFF — often non-manifold,
mis-wound, or full of duplicate/degenerate faces) and turns them into a clean,
manifold half-edge structure you can navigate, edit, measure, simplify, subdivide,
smooth, validate, and save. The C++ is small and readable; the Python module mirrors
it exactly.

```python
import simplemesh as sm

mesh = sm.load_and_repair("bunny.stl")     # load + auto-repair to a clean manifold
print(sm.MeshChecker(mesh).check())        # full validity report
sm.quadric_decimate(mesh, target_faces=5000)
mesh.write_obj("bunny_clean.obj")
```

---

## Why SimpleMesh

- **It ingests broken meshes.** A dedicated repair layer fixes topology *before* it
  reaches the manifold kernel — welding coincident vertices, fixing winding, dropping
  degenerate/duplicate faces, and splitting non-manifold edges & bow-tie vertices —
  so nothing is silently dropped.
- **It tells you what's wrong.** A thorough checker reports manifoldness, Euler
  characteristic & genus, components, boundaries, duplicates, and self-intersections.
- **It's a complete processing toolkit** — subdivision, decimation, smoothing, hole
  filling, and the core edit operators, all on one data structure.
- **It's Pythonic and fast where it counts** — zero-copy NumPy views and bulk
  transfer for vertices, faces, properties, and normals.
- **It's dependency-light** — header + source C++17, pybind11 for Python, NumPy
  optional at runtime. `pip install .` or drop the headers into your CMake project.

## Features

### Data structure & navigation
- Half-edge connectivity for arbitrary **polygon** meshes (triangles, quads, n-gons).
- Topologically-correct `add_vertex` / `add_face` / `add_triangle`.
- Full navigation: `next` / `prev` / `opposite` half-edges, `to`/`from` vertex,
  incident face/edge, boundary tests.
- Circulators (vertex 1-ring, incident faces/edges/half-edges; face vertices/edges/
  neighbours) as both eager vectors and **lazy zero-allocation iterators**.
- **Smart handles** for fluent chaining: `mesh.smart(h).next().opp().to()`.

### Editing
- `flip`, `split` (edge & face), `collapse`, and **`vertex_split`** (the inverse of
  collapse) — each guarded by an `is_*_ok` predicate.
- `triangulate()`, lazy **deletion** + **garbage collection** with handle remapping.

### Repair & validation
- **`MeshChecker`** — connectivity integrity, non-manifold vertices/edges, isolated/
  degenerate elements, duplicate vertices/faces, connected components, Euler/genus,
  watertightness, and optional triangle/triangle self-intersections.
- **`repair_mesh`** — weld, coherently reorient, drop degenerate/duplicate faces &
  unused vertices.
- **`PreMesh` / `load_and_repair`** — a pre-kernel staging layer (a static
  *radial-edge* structure) that makes raw "triangle soup" manifold-legal: it welds
  colocated vertices (neighbour-aware spatial hash), reorients, drops degenerate/
  duplicate faces, splits non-manifold edges, dissociates bow-tie vertices, and
  optionally removes **slivers** (`min_face_area`) and **tiny islands**
  (`min_component_faces`) — so the half-edge kernel never has to reject a face.
- **`remove_small_components`** — delete stray islands / scan-print debris.
- **`fill_holes`** — close boundary loops by triangulation, with a `max_edges`
  guard so a genuine open boundary isn't wrongly capped.

### Geometry & properties
- Face/vertex normals, areas, centroids, edge lengths/vectors/midpoints, dihedral &
  sector angles; whole-mesh surface area, center of mass, bounding box.
- Runtime **typed properties** (`float`/`int`/`Vec3`) on any element kind, plus
  per-element **status bits** (selected/tagged/locked/feature/hidden); both grow with
  the mesh and survive garbage collection.

### Algorithms
- **Subdivision:** Loop, Catmull–Clark, √3, Midpoint, Butterfly, and adaptive
  Longest-edge (Rivara) refinement.
- **Decimation:** Garland–Heckbert QEM, plus a modular `Decimator` with pluggable
  priority (quadric) and binary veto modules (edge length, normal flipping, aspect
  ratio).
- **Smoothing:** uniform Laplacian and cotangent-weighted (Laplace–Beltrami).

### I/O & NumPy
- Read/write **OBJ, PLY (ascii + binary), STL (ascii + binary), OFF**.
- Zero-copy `points_view`, bulk `add_vertices` / `get_vertices` / `add_triangles` /
  `get_triangles`, per-property NumPy transfer, and bulk normals.

## Install

**Python package:**
```bash
pip install .          # builds the C++ extension via scikit-build-core + CMake
python -c "import simplemesh as sm; print(sm.Mesh())"
```

**C++ library (standalone):**
```bash
cmake -S . -B build            # add -G Ninja if you like
cmake --build build --config Release
./build/bin/simplemesh_tests   # run the test suite
```
On Windows, configure from a Developer prompt (VS 2022/2026 C++ tools) and run
`build\bin\simplemesh_tests.exe`.

Requirements: a C++17 compiler and CMake ≥ 3.18. pybind11 is fetched automatically
(or supplied by `pip`). NumPy is optional and only needed for the array helpers.

## Quick start

```python
import simplemesh as sm

# build a mesh by hand
m = sm.Mesh()
a = m.add_vertex(sm.Vec3(0, 0, 0))
b = m.add_vertex(sm.Vec3(1, 0, 0))
c = m.add_vertex(sm.Vec3(0, 1, 0))
f = m.add_triangle(a, b, c)
print(m.calc_face_normal(f))               # Vec3(0, 0, 1)

# or load + repair an arbitrary file, then process it
mesh = sm.load_and_repair("model.stl")
report = sm.MeshChecker(mesh).check()
print("manifold:", report.is_manifold, "components:", report.n_components)

fine = sm.loop_subdivide(mesh, iterations=2)
fine.write_ply("subdivided.ply")
```

## Examples

Each topic in [`examples/`](examples/) ships as **both** a Python (`.py`) and a C++
(`.cpp`) file, demonstrating **one** capability so you can learn it in isolation:

| # | Shows |
|---|-------|
| `01_build_mesh`        | building meshes (`add_vertex`/`add_face`), counts |
| `02_io`                | read/write OBJ · PLY · STL · OFF round-trips |
| `03_navigation`        | half-edge navigation, circulators, lazy iterators |
| `04_smart_handles`     | fluent `smart(...)` chaining |
| `05_editing`           | flip · split · collapse · delete + GC |
| `06_properties`        | typed per-element properties & status bits |
| `07_numpy` / `07_buffers` | NumPy zero-copy views (Py) · contiguous buffer access (C++) |
| `08_geometry`          | normals, areas, angles, bounding box |
| `09_subdivision`       | Loop · Catmull–Clark · √3 · midpoint · butterfly · longest-edge |
| `10_decimation`        | QEM + the modular `Decimator` |
| `11_smoothing`         | uniform & cotangent Laplacian |
| `12_check_and_repair`  | `MeshChecker`, `PreMesh`/`load_and_repair`, `fill_holes` |

Run the Python version (after building the module):
```bash
PYTHONPATH=build/bin python examples/09_subdivision.py
```
The C++ versions are built as standalone executables (one per file):
```bash
cmake --build build
./build/bin/09_subdivision        # build\bin\09_subdivision.exe on Windows
```

## License

BSD-3-Clause — see [LICENSE](LICENSE).
