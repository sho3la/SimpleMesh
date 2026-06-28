"""05 - Topological editing: flip, split, collapse, vertex_split, delete + GC.

Each operator is guarded by an `is_*_ok` predicate; deletion is lazy until
garbage_collection compacts the arrays.
Run:  PYTHONPATH=build/bin python examples/05_editing.py
"""
import simplemesh as sm


def quad():
    m = sm.Mesh()
    a = m.add_vertex(sm.Vec3(0, 0, 0)); b = m.add_vertex(sm.Vec3(1, 0, 0))
    c = m.add_vertex(sm.Vec3(1, 1, 0)); d = m.add_vertex(sm.Vec3(0, 1, 0))
    m.add_triangle(a, b, c); m.add_triangle(a, c, d)
    return m, (a, b, c, d)

# --- edge flip ---
m, (a, b, c, d) = quad()
e = m.edge(m.find_halfedge(a, c))
print("flip ok?", m.is_flip_ok(e)); m.flip(e)
print("  after flip: edge a-c gone, b-d present:",
      not m.find_halfedge(a, c).is_valid(), m.find_halfedge(b, d).is_valid())

# --- edge split ---
m, (a, b, c, d) = quad()
mid = m.add_vertex(sm.Vec3(0.5, 0.5, 0))     # will sit on diagonal a-c
m.split_edge(m.edge(m.find_halfedge(a, c)), mid)
print("split: faces", m.n_faces, "(was 2)")

# --- collapse ---
m, (a, b, c, d) = quad()
h = m.find_halfedge(a, c)
print("collapse ok?", m.is_collapse_ok(h))
if m.is_collapse_ok(h):
    m.collapse(h); m.garbage_collection()
print("  after collapse+GC: V", m.n_vertices, "F", m.n_faces)

# --- delete + garbage collection ---
m, _ = quad()
m.delete_face(sm.FaceHandle(0))
print("deleted face flagged:", m.is_deleted_face(sm.FaceHandle(0)))
m.garbage_collection()
print("  after GC: F", m.n_faces)
