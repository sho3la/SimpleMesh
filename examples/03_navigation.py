"""03 - Navigating the half-edge structure.

The eight primitives (next/prev/opposite, to/from vertex, face, edge) plus the
circulators and lazy iterators that traverse a neighbourhood.
Run:  PYTHONPATH=build/bin python examples/03_navigation.py
"""
import simplemesh as sm

# two triangles sharing the diagonal v0-v2 (a split quad)
m = sm.Mesh()
v0 = m.add_vertex(sm.Vec3(0, 0, 0))
v1 = m.add_vertex(sm.Vec3(1, 0, 0))
v2 = m.add_vertex(sm.Vec3(1, 1, 0))
v3 = m.add_vertex(sm.Vec3(0, 1, 0))
m.add_triangle(v0, v1, v2)
m.add_triangle(v0, v2, v3)

h = m.find_halfedge(v0, v2)                 # the shared diagonal
print("halfedge v0->v2:", h.idx())
print("  to_vertex      :", m.to_vertex(h).idx())
print("  from_vertex    :", m.from_vertex(h).idx())
print("  opposite       :", m.opposite_halfedge(h).idx())
print("  next / prev    :", m.next_halfedge(h).idx(), "/", m.prev_halfedge(h).idx())
print("  interior edge? :", not m.is_boundary_edge(m.edge(h)))

print("1-ring of v0     :", [w.idx() for w in m.vertex_vertices(v0)])
print("faces around v0  :", [f.idx() for f in m.vertex_faces(v0)])
print("valence of v0    :", m.valence(v0))

# lazy iterators (no allocation) work in for-loops:
print("all vertices     :", [w.idx() for w in m.vertices()])
print("neighbours of v0 :", [w.idx() for w in m.vv(v0)])
