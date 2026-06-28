"""04 - Smart handles: fluent navigation.

A smart handle bundles an element with its mesh so you can chain navigation
instead of threading the mesh through every call.
Run:  PYTHONPATH=build/bin python examples/04_smart_handles.py
"""
import simplemesh as sm

m = sm.Mesh()
v0 = m.add_vertex(sm.Vec3(0, 0, 0))
v1 = m.add_vertex(sm.Vec3(1, 0, 0))
v2 = m.add_vertex(sm.Vec3(1, 1, 0))
v3 = m.add_vertex(sm.Vec3(0, 1, 0))
f0 = m.add_triangle(v0, v1, v2)
m.add_triangle(v0, v2, v3)

h = m.smart_halfedge(m.face_halfedge(f0))
# chain: walk to next, flip to the opposite half-edge, read its target vertex
print("h.next().opp().to() =", h.next().opp().to().idx())
print("opp().opp() == h    :", h.opp().opp().idx() == h.idx())

sv = m.smart_vertex(v0)
print("v0 1-ring   :", [w.idx() for w in sv.vertices()])
print("v0 valence  :", sv.valence(), " boundary:", sv.is_boundary())

sf = m.smart_face(f0)
print("face area   :", sf.area())
print("face normal :", sf.normal())
print("face verts  :", [w.idx() for w in sf.vertices()])
