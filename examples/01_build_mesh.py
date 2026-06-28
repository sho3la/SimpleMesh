"""01 - Build a mesh by hand.

Create vertices and faces, then read back the basic counts and connectivity.
Run:  PYTHONPATH=build/bin python examples/01_build_mesh.py
"""
import simplemesh as sm

m = sm.Mesh()

# a tetrahedron: 4 vertices, 4 triangular faces
v = [m.add_vertex(sm.Vec3(*p)) for p in
     [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]]
m.add_triangle(v[0], v[2], v[1])
m.add_triangle(v[0], v[1], v[3])
m.add_triangle(v[1], v[2], v[3])
m.add_triangle(v[2], v[0], v[3])

# add_face also accepts polygons of any size:
#   m.add_face([v0, v1, v2, v3])   # a quad

print(m)                                   # <simplemesh.Mesh V=4 E=6 F=4>
print("vertices :", m.n_vertices)
print("edges    :", m.n_edges)
print("faces    :", m.n_faces)
print("closed   :", not m.is_boundary_vertex(v[0]))   # tetra is watertight

f0 = sm.FaceHandle(0)
print("face 0 vertices:", [h.idx() for h in m.face_vertices(f0)])
