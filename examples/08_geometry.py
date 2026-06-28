"""08 - Geometry queries.

Per-element measures (normals, area, centroid, edge length, angles) and
whole-mesh measures (surface area, center of mass, bounding box).
Run:  PYTHONPATH=build/bin python examples/08_geometry.py
"""
import simplemesh as sm

m = sm.Mesh()
a = m.add_vertex(sm.Vec3(0, 0, 0))
b = m.add_vertex(sm.Vec3(2, 0, 0))
c = m.add_vertex(sm.Vec3(0, 2, 0))
f = m.add_triangle(a, b, c)

print("face normal   :", m.calc_face_normal(f))
print("face area     :", m.calc_face_area(f))         # 2.0
print("face centroid :", m.calc_face_centroid(f))
e = m.edge(m.find_halfedge(a, b))
print("edge length   :", m.calc_edge_length(e))       # 2.0
print("edge midpoint :", m.calc_edge_midpoint(e))
print("vertex normal :", m.calc_vertex_normal(a))

print("surface area  :", m.surface_area())
print("center of mass:", m.center_of_mass())
lo, hi = m.bounding_box()
print("bounding box  :", lo, "->", hi)
