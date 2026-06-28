"""06 - Custom properties and status bits.

Attach named typed data (float / int / Vec3) to any element kind, and toggle the
built-in status flags. Both grow with the mesh automatically.
Run:  PYTHONPATH=build/bin python examples/06_properties.py
"""
import simplemesh as sm

m = sm.Mesh()
a = m.add_vertex(sm.Vec3(0, 0, 0))
b = m.add_vertex(sm.Vec3(1, 0, 0))
c = m.add_vertex(sm.Vec3(0, 1, 0))
f = m.add_triangle(a, b, c)

# a per-vertex scalar "quality"
quality = m.add_vertex_float_property("quality", 0.0)
m.set_vertex_float(quality, a, 0.9)
m.set_vertex_float(quality, b, 0.5)
print("quality[a] =", m.get_vertex_float(quality, a))
print("n vertex properties:", m.n_vertex_properties)

# a per-face Vec3 "color"
color = m.add_face_vec3_property("color", sm.Vec3(0, 0, 0))
m.set_face_vec3(color, f, sm.Vec3(1, 0, 0))
print("color[f]   =", m.get_face_vec3(color, f))

# status bits
m.set_vertex_selected(a, True)
m.set_vertex_locked(b, True)
print("a selected :", m.is_vertex_selected(a))
print("b locked   :", m.is_vertex_locked(b))
print("a status word:", m.vertex_status(a), "(SELECTED =", sm.SELECTED, ")")
