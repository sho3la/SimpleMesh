"""07 - NumPy interchange: bulk transfer and zero-copy views.

Move whole vertex/face arrays in and out, and edit positions in place through a
view that aliases the mesh's own buffer. Requires NumPy.
Run:  PYTHONPATH=build/bin python examples/07_numpy.py
"""
import numpy as np
import simplemesh as sm

# build from arrays (one quad as two triangles)
verts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
tris = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int32)

m = sm.Mesh()
m.add_vertices(verts)
m.add_triangles(tris)
print("loaded:", m.n_vertices, "V", m.n_faces, "F")

# copy data back out
print("get_vertices:\n", m.get_vertices())
print("get_triangles:\n", m.get_triangles())

# zero-copy view: editing the array edits the mesh in place
view = m.points_view()          # shape (N, 3), aliases the mesh buffer
view[:, 2] += 5.0               # lift every vertex in Z
print("vertex 0 after in-place edit:", m.point(sm.VertexHandle(0)))

# bulk per-vertex normals as an (N,3) array
print("vertex normals:\n", m.vertex_normals_to_numpy())
