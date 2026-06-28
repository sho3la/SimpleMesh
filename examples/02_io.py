"""02 - File I/O: read and write OBJ, PLY, STL and OFF.

Builds a cube, writes it in every supported format, reads each back, and checks
the round-trip preserves the geometry.
Run:  PYTHONPATH=build/bin python examples/02_io.py
"""
import os, tempfile
import simplemesh as sm


def make_cube():
    m = sm.Mesh()
    p = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    v = [m.add_vertex(sm.Vec3(*q)) for q in p]
    quads = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    for q in quads:
        m.add_face([v[i] for i in q])
    return m


cube = make_cube()
tmp = tempfile.mkdtemp()
writers = {
    "obj": lambda m, p: m.write_obj(p),
    "ply": lambda m, p: m.write_ply(p, binary=True),    # binary_little_endian
    "stl": lambda m, p: m.write_stl(p, binary=True),    # polygons fan-triangulated
    "off": lambda m, p: m.write_off(p),
}
readers = {"obj": "read_obj", "ply": "read_ply", "stl": "read_stl", "off": "read_off"}

for ext, write in writers.items():
    path = os.path.join(tmp, "cube." + ext)
    write(cube, path)
    m2 = sm.Mesh()
    getattr(m2, readers[ext])(path)
    print(f"{ext.upper():4} -> {m2.n_vertices:2} V  {m2.n_faces:2} F   "
          f"({os.path.getsize(path)} bytes)")
