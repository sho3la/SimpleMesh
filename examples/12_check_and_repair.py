"""12 - Validate and repair.

MeshChecker reports what is wrong; repair_mesh / load_and_repair fix it; fill_holes
closes boundary loops. The repair layer turns non-manifold "soup" into a clean
manifold the half-edge kernel accepts.
Run:  PYTHONPATH=build/bin python examples/12_check_and_repair.py
"""
import simplemesh as sm

# --- check a clean mesh ---
m = sm.Mesh()
v = [m.add_vertex(sm.Vec3(*p)) for p in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]]
for t in [(0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)]:
    m.add_triangle(v[t[0]], v[t[1]], v[t[2]])
r = sm.MeshChecker(m).check()
print("tetra: manifold=%s closed=%s euler=%d genus=%d" %
      (r.is_manifold, r.is_closed, r.euler, r.genus))

# --- repair a non-manifold "soup": 3 triangles sharing one edge (a book spine) ---
# Write it as a tiny OBJ, then load_and_repair makes it halfedge-legal.
import tempfile, os
pos = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1)]
faces = [[0, 1, 2], [0, 1, 3], [0, 1, 4]]            # edge 0-1 shared by 3 faces
objp = os.path.join(tempfile.mkdtemp(), "spine.obj")
with open(objp, "w") as fh:
    for p in pos: fh.write("v %g %g %g\n" % p)
    for f in faces: fh.write("f %d %d %d\n" % (f[0] + 1, f[1] + 1, f[2] + 1))
clean = sm.load_and_repair(objp)
cr = sm.MeshChecker(clean).check()
print("spine soup -> manifold=%s, %d faces (all recovered, none dropped)" %
      (cr.is_manifold, clean.n_faces))

# --- fill a hole: an open box (cube missing its top) ---
box = sm.Mesh()
p = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
     (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
b = [box.add_vertex(sm.Vec3(*q)) for q in p]
for q in [(0, 3, 2, 1), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]:
    box.add_face([b[i] for i in q])                  # 5 faces: top is open
before = sm.MeshChecker(box).check()
n = sm.fill_holes(box)
after = sm.MeshChecker(box).check()
print("open box: filled %d hole(s) -> closed=%s, %d faces" %
      (n, after.is_closed, box.n_faces))
