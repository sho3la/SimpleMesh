"""10 - Decimation: QEM and the modular Decimater.

quadric_decimate is the one-call Garland-Heckbert simplifier; Decimator lets you
combine a priority module with binary veto modules (no foldovers, length caps...).
Run:  PYTHONPATH=build/bin python examples/10_decimation.py
"""
import simplemesh as sm


def sphere(iters=3):
    m = sm.Mesh()
    v = [m.add_vertex(sm.Vec3(*p)) for p in
         [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]]
    for t in [(0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)]:
        m.add_triangle(v[t[0]], v[t[1]], v[t[2]])
    return sm.loop_subdivide(m, iters)


# --- one-call QEM ---
m = sphere()
print("before        :", m.n_faces, "faces")
sm.quadric_decimate(m, target_faces=40)
print("quadric -> 40 :", m.n_faces, "faces")

# --- modular: quadric priority + a no-foldover guard ---
m = sphere()
dec = sm.Decimator(m)
dec.add_quadric_module()                 # priority: shape error
dec.add_normal_flipping_module(1.0)      # binary: reject collapses that fold > ~57 deg
n = dec.decimate_to_faces(40)
print("modular -> 40 :", m.n_faces, "faces in", n, "collapses")
