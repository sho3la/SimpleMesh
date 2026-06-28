"""09 - Subdivision schemes.

Six schemes on one tetrahedron. Approximating (Loop, Catmull-Clark, sqrt3) move
the surface toward a smooth limit; interpolating (midpoint, butterfly) keep the
original vertices; longest-edge refines adaptively in place.
Run:  PYTHONPATH=build/bin python examples/09_subdivision.py
"""
import simplemesh as sm


def tetra():
    m = sm.Mesh()
    v = [m.add_vertex(sm.Vec3(*p)) for p in
         [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]]
    for t in [(0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)]:
        m.add_triangle(v[t[0]], v[t[1]], v[t[2]])
    return m


base = tetra()
print("base               :", base.n_faces, "faces")
print("loop x1            :", sm.loop_subdivide(base, 1).n_faces, "faces (x4)")
print("sqrt3 x1           :", sm.sqrt3_subdivide(base, 1).n_faces, "faces (x3)")
print("catmull-clark x1   :", sm.catmull_clark(base, 1).n_faces, "quad faces")
print("midpoint x1        :", sm.midpoint_subdivide(base, 1).n_faces, "faces (interpolating)")
print("butterfly x1       :", sm.butterfly_subdivide(base, 1).n_faces, "faces (interpolating)")

# adaptive: bisect every edge longer than the bound, in place
big = sm.Mesh()
a = big.add_vertex(sm.Vec3(0, 0, 0)); b = big.add_vertex(sm.Vec3(4, 0, 0))
c = big.add_vertex(sm.Vec3(0, 4, 0))
big.add_triangle(a, b, c)
sm.longest_edge_subdivide(big, 1.5)
print("longest-edge<=1.5  :", big.n_faces, "faces (adaptive, was 1)")
