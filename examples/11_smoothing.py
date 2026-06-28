"""11 - Smoothing: uniform vs cotangent Laplacian.

Both pull a perturbed vertex back toward its neighbourhood; cotangent weights
respect the surface metric so the result is independent of triangulation.
Boundary vertices stay fixed. Run:
    PYTHONPATH=build/bin python examples/11_smoothing.py
"""
import simplemesh as sm


def grid(n=6):
    m = sm.Mesh()
    v = [[m.add_vertex(sm.Vec3(x, y, 0)) for x in range(n)] for y in range(n)]
    for y in range(n - 1):
        for x in range(n - 1):
            m.add_triangle(v[y][x], v[y][x + 1], v[y + 1][x + 1])
            m.add_triangle(v[y][x], v[y + 1][x + 1], v[y + 1][x])
    return m, v


def spike(m, v):
    p = m.point(v); m.set_point(v, sm.Vec3(p.x, p.y, 1.0))   # lift one vertex

m, v = grid(); mid = v[3][3]
spike(m, mid)
sm.laplacian_smooth(m, iterations=10, lambda_=0.5)
print("uniform  Laplacian: spike z =", round(m.point(mid).z, 4))

m, v = grid(); mid = v[3][3]
spike(m, mid)
sm.cotan_smooth(m, iterations=10, lambda_=0.5)
print("cotangent Laplacian: spike z =", round(m.point(mid).z, 4))
