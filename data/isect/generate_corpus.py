#!/usr/bin/env python3
"""
generate_corpus.py - reproducibly emit the self-intersection test corpus.

Each mesh is small and hand-designed so its *expected* intersection behaviour is
known by construction (documented per function). These OBJs are the fixed inputs
for every stage of the exact-arrangement pipeline: the same files feed the
numeric golden-file tests, the per-facet CDT tests, and the end-to-end
arrangement test.

Run from anywhere:  python tests/data/isect/generate_corpus.py
Writes the .obj files next to this script.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def write_obj(name, verts, faces, doc):
    """verts: list of (x,y,z); faces: list of 1-based index tuples."""
    path = os.path.join(HERE, name)
    with open(path, "w") as f:
        f.write(f"# {name}\n")
        for line in doc.strip().splitlines():
            f.write(f"# {line.strip()}\n")
        for x, y, z in verts:
            f.write(f"v {x:g} {y:g} {z:g}\n")
        for face in faces:
            f.write("f " + " ".join(str(i) for i in face) + "\n")
    print(f"wrote {name}: {len(verts)} verts, {len(faces)} faces")


def cube(ox, oy, oz, s=1.0):
    """8 corners + 12 CCW (outward) triangles for an axis-aligned cube at origin
    (ox,oy,oz) with side s. Returns (verts, faces) with faces 0-based here."""
    v = [
        (ox,     oy,     oz),      # 0
        (ox + s, oy,     oz),      # 1
        (ox + s, oy + s, oz),      # 2
        (ox,     oy + s, oz),      # 3
        (ox,     oy,     oz + s),  # 4
        (ox + s, oy,     oz + s),  # 5
        (ox + s, oy + s, oz + s),  # 6
        (ox,     oy + s, oz + s),  # 7
    ]
    f = [
        (0, 3, 2), (0, 2, 1),   # bottom  z=oz   (-z out)
        (4, 5, 6), (4, 6, 7),   # top     z=oz+s (+z out)
        (0, 1, 5), (0, 5, 4),   # front   y=oy   (-y out)
        (3, 7, 6), (3, 6, 2),   # back    y=oy+s (+y out)
        (0, 4, 7), (0, 7, 3),   # left    x=ox   (-x out)
        (1, 2, 6), (1, 6, 5),   # right   x=ox+s (+x out)
    ]
    return v, f


def main():
    # 1. two_tris_cross: triangle A in z=0 plane, triangle B vertical (plane
    #    y=0.5) passing straight through it. EXPECT: 1 intersecting pair, the
    #    intersection is an interior segment of both triangles.
    write_obj(
        "two_tris_cross.obj",
        [(0, 0, 0), (2, 0, 0), (0, 2, 0),
         (0.5, 0.5, -1), (0.5, 0.5, 1), (1.5, 0.5, 0)],
        [(1, 2, 3), (4, 5, 6)],
        "Two triangles crossing in an X. EXPECT exactly one intersecting pair.",
    )

    # 2. two_quads_overlap: quad A in z=0, quad B vertical (plane y=1) crossing
    #    through. Each quad is two triangles. EXPECT: B's triangles cross A's.
    write_obj(
        "two_quads_overlap.obj",
        [(0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0),
         (0, 1, -1), (2, 1, -1), (2, 1, 1), (0, 1, 1)],
        [(1, 2, 3), (1, 3, 4), (5, 6, 7), (5, 7, 8)],
        "Two quads (4 triangles) crossing. EXPECT crossing pairs between sheets.",
    )

    # 3. cube_x_cube: unit cube [0,1]^3 and cube [0.5,1.5]^3. Canonical boolean
    #    test (union/intersection/difference). EXPECT many intersecting pairs.
    va, fa = cube(0.0, 0.0, 0.0, 1.0)
    vb, fb = cube(0.5, 0.5, 0.5, 1.0)
    verts = va + vb
    faces = [tuple(i + 1 for i in t) for t in fa] + \
            [tuple(i + 1 + len(va) for i in t) for t in fb]
    write_obj("cube_x_cube.obj", verts, faces,
              "Two overlapping unit cubes. Canonical boolean op input.")

    # 4. self_fold: a connected ribbon of 5 quads that goes out flat, climbs,
    #    arcs back over the top, then plunges DOWN through the very first (flat)
    #    quad. The piercing quad (cols 4-5, verts 9..12) shares no vertices with
    #    the pierced flat quad (cols 0-1, verts 1..4), so it is a genuine
    #    non-adjacent self-intersection inside ONE connected component.
    #    Columns (each is a y=0 / y=1 vertex pair):
    #      col0 x=0 z=0   col1 x=1 z=0   col2 x=2 z=0
    #      col3 x=2 z=2   col4 x=1 z=2   col5 x=0 z=-1
    #    Segment col4->col5 runs z:2->-1 over x:1->0, crossing z=0 near x=0.33,
    #    inside the flat quad's footprint [0,1]x[0,1]. EXPECT >=1 crossing pair.
    write_obj(
        "self_fold.obj",
        [(0, 0, 0), (0, 1, 0),      # col0  v1,v2
         (1, 0, 0), (1, 1, 0),      # col1  v3,v4
         (2, 0, 0), (2, 1, 0),      # col2  v5,v6
         (2, 0, 2), (2, 1, 2),      # col3  v7,v8
         (1, 0, 2), (1, 1, 2),      # col4  v9,v10
         (0, 0, -1), (0, 1, -1)],   # col5  v11,v12
        [(1, 3, 4), (1, 4, 2),      # Q0 flat (the pierced quad)
         (3, 5, 6), (3, 6, 4),      # Q1
         (5, 7, 8), (5, 8, 6),      # Q2 up
         (7, 9, 10), (7, 10, 8),    # Q3 over the top
         (9, 11, 12), (9, 12, 10)], # Q4 plunges through Q0
        "Folded ribbon; final quad pierces the first quad. Self-intersection.",
    )

    # 5. coplanar_overlap: two coplanar triangles in z=0 with overlapping area.
    #    Degenerate (the intersection is 2D area, not a clean segment). Stress
    #    case for the exact coplanar handling.
    write_obj(
        "coplanar_overlap.obj",
        [(0, 0, 0), (2, 0, 0), (0, 2, 0),
         (1.5, -0.5, 0), (1.5, 1.5, 0), (-0.5, 1.5, 0)],
        [(1, 2, 3), (4, 5, 6)],
        "Two coplanar overlapping triangles. Degenerate area overlap.",
    )

    # 6. tangent_touch: two triangles meeting at a single shared point (the
    #    origin) but otherwise disjoint and non-coplanar. EXPECT: a boundary
    #    touch, NOT an interior crossing - the exact path treats this symbolically
    #    (vertex region), and no shell removal should be needed.
    write_obj(
        "tangent_touch.obj",
        [(0, 0, 0), (1, 0, 0), (0, 1, 0),
         (0, 0, 0), (-1, 0, 0.0001), (0, -1, 1)],
        [(1, 2, 3), (4, 5, 6)],
        "Two triangles touching at one point. Boundary (vertex) contact only.",
    )


if __name__ == "__main__":
    main()
