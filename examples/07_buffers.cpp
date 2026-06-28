// 07 - Direct buffer access (C++ analogue of the NumPy example).
//   In C++ the position buffer is a contiguous std::vector<Vec3> you can read
//   and edit in place; the same memory is what the Python `points_view` aliases.
#include "simplemesh/Mesh.h"
#include <iostream>

int main() {
    sm::Mesh m;
    auto v0 = m.add_vertex({0, 0, 0});
    auto v1 = m.add_vertex({1, 0, 0});
    auto v2 = m.add_vertex({1, 1, 0});
    auto v3 = m.add_vertex({0, 1, 0});
    m.add_triangle(v0, v1, v2);
    m.add_triangle(v0, v2, v3);

    // contiguous position buffer (3*N packed doubles) - edit in place
    std::vector<sm::Vec3>& pts = m.points();
    for (auto& p : pts) p.z += 5.0;                 // lift every vertex in Z
    std::cout << "vertex 0 after in-place edit: " << m.point(v0) << "\n";

    std::cout << "per-vertex normals:\n";
    for (auto v : m.all_vertices())
        std::cout << "  v" << v.idx() << " " << m.calc_vertex_normal(v) << "\n";
    return 0;
}
