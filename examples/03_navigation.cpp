// 03 - Navigating the half-edge structure (C++).
//   The navigation primitives, the circulators, and the lazy ranges.
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

    auto h = m.find_halfedge(v0, v2);                 // the shared diagonal
    std::cout << "halfedge v0->v2 : " << h.idx() << "\n";
    std::cout << "  to / from     : " << m.to_vertex(h).idx() << " / " << m.from_vertex(h).idx() << "\n";
    std::cout << "  opposite      : " << m.opposite_halfedge(h).idx() << "\n";
    std::cout << "  next / prev   : " << m.next_halfedge(h).idx() << " / " << m.prev_halfedge(h).idx() << "\n";
    std::cout << "  interior edge : " << (m.is_boundary(m.edge(h)) ? "no" : "yes") << "\n";

    std::cout << "1-ring of v0    :";
    for (auto w : m.vertex_vertices(v0)) std::cout << " " << w.idx();
    std::cout << "\nvalence of v0   : " << m.valence(v0) << "\n";

    // lazy ranges (zero allocation) work in range-for:
    std::cout << "neighbours (lazy):";
    for (auto w : m.vv_range(v0)) std::cout << " " << w.idx();
    std::cout << "\n";
    return 0;
}
