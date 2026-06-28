// 04 - Smart handles: fluent navigation (C++).
//   A smart handle bundles an element with its mesh so navigation chains.
#include "simplemesh/Mesh.h"
#include <iostream>

int main() {
    sm::Mesh m;
    auto v0 = m.add_vertex({0, 0, 0});
    auto v1 = m.add_vertex({1, 0, 0});
    auto v2 = m.add_vertex({1, 1, 0});
    auto v3 = m.add_vertex({0, 1, 0});
    auto f0 = m.add_triangle(v0, v1, v2);
    m.add_triangle(v0, v2, v3);

    auto h = m.smart(m.halfedge(f0));
    std::cout << "h.next().opp().to() = " << h.next().opp().to().idx() << "\n";
    std::cout << "opp().opp() == h    : " << (h.opp().opp().idx() == h.idx()) << "\n";

    auto sv = m.smart(v0);
    std::cout << "v0 1-ring  :";
    for (auto w : sv.vertices()) std::cout << " " << w.idx();
    std::cout << "\nv0 valence : " << sv.valence() << "  boundary: " << sv.is_boundary() << "\n";

    auto sf = m.smart(f0);
    std::cout << "face area   : " << sf.area() << "\n";
    std::cout << "face normal : " << sf.normal() << "\n";
    return 0;
}
