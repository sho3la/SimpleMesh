// 09 - Subdivision schemes (C++).
//   Six schemes on one tetrahedron, plus adaptive longest-edge refinement.
#include "simplemesh/Mesh.h"
#include "simplemesh/Algorithms.h"
#include <iostream>

static sm::Mesh tetra() {
    sm::Mesh m;
    sm::VertexHandle v[4] = {
        m.add_vertex({0, 0, 0}), m.add_vertex({1, 0, 0}),
        m.add_vertex({0, 1, 0}), m.add_vertex({0, 0, 1})
    };
    int t[4][3] = {{0,2,1},{0,1,3},{1,2,3},{2,0,3}};
    for (auto& f : t) m.add_triangle(v[f[0]], v[f[1]], v[f[2]]);
    return m;
}

int main() {
    sm::Mesh base = tetra();
    std::cout << "base             : " << base.n_faces() << " faces\n";
    std::cout << "loop x1          : " << sm::loop_subdivide(base, 1).n_faces() << " faces (x4)\n";
    std::cout << "sqrt3 x1         : " << sm::sqrt3_subdivide(base, 1).n_faces() << " faces (x3)\n";
    std::cout << "catmull-clark x1 : " << sm::catmull_clark(base, 1).n_faces() << " quad faces\n";
    std::cout << "midpoint x1      : " << sm::midpoint_subdivide(base, 1).n_faces() << " faces\n";
    std::cout << "butterfly x1     : " << sm::butterfly_subdivide(base, 1).n_faces() << " faces\n";

    sm::Mesh big;
    auto a = big.add_vertex({0, 0, 0}), b = big.add_vertex({4, 0, 0}), c = big.add_vertex({0, 4, 0});
    big.add_triangle(a, b, c);
    sm::longest_edge_subdivide(big, 1.5);
    std::cout << "longest-edge<=1.5: " << big.n_faces() << " faces (adaptive, was 1)\n";
    return 0;
}
