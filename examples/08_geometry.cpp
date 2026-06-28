// 08 - Geometry queries (C++).
#include "simplemesh/Mesh.h"
#include <iostream>

int main() {
    sm::Mesh m;
    auto a = m.add_vertex({0, 0, 0});
    auto b = m.add_vertex({2, 0, 0});
    auto c = m.add_vertex({0, 2, 0});
    auto f = m.add_triangle(a, b, c);

    std::cout << "face normal   : " << m.calc_face_normal(f)   << "\n";
    std::cout << "face area      : " << m.calc_face_area(f)    << "\n";   // 2
    std::cout << "face centroid : " << m.calc_face_centroid(f) << "\n";
    auto e = m.edge(m.find_halfedge(a, b));
    std::cout << "edge length   : " << m.calc_edge_length(e)   << "\n";   // 2
    std::cout << "edge midpoint : " << m.calc_edge_midpoint(e) << "\n";
    std::cout << "vertex normal : " << m.calc_vertex_normal(a) << "\n";

    std::cout << "surface area  : " << m.surface_area()   << "\n";
    std::cout << "center of mass: " << m.center_of_mass() << "\n";
    auto bb = m.bounding_box();
    std::cout << "bounding box  : " << bb.first << " -> " << bb.second << "\n";
    return 0;
}
