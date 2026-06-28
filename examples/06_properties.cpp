// 06 - Custom properties and status bits (C++).
//   Attach named typed data to any element kind; toggle status flags.
#include "simplemesh/Mesh.h"
#include <iostream>

int main() {
    sm::Mesh m;
    auto a = m.add_vertex({0, 0, 0});
    auto b = m.add_vertex({1, 0, 0});
    auto c = m.add_vertex({0, 1, 0});
    auto f = m.add_triangle(a, b, c);

    // a per-vertex scalar "quality"
    auto quality = m.add_vertex_property<double>("quality", 0.0);
    m.property(quality, a) = 0.9;
    m.property(quality, b) = 0.5;
    std::cout << "quality[a] = " << m.property(quality, a) << "\n";
    std::cout << "n vertex properties: " << m.n_vertex_properties() << "\n";

    // a per-face Vec3 "color"
    auto color = m.add_face_property<sm::Vec3>("color", {0, 0, 0});
    m.property(color, f) = {1, 0, 0};
    std::cout << "color[f]   = " << m.property(color, f) << "\n";

    // status bits
    m.set_selected(a);
    m.set_locked(b);
    std::cout << "a selected : " << m.is_selected(a) << "\n";
    std::cout << "b locked   : " << m.is_locked(b) << "\n";
    std::cout << "a status word: " << m.get_status(a)
              << " (SELECTED = " << (int)sm::Mesh::SELECTED << ")\n";
    return 0;
}
