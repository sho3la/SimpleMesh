// 01 - Build a mesh by hand (C++).
//   Create vertices and faces, then read back basic counts and connectivity.
#include "simplemesh/Mesh.h"
#include <iostream>

int main() {
    sm::Mesh m;

    // a tetrahedron: 4 vertices, 4 triangular faces
    sm::VertexHandle v[4] = {
        m.add_vertex({0, 0, 0}), m.add_vertex({1, 0, 0}),
        m.add_vertex({0, 1, 0}), m.add_vertex({0, 0, 1})
    };
    m.add_triangle(v[0], v[2], v[1]);
    m.add_triangle(v[0], v[1], v[3]);
    m.add_triangle(v[1], v[2], v[3]);
    m.add_triangle(v[2], v[0], v[3]);

    // add_face also accepts polygons of any size, e.g. m.add_face({a,b,c,d});

    std::cout << "vertices : " << m.n_vertices() << "\n";
    std::cout << "edges    : " << m.n_edges()    << "\n";
    std::cout << "faces    : " << m.n_faces()    << "\n";
    std::cout << "closed   : " << (m.is_boundary(v[0]) ? "no" : "yes") << "\n";

    std::cout << "face 0 vertices:";
    for (auto vh : m.face_vertices(sm::FaceHandle(0))) std::cout << " " << vh.idx();
    std::cout << "\n";
    return 0;
}
