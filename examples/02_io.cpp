// 02 - File I/O: read and write OBJ, PLY, STL and OFF (C++).
//   Build a cube, write every format, read each back, check the round-trip.
#include "simplemesh/Mesh.h"
#include <iostream>

static sm::Mesh make_cube() {
    sm::Mesh m;
    sm::VertexHandle v[8];
    double p[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for (int i = 0; i < 8; ++i) v[i] = m.add_vertex({p[i][0], p[i][1], p[i][2]});
    int q[6][4] = {{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}};
    for (auto& f : q) m.add_face({v[f[0]], v[f[1]], v[f[2]], v[f[3]]});
    return m;
}

int main() {
    sm::Mesh cube = make_cube();

    cube.write_obj("cube.obj");
    cube.write_ply("cube.ply", /*binary=*/true);
    cube.write_stl("cube.stl", /*binary=*/true);   // quads fan-triangulated
    cube.write_off("cube.off");

    sm::Mesh a; a.read_obj("cube.obj");
    sm::Mesh b; b.read_ply("cube.ply");
    sm::Mesh c; c.read_stl("cube.stl");
    sm::Mesh d; d.read_off("cube.off");

    std::cout << "OBJ -> " << a.n_vertices() << " V " << a.n_faces() << " F\n";
    std::cout << "PLY -> " << b.n_vertices() << " V " << b.n_faces() << " F\n";
    std::cout << "STL -> " << c.n_vertices() << " V " << c.n_faces() << " F (welded, triangulated)\n";
    std::cout << "OFF -> " << d.n_vertices() << " V " << d.n_faces() << " F\n";
    return 0;
}
