// 10 - Decimation: QEM and the modular Decimater (C++).
#include "simplemesh/Mesh.h"
#include "simplemesh/Algorithms.h"
#include "simplemesh/Decimater.h"
#include <iostream>

static sm::Mesh sphere(int iters) {
    sm::Mesh m;
    sm::VertexHandle v[4] = {
        m.add_vertex({0, 0, 0}), m.add_vertex({1, 0, 0}),
        m.add_vertex({0, 1, 0}), m.add_vertex({0, 0, 1})
    };
    int t[4][3] = {{0,2,1},{0,1,3},{1,2,3},{2,0,3}};
    for (auto& f : t) m.add_triangle(v[f[0]], v[f[1]], v[f[2]]);
    return sm::loop_subdivide(m, iters);
}

int main() {
    {   // one-call QEM
        sm::Mesh m = sphere(3);
        std::cout << "before        : " << m.n_faces() << " faces\n";
        sm::quadric_decimate(m, 40);
        std::cout << "quadric -> 40 : " << m.n_faces() << " faces\n";
    }
    {   // modular: quadric priority + a no-foldover guard
        sm::Mesh m = sphere(3);
        sm::Decimator dec(m);
        dec.add_module<sm::ModQuadric>();             // priority: shape error
        dec.add_module<sm::ModNormalFlipping>(1.0);   // binary: no big foldovers
        size_t n = dec.decimate_to_faces(40);
        std::cout << "modular -> 40 : " << m.n_faces() << " faces in " << n << " collapses\n";
    }
    return 0;
}
