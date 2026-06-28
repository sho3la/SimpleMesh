// 11 - Smoothing: uniform vs cotangent Laplacian (C++).
//   Both pull a perturbed vertex back toward its neighbourhood; boundary fixed.
#include "simplemesh/Mesh.h"
#include "simplemesh/Algorithms.h"
#include <iostream>
#include <vector>

static sm::Mesh grid(int n, sm::VertexHandle& mid) {
    sm::Mesh m;
    std::vector<std::vector<sm::VertexHandle>> v(n, std::vector<sm::VertexHandle>(n));
    for (int y = 0; y < n; ++y)
        for (int x = 0; x < n; ++x)
            v[y][x] = m.add_vertex({double(x), double(y), 0});
    for (int y = 0; y < n - 1; ++y)
        for (int x = 0; x < n - 1; ++x) {
            m.add_triangle(v[y][x], v[y][x+1], v[y+1][x+1]);
            m.add_triangle(v[y][x], v[y+1][x+1], v[y+1][x]);
        }
    mid = v[n/2][n/2];
    return m;
}

int main() {
    sm::VertexHandle mid;
    {
        sm::Mesh m = grid(6, mid);
        sm::Vec3 p = m.point(mid); p.z = 1.0; m.set_point(mid, p);  // spike
        sm::laplacian_smooth(m, 10, 0.5);
        std::cout << "uniform   Laplacian: spike z = " << m.point(mid).z << "\n";
    }
    {
        sm::Mesh m = grid(6, mid);
        sm::Vec3 p = m.point(mid); p.z = 1.0; m.set_point(mid, p);
        sm::cotan_smooth(m, 10, 0.5);
        std::cout << "cotangent Laplacian: spike z = " << m.point(mid).z << "\n";
    }
    return 0;
}
