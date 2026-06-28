// 12 - Validate and repair (C++).
//   MeshChecker reports defects; PreMesh / repair fix them; fill_holes closes loops.
#include "simplemesh/Mesh.h"
#include "simplemesh/Algorithms.h"
#include "simplemesh/MeshChecker.h"
#include "simplemesh/PreMesh.h"
#include <iostream>

int main() {
    // --- check a clean mesh ---
    {
        sm::Mesh m;
        sm::VertexHandle v[4] = {
            m.add_vertex({0,0,0}), m.add_vertex({1,0,0}),
            m.add_vertex({0,1,0}), m.add_vertex({0,0,1})
        };
        int t[4][3] = {{0,2,1},{0,1,3},{1,2,3},{2,0,3}};
        for (auto& f : t) m.add_triangle(v[f[0]], v[f[1]], v[f[2]]);
        sm::MeshCheckReport r = sm::MeshChecker(m).check();
        std::cout << "tetra: manifold=" << r.is_manifold << " closed=" << r.is_closed
                  << " euler=" << r.euler << " genus=" << r.genus << "\n";
    }

    // --- repair a non-manifold "soup": 3 triangles sharing one edge (book spine) ---
    {
        std::vector<sm::Vec3> pos = {{0,0,0},{1,0,0},{0,1,0},{0,-1,0},{0,0,1}};
        std::vector<std::vector<int>> faces = {{0,1,2},{0,1,3},{0,1,4}};  // edge 0-1 x3
        sm::PreMesh p = sm::PreMesh::from_soup(pos, faces);
        sm::RepairReport rr = p.repair();
        sm::Mesh clean = p.to_mesh();
        std::cout << "spine soup -> manifold=" << sm::MeshChecker(clean).check().is_manifold
                  << ", " << clean.n_faces() << " faces (nm_edges found=" << rr.nm_edges << ")\n";
    }

    // --- fill a hole: an open box (cube missing its top) ---
    {
        sm::Mesh box;
        sm::VertexHandle b[8];
        double p[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
        for (int i = 0; i < 8; ++i) b[i] = box.add_vertex({p[i][0], p[i][1], p[i][2]});
        int q[5][4] = {{0,3,2,1},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}};  // top open
        for (auto& f : q) box.add_face({b[f[0]], b[f[1]], b[f[2]], b[f[3]]});
        int n = sm::fill_holes(box);
        std::cout << "open box: filled " << n << " hole(s) -> closed="
                  << sm::MeshChecker(box).check().is_closed << ", " << box.n_faces() << " faces\n";
    }
    return 0;
}
