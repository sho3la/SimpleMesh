// 05 - Topological editing: flip, split, collapse, delete + GC (C++).
//   Each operator is guarded by an is_*_ok predicate; deletion is lazy.
#include "simplemesh/Mesh.h"
#include <iostream>

static sm::Mesh quad(sm::VertexHandle& a, sm::VertexHandle& b,
                     sm::VertexHandle& c, sm::VertexHandle& d) {
    sm::Mesh m;
    a = m.add_vertex({0, 0, 0}); b = m.add_vertex({1, 0, 0});
    c = m.add_vertex({1, 1, 0}); d = m.add_vertex({0, 1, 0});
    m.add_triangle(a, b, c); m.add_triangle(a, c, d);
    return m;
}

int main() {
    sm::VertexHandle a, b, c, d;

    {   // edge flip
        sm::Mesh m = quad(a, b, c, d);
        auto e = m.edge(m.find_halfedge(a, c));
        std::cout << "flip ok? " << m.is_flip_ok(e) << "\n";
        m.flip(e);
        std::cout << "  a-c gone: " << !m.find_halfedge(a, c).is_valid()
                  << ", b-d present: " << m.find_halfedge(b, d).is_valid() << "\n";
    }
    {   // edge split
        sm::Mesh m = quad(a, b, c, d);
        auto mid = m.add_vertex({0.5, 0.5, 0});
        m.split(m.edge(m.find_halfedge(a, c)), mid);
        std::cout << "split: faces " << m.n_faces() << " (was 2)\n";
    }
    {   // collapse
        sm::Mesh m = quad(a, b, c, d);
        auto h = m.find_halfedge(a, c);
        std::cout << "collapse ok? " << m.is_collapse_ok(h) << "\n";
        if (m.is_collapse_ok(h)) { m.collapse(h); m.garbage_collection(); }
        std::cout << "  after collapse+GC: V " << m.n_vertices() << " F " << m.n_faces() << "\n";
    }
    {   // delete + garbage collection
        sm::Mesh m = quad(a, b, c, d);
        m.delete_face(sm::FaceHandle(0));
        std::cout << "deleted flagged: " << m.is_deleted(sm::FaceHandle(0)) << "\n";
        m.garbage_collection();
        std::cout << "  after GC: F " << m.n_faces() << "\n";
    }
    return 0;
}
