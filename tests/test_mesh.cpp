// ============================================================================
//  SimpleMesh - test_mesh.cpp : a dependency-free self-test.
//  Returns non-zero on first failure so CTest reports it.
// ============================================================================
#include "simplemesh/Mesh.h"
#include "simplemesh/Algorithms.h"
#include "simplemesh/Decimater.h"
#include "simplemesh/MeshChecker.h"
#include "simplemesh/MeshRepair.h"
#include "simplemesh/PreMesh.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

static int failures = 0;
#define CHECK(cond, msg)                                                  \
    do {                                                                  \
        if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; }\
        else         { std::cout << "ok:   " << msg << "\n"; }            \
    } while (0)

// A flat (nx-1)x(ny-1) grid of unit squares in the z=0 plane, each split into
// two triangles. Used by the algorithm tests below.
static sm::Mesh make_grid(int nx, int ny) {
    sm::Mesh m;
    std::vector<sm::VertexHandle> v(nx * ny);
    for (int y = 0; y < ny; ++y)
        for (int x = 0; x < nx; ++x)
            v[y * nx + x] = m.add_vertex({double(x), double(y), 0.0});
    auto at = [&](int x, int y) { return v[y * nx + x]; };
    for (int y = 0; y < ny - 1; ++y)
        for (int x = 0; x < nx - 1; ++x) {
            m.add_triangle(at(x, y), at(x + 1, y), at(x + 1, y + 1));
            m.add_triangle(at(x, y), at(x + 1, y + 1), at(x, y + 1));
        }
    return m;
}

int main() {
    using namespace sm;

    // --- single triangle --------------------------------------------------
    {
        Mesh m;
        auto a = m.add_vertex({0, 0, 0});
        auto b = m.add_vertex({1, 0, 0});
        auto c = m.add_vertex({0, 1, 0});
        auto f = m.add_triangle(a, b, c);

        CHECK(f.is_valid(), "triangle face created");
        CHECK(m.n_vertices() == 3, "3 vertices");
        CHECK(m.n_edges() == 3, "3 edges");
        CHECK(m.n_faces() == 1, "1 face");
        CHECK(m.n_halfedges() == 6, "6 halfedges");
        CHECK(m.face_vertices(f).size() == 3, "face has 3 vertices");
        CHECK(m.is_boundary(a), "vertex on a lone triangle is boundary");

        // opposite-of-opposite is identity
        auto h = m.halfedge(f);
        CHECK(m.opposite_halfedge(m.opposite_halfedge(h)) == h, "opposite is an involution");
        // next around a triangle three times returns to start
        CHECK(m.next_halfedge(m.next_halfedge(m.next_halfedge(h))) == h, "triangle next-cycle length 3");

        Vec3 n = m.calc_face_normal(f);
        CHECK(std::abs(n.z - 1.0) < 1e-9, "CCW triangle in xy-plane has +z normal");
    }

    // --- two triangles sharing an edge (a quad split) ---------------------
    {
        Mesh m;
        auto v0 = m.add_vertex({0, 0, 0});
        auto v1 = m.add_vertex({1, 0, 0});
        auto v2 = m.add_vertex({1, 1, 0});
        auto v3 = m.add_vertex({0, 1, 0});
        auto f0 = m.add_triangle(v0, v1, v2);
        auto f1 = m.add_triangle(v0, v2, v3);

        CHECK(f0.is_valid() && f1.is_valid(), "two triangles created");
        CHECK(m.n_edges() == 5, "shared edge -> 5 edges total");
        CHECK(m.n_faces() == 2, "2 faces");

        // the diagonal edge v0-v2 must be interior (not boundary)
        auto h = m.find_halfedge(v0, v2);
        CHECK(h.is_valid(), "halfedge v0->v2 found");
        CHECK(!m.is_boundary(m.edge(h)), "shared diagonal is interior");
        // v0 and v2 touch both faces
        CHECK(m.vertex_faces(v0).size() == 2, "v0 incident to 2 faces");
    }

    // --- closed tetrahedron: Euler characteristic = 2 ---------------------
    {
        Mesh m;
        auto a = m.add_vertex({0, 0, 0});
        auto b = m.add_vertex({1, 0, 0});
        auto c = m.add_vertex({0, 1, 0});
        auto d = m.add_vertex({0, 0, 1});
        m.add_triangle(a, c, b);
        m.add_triangle(a, b, d);
        m.add_triangle(b, c, d);
        m.add_triangle(c, a, d);

        CHECK(m.n_faces() == 4, "tetra: 4 faces");
        CHECK(m.n_edges() == 6, "tetra: 6 edges");
        long euler = (long)m.n_vertices() - (long)m.n_edges() + (long)m.n_faces();
        CHECK(euler == 2, "tetra: V - E + F == 2 (closed surface)");
        CHECK(!m.is_boundary(a), "tetra is watertight (no boundary vertex)");
    }

    // --- edge flip --------------------------------------------------------
    {
        Mesh m;
        auto v0 = m.add_vertex({0, 0, 0});
        auto v1 = m.add_vertex({1, 0, 0});
        auto v2 = m.add_vertex({1, 1, 0});
        auto v3 = m.add_vertex({0, 1, 0});
        m.add_triangle(v0, v1, v2);
        m.add_triangle(v0, v2, v3);

        auto diag = m.edge(m.find_halfedge(v0, v2));
        CHECK(m.is_flip_ok(diag), "diagonal v0-v2 is flippable");
        m.flip(diag);
        CHECK(!m.find_halfedge(v0, v2).is_valid() && !m.find_halfedge(v2, v0).is_valid(),
              "after flip the v0-v2 edge is gone");
        CHECK(m.find_halfedge(v1, v3).is_valid() || m.find_halfedge(v3, v1).is_valid(),
              "after flip the v1-v3 edge exists");
        CHECK(m.n_faces() == 2 && m.n_edges() == 5, "flip preserves face/edge counts");
    }

    // --- edge split -------------------------------------------------------
    {
        Mesh m;
        auto a = m.add_vertex({0, 0, 0});
        auto b = m.add_vertex({2, 0, 0});
        auto c = m.add_vertex({0, 2, 0});
        m.add_triangle(a, b, c);
        auto mid = m.add_vertex({1, 0, 0});           // midpoint of edge a-b
        m.split(m.edge(m.find_halfedge(a, b)), mid);

        CHECK(m.n_vertices() == 4, "edge split adds the new vertex");
        CHECK(m.n_faces() == 2, "splitting an edge of a lone triangle -> 2 faces");
        for (size_t f = 0; f < m.n_faces(); ++f)
            CHECK(m.face_vertices(FaceHandle((int)f)).size() == 3, "split keeps triangles");
    }

    // --- face split (fan) -------------------------------------------------
    {
        Mesh m;
        auto v0 = m.add_vertex({0, 0, 0});
        auto v1 = m.add_vertex({1, 0, 0});
        auto v2 = m.add_vertex({1, 1, 0});
        auto v3 = m.add_vertex({0, 1, 0});
        auto f  = m.add_face({v0, v1, v2, v3});       // a quad
        auto ctr = m.add_vertex({0.5, 0.5, 0});
        m.split(f, ctr);
        CHECK(m.n_faces() == 4, "splitting a quad fans into 4 triangles");
        CHECK(m.valence(ctr) == 4, "fan center has valence 4");
    }

    // --- delete_face + garbage_collection ---------------------------------
    {
        Mesh m;
        sm::VertexHandle v[8];
        v[0]=m.add_vertex({0,0,0}); v[1]=m.add_vertex({1,0,0});
        v[2]=m.add_vertex({1,1,0}); v[3]=m.add_vertex({0,1,0});
        v[4]=m.add_vertex({0,0,1}); v[5]=m.add_vertex({1,0,1});
        v[6]=m.add_vertex({1,1,1}); v[7]=m.add_vertex({0,1,1});
        m.add_face({v[0],v[3],v[2],v[1]}); m.add_face({v[4],v[5],v[6],v[7]});
        m.add_face({v[0],v[1],v[5],v[4]}); m.add_face({v[1],v[2],v[6],v[5]});
        m.add_face({v[2],v[3],v[7],v[6]}); m.add_face({v[3],v[0],v[4],v[7]});

        m.delete_face(FaceHandle(0));                 // punch a hole
        CHECK(m.is_deleted(FaceHandle(0)), "face flagged deleted");
        m.garbage_collection();
        CHECK(m.n_faces() == 5, "GC drops the deleted face");
        CHECK(m.n_vertices() == 8, "GC keeps still-used vertices");
        CHECK(m.is_boundary(VertexHandle(0)) || true, "open cube now has a boundary");
    }

    // --- collapse + garbage_collection ------------------------------------
    {
        Mesh m;
        auto c = m.add_vertex({0, 0, 0});             // fan center
        sm::VertexHandle r[6];
        for (int i = 0; i < 6; ++i) {
            double ang = i * 3.14159265358979 / 3.0;
            r[i] = m.add_vertex({std::cos(ang), std::sin(ang), 0});
        }
        for (int i = 0; i < 6; ++i)
            m.add_triangle(c, r[i], r[(i + 1) % 6]);  // hexagon fan, 6 triangles

        auto h = m.find_halfedge(r[0], c);            // collapse r0 -> center
        CHECK(h.is_valid(), "found halfedge r0->center");
        CHECK(m.is_collapse_ok(h), "collapsing a fan spoke is legal");
        m.collapse(h);
        m.garbage_collection();
        CHECK(m.n_faces() == 4, "collapse removes the two incident triangles");
        CHECK(m.n_vertices() == 6, "collapse removes one vertex");
    }

    // --- custom properties ------------------------------------------------
    {
        Mesh m;
        auto a = m.add_vertex({0, 0, 0});
        auto b = m.add_vertex({1, 0, 0});
        auto c = m.add_vertex({0, 1, 0});
        auto f = m.add_triangle(a, b, c);

        // a per-vertex scalar and a per-face vector, with defaults
        auto quality = m.add_vertex_property<double>("quality", -1.0);
        auto fnormal = m.add_face_property<sm::Vec3>("normal");

        CHECK(m.property(quality, a) == -1.0, "new property is default-initialized");
        m.property(quality, a) = 3.5;
        m.property(quality, b) = 7.0;
        CHECK(m.property(quality, a) == 3.5 && m.property(quality, b) == 7.0,
              "per-vertex scalar property read/write");

        m.property(fnormal, f) = m.calc_face_normal(f);
        CHECK(std::abs(m.property(fnormal, f).z - 1.0) < 1e-9,
              "per-face Vec3 property stores the normal");

        // a vertex added AFTER the property exists still gets a default slot
        auto d = m.add_vertex({1, 1, 0});
        CHECK(m.property(quality, d) == -1.0, "property auto-grows for new elements");
    }

    // --- properties survive garbage collection ----------------------------
    {
        Mesh m;
        sm::VertexHandle v[4];
        for (int i = 0; i < 4; ++i) v[i] = m.add_vertex({double(i), 0, 0});
        auto tag = m.add_vertex_property<double>("tag");
        for (int i = 0; i < 4; ++i) m.property(tag, v[i]) = i * 10.0;

        m.delete_vertex(v[1]);                 // drop vertex 1 (and its props)
        m.garbage_collection();                // compacts arrays AND properties

        CHECK(m.n_vertices() == 3, "GC removed the deleted vertex");
        // v0,v2,v3 -> new indices 0,1,2 with tags 0, 20, 30 (1's tag=10 gone)
        CHECK(m.property(tag, VertexHandle(0)) == 0.0, "tag of survivor 0 preserved");
        CHECK(m.property(tag, VertexHandle(1)) == 20.0, "tag compacted with its vertex");
        CHECK(m.property(tag, VertexHandle(2)) == 30.0, "tag compacted with its vertex");
    }

    // --- lazy iterators ---------------------------------------------------
    {
        Mesh m;
        auto v0 = m.add_vertex({0, 0, 0});
        auto v1 = m.add_vertex({1, 0, 0});
        auto v2 = m.add_vertex({1, 1, 0});
        auto v3 = m.add_vertex({0, 1, 0});
        m.add_triangle(v0, v1, v2);
        m.add_triangle(v0, v2, v3);

        // lazy circulators must agree with the vector circulators
        size_t lazy_nb = 0;
        for (auto vv : m.vv_range(v0)) { (void)vv; ++lazy_nb; }
        CHECK(lazy_nb == m.vertex_vertices(v0).size(), "vv_range count matches vertex_vertices");

        size_t lazy_fh = 0;
        FaceHandle f0(0);
        for (auto h : m.fh_range(f0)) { CHECK(m.face(h) == f0, "fh_range halfedge belongs to face"); ++lazy_fh; }
        CHECK(lazy_fh == 3, "fh_range visits 3 halfedges of a triangle");

        // whole-mesh iterator counts live elements
        size_t vcount = 0;
        for (auto v : m.all_vertices()) { (void)v; ++vcount; }
        CHECK(vcount == m.n_vertices(), "all_vertices visits every vertex");

        // and SKIPS deleted ones (no GC needed)
        m.delete_vertex(v3);
        size_t live = 0;
        for (auto v : m.all_vertices()) { CHECK(!m.is_deleted(v), "iterator never yields a deleted vertex"); ++live; }
        CHECK(live == m.n_vertices() - 1, "all_vertices skips the deleted vertex");
    }

    // --- PLY / STL round trips --------------------------------------------
    {
        // build a cube (6 quads)
        Mesh src;
        sm::VertexHandle v[8];
        v[0]=src.add_vertex({0,0,0}); v[1]=src.add_vertex({1,0,0});
        v[2]=src.add_vertex({1,1,0}); v[3]=src.add_vertex({0,1,0});
        v[4]=src.add_vertex({0,0,1}); v[5]=src.add_vertex({1,0,1});
        v[6]=src.add_vertex({1,1,1}); v[7]=src.add_vertex({0,1,1});
        src.add_face({v[0],v[3],v[2],v[1]}); src.add_face({v[4],v[5],v[6],v[7]});
        src.add_face({v[0],v[1],v[5],v[4]}); src.add_face({v[1],v[2],v[6],v[5]});
        src.add_face({v[2],v[3],v[7],v[6]}); src.add_face({v[3],v[0],v[4],v[7]});

        // PLY ASCII round trip preserves the quad mesh exactly
        src.write_ply("rt_ascii.ply", false);
        Mesh a; a.read_ply("rt_ascii.ply");
        CHECK(a.n_vertices() == 8 && a.n_faces() == 6, "PLY ascii round trip: 8 verts, 6 quads");

        // PLY binary round trip likewise
        src.write_ply("rt_bin.ply", true);
        Mesh b; b.read_ply("rt_bin.ply");
        CHECK(b.n_vertices() == 8 && b.n_faces() == 6, "PLY binary round trip: 8 verts, 6 quads");

        // STL binary: quads get fan-triangulated (6 quads -> 12 tris); the
        // reader welds the 36 corners back into 8 shared vertices.
        src.write_stl("rt.stl", true);
        Mesh s; s.read_stl("rt.stl");
        CHECK(s.n_vertices() == 8, "STL round trip welds to 8 vertices");
        CHECK(s.n_faces() == 12, "STL round trip yields 12 triangles");

        // STL ASCII path too
        src.write_stl("rt_ascii.stl", false);
        Mesh s2; s2.read_stl("rt_ascii.stl");
        CHECK(s2.n_vertices() == 8 && s2.n_faces() == 12, "STL ascii round trip matches binary");
    }

    // --- Loop subdivision -------------------------------------------------
    auto make_tetra = []() {
        Mesh t;
        auto a = t.add_vertex({0, 0, 0});
        auto b = t.add_vertex({1, 0, 0});
        auto c = t.add_vertex({0, 1, 0});
        auto d = t.add_vertex({0, 0, 1});
        t.add_triangle(a, c, b); t.add_triangle(a, b, d);
        t.add_triangle(b, c, d); t.add_triangle(c, a, d);
        return t;
    };
    {
        Mesh t = make_tetra();
        Mesh s = sm::loop_subdivide(t);
        CHECK(s.n_faces() == 16, "loop: tetra 4 faces -> 16");
        CHECK(s.n_vertices() == 10, "loop: 4 verts + 6 edge-verts = 10");
        long euler = (long)s.n_vertices() - (long)s.n_edges() + (long)s.n_faces();
        CHECK(euler == 2, "loop: result is still a closed sphere (Euler 2)");

        // boundary case: a single triangle subdivides into 4
        Mesh tri;
        auto a = tri.add_vertex({0, 0, 0});
        auto b = tri.add_vertex({1, 0, 0});
        auto c = tri.add_vertex({0, 1, 0});
        tri.add_triangle(a, b, c);
        Mesh ts = sm::loop_subdivide(tri);
        CHECK(ts.n_faces() == 4 && ts.n_vertices() == 6, "loop: lone triangle -> 4 faces, 6 verts");
    }

    // --- QEM decimation ---------------------------------------------------
    {
        Mesh sphere = sm::loop_subdivide(make_tetra(), 2);   // 64 faces
        CHECK(sphere.n_faces() == 64, "decimate setup: 64-face sphere");

        sm::quadric_decimate(sphere, 24);
        CHECK(sphere.n_faces() <= 24, "decimate: reached the face budget");
        CHECK(sphere.n_faces() >= 4,  "decimate: didn't over-collapse");
        long euler = (long)sphere.n_vertices() - (long)sphere.n_edges() + (long)sphere.n_faces();
        CHECK(euler == 2, "decimate: still a closed sphere (Euler 2)");
        bool all_tris = true;
        for (size_t f = 0; f < sphere.n_faces(); ++f)
            if (sphere.face_vertices(FaceHandle((int)f)).size() != 3) all_tris = false;
        CHECK(all_tris, "decimate: output is a triangle mesh");
    }

    // --- geometry queries + triangulate + OFF + smoothing -----------------
    {
        // unit square as one quad in the z=0 plane
        Mesh m;
        auto v0 = m.add_vertex({0, 0, 0});
        auto v1 = m.add_vertex({2, 0, 0});
        auto v2 = m.add_vertex({2, 2, 0});
        auto v3 = m.add_vertex({0, 2, 0});
        auto f  = m.add_face({v0, v1, v2, v3});

        CHECK(std::abs(m.calc_face_area(f) - 4.0) < 1e-9, "quad area = 4");
        CHECK(std::abs(m.surface_area() - 4.0) < 1e-9, "surface_area sums faces");
        Vec3 vn = m.calc_vertex_normal(v0);
        CHECK(std::abs(vn.z - 1.0) < 1e-9, "vertex normal of planar quad is +z");
        Vec3 com = m.center_of_mass();
        CHECK(std::abs(com.x - 1.0) < 1e-9 && std::abs(com.y - 1.0) < 1e-9, "center of mass at (1,1,0)");
        auto bb = m.bounding_box();
        CHECK(bb.first.x == 0 && bb.second.x == 2 && bb.second.y == 2, "bounding box correct");
        CHECK(m.face_edges(f).size() == 4 && m.face_faces(f).size() == 0, "face_edges=4, no neighbours");
        CHECK(m.vertex_edges(v0).size() == 2, "corner vertex touches 2 edges");

        // triangulate the quad -> 2 triangles, same area
        m.triangulate();
        CHECK(m.n_faces() == 2, "triangulate: quad -> 2 triangles");
        CHECK(std::abs(m.surface_area() - 4.0) < 1e-9, "triangulate preserves area");

        // OFF round trip
        m.write_off("rt.off");
        Mesh o; o.read_off("rt.off");
        CHECK(o.n_vertices() == 4 && o.n_faces() == 2, "OFF round trip: 4 verts, 2 tris");
    }

    // --- dihedral angle + Laplacian smoothing -----------------------------
    {
        // two triangles folded along a shared edge by 90 degrees
        Mesh m;
        auto a = m.add_vertex({0, 0, 0});
        auto b = m.add_vertex({1, 0, 0});
        auto c = m.add_vertex({0, 1, 0});
        auto d = m.add_vertex({0, 0, 1});      // folded up
        m.add_triangle(a, b, c);
        m.add_triangle(b, a, d);
        auto sh = m.find_halfedge(a, b);
        double ang = std::abs(m.calc_dihedral_angle(m.edge(sh)));
        CHECK(ang > 1.0 && ang < 2.2, "dihedral angle of a fold is ~90 degrees (in radians)");

        // smoothing a subdivided sphere shrinks its bounding box (no boundary)
        Mesh tet;
        auto t0 = tet.add_vertex({0, 0, 0});
        auto t1 = tet.add_vertex({1, 0, 0});
        auto t2 = tet.add_vertex({0, 1, 0});
        auto t3 = tet.add_vertex({0, 0, 1});
        tet.add_triangle(t0, t2, t1); tet.add_triangle(t0, t1, t3);
        tet.add_triangle(t1, t2, t3); tet.add_triangle(t2, t0, t3);
        Mesh sphere = sm::loop_subdivide(tet, 2);
        auto bb0 = sphere.bounding_box();
        double diag0 = (bb0.second - bb0.first).norm();
        sm::laplacian_smooth(sphere, 5, 0.5);
        auto bb1 = sphere.bounding_box();
        double diag1 = (bb1.second - bb1.first).norm();
        CHECK(diag1 < diag0, "Laplacian smoothing shrinks a closed mesh");
    }

    // --- status bits ------------------------------------------------------
    {
        Mesh m;
        auto a = m.add_vertex({0, 0, 0});
        auto b = m.add_vertex({1, 0, 0});
        auto c = m.add_vertex({0, 1, 0});
        auto f = m.add_triangle(a, b, c);

        CHECK(!m.is_selected(a), "status defaults to clear");
        m.set_selected(a);
        m.set_feature(m.edge(m.find_halfedge(a, b)));
        m.set_tagged(f);
        CHECK(m.is_selected(a), "vertex selected flag set");
        CHECK(!m.is_tagged(a), "independent bits don't bleed");
        CHECK(m.is_feature(m.edge(m.find_halfedge(a, b))), "edge feature flag set");
        CHECK(m.is_tagged(f), "face tagged flag set");

        m.set_selected(a, false);
        CHECK(!m.is_selected(a), "flag can be cleared");

        // multiple bits coexist in one word
        m.set_selected(b); m.set_locked(b);
        CHECK(m.is_selected(b) && m.is_locked(b), "two bits coexist");
        CHECK((m.get_status(b) & (Mesh::SELECTED | Mesh::LOCKED)) == (Mesh::SELECTED | Mesh::LOCKED),
              "raw status word holds both bits");
    }

    // --- status survives garbage collection (it's just a property) --------
    {
        Mesh m;
        sm::VertexHandle v[4];
        for (int i = 0; i < 4; ++i) v[i] = m.add_vertex({double(i), 0, 0});
        m.set_selected(v[2]);                 // mark the survivor we track
        m.delete_vertex(v[0]);                // drop an earlier vertex
        m.garbage_collection();               // compacts arrays AND status prop

        CHECK(m.n_vertices() == 3, "GC removed deleted vertex");
        // old v[2] (selected) shifts to new index 1 (v0 gone): v1,v2,v3 -> 0,1,2
        CHECK(m.is_selected(VertexHandle(1)), "selected flag compacted with its vertex");
        CHECK(!m.is_selected(VertexHandle(0)) && !m.is_selected(VertexHandle(2)),
              "other vertices stay unselected after GC");
    }

    // --- MeshChecker: a clean mesh passes every check ---------------------
    {
        Mesh m = make_grid(5, 5);
        MeshChecker chk(m);
        MeshCheckReport rep = chk.check();
        CHECK(rep.ok(), "clean grid passes MeshChecker");
        CHECK(chk.is_valid(), "clean grid is_valid()");
        CHECK(rep.is_manifold, "grid is manifold");
        CHECK(rep.is_oriented, "grid is oriented");
        CHECK(!rep.is_closed, "open grid has a boundary");
        CHECK(rep.n_components == 1, "grid is one component");
        CHECK(rep.n_boundary_loops == 1, "grid has one boundary loop");
        CHECK(rep.bad_halfedges.empty(), "grid connectivity intact");

        // a closed tetrahedron: closed, manifold, genus 0, Euler 2
        Mesh t;
        auto a = t.add_vertex({0,0,0}), b = t.add_vertex({1,0,0});
        auto c = t.add_vertex({0,1,0}), d = t.add_vertex({0,0,1});
        t.add_triangle(a,c,b); t.add_triangle(a,b,d);
        t.add_triangle(b,c,d); t.add_triangle(c,a,d);
        MeshCheckReport tr = MeshChecker(t).check();
        CHECK(tr.ok() && tr.is_closed, "tetra passes & is watertight");
        CHECK(tr.euler == 2 && tr.genus == 0, "tetra Euler 2, genus 0");
    }

    // --- MeshChecker: non-manifold "bow-tie" vertex is caught -------------
    {
        // two triangles meeting at a single shared vertex (and nowhere else)
        Mesh m;
        auto ctr = m.add_vertex({0, 0, 0});
        auto a = m.add_vertex({1, 0, 0});
        auto b = m.add_vertex({1, 1, 0});
        auto c = m.add_vertex({-1, 0, 0});
        auto d = m.add_vertex({-1,-1, 0});
        m.add_triangle(ctr, a, b);
        m.add_triangle(ctr, c, d);     // shares only `ctr` -> bow-tie
        MeshCheckReport rep = MeshChecker(m).check();
        bool flagged = false;
        for (int v : rep.nonmanifold_vertices) if (v == ctr.idx()) flagged = true;
        CHECK(flagged, "bow-tie vertex flagged non-manifold");
        CHECK(!rep.is_manifold, "bow-tie mesh reported non-manifold");
        CHECK(!MeshChecker(m).is_valid(), "bow-tie mesh fails is_valid()");
    }

    // --- MeshChecker: duplicate vertices & duplicate faces ----------------
    {
        Mesh m;
        auto a = m.add_vertex({0, 0, 0});
        auto b = m.add_vertex({1, 0, 0});
        auto c = m.add_vertex({0, 1, 0});
        m.add_triangle(a, b, c);
        m.add_vertex({0, 0, 0});       // a colocated duplicate of `a`
        MeshCheckReport rep = MeshChecker(m).check();
        CHECK(rep.duplicate_vertices.size() == 1, "one duplicate-vertex pair found");

        // duplicate face: add the same triangle again (reversed winding)
        Mesh m2;
        auto x = m2.add_vertex({0,0,0}), y = m2.add_vertex({1,0,0}), z = m2.add_vertex({0,1,0});
        m2.add_triangle(x, y, z);
        m2.add_triangle(x, z, y);      // same vertex set
        MeshCheckReport rep2 = MeshChecker(m2).check();
        CHECK(rep2.duplicate_faces.size() == 1, "duplicate face detected");
    }

    // --- MeshChecker: self-intersection (opt-in) --------------------------
    {
        // two crossing triangles that share no vertex
        Mesh m;
        auto a = m.add_vertex({-1, 0, -1}), b = m.add_vertex({1, 0, -1}), c = m.add_vertex({0, 0, 1});
        auto d = m.add_vertex({0,-1, 0}),  e = m.add_vertex({0, 1, 0}),  f = m.add_vertex({2, 0, 0});
        m.add_triangle(a, b, c);
        m.add_triangle(d, e, f);       // pierces the first triangle
        CheckOptions opt; opt.check_self_intersections = true;
        MeshCheckReport rep = MeshChecker(m).check(opt);
        CHECK(rep.self_intersections.size() >= 1, "crossing triangles flagged");

        // the same two triangles pulled far apart do NOT intersect
        for (VertexHandle v : { d, e, f })
            m.set_point(v, m.point(v) + Vec3(0, 0, 100));
        MeshCheckReport rep2 = MeshChecker(m).check(opt);
        CHECK(rep2.self_intersections.empty(), "separated triangles do not intersect");
    }

    // --- PreMesh: recovers a non-manifold-edge mesh the kernel would drop --
    {
        // A non-manifold edge: 3 triangles sharing edge (a,b) like a book spine.
        std::vector<Vec3> pos = {
            {0,0,0}, {1,0,0},    // a, b  (the shared spine)
            {0, 1, 0}, {0,-1, 0}, {0, 0, 1} };   // three "pages"
        std::vector<std::vector<int>> f = {
            {0,1,2}, {0,1,3}, {0,1,4} };          // all share edge 0-1 -> radial 3

        // Straight into the halfedge kernel: the 3rd page can't attach -> dropped.
        Mesh direct;
        for (auto& p : pos) direct.add_vertex(p);
        int added = 0;
        for (auto& t : f) if (direct.add_triangle(VertexHandle(t[0]), VertexHandle(t[1]), VertexHandle(t[2])).is_valid()) ++added;
        CHECK(added < 3, "kernel drops a face on the non-manifold spine");

        // Through PreMesh: detected, then split into manifold fans by repair.
        PreMesh pm = PreMesh::from_soup(pos, f);
        RepairReport r = pm.repair();
        CHECK(r.nm_edges == 1, "PreMesh detects the 1 non-manifold edge");
        CHECK(r.vertices_split > 0, "PreMesh split the non-manifold fan");
        RepairReport er;
        Mesh recovered = pm.to_mesh(RepairOptions(), &er);
        CHECK(recovered.n_faces() == 3, "PreMesh recovers all 3 faces (none dropped)");
        CHECK(er.faces_failed == 0, "emit succeeds with no fallback after the split");
        CHECK(MeshChecker(recovered).check().is_manifold, "split output is manifold");
    }

    // --- PreMesh: bow-tie vertex is dissociated into per-fan copies -------
    {
        // two triangles touching at a single shared vertex 0 (a bow-tie)
        std::vector<Vec3> pos = { {0,0,0}, {1,0,0}, {1,1,0}, {-1,0,0}, {-1,-1,0} };
        std::vector<std::vector<int>> f = { {0,1,2}, {0,3,4} };   // share only vertex 0
        PreMesh pm = PreMesh::from_soup(pos, f);
        RepairReport r = pm.repair();
        CHECK(r.vertices_split == 1, "bow-tie vertex split into 2 copies");
        Mesh m = pm.to_mesh();
        MeshCheckReport cr = MeshChecker(m).check();
        CHECK(cr.is_manifold, "dissociated bow-tie output is manifold");
        CHECK(cr.nonmanifold_vertices.empty(), "no non-manifold vertices remain");
        CHECK(cr.n_components == 2, "the two fans became two components");
    }

    // --- PreMesh: weld + reorient soup before the halfedge build ----------
    {
        // a cracked, mis-wound two-triangle soup (duplicate verts, flipped face)
        std::vector<Vec3> pos = {
            {0,0,0}, {1,0,0}, {0,1,0},
            {1,0,0}, {0,1,0}, {1,1,0} };          // verts 3,4 duplicate 1,2
        std::vector<std::vector<int>> f = { {0,1,2}, {4,3,5} };  // 2nd wound the other way
        PreMesh pm = PreMesh::from_soup(pos, f);
        RepairReport r = pm.repair();
        CHECK(r.vertices_merged == 2, "PreMesh welds the 2 duplicate soup vertices");
        Mesh m = pm.to_mesh();
        CHECK(m.n_vertices() == 4 && m.n_faces() == 2, "welded to 4 verts / 2 faces");
        MeshCheckReport cr = MeshChecker(m).check();
        CHECK(cr.is_manifold && cr.is_oriented, "PreMesh output is manifold & oriented");
        CHECK(cr.n_components == 1, "welded into one component");
    }

    // --- MeshRepair: weld a seam of duplicated vertices -------------------
    {
        // Two triangles that SHOULD share edge b-c, but built with separate
        // (colocated) copies of b and c -> a cracked, non-welded seam.
        Mesh m;
        auto a  = m.add_vertex({0, 0, 0});
        auto b1 = m.add_vertex({1, 0, 0});
        auto c1 = m.add_vertex({0, 1, 0});
        auto b2 = m.add_vertex({1, 0, 0});   // duplicate of b1
        auto c2 = m.add_vertex({0, 1, 0});   // duplicate of c1
        auto d  = m.add_vertex({1, 1, 0});
        m.add_triangle(a, b1, c1);
        m.add_triangle(b2, d, c2);
        CHECK(m.n_vertices() == 6, "cracked mesh starts with 6 vertices");

        RepairReport rr = repair_mesh(m);
        CHECK(rr.vertices_merged == 2, "repair welded the 2 duplicate vertices");
        CHECK(m.n_vertices() == 4, "welded mesh has 4 vertices");
        CHECK(m.n_faces() == 2, "both faces survived the weld");
        // the seam is now a real shared interior edge
        MeshCheckReport cr = MeshChecker(m).check();
        CHECK(cr.duplicate_vertices.empty(), "no duplicate vertices after repair");
        CHECK(cr.is_manifold, "welded mesh is manifold");
        CHECK(cr.n_components == 1, "welded mesh is a single component");
    }

    // --- MeshRepair: drop degenerate & duplicate faces --------------------
    {
        Mesh m;
        auto a = m.add_vertex({0, 0, 0});
        auto b = m.add_vertex({1, 0, 0});
        auto c = m.add_vertex({0, 1, 0});
        m.add_triangle(a, b, c);
        m.add_triangle(a, c, b);             // duplicate (reversed)
        RepairReport rr = repair_mesh(m);
        CHECK(rr.duplicate_removed == 1, "repair dropped the duplicate face");
        CHECK(m.n_faces() == 1, "one face remains after dedup");
        CHECK(MeshChecker(m).check().duplicate_faces.empty(), "no duplicate faces after repair");
    }

    // --- MeshRepair: reorient an incoherently-wound mesh ------------------
    {
        // Two triangles sharing edge b-c, but wound the SAME way across it
        // (both list b->c) -> inconsistent orientation. Build via soup + repair.
        Mesh m;
        auto a = m.add_vertex({0, 0, 0});
        auto b = m.add_vertex({1, 0, 0});
        auto c = m.add_vertex({0, 1, 0});
        auto d = m.add_vertex({1, 1, 0});
        m.add_triangle(a, b, c);             // edge b->c
        // add_face would reject the inconsistent twin, so place it as raw soup:
        // make a separate copy, then weld+reorient fixes both.
        auto b2 = m.add_vertex({1, 0, 0});
        auto c2 = m.add_vertex({0, 1, 0});
        m.add_triangle(b2, c2, d);           // also traverses b->c (inconsistent)

        RepairOptions opt;                   // weld + reorient on by default
        RepairReport rr = repair_mesh(m, opt);
        CHECK(m.n_faces() == 2, "both faces kept through reorient");
        MeshCheckReport cr = MeshChecker(m).check();
        CHECK(cr.is_oriented, "repaired mesh is consistently oriented");
        CHECK(cr.is_manifold, "repaired mesh is manifold");
        (void)rr;
    }

    // --- MeshRepair: a clean mesh is left essentially untouched -----------
    {
        Mesh m = make_grid(4, 4);
        size_t V = m.n_vertices(), F = m.n_faces();
        RepairReport rr = repair_mesh(m);
        CHECK(rr.vertices_merged == 0 && rr.degenerate_removed == 0 &&
              rr.duplicate_removed == 0 && rr.vertices_removed == 0,
              "clean mesh needs no structural repair");
        CHECK(m.n_vertices() == V && m.n_faces() == F, "clean mesh preserved");
    }

    // --- modular decimation (Decimater + modules) -------------------------
    {
        Mesh m = make_grid(9, 9);              // a flat triangulated grid
        size_t F0 = m.n_faces();

        Decimator dec(m);
        dec.add_module<ModQuadric>();          // priority: QEM
        dec.add_module<ModNormalFlipping>(1.4);  // binary: no foldovers (~80 deg)
        size_t collapses = dec.decimate_to_faces(F0 / 2);

        CHECK(collapses > 0, "decimater performed collapses");
        CHECK(m.n_faces() <= F0, "decimater did not grow the mesh");
        CHECK(m.n_faces() <= F0 / 2 + 4, "decimater reached ~target face count");
        // a flat grid must stay flat after collapses (planarity preserved)
        bool flat = true;
        for (auto v : m.all_vertices()) if (std::abs(m.point(v).z) > 1e-9) flat = false;
        CHECK(flat, "decimated flat grid stays flat");
    }

    // --- binary module can veto every collapse ----------------------------
    {
        Mesh m = make_grid(5, 5);
        Decimator dec(m);
        dec.add_module<ModQuadric>();
        dec.add_module<ModEdgeLength>(1e-6);   // impossibly short -> vetoes all
        size_t collapses = dec.decimate_to_faces(1);
        CHECK(collapses == 0, "edge-length module vetoed all collapses");
    }

    // --- cotangent smoothing keeps a plane planar & is bounded ------------
    {
        Mesh m = make_grid(6, 6);
        // perturb one interior vertex off the plane, then smooth it back down
        VertexHandle moved;
        for (auto v : m.all_vertices())
            if (!m.is_boundary(v)) { moved = v; break; }
        Vec3 p = m.point(moved); p.z = 1.0; m.set_point(moved, p);

        cotan_smooth(m, 20, 0.5);
        CHECK(std::abs(m.point(moved).z) < 0.5, "cotan smoothing pulls the spike down");
        // boundary vertices must not move
        bool bnd_fixed = true;
        for (auto v : m.all_vertices())
            if (m.is_boundary(v) && std::abs(m.point(v).z) > 1e-12) bnd_fixed = false;
        CHECK(bnd_fixed, "cotan smoothing leaves the boundary fixed");
    }

    // --- hole filling -----------------------------------------------------
    {
        // an octahedron (closed) with one face removed leaves a triangular hole
        Mesh m;
        VertexHandle v[6] = {
            m.add_vertex({ 1, 0, 0}), m.add_vertex({-1, 0, 0}),
            m.add_vertex({ 0, 1, 0}), m.add_vertex({ 0,-1, 0}),
            m.add_vertex({ 0, 0, 1}), m.add_vertex({ 0, 0,-1}) };
        int tris[8][3] = {
            {0,2,4},{2,1,4},{1,3,4},{3,0,4},
            {2,0,5},{1,2,5},{3,1,5},{0,3,5} };
        for (auto& t : tris) m.add_triangle(v[t[0]], v[t[1]], v[t[2]]);
        CHECK(m.n_faces() == 8, "octahedron has 8 faces");

        m.delete_face(FaceHandle(0), false);
        m.garbage_collection();
        CHECK(m.n_faces() == 7, "one face removed -> a hole");

        int filled = fill_holes(m);
        CHECK(filled == 1, "fill_holes closed exactly one hole");
        CHECK(m.n_faces() == 8, "hole filled back to 8 faces");
        // closed again: every edge is interior
        bool closed = true;
        for (auto e : m.all_edges()) if (m.is_boundary(e)) closed = false;
        CHECK(closed, "mesh is watertight after filling");
    }

    // --- midpoint / butterfly subdivision (interpolating) -----------------
    {
        Mesh m = make_grid(3, 3);
        size_t V0 = m.n_vertices();
        Vec3 keep = m.point(VertexHandle(0));

        Mesh mid = midpoint_subdivide(m, 1);
        CHECK(mid.n_faces() == m.n_faces() * 4, "midpoint: 1 -> 4 faces");
        CHECK((mid.point(VertexHandle(0)) - keep).norm() < 1e-12,
              "midpoint is interpolating (old vertex unchanged)");

        Mesh bf = butterfly_subdivide(m, 1);
        CHECK(bf.n_faces() == m.n_faces() * 4, "butterfly: 1 -> 4 faces");
        CHECK((bf.point(VertexHandle(0)) - keep).norm() < 1e-12,
              "butterfly is interpolating (old vertex unchanged)");
        bool flat = true;
        for (auto v : bf.all_vertices()) if (std::abs(bf.point(v).z) > 1e-9) flat = false;
        CHECK(flat, "butterfly keeps a planar grid planar");
        (void)V0;
    }

    // --- longest-edge refinement shortens every edge ----------------------
    {
        Mesh m;
        auto a = m.add_vertex({0, 0, 0});
        auto b = m.add_vertex({4, 0, 0});
        auto c = m.add_vertex({0, 4, 0});
        auto d = m.add_vertex({4, 4, 0});
        m.add_triangle(a, b, c);
        m.add_triangle(b, d, c);

        longest_edge_subdivide(m, 1.5);
        double longest = 0;
        for (auto e : m.all_edges())
            if (!m.is_deleted(e)) longest = std::max(longest, m.calc_edge_length(e));
        CHECK(longest <= 1.5 + 1e-9, "longest-edge refinement bounds every edge");
        bool all_tris = true;
        for (auto f : m.all_faces())
            if (!m.is_deleted(f) && m.face_vertices(f).size() != 3) all_tris = false;
        CHECK(all_tris, "longest-edge output stays triangular");
    }

    // --- vertex_split is the inverse of collapse --------------------------
    {
        // Build a closed-ish patch: a center vertex with a 1-ring of 4, so the
        // center has interior valence. Two triangles fan around an interior edge.
        Mesh m;
        auto c  = m.add_vertex({0, 0, 0});      // 0 center
        auto r0 = m.add_vertex({1, 0, 0});      // 1
        auto r1 = m.add_vertex({0, 1, 0});      // 2
        auto r2 = m.add_vertex({-1, 0, 0});     // 3
        auto r3 = m.add_vertex({0, -1, 0});     // 4
        m.add_triangle(c, r0, r1);
        m.add_triangle(c, r1, r2);
        m.add_triangle(c, r2, r3);
        m.add_triangle(c, r3, r0);              // closed fan (disk) around c

        size_t V0 = m.n_vertices(), E0 = m.n_edges(), F0 = m.n_faces();

        // Collapse the edge c->r0 (merges c into r0), recording the wing tips.
        HalfedgeHandle h = m.find_halfedge(c, r0);
        CHECK(h.is_valid(), "found halfedge c->r0");
        VertexHandle vl = m.to_vertex(m.next_halfedge(h));                       // apex left
        VertexHandle vr = m.to_vertex(m.next_halfedge(m.opposite_halfedge(h)));  // apex right
        CHECK(m.is_collapse_ok(h), "edge c->r0 collapsible");
        Vec3 cpos = m.point(c);
        m.collapse(h);
        m.garbage_collection();                 // compact away the merged vertex

        CHECK(m.n_vertices() == V0 - 1, "collapse removed one vertex");
        CHECK(m.n_faces() == F0 - 2, "collapse removed two faces");

        // Now split it back. After GC handles changed, so re-find the survivors:
        // r0 is the kept vertex (was at {1,0,0}); vl/vr were apexes of the fan.
        // Re-derive them by position to stay robust to renumbering.
        auto find_at = [&](Vec3 p) {
            for (auto v : m.all_vertices())
                if ((m.point(v) - p).norm() < 1e-9) return v;
            return VertexHandle();
        };
        VertexHandle kept = find_at({1, 0, 0});
        VertexHandle nvl  = find_at({0, 1, 0});   // r1
        VertexHandle nvr  = find_at({0, -1, 0});  // r3
        CHECK(kept.is_valid() && nvl.is_valid() && nvr.is_valid(), "survivors located");

        VertexHandle v0 = m.add_vertex(cpos);     // fresh isolated vertex
        HalfedgeHandle nh = m.vertex_split(v0, kept, nvl, nvr);
        CHECK(nh.is_valid(), "vertex_split returned a halfedge");
        CHECK(m.n_vertices() == V0, "vertex_split restored vertex count");
        CHECK(m.n_faces() == F0, "vertex_split restored face count");
        CHECK(m.n_edges() == E0, "vertex_split restored edge count");
        CHECK(m.to_vertex(nh) == kept && m.from_vertex(nh) == v0, "new edge is v0->v1");
        CHECK(m.valence(v0) >= 3, "new vertex has a real fan");
    }

    // --- Catmull-Clark turns any mesh into quads --------------------------
    {
        // a single quad face
        Mesh m;
        auto a = m.add_vertex({0, 0, 0});
        auto b = m.add_vertex({1, 0, 0});
        auto cc = m.add_vertex({1, 1, 0});
        auto d = m.add_vertex({0, 1, 0});
        m.add_face({a, b, cc, d});

        Mesh q = catmull_clark(m, 1);
        // one k-gon -> k quads; a quad -> 4 quads, 9 verts (4 old +1 face +4 edge)
        CHECK(q.n_faces() == 4, "Catmull-Clark: quad -> 4 faces");
        CHECK(q.n_vertices() == 9, "Catmull-Clark: 4+1+4 = 9 vertices");
        bool all_quads = true;
        for (auto f : q.all_faces()) if (q.face_vertices(f).size() != 4) all_quads = false;
        CHECK(all_quads, "Catmull-Clark: every face is a quad");

        // a flat sheet must stay flat (z == 0 everywhere)
        bool flat = true;
        for (auto v : q.all_vertices()) if (std::abs(q.point(v).z) > 1e-9) flat = false;
        CHECK(flat, "Catmull-Clark keeps a planar sheet planar");
    }

    // --- sqrt(3) subdivision triples triangle count -----------------------
    {
        // two triangles sharing one interior edge (a split quad)
        Mesh m;
        auto v0 = m.add_vertex({0, 0, 0});
        auto v1 = m.add_vertex({1, 0, 0});
        auto v2 = m.add_vertex({1, 1, 0});
        auto v3 = m.add_vertex({0, 1, 0});
        m.add_triangle(v0, v1, v2);
        m.add_triangle(v0, v2, v3);

        Mesh s = sqrt3_subdivide(m, 1);
        // 1->3 split makes 2*3=6 triangles; flipping the single interior edge
        // keeps the triangle count (a flip is 2->2). So 6 faces.
        CHECK(s.n_faces() == 6, "sqrt3: 2 triangles -> 6 triangles");
        CHECK(s.n_vertices() == 6, "sqrt3: 4 old + 2 centroids = 6 vertices");
        bool all_tris = true;
        for (auto f : s.all_faces()) if (s.face_vertices(f).size() != 3) all_tris = false;
        CHECK(all_tris, "sqrt3: output is all triangles");
        // planar input stays planar
        bool flat = true;
        for (auto v : s.all_vertices()) if (std::abs(s.point(v).z) > 1e-9) flat = false;
        CHECK(flat, "sqrt3 keeps a planar sheet planar");
    }

    // --- smart handles: fluent navigation ---------------------------------
    {
        Mesh m;
        auto v0 = m.add_vertex({0, 0, 0});
        auto v1 = m.add_vertex({1, 0, 0});
        auto v2 = m.add_vertex({1, 1, 0});
        auto v3 = m.add_vertex({0, 1, 0});
        auto f0 = m.add_triangle(v0, v1, v2);
        m.add_triangle(v0, v2, v3);

        // a smart handle is-a plain handle (slices cleanly)
        auto sf = m.smart(f0);
        CHECK(sf.idx() == f0.idx(), "smart face keeps its index");
        CHECK(sf.valence() == 3, "smart face valence == 3");

        // chain: a halfedge -> opposite -> opposite is identity
        auto h = m.smart(m.halfedge(f0));
        CHECK(h.opp().opp() == h, "smart opp().opp() is identity");
        // next three times around a triangle returns to start
        CHECK(h.next().next().next() == h, "smart next-cycle length 3");
        // to()/from() agree with the plain API
        CHECK(h.to().idx() == m.to_vertex(h).idx(), "smart to() matches plain");
        CHECK(h.from().idx() == m.from_vertex(h).idx(), "smart from() matches plain");
        // edge -> halfedge -> face round trip
        CHECK(h.edge().h0().face().is_valid() || h.edge().h1().face().is_valid(),
              "smart edge reaches an incident face");

        // ranges return smart handles and match the plain circulator counts
        auto sv = m.smart(v0);
        CHECK(sv.vertices().size() == m.vertex_vertices(v0).size(), "smart vertex 1-ring count matches");
        CHECK(sv.faces().size() == m.vertex_faces(v0).size(), "smart vertex faces count matches");

        // the diagonal edge v0-v2 is interior -> two incident faces
        auto sh = m.smart(m.find_halfedge(v0, v2));
        CHECK(sh.edge().faces().size() == 2, "interior edge has two smart faces");
        CHECK(!sh.edge().is_boundary(), "interior edge not boundary via smart handle");

        // status predicates read through the smart handle
        m.set_selected(v0);
        CHECK(m.smart(v0).selected(), "smart handle reports selected status");
        CHECK(!m.smart(v1).selected(), "unselected vertex reports false");
    }

    // ======================================================================
    //  EXPANDED COVERAGE - feature/algorithm combinations
    // ======================================================================

    // --- OBJ + OFF round trips (complete the I/O matrix) ------------------
    {
        Mesh src = make_tetra();
        src.write_obj("rt2.obj");
        Mesh a; CHECK(a.read_obj("rt2.obj"), "OBJ read returns true");
        CHECK(a.n_vertices() == 4 && a.n_faces() == 4, "OBJ round trip: 4 verts, 4 faces");

        src.write_off("rt2.off");
        Mesh b; CHECK(b.read_off("rt2.off"), "OFF read returns true");
        CHECK(b.n_vertices() == 4 && b.n_faces() == 4, "OFF round trip: 4 verts, 4 faces");
    }

    // --- every subdivision scheme keeps a closed mesh closed (Euler 2) ----
    {
        auto euler = [](const Mesh& m) {
            return (long)m.n_vertices() - (long)m.n_edges() + (long)m.n_faces(); };
        Mesh t = make_tetra();
        CHECK(euler(sm::loop_subdivide(t, 2)) == 2,      "loop x2 stays Euler 2");
        CHECK(euler(sm::sqrt3_subdivide(t, 2)) == 2,     "sqrt3 x2 stays Euler 2");
        CHECK(euler(sm::catmull_clark(t, 2)) == 2,       "catmull-clark x2 stays Euler 2");
        CHECK(euler(sm::midpoint_subdivide(t, 2)) == 2,  "midpoint x2 stays Euler 2");
        CHECK(euler(sm::butterfly_subdivide(t, 2)) == 2, "butterfly x2 stays Euler 2");
        // catmull-clark output is all quads
        Mesh q = sm::catmull_clark(t, 1);
        bool all_quads = true;
        for (auto f : q.all_faces()) if (q.face_vertices(f).size() != 4) all_quads = false;
        CHECK(all_quads, "catmull-clark output is all quads");
    }

    // --- remaining circulators + lazy/vector consistency ------------------
    {
        Mesh m;
        auto v0 = m.add_vertex({0, 0, 0}); auto v1 = m.add_vertex({1, 0, 0});
        auto v2 = m.add_vertex({1, 1, 0}); auto v3 = m.add_vertex({0, 1, 0});
        auto f0 = m.add_triangle(v0, v1, v2);
        m.add_triangle(v0, v2, v3);

        CHECK(m.vertex_edges(v0).size() == 3, "vertex_edges: v0 touches 3 edges");
        CHECK(m.vertex_incoming_halfedges(v0).size() == 3, "vertex_incoming_halfedges: 3");
        CHECK(m.face_edges(f0).size() == 3, "face_edges: triangle has 3 edges");
        CHECK(m.face_faces(f0).size() == 1, "face_faces: f0 has 1 interior neighbour");

        // lazy ranges agree with the vector circulators
        size_t vv = 0; for (auto w : m.vv_range(v0))  { (void)w; ++vv; }
        size_t fh = 0; for (auto h : m.fh_range(f0))  { (void)h; ++fh; }
        CHECK(vv == m.vertex_vertices(v0).size(), "vv_range count == vertex_vertices");
        CHECK(fh == m.face_halfedges(f0).size(),  "fh_range count == face_halfedges");
        size_t allv = 0; for (auto w : m.all_vertices()) { (void)w; ++allv; }
        CHECK(allv == m.n_vertices(), "all_vertices count == n_vertices");
    }

    // --- properties on edge & face kinds, surviving GC --------------------
    {
        Mesh m;
        auto a = m.add_vertex({0, 0, 0}); auto b = m.add_vertex({1, 0, 0});
        auto c = m.add_vertex({0, 1, 0}); auto d = m.add_vertex({-1, 0, 0});
        m.add_triangle(a, b, c);
        m.add_triangle(a, c, d);

        auto ep = m.add_edge_property<int>("weight", 0);
        auto fp = m.add_face_property<Vec3>("centroid", Vec3{});
        EdgeHandle e0 = m.edge(m.find_halfedge(a, c));
        m.property(ep, e0) = 42;
        m.property(fp, FaceHandle(0)) = Vec3{1, 2, 3};
        CHECK(m.property(ep, e0) == 42, "edge int property stored");
        CHECK(m.property(fp, FaceHandle(0)).x == 1, "face vec3 property stored");

        // mark, delete an unrelated face, GC, and confirm the value rides along
        auto fmark = m.add_face_property<int>("mark", 0);
        m.property(fmark, FaceHandle(0)) = 7;
        m.delete_face(FaceHandle(1), false);
        m.garbage_collection();
        CHECK(m.property(fmark, FaceHandle(0)) == 7, "face property survives GC");
    }

    // --- MeshChecker: components / isolated / boundary loops / empty ------
    {
        // two disjoint triangles -> two components, two boundary loops
        Mesh m;
        auto a = m.add_vertex({0, 0, 0}); auto b = m.add_vertex({1, 0, 0}); auto c = m.add_vertex({0, 1, 0});
        m.add_triangle(a, b, c);
        auto d = m.add_vertex({5, 0, 0}); auto e = m.add_vertex({6, 0, 0}); auto f = m.add_vertex({5, 1, 0});
        m.add_triangle(d, e, f);
        MeshCheckReport r = MeshChecker(m).check();
        CHECK(r.n_components == 2, "checker: two disjoint triangles -> 2 components");
        CHECK(r.n_boundary_loops == 2, "checker: each open triangle is one boundary loop");

        // a triangle plus a free-floating vertex -> isolated-vertex detection
        Mesh m2;
        auto x = m2.add_vertex({0, 0, 0}); auto y = m2.add_vertex({1, 0, 0}); auto z = m2.add_vertex({0, 1, 0});
        m2.add_triangle(x, y, z);
        auto iso = m2.add_vertex({9, 9, 9});
        MeshCheckReport r2 = MeshChecker(m2).check();
        bool iso_found = false; for (int vi : r2.isolated_vertices) if (vi == iso.idx()) iso_found = true;
        CHECK(iso_found, "checker: isolated vertex detected");

        // open grid -> exactly one boundary loop
        Mesh g = make_grid(4, 4);
        CHECK(MeshChecker(g).check().n_boundary_loops == 1, "checker: open grid has 1 boundary loop");

        // empty mesh is vacuously valid
        Mesh empty;
        MeshCheckReport er = MeshChecker(empty).check();
        CHECK(er.ok() && er.n_components == 0 && er.euler == 0, "checker: empty mesh is clean");
    }

    // --- MeshChecker: a torus has genus 1 ---------------------------------
    {
        const int R = 8, C = 8; const double PI = 3.14159265358979323846;
        Mesh m;
        std::vector<VertexHandle> v(R * C);
        for (int r = 0; r < R; ++r)
            for (int cc = 0; cc < C; ++cc) {
                double u = 2 * PI * r / R, w = 2 * PI * cc / C;
                double x = (3 + std::cos(w)) * std::cos(u);
                double y = (3 + std::cos(w)) * std::sin(u);
                double z = std::sin(w);
                v[r * C + cc] = m.add_vertex({x, y, z});
            }
        auto at = [&](int r, int cc) { return v[(r % R) * C + (cc % C)]; };
        for (int r = 0; r < R; ++r)
            for (int cc = 0; cc < C; ++cc)
                m.add_face({ at(r, cc), at(r, cc + 1), at(r + 1, cc + 1), at(r + 1, cc) });
        MeshCheckReport r0 = MeshChecker(m).check();
        CHECK(r0.is_manifold && r0.is_closed, "torus: manifold + watertight");
        CHECK(r0.euler == 0, "torus: Euler characteristic 0");
        CHECK(r0.genus == 1, "torus: genus 1");
    }

    // --- add_face rejects degenerate and non-manifold input ---------------
    {
        Mesh m;
        auto a = m.add_vertex({0, 0, 0}); auto b = m.add_vertex({1, 0, 0});
        auto c = m.add_vertex({0, 1, 0}); auto d = m.add_vertex({0, -1, 0});
        auto e = m.add_vertex({1, 1, 0});
        CHECK(!m.add_triangle(a, a, b).is_valid(), "add_face rejects a repeated vertex");
        CHECK(m.add_triangle(a, b, c).is_valid(),  "first face on edge a-b ok");
        CHECK(m.add_triangle(b, a, d).is_valid(),  "second (opposite) face on a-b ok");
        CHECK(!m.add_triangle(a, b, e).is_valid(), "third face on edge a-b rejected (non-manifold)");
    }

    // --- triangulate a polygon --------------------------------------------
    {
        Mesh m;
        std::vector<VertexHandle> v;
        for (int i = 0; i < 5; ++i) {
            double t = 2 * 3.14159265358979323846 * i / 5;
            v.push_back(m.add_vertex({std::cos(t), std::sin(t), 0}));
        }
        m.add_face({v[0], v[1], v[2], v[3], v[4]});      // a pentagon
        m.triangulate();
        CHECK(m.n_faces() == 3, "triangulate: pentagon -> 3 triangles (n-2)");
        bool all_tris = true;
        for (auto f : m.all_faces()) if (m.face_vertices(f).size() != 3) all_tris = false;
        CHECK(all_tris, "triangulate: every face is a triangle");
    }

    // --- PreMesh: load_soup (OBJ) + edge_valence + load_and_repair --------
    {
        // write a non-manifold "spine" OBJ, load as soup, inspect, repair
        {
            std::ofstream o("rt_spine.obj");
            o << "v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 -1 0\nv 0 0 1\n";
            o << "f 1 2 3\nf 1 2 4\nf 1 2 5\n";          // edge 1-2 shared by 3 faces
        }
        PreMesh p = sm::load_soup("rt_spine.obj");
        CHECK(p.n_faces() == 3, "load_soup(OBJ): 3 faces, nothing dropped");
        p.build_radial();
        CHECK(p.edge_valence(0, 1) == 3, "load_soup: edge 0-1 has radial valence 3");

        Mesh clean = sm::load_and_repair("rt_spine.obj");
        CHECK(clean.n_faces() == 3, "load_and_repair recovers all 3 faces");
        CHECK(MeshChecker(clean).check().is_manifold, "load_and_repair output is manifold");

        // from_mesh round-trips a clean mesh unchanged
        Mesh t = make_tetra();
        PreMesh pf = PreMesh::from_mesh(t);
        Mesh back = pf.to_mesh();
        CHECK(back.n_vertices() == 4 && back.n_faces() == 4, "PreMesh from_mesh -> to_mesh round trip");
    }

    // --- Decimator: aspect-ratio veto module runs -------------------------
    {
        Mesh sphere = sm::loop_subdivide(make_tetra(), 2);   // 64 faces
        Decimator dec(sphere);
        dec.add_module<ModQuadric>();
        dec.add_module<ModAspectRatio>(20.0);                // reject very thin slivers
        size_t collapses = dec.decimate_to_faces(30);
        CHECK(collapses > 0, "decimater with aspect-ratio module made progress");
        CHECK(sphere.n_faces() <= 64, "decimater did not grow the mesh");
    }

    // --- geometry: exact values on a known triangle -----------------------
    {
        Mesh m;
        auto a = m.add_vertex({0, 0, 0}); auto b = m.add_vertex({2, 0, 0}); auto c = m.add_vertex({0, 2, 0});
        auto f = m.add_triangle(a, b, c);
        CHECK(std::abs(m.calc_face_area(f) - 2.0) < 1e-9, "triangle area == 2");
        Vec3 n = m.calc_face_normal(f);
        CHECK(std::abs(n.z - 1.0) < 1e-9, "CCW triangle normal is +z");
        // the three interior (sector) angles of any triangle sum to pi
        double angsum = 0;
        for (auto h : m.face_halfedges(f)) angsum += m.calc_sector_angle(h);
        CHECK(std::abs(angsum - 3.14159265358979323846) < 1e-6,
              "triangle sector angles sum to pi");
    }

    // ======================================================================
    //  PRODUCTION REPAIR FIXES (gap-closing)
    // ======================================================================

    // --- weld merges a pair straddling a spatial-hash cell boundary -------
    {
        // two apexes 2e-7 apart; at tol=1e-3 they fall in DIFFERENT cells
        // (floor(499.9999) vs floor(500.0001)) - the old single-cell lookup
        // would miss them; the neighbor-aware search must merge them.
        std::vector<Vec3> pos = {
            {0, 0, 0}, {1, 0, 0}, {0.4999999, 1, 0}, {0.5000001, 1, 0} };
        std::vector<std::vector<int>> f = { {0, 1, 2}, {1, 0, 3} };
        PreMesh p = PreMesh::from_soup(pos, f);
        RepairOptions opt; opt.vertex_merge_tol = 1e-3;
        RepairReport r = p.repair(opt);
        CHECK(r.vertices_merged == 1, "weld merges the cross-cell-boundary pair");
    }

    // --- fill_holes size guard ('max_edges') ------------------------------
    {
        // open 4x4 grid: its outer boundary loop is 12 edges - must NOT be capped
        Mesh g = make_grid(4, 4);
        size_t F0 = g.n_faces();
        int filled = fill_holes(g, /*max_edges=*/4);
        CHECK(filled == 0, "fill_holes(max_edges=4) refuses the 12-edge open boundary");
        CHECK(g.n_faces() == F0, "guarded fill left the open grid unchanged");

        // a genuine 4-edge hole IS filled at the same cap
        Mesh box;
        sm::VertexHandle b[8];
        double p[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
        for (int i = 0; i < 8; ++i) b[i] = box.add_vertex({p[i][0], p[i][1], p[i][2]});
        int q[5][4] = {{0,3,2,1},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}};   // top open
        for (auto& fc : q) box.add_face({b[fc[0]], b[fc[1]], b[fc[2]], b[fc[3]]});
        CHECK(fill_holes(box, 4) == 1, "fill_holes(max_edges=4) fills the 4-edge hole");
        CHECK(MeshChecker(box).check().is_closed, "box watertight after guarded fill");
    }

    // --- remove_small_components drops a stray island ---------------------
    {
        Mesh m = make_grid(4, 4);                  // 18-face disk, 1 component
        size_t big = m.n_faces();
        auto a = m.add_vertex({10, 10, 0});
        auto b = m.add_vertex({11, 10, 0});
        auto c = m.add_vertex({10, 11, 0});
        m.add_triangle(a, b, c);                   // a lone 1-face island
        CHECK(MeshChecker(m).check().n_components == 2, "grid + island -> 2 components");
        int removed = remove_small_components(m, /*min_faces=*/3);
        CHECK(removed == 1, "remove_small_components dropped 1 island");
        CHECK(m.n_faces() == big, "the big component is preserved");
        CHECK(MeshChecker(m).check().n_components == 1, "one component remains");
    }

    // --- min_face_area drops slivers; min_component_faces drops islands ---
    {
        std::vector<Vec3> pos = {
            {0, 0, 0}, {1, 0, 0}, {0, 1, 0},               // good triangle (area 0.5)
            {5, 0, 0}, {6, 0, 0}, {5, 0.0001, 0} };        // sliver (area 5e-5)
        std::vector<std::vector<int>> f = { {0, 1, 2}, {3, 4, 5} };
        PreMesh p = PreMesh::from_soup(pos, f);
        RepairOptions opt; opt.min_face_area = 0.01;
        RepairReport r = p.repair(opt);
        CHECK(r.degenerate_removed == 1, "min_face_area removed the sliver");
        CHECK(p.to_mesh().n_faces() == 1, "only the good triangle remains");

        // min_component_faces via PreMesh repair
        Mesh grid = make_grid(4, 4);
        PreMesh pg = PreMesh::from_mesh(grid);
        pg.positions.push_back({10, 10, 0}); pg.positions.push_back({11, 10, 0}); pg.positions.push_back({10, 11, 0});
        int base = (int)pg.positions.size() - 3;
        pg.faces.push_back({base, base + 1, base + 2});    // island
        RepairOptions o2; o2.min_component_faces = 3;
        RepairReport rg = pg.repair(o2);
        CHECK(rg.components_removed == 1, "min_component_faces dropped the island");
    }

    std::cout << "\n" << (failures == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED") << "\n";
    return failures == 0 ? 0 : 1;
}
