// ============================================================================
//  SimpleMesh - MeshIO.cpp : PLY and STL readers/writers
// ----------------------------------------------------------------------------
//  Kept separate from Mesh.cpp because file I/O is a self-contained concern with
//  its own helpers (binary scalar reads, vertex welding, ASCII/binary sniffing).
//
//  Supported:
//    * PLY  - read/write, both `ascii` and `binary_little_endian` 1.0, for the
//             standard layout: vertex x,y,z (+ extra scalar props skipped) and
//             a face `vertex_indices` list. Arbitrary polygon faces.
//    * STL  - read (auto-detects ASCII vs binary) and write (ascii or binary).
//             STL is a triangle soup with no shared vertices, so the reader
//             WELDS coincident corners back into shared vertices, and the
//             writer fan-triangulates any non-triangle face.
// ============================================================================
#include "simplemesh/Mesh.h"

#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace sm {

namespace {

// ---- little-endian binary scalar helpers -----------------------------------
template <class T>
void write_le(std::ostream& os, T v) {
    // x86/x64 are little-endian, so a raw write is already LE. Kept explicit for
    // clarity about what the format requires.
    os.write(reinterpret_cast<const char*>(&v), sizeof(T));
}
template <class T>
T read_le(std::istream& is) {
    T v{};
    is.read(reinterpret_cast<char*>(&v), sizeof(T));
    return v;
}

// Read one PLY scalar of the named type and return it as double (for positions)
// or as long (for indices) - we just need a numeric value and the right width.
double read_ply_scalar_ascii(std::istream& is, const std::string& /*type*/) {
    double d; is >> d; return d;
}
double read_ply_scalar_binary(std::istream& is, const std::string& type) {
    if (type == "float"  || type == "float32") return read_le<float>(is);
    if (type == "double" || type == "float64") return read_le<double>(is);
    if (type == "char"   || type == "int8")    return read_le<int8_t>(is);
    if (type == "uchar"  || type == "uint8")   return read_le<uint8_t>(is);
    if (type == "short"  || type == "int16")   return read_le<int16_t>(is);
    if (type == "ushort" || type == "uint16")  return read_le<uint16_t>(is);
    if (type == "int"    || type == "int32")   return read_le<int32_t>(is);
    if (type == "uint"   || type == "uint32")  return read_le<uint32_t>(is);
    return read_le<float>(is);  // sensible fallback
}

} // anonymous namespace

// ============================================================================
//  PLY
// ============================================================================

bool Mesh::write_ply(const std::string& path, bool binary) const {
    std::ofstream out(path, binary ? (std::ios::binary | std::ios::out) : std::ios::out);
    if (!out) { std::cerr << "write_ply: cannot open " << path << "\n"; return false; }

    // Compact index map over live vertices (handles may be sparse pre-GC).
    std::vector<int> idx(points_.size(), -1);
    int nv = 0;
    for (size_t i = 0; i < points_.size(); ++i) if (!v_deleted_[i]) idx[i] = nv++;
    int nf = 0;
    for (size_t f = 0; f < faces_.size(); ++f) if (!f_deleted_[f]) ++nf;

    out << "ply\n";
    out << "format " << (binary ? "binary_little_endian" : "ascii") << " 1.0\n";
    out << "comment written by SimpleMesh\n";
    out << "element vertex " << nv << "\n";
    out << "property float x\nproperty float y\nproperty float z\n";
    out << "element face " << nf << "\n";
    out << "property list uchar int vertex_indices\n";
    out << "end_header\n";

    if (binary) {
        for (size_t i = 0; i < points_.size(); ++i) {
            if (v_deleted_[i]) continue;
            write_le<float>(out, (float)points_[i].x);
            write_le<float>(out, (float)points_[i].y);
            write_le<float>(out, (float)points_[i].z);
        }
        for (size_t f = 0; f < faces_.size(); ++f) {
            if (f_deleted_[f]) continue;
            auto vs = face_vertices(FaceHandle((int)f));
            write_le<uint8_t>(out, (uint8_t)vs.size());
            for (VertexHandle v : vs) write_le<int32_t>(out, idx[v.idx()]);
        }
    } else {
        for (size_t i = 0; i < points_.size(); ++i) {
            if (v_deleted_[i]) continue;
            out << points_[i].x << ' ' << points_[i].y << ' ' << points_[i].z << '\n';
        }
        for (size_t f = 0; f < faces_.size(); ++f) {
            if (f_deleted_[f]) continue;
            auto vs = face_vertices(FaceHandle((int)f));
            out << vs.size();
            for (VertexHandle v : vs) out << ' ' << idx[v.idx()];
            out << '\n';
        }
    }
    return true;
}

bool Mesh::read_ply(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) { std::cerr << "read_ply: cannot open " << path << "\n"; return false; }

    std::string line;
    std::getline(in, line);
    if (line.rfind("ply", 0) != 0) { std::cerr << "read_ply: not a PLY file\n"; return false; }

    bool binary = false;
    int  n_vertices = 0, n_faces = 0;
    // vertex scalar property types, in order; we record which index is x/y/z.
    std::vector<std::string> vprop_types;
    int xi = -1, yi = -1, zi = -1;
    std::string face_count_type = "uchar", face_index_type = "int";
    enum { NONE, VERTEX, FACE } cur = NONE;

    while (std::getline(in, line)) {
        std::istringstream ls(line);
        std::string kw; ls >> kw;
        if (kw == "format") {
            std::string fmt; ls >> fmt;
            binary = (fmt.rfind("binary", 0) == 0);   // we assume little-endian
        } else if (kw == "comment") {
            // ignore
        } else if (kw == "element") {
            std::string what; ls >> what;
            if (what == "vertex") { ls >> n_vertices; cur = VERTEX; }
            else if (what == "face") { ls >> n_faces; cur = FACE; }
            else cur = NONE;
        } else if (kw == "property") {
            if (cur == VERTEX) {
                std::string type, name; ls >> type >> name;
                if (name == "x") xi = (int)vprop_types.size();
                else if (name == "y") yi = (int)vprop_types.size();
                else if (name == "z") zi = (int)vprop_types.size();
                vprop_types.push_back(type);
            } else if (cur == FACE) {
                std::string kind; ls >> kind;     // expect "list"
                if (kind == "list") ls >> face_count_type >> face_index_type;
            }
        } else if (kw == "end_header") {
            break;
        }
    }
    if (xi < 0 || yi < 0 || zi < 0) { std::cerr << "read_ply: no x/y/z properties\n"; return false; }

    std::vector<VertexHandle> vmap(n_vertices);
    for (int i = 0; i < n_vertices; ++i) {
        std::vector<double> vals(vprop_types.size());
        if (binary) {
            for (size_t p = 0; p < vprop_types.size(); ++p)
                vals[p] = read_ply_scalar_binary(in, vprop_types[p]);
        } else {
            std::string row; std::getline(in, row); std::istringstream rs(row);
            for (size_t p = 0; p < vprop_types.size(); ++p) rs >> vals[p];
        }
        vmap[i] = add_vertex(Vec3(vals[xi], vals[yi], vals[zi]));
    }

    for (int f = 0; f < n_faces; ++f) {
        int count;
        std::vector<VertexHandle> fv;
        if (binary) {
            count = (int)read_ply_scalar_binary(in, face_count_type);
            for (int k = 0; k < count; ++k)
                fv.push_back(vmap[(int)read_ply_scalar_binary(in, face_index_type)]);
        } else {
            std::string row; std::getline(in, row); std::istringstream rs(row);
            rs >> count;
            for (int k = 0; k < count; ++k) { int vi; rs >> vi; fv.push_back(vmap[vi]); }
        }
        if (fv.size() >= 3) add_face(fv);
    }
    return true;
}

// ============================================================================
//  STL
// ============================================================================

bool Mesh::write_stl(const std::string& path, bool binary) const {
    // Gather triangles (fan-triangulating polygons) so we can write a count.
    std::vector<std::array<Vec3, 3>> tris;
    for (size_t f = 0; f < faces_.size(); ++f) {
        if (f_deleted_[f]) continue;
        auto vs = face_vertices(FaceHandle((int)f));
        for (size_t k = 1; k + 1 < vs.size(); ++k)            // fan: (0, k, k+1)
            tris.push_back({ point(vs[0]), point(vs[k]), point(vs[k + 1]) });
    }

    auto normal = [](const std::array<Vec3, 3>& t) {
        return (t[1] - t[0]).cross(t[2] - t[0]).normalized();
    };

    if (binary) {
        std::ofstream out(path, std::ios::binary);
        if (!out) { std::cerr << "write_stl: cannot open " << path << "\n"; return false; }
        char header[80] = {0};
        std::strncpy(header, "SimpleMesh binary STL", sizeof(header) - 1);
        out.write(header, 80);
        write_le<uint32_t>(out, (uint32_t)tris.size());
        for (auto& t : tris) {
            Vec3 n = normal(t);
            for (double c : { n.x, n.y, n.z }) write_le<float>(out, (float)c);
            for (auto& v : t) for (double c : { v.x, v.y, v.z }) write_le<float>(out, (float)c);
            write_le<uint16_t>(out, 0);   // attribute byte count
        }
    } else {
        std::ofstream out(path);
        if (!out) { std::cerr << "write_stl: cannot open " << path << "\n"; return false; }
        out << "solid SimpleMesh\n";
        for (auto& t : tris) {
            Vec3 n = normal(t);
            out << "  facet normal " << n.x << ' ' << n.y << ' ' << n.z << "\n    outer loop\n";
            for (auto& v : t) out << "      vertex " << v.x << ' ' << v.y << ' ' << v.z << "\n";
            out << "    endloop\n  endfacet\n";
        }
        out << "endsolid SimpleMesh\n";
    }
    return true;
}

bool Mesh::read_stl(const std::string& path) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in) { std::cerr << "read_stl: cannot open " << path << "\n"; return false; }
    const std::streamsize size = in.tellg();
    in.seekg(0);

    // Reliable ASCII/binary sniff: a binary STL is exactly 84 + 50*count bytes.
    bool binary = false;
    if (size >= 84) {
        in.seekg(80);
        uint32_t count = read_le<uint32_t>(in);
        if ((std::streamsize)(84 + 50ull * count) == size) binary = true;
        in.seekg(0);
    }

    // Weld coincident corners (stored as floats) into shared vertices.
    std::map<std::array<float, 3>, VertexHandle> weld;
    auto vertex_for = [&](const Vec3& p) {
        std::array<float, 3> key{ (float)p.x, (float)p.y, (float)p.z };
        auto it = weld.find(key);
        if (it != weld.end()) return it->second;
        VertexHandle v = add_vertex(Vec3(key[0], key[1], key[2]));
        weld.emplace(key, v);
        return v;
    };

    if (binary) {
        in.seekg(80);
        uint32_t count = read_le<uint32_t>(in);
        for (uint32_t i = 0; i < count; ++i) {
            read_le<float>(in); read_le<float>(in); read_le<float>(in);   // skip normal
            Vec3 p[3];
            for (auto& q : p) { q.x = read_le<float>(in); q.y = read_le<float>(in); q.z = read_le<float>(in); }
            read_le<uint16_t>(in);                                        // skip attribute
            add_triangle(vertex_for(p[0]), vertex_for(p[1]), vertex_for(p[2]));
        }
    } else {
        in.close();
        std::ifstream tin(path);
        std::string tok;
        std::vector<Vec3> corners;
        while (tin >> tok) {
            if (tok == "vertex") {
                Vec3 q; tin >> q.x >> q.y >> q.z;
                corners.push_back(q);
                if (corners.size() == 3) {
                    add_triangle(vertex_for(corners[0]), vertex_for(corners[1]), vertex_for(corners[2]));
                    corners.clear();
                }
            }
        }
    }
    return true;
}

// ============================================================================
//  OFF (Object File Format) - the simplest of the classic formats.
//  Header "OFF", then a line "nV nF nE", then nV vertex lines "x y z", then
//  nF face lines "k i0 i1 ... i(k-1)". ASCII only here.
// ============================================================================

bool Mesh::write_off(const std::string& path) const {
    std::ofstream out(path);
    if (!out) { std::cerr << "write_off: cannot open " << path << "\n"; return false; }

    std::vector<int> idx(points_.size(), -1);
    int nv = 0;
    for (size_t i = 0; i < points_.size(); ++i) if (!v_deleted_[i]) idx[i] = nv++;
    int nf = 0;
    for (size_t f = 0; f < faces_.size(); ++f) if (!f_deleted_[f]) ++nf;

    out << "OFF\n" << nv << ' ' << nf << " 0\n";
    for (size_t i = 0; i < points_.size(); ++i) {
        if (v_deleted_[i]) continue;
        out << points_[i].x << ' ' << points_[i].y << ' ' << points_[i].z << '\n';
    }
    for (size_t f = 0; f < faces_.size(); ++f) {
        if (f_deleted_[f]) continue;
        auto vs = face_vertices(FaceHandle((int)f));
        out << vs.size();
        for (VertexHandle v : vs) out << ' ' << idx[v.idx()];
        out << '\n';
    }
    return true;
}

bool Mesh::read_off(const std::string& path) {
    std::ifstream in(path);
    if (!in) { std::cerr << "read_off: cannot open " << path << "\n"; return false; }

    auto next_data_line = [&](std::string& line) {
        while (std::getline(in, line)) {
            // skip blank lines and comments (#)
            size_t s = line.find_first_not_of(" \t\r");
            if (s == std::string::npos || line[s] == '#') continue;
            return true;
        }
        return false;
    };

    std::string line;
    if (!next_data_line(line)) return false;
    // first token should be OFF (possibly with a count on the same line in some
    // variants, but the standard puts counts on the next line)
    if (line.rfind("OFF", 0) != 0 && line.rfind("off", 0) != 0) {
        std::cerr << "read_off: not an OFF file\n"; return false;
    }
    int nv = 0, nf = 0, ne = 0;
    // counts may follow "OFF" on the same line; try that first
    {
        std::istringstream hs(line.substr(line.find("OFF") + 3));
        if (!(hs >> nv >> nf >> ne)) {
            if (!next_data_line(line)) return false;
            std::istringstream cs(line);
            cs >> nv >> nf >> ne;
        }
    }

    std::vector<VertexHandle> vmap(nv);
    for (int i = 0; i < nv; ++i) {
        if (!next_data_line(line)) return false;
        std::istringstream vs(line);
        double x, y, z; vs >> x >> y >> z;
        vmap[i] = add_vertex({x, y, z});
    }
    for (int f = 0; f < nf; ++f) {
        if (!next_data_line(line)) return false;
        std::istringstream fs(line);
        int k; fs >> k;
        std::vector<VertexHandle> fv;
        for (int j = 0; j < k; ++j) { int vi; fs >> vi; fv.push_back(vmap[vi]); }
        if (fv.size() >= 3) add_face(fv);
    }
    return true;
}

} // namespace sm
