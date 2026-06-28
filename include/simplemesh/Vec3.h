// ============================================================================
//  SimpleMesh - a teaching-oriented, simplified halfedge mesh library
//  Vec3.h : a minimal 3D vector type
// ----------------------------------------------------------------------------
//  A general-purpose vector template could support any dimension and scalar
//  type. For learning purposes we only need a concrete 3D vector of doubles, so
//  we implement a small, self-contained struct. Everything is header-only and
//  constexpr-friendly.
// ============================================================================
#pragma once

#include <cmath>
#include <ostream>

namespace sm {

/// A minimal 3-component vector of doubles.
///
/// This is intentionally tiny: just enough arithmetic to compute face normals,
/// centroids and edge lengths.
struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    // --- element access (handy for generic / Python-style code) ------------
    double  operator[](int i) const { return (&x)[i]; }
    double& operator[](int i)       { return (&x)[i]; }

    // --- vector arithmetic -------------------------------------------------
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s)      const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s)      const { return {x / s, y / s, z / s}; }

    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(double s)      { x *= s; y *= s; z *= s; return *this; }

    bool operator==(const Vec3& o) const { return x == o.x && y == o.y && z == o.z; }

    // --- geometric helpers -------------------------------------------------
    double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }

    Vec3 cross(const Vec3& o) const {
        return {y * o.z - z * o.y,
                z * o.x - x * o.z,
                x * o.y - y * o.x};
    }

    double sqrnorm() const { return dot(*this); }
    double norm()    const { return std::sqrt(sqrnorm()); }

    /// Returns a unit-length copy. If the vector is (near) zero it is returned
    /// unchanged to avoid division by zero.
    Vec3 normalized() const {
        const double n = norm();
        return (n > 1e-30) ? (*this) / n : *this;
    }
};

inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    return os << '(' << v.x << ", " << v.y << ", " << v.z << ')';
}

} // namespace sm
