// ============================================================================
//  SimpleMesh - Handles.h : typed index handles
// ----------------------------------------------------------------------------
//  A "handle" is just an integer index into one of the mesh's element arrays,
//  wrapped in a distinct type so the compiler stops you from passing a
//  VertexHandle where a FaceHandle is expected. It is the single most
//  important idea in the whole data structure:
//
//      * Elements are stored contiguously in std::vector arrays.
//      * We never hand out raw pointers or iterators that could dangle when a
//        vector reallocates. We hand out indices.
//      * An index of -1 means "invalid / no element".
//
//  Because handles are tiny value types they are cheap to copy, store and
//  compare, and they survive reallocation of the underlying arrays.
// ============================================================================
#pragma once

#include <cstddef>
#include <functional>
#include <ostream>

namespace sm {

/// Base class for all handles: a single signed integer index.
class BaseHandle {
public:
    explicit BaseHandle(int idx = -1) : idx_(idx) {}

    /// The raw array index this handle refers to.
    int idx() const { return idx_; }

    /// A handle is valid iff its index is non-negative.
    bool is_valid() const { return idx_ >= 0; }

    /// Make this handle refer to nothing.
    void invalidate() { idx_ = -1; }

    bool operator==(const BaseHandle& o) const { return idx_ == o.idx_; }
    bool operator!=(const BaseHandle& o) const { return idx_ != o.idx_; }
    bool operator< (const BaseHandle& o) const { return idx_ <  o.idx_; }

private:
    int idx_;
};

// Each entity gets its own strong type. They add no data, only a name, so the
// type system can distinguish them. This is the C++ "strong typedef" idiom.
struct VertexHandle   : public BaseHandle { explicit VertexHandle  (int i = -1) : BaseHandle(i) {} };
struct HalfedgeHandle : public BaseHandle { explicit HalfedgeHandle(int i = -1) : BaseHandle(i) {} };
struct EdgeHandle     : public BaseHandle { explicit EdgeHandle    (int i = -1) : BaseHandle(i) {} };
struct FaceHandle     : public BaseHandle { explicit FaceHandle    (int i = -1) : BaseHandle(i) {} };

inline std::ostream& operator<<(std::ostream& os, const BaseHandle& h) {
    return os << h.idx();
}

} // namespace sm

// Allow handles to be used as keys in std::unordered_map / unordered_set.
namespace std {
template <> struct hash<sm::VertexHandle>   { size_t operator()(const sm::VertexHandle&   h) const { return static_cast<size_t>(h.idx()); } };
template <> struct hash<sm::HalfedgeHandle> { size_t operator()(const sm::HalfedgeHandle& h) const { return static_cast<size_t>(h.idx()); } };
template <> struct hash<sm::EdgeHandle>     { size_t operator()(const sm::EdgeHandle&     h) const { return static_cast<size_t>(h.idx()); } };
template <> struct hash<sm::FaceHandle>     { size_t operator()(const sm::FaceHandle&     h) const { return static_cast<size_t>(h.idx()); } };
} // namespace std
