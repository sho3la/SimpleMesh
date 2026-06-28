// ============================================================================
//  SimpleMesh - Properties.h : a custom per-element property system
// ----------------------------------------------------------------------------
//  A simplified per-element property system. It lets you attach an
//  arbitrary number of named, typed data arrays to the mesh's vertices, edges,
//  halfedges or faces - e.g. a per-vertex scalar "quality", a per-face normal,
//  per-edge weights, and so on - without ever changing the Mesh class itself.
//
//  Two ideas make this work:
//
//    1. TYPE ERASURE. We want to store properties of *different* element types
//       (double, Vec3, int, ...) side by side in one container. C++ can't put
//       differently-typed objects in one std::vector directly, so we hide the
//       type behind an abstract base class `BasePropertyArray` and keep
//       `unique_ptr<BasePropertyArray>`s. The concrete `PropertyArray<T>` knows
//       its real type; the container only ever calls the virtual interface.
//
//    2. PARALLEL ARRAYS. Every property array has exactly one slot per element,
//       indexed by the element's handle index. So `quality[v.idx()]` is the
//       quality of vertex `v`. The Mesh keeps the arrays in lock-step with the
//       elements by forwarding push_back / resize / compaction to the container.
// ============================================================================
#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace sm {

// ----------------------------------------------------------------------------
//  The type-erased interface. The container manipulates properties only through
//  these virtual calls, so it never needs to know T.
// ----------------------------------------------------------------------------
class BasePropertyArray {
public:
    explicit BasePropertyArray(std::string name) : name_(std::move(name)) {}
    virtual ~BasePropertyArray() = default;

    virtual void reserve(size_t n) = 0;          // grow capacity
    virtual void resize(size_t n)  = 0;          // set element count
    virtual void push_back()       = 0;          // append one default-valued slot
    virtual void swap(size_t a, size_t b) = 0;   // swap two slots
    /// Rebuild the array keeping only mapped slots: new[old_to_new[i]] = old[i].
    virtual void compact(const std::vector<int>& old_to_new, size_t new_size) = 0;
    virtual std::unique_ptr<BasePropertyArray> clone() const = 0;

    const std::string& name() const { return name_; }

protected:
    std::string name_;
};

// ----------------------------------------------------------------------------
//  The concrete, typed property array. Just a std::vector<T> plus a default
//  value used when new elements appear.
// ----------------------------------------------------------------------------
template <class T>
class PropertyArray : public BasePropertyArray {
public:
    explicit PropertyArray(std::string name, T def = T())
        : BasePropertyArray(std::move(name)), default_(std::move(def)) {}

    void reserve(size_t n) override { data_.reserve(n); }
    void resize(size_t n)  override { data_.resize(n, default_); }
    void push_back()       override { data_.push_back(default_); }
    void swap(size_t a, size_t b) override { std::swap(data_[a], data_[b]); }

    void compact(const std::vector<int>& old_to_new, size_t new_size) override {
        std::vector<T> nd(new_size, default_);
        for (size_t i = 0; i < old_to_new.size(); ++i)
            if (old_to_new[i] >= 0)
                nd[old_to_new[i]] = std::move(data_[i]);
        data_.swap(nd);
    }

    std::unique_ptr<BasePropertyArray> clone() const override {
        return std::unique_ptr<BasePropertyArray>(new PropertyArray<T>(*this));
    }

    // typed element access
    T&       operator[](size_t i)       { return data_[i]; }
    const T& operator[](size_t i) const { return data_[i]; }

private:
    std::vector<T> data_;
    T              default_;
};

// ----------------------------------------------------------------------------
//  A bag of properties for ONE element kind (all the vertex properties, say).
//  A property is referenced by its integer slot index. Removed slots become
//  null but keep their index so other handles stay valid (same philosophy as
//  the mesh's lazy element deletion).
// ----------------------------------------------------------------------------
class PropertyContainer {
public:
    // A container owns its properties via unique_ptr, so it is move-only.
    // Declaring these explicitly guarantees the move operations exist (and that
    // copies are a clean compile error), which in turn lets the Mesh - and
    // pybind11's return-by-value - move a whole mesh cheaply.
    PropertyContainer() = default;
    PropertyContainer(PropertyContainer&&) = default;
    PropertyContainer& operator=(PropertyContainer&&) = default;
    PropertyContainer(const PropertyContainer&) = delete;
    PropertyContainer& operator=(const PropertyContainer&) = delete;

    size_t size() const { return size_; }
    size_t n_properties() const { return arrays_.size(); }

    /// Create a new property of type T; returns its slot index.
    template <class T>
    int add(const std::string& name, const T& def) {
        auto p = std::unique_ptr<PropertyArray<T>>(new PropertyArray<T>(name, def));
        p->resize(size_);                       // line it up with current elements
        arrays_.push_back(std::move(p));
        return static_cast<int>(arrays_.size()) - 1;
    }

    /// Typed access to a property array. The caller guarantees T matches what
    /// was used in add() (the handle carries T, so this is type-safe in practice).
    template <class T>
    PropertyArray<T>* get(int id) {
        return static_cast<PropertyArray<T>*>(arrays_[id].get());
    }
    template <class T>
    const PropertyArray<T>* get(int id) const {
        return static_cast<const PropertyArray<T>*>(arrays_[id].get());
    }

    int find(const std::string& name) const {
        for (size_t i = 0; i < arrays_.size(); ++i)
            if (arrays_[i] && arrays_[i]->name() == name) return static_cast<int>(i);
        return -1;
    }

    void remove(int id) { arrays_[id].reset(); }   // free, keep the slot

    // --- forwarded lifetime operations (driven by the Mesh) ---------------
    void push_back() { ++size_; for (auto& a : arrays_) if (a) a->push_back(); }
    void reserve(size_t n) { for (auto& a : arrays_) if (a) a->reserve(n); }
    void resize(size_t n)  { size_ = n; for (auto& a : arrays_) if (a) a->resize(n); }
    void compact(const std::vector<int>& old_to_new, size_t new_size) {
        size_ = new_size;
        for (auto& a : arrays_) if (a) a->compact(old_to_new, new_size);
    }

private:
    std::vector<std::unique_ptr<BasePropertyArray>> arrays_;
    size_t size_ = 0;
};

// ----------------------------------------------------------------------------
//  Typed property handles. Each carries the value type T *and* the element kind
//  (via the distinct struct), so the compiler stops you from reading a vertex
//  property with a face handle, or as the wrong type. Mirrors the handle idiom
//  from Handles.h, one level up.
// ----------------------------------------------------------------------------
template <class T> struct VertexPropHandle   { int id = -1; bool is_valid() const { return id >= 0; } };
template <class T> struct HalfedgePropHandle { int id = -1; bool is_valid() const { return id >= 0; } };
template <class T> struct EdgePropHandle     { int id = -1; bool is_valid() const { return id >= 0; } };
template <class T> struct FacePropHandle     { int id = -1; bool is_valid() const { return id >= 0; } };

} // namespace sm
