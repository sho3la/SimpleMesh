// ============================================================================
//  SimpleMesh - Decimater.h : a modular, pluggable decimation framework
// ----------------------------------------------------------------------------
//  The iteration-7 `quadric_decimate` is one hard-wired algorithm. Production
//  decimators instead build decimation from swappable MODULES, so you can mix
//  criteria without rewriting the engine:
//
//    * a single PRIORITY module scores each candidate collapse (lower = better);
//    * any number of BINARY modules veto collapses that violate a constraint
//      (edge too long, would flip a normal, would make a sliver, ...).
//
//  The `Decimator` engine owns the greedy priority queue, the legality checks
//  and the collapse loop; the modules only answer "how good?" / "is it allowed?".
//  This is the classic Strategy pattern, and the single most important idea in a
//  production decimator.
// ============================================================================
#pragma once

#include "Mesh.h"
#include "Algorithms.h"   // Quadric

#include <memory>
#include <utility>
#include <vector>

namespace sm {

/// Everything a module needs to judge one candidate half-edge collapse that
/// would merge `v0` (removed) into `v1` (kept).
struct CollapseInfo {
    HalfedgeHandle heh;    // the half-edge collapsed: v0 = from, v1 = to
    VertexHandle   v0, v1; // v0 disappears, v1 survives
    Vec3           p0, p1; // their positions
};

// Sentinel return values for collapse_priority (floats so they share one API).
constexpr float ILLEGAL_COLLAPSE = -1.0f;   // veto this collapse
constexpr float LEGAL_COLLAPSE   =  0.0f;   // binary module: allowed

// ----------------------------------------------------------------------------
//  Module interface. A module is either binary (a yes/no filter) or the one
//  priority module (a cost function).
// ----------------------------------------------------------------------------
class DecimationModule {
public:
    virtual ~DecimationModule() = default;

    /// Binary modules only veto; the (single) non-binary module scores.
    virtual bool is_binary() const = 0;

    /// Called once before the run (e.g. to precompute per-vertex quadrics).
    virtual void initialize(Mesh& /*m*/) {}

    /// Binary  : return LEGAL_COLLAPSE or ILLEGAL_COLLAPSE.
    /// Priority: return a cost >= 0, or ILLEGAL_COLLAPSE to refuse.
    virtual float collapse_priority(const Mesh& m, const CollapseInfo& ci) = 0;

    /// Called after a collapse is committed (e.g. to merge quadrics).
    virtual void postprocess_collapse(Mesh& /*m*/, const CollapseInfo& /*ci*/) {}
};

// ----------------------------------------------------------------------------
//  Priority module: Garland-Heckbert quadric error (the gold standard).
// ----------------------------------------------------------------------------
class ModQuadric : public DecimationModule {
public:
    bool  is_binary() const override { return false; }
    void  initialize(Mesh& m) override;
    float collapse_priority(const Mesh& m, const CollapseInfo& ci) override;
    void  postprocess_collapse(Mesh& m, const CollapseInfo& ci) override;
private:
    std::vector<Quadric> Q_;   // one accumulated quadric per vertex index
};

// ----------------------------------------------------------------------------
//  Binary module: refuse collapses whose surviving edge would be too long.
// ----------------------------------------------------------------------------
class ModEdgeLength : public DecimationModule {
public:
    explicit ModEdgeLength(double max_length) : max_length_(max_length) {}
    bool  is_binary() const override { return true; }
    float collapse_priority(const Mesh& m, const CollapseInfo& ci) override;
private:
    double max_length_;
};

// ----------------------------------------------------------------------------
//  Binary module: refuse collapses that flip an incident face normal by more
//  than `max_angle` radians (prevents foldovers).
// ----------------------------------------------------------------------------
class ModNormalFlipping : public DecimationModule {
public:
    explicit ModNormalFlipping(double max_angle_rad) : max_angle_(max_angle_rad) {}
    bool  is_binary() const override { return true; }
    float collapse_priority(const Mesh& m, const CollapseInfo& ci) override;
private:
    double max_angle_;
};

// ----------------------------------------------------------------------------
//  Binary module: refuse collapses that would create a sliver triangle whose
//  aspect ratio (longest edge / shortest altitude) exceeds `max_aspect`.
// ----------------------------------------------------------------------------
class ModAspectRatio : public DecimationModule {
public:
    explicit ModAspectRatio(double max_aspect) : max_aspect_(max_aspect) {}
    bool  is_binary() const override { return true; }
    float collapse_priority(const Mesh& m, const CollapseInfo& ci) override;
private:
    double max_aspect_;
};

// ----------------------------------------------------------------------------
//  The engine. Register modules, then decimate.
// ----------------------------------------------------------------------------
class Decimator {
public:
    explicit Decimator(Mesh& m) : mesh_(m) {}

    /// Construct a module in place and register it; returns a borrowed pointer.
    ///     dec.add_module<ModQuadric>();
    ///     dec.add_module<ModNormalFlipping>(M_PI / 4);
    template <class Mod, class... Args>
    Mod* add_module(Args&&... args) {
        auto* p = new Mod(std::forward<Args>(args)...);
        modules_.emplace_back(p);
        return p;
    }

    /// Greedily collapse until the mesh has at most `target_faces` faces (or no
    /// legal collapse remains). Runs garbage_collection at the end, so all
    /// handles into the mesh are invalidated. Returns the number of collapses.
    size_t decimate_to_faces(size_t target_faces);

private:
    CollapseInfo make_info(HalfedgeHandle h) const;
    bool         passes_binary(const CollapseInfo& ci) const;
    // cheapest legal orientation of edge e: (ok, cost, heh)
    bool best_collapse(EdgeHandle e, float& cost_out, HalfedgeHandle& heh_out);

    Mesh& mesh_;
    std::vector<std::unique_ptr<DecimationModule>> modules_;
    DecimationModule* priority_ = nullptr;   // the one non-binary module
};

} // namespace sm
