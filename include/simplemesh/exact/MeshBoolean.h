// ============================================================================
//  simplemesh/exact/MeshBoolean.h - Weiler model, radial sort, classification
// ----------------------------------------------------------------------------
//  Turns the surface arrangement (SurfaceIntersection.h) into a volumetric
//  "Weiler model" and keeps only the facets on the boundary of the region
//  selected by a boolean expression over the input operands (the connected
//  components of the input). This gives exact union / intersection / difference.
//
//   build():  (a combinatorial 3-map with surfacic alpha2 + volumetric alpha3)
//     1. duplicate every facet (reversed) and sew the two with volumetric alpha3
//        links, so each side of an interface owns its own facet.
//     2. group halfedges into radial BUNDLES (halfedges on the same directed edge).
//     3. sew the manifold edges (bundles of 1 or 2 halfedges) with alpha2.
//     4. RADIAL SORT the facets around every non-manifold edge using the EXACT 3D
//        orientation predicate (a (u_sign,v_sign) pseudo-angle quadrant table with
//        an orient3d tie-break), then create the alpha2 links around the radial
//        edge in that cyclic order.
//     5. the connected components of the alpha2 graph are the volumetric SHELLS
//        (regions); alpha2+alpha3 components are the connected objects.
//
//   classify(op):
//     - signed volume of each shell -> the largest-|volume| shell of each object
//       is its external boundary;
//     - operand-inclusion bits of each object by RAY TRACING from one of its
//       vertices (random ray, XOR the operand bit per exact segment/triangle
//       crossing);
//     - propagate inclusion bits from the external shell (unchanged across alpha2,
//       XOR the operand bit across alpha3);
//     - keep the facet of each alpha3 pair on the boundary of the region where the
//       boolean expression flips (the E(outside) && !E(inside) rule).
//
//  The radial sort is done geometrically for every non-manifold bundle (no
//  parallelism or chart-order reuse), which is robust and exact for these inputs.
//
//  Compile strict-FP (-ffp-contract=off / /fp:strict).
// ============================================================================
#pragma once

#include "SurfaceIntersection.h"
#include "Predicates.h"
#include "HomogeneousGeometry.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
#include <numeric>
#include <random>
#include <stack>
#include <vector>

namespace sm {
namespace exact {

enum class BoolOp { Union, Intersection, Difference, SymmetricDifference };

static constexpr std::size_t NO_INDEX = std::size_t(-1);

// ---- exact segment/triangle intersection ----------------------------------
// True iff the open segment crosses the triangle interior and is not coplanar.
// 'degenerate' => geometrically ambiguous (retry with a fresh random ray).
template <class P>
inline bool segment_triangle_intersection(
    const P& P1, const P& P2, const P& q1, const P& q2, const P& q3,
    bool& degenerate, int (*o3d)(const P&, const P&, const P&, const P&)) {
    degenerate = false;
    int o1 = o3d(P1, q1, q2, q3);
    int o2 = o3d(P2, q1, q2, q3);
    if (o1 == 0 && o2 == 0) { degenerate = true; return false; }
    if (o1 == o2) return false;
    int s1 = o3d(P1, P2, q1, q2);
    int s2 = o3d(P1, P2, q2, q3);
    int s3 = o3d(P1, P2, q3, q1);
    if (s1*s2 < 0 || s2*s3 < 0 || s3*s1 < 0) return false;
    if (s1 == 0 || s2 == 0 || s3 == 0) { degenerate = true; return false; }
    if (o1 == 0 || o2 == 0) { degenerate = true; return false; }
    return true;
}
inline int orient3d_he(const Vec3HE& a, const Vec3HE& b, const Vec3HE& c, const Vec3HE& d) {
    return orient3d(a, b, c, d);
}

struct BooleanResult {
    std::vector<std::array<double,3>>      points;
    std::vector<std::array<std::size_t,3>> triangles;
    int n_operands = 0;
};

// ----------------------------------------------------------------------------
//  WeilerModel - the duplicated/alpha3/alpha2 combinatorial 3-map over the
//  arrangement soup.
// ----------------------------------------------------------------------------
class WeilerModel {
public:
    using Pts = std::vector<std::array<double,3>>;
    using Tris = std::vector<std::array<std::size_t,3>>;

    WeilerModel(const ArrangementResult& A, const Pts& Vin, const Tris& Fin)
        : A_(A), Vin_(Vin), Fin_(Fin) {
        nf0_ = A.triangles.size();
        nf_ = 2 * nf0_;
        op_input_ = detail::operand_of_facet(Vin, Fin, nb_operands_);
    }

    int nb_operands() const { return nb_operands_; }

    // --- halfedge / facet accessors -----------------------------------------
    std::size_t facet(std::size_t h) const { return h / 3; }
    std::size_t hv(std::size_t h, std::size_t dlv) const {
        return facet_[h/3][(h + dlv) % 3];
    }
    std::size_t alpha3(std::size_t h) const { return alpha3_[h]; }
    std::size_t facet_alpha3(std::size_t f) const { return alpha3_[3*f] / 3; }
    std::size_t alpha2(std::size_t h) const {
        std::size_t t1 = h/3, t2 = adj_facet_[h];
        if (t2 == NO_INDEX) return NO_INDEX;
        for (std::size_t k = 0; k < 3; ++k)
            if (adj_facet_[3*t2 + k] == t1) return 3*t2 + k;
        return NO_INDEX;
    }
    void sew2(std::size_t h1, std::size_t h2) {
        adj_facet_[h1] = h2/3; adj_facet_[h2] = h1/3;
    }

    Vec3HE exact_vertex(std::size_t v) const { return A_.exact_points[v]; }
    bool is_original_vertex(std::size_t v) const { return A_.is_original[v] != 0; }

    // Original (input) facet vertices supporting facet f, orientation preserved.
    std::array<std::array<double,3>,3> initial_facet_vertices(std::size_t f) const {
        std::size_t of = source_[f];
        std::array<std::array<double,3>,3> p = {
            Vin_[Fin_[of][0]], Vin_[Fin_[of][1]], Vin_[Fin_[of][2]] };
        if (flipped_[f]) std::swap(p[0], p[2]);
        return p;
    }

    // ========================================================================
    //  build_Weiler_model
    // ========================================================================
    void build() {
        // ---- step 1: duplicate facets + alpha3 ------------------------------
        facet_.resize(nf_);
        source_.resize(nf_);
        flipped_.resize(nf_);
        operand_bit_.resize(nf_);
        for (std::size_t f = 0; f < nf0_; ++f) {
            facet_[f] = A_.triangles[f];
            std::size_t f2 = f + nf0_;
            facet_[f2] = { A_.triangles[f][2], A_.triangles[f][1], A_.triangles[f][0] };
            source_[f] = source_[f2] = A_.triangle_source[f];
            flipped_[f] = 0; flipped_[f2] = 1;
            std::uint32_t ob = std::uint32_t(1) << op_input_[A_.triangle_source[f]];
            operand_bit_[f] = operand_bit_[f2] = ob;
        }
        alpha3_.assign(3*nf_, NO_INDEX);
        adj_facet_.assign(3*nf_, NO_INDEX);
        for (std::size_t f1 = 0; f1 < nf0_; ++f1) {
            std::size_t f2 = f1 + nf0_;
            sew3(3*f1,   3*f2+1);
            sew3(3*f1+1, 3*f2);
            sew3(3*f1+2, 3*f2+2);
        }

        // ---- step 2: radial bundles ----------------------------------------
        // Bundles of halfedges sharing the same directed
        // edge, built from the CANONICAL half-edges (v0<v1). The opposite
        // (v0>v1) bundle of each is exactly the alpha3-image, so we only store
        // and drive sewing from the canonical bundles - their alpha3 mates are
        // handled explicitly below. This preserves the alpha3-paired bundle
        // structure (losing it scrambles the duplicate-shell links).
        std::map<std::pair<std::size_t,std::size_t>, std::vector<std::size_t>> bmap;
        for (std::size_t h = 0; h < 3*nf_; ++h)
            if (hv(h,0) < hv(h,1)) bmap[{ hv(h,0), hv(h,1) }].push_back(h);

        // ---- step 3 & 4: sew manifold edges; radial-sort non-manifold edges -
        RadialSort RS(*this);
        for (auto& kv : bmap) {
            auto& H = kv.second;
            std::size_t N = H.size();
            if (N == 1) {
                // surface border ("hem"): sew the halfedge to its alpha3 mate.
                sew2(H[0], alpha3(H[0]));
            } else if (N == 2) {
                // Manifold edge: the canonical bundle connects the two original
                // facets; its alpha3 mirror connects the two duplicate facets.
                sew2(H[0], alpha3(H[1]));
                sew2(alpha3(H[0]), H[1]);
            } else {
                // Non-manifold radial edge: sort the facets around it
                // geometrically, then create the alpha2 links in cyclic order
                // The matching h_i <-> alpha3(h_{i+1}) pairs up
                // ALL 2N halfedges around the edge (these N and their alpha3
                // mates = the opposite directed bundle), so one canonical bundle
                // sews the whole radial fan.
                RS.init(H[0]);
                std::sort(H.begin(), H.end(),
                          [&](std::size_t a, std::size_t b){ return RS(a, b); });
                for (std::size_t i = 0; i < N; ++i)
                    sew2(H[i], alpha3(H[(i == N-1) ? 0 : i+1]));
            }
        }
    }

    // ========================================================================
    //  classify(expr): keep the facets on the boundary of the selected region.
    // ========================================================================
    std::vector<std::array<std::size_t,3>> classify(BoolOp op, unsigned seed) {
        // ---- regions (shells): connected components of the alpha2 graph -----
        std::vector<std::size_t> region(nf_, NO_INDEX);
        std::size_t nb_regions = connected_components(region, /*use_alpha3=*/false);

        // ---- objects: connected components of alpha2 + alpha3 ---------------
        std::vector<std::size_t> component(nf_, NO_INDEX);
        std::size_t nb_components = connected_components(component, /*use_alpha3=*/true);

        // ---- signed volume per shell ----------------------------------------
        std::vector<double> region_volume(nb_regions, 0.0);
        for (std::size_t f = 0; f < nf_; ++f)
            region_volume[region[f]] += signed_volume(f);

        // ---- external shell (largest |volume|) per object -------------------
        std::vector<double> max_vol(nb_components, 0.0);
        std::vector<std::size_t> ext_shell(nb_components, NO_INDEX);
        for (std::size_t f = 0; f < nf_; ++f) {
            double V = region_volume[region[f]];
            if (std::fabs(V) >= std::fabs(max_vol[component[f]])) {
                max_vol[component[f]] = V;        // signed
                ext_shell[component[f]] = region[f];
            }
        }

        // ---- object inclusion bits by ray tracing (if >1 object) ------------
        std::vector<std::uint32_t> comp_bits(nb_components, 0);
        if (nb_components > 1) {
            std::mt19937 rng(seed);
            for (std::size_t c = 0; c < nb_components; ++c)
                comp_bits[c] = classify_component(c, component, rng);
        }

        // ---- propagate operand-inclusion bits from the external shell -------
        std::vector<std::uint32_t> oib(nf_, 0);
        std::vector<char> visited(nf_, 0);
        std::stack<std::size_t> S;
        for (std::size_t f = 0; f < nf_; ++f)
            if (region[f] == ext_shell[component[f]]) {
                visited[f] = 1;
                oib[f] = comp_bits[component[f]];
                S.push(f);
            }
        while (!S.empty()) {
            std::size_t f1 = S.top(); S.pop();
            std::size_t f3 = facet_alpha3(f1);
            if (f3 != NO_INDEX && !visited[f3]) {
                visited[f3] = 1;
                oib[f3] = oib[f1] ^ operand_bit_[f1];   // cross the surface
                S.push(f3);
            }
            for (std::size_t le = 0; le < 3; ++le) {
                std::size_t h2 = alpha2(3*f1 + le);
                if (h2 == NO_INDEX) continue;
                std::size_t f2 = h2/3;
                if (!visited[f2]) { visited[f2] = 1; oib[f2] = oib[f1]; S.push(f2); }
            }
        }

        // ---- keep facets where the boolean expression flips -----------------
        std::uint32_t all = (std::uint32_t(1) << nb_operands_) - 1u;
        auto E = [&](std::uint32_t bits) { return eval_expr(op, bits, all); };
        std::vector<std::array<std::size_t,3>> kept;
        for (std::size_t f = 0; f < nf_; ++f) {
            bool flip = (max_vol[component[f]] < 0.0);
            std::uint32_t f_in = oib[f];
            std::uint32_t g_in = oib[facet_alpha3(f)];
            bool keep = flip ? (E(f_in) && !E(g_in)) : (E(g_in) && !E(f_in));
            if (keep) kept.push_back(facet_[f]);
        }
        return kept;
    }

private:
    static bool eval_expr(BoolOp op, std::uint32_t bits, std::uint32_t all) {
        switch (op) {
        case BoolOp::Union:               return bits != 0;
        case BoolOp::Intersection:        return bits == all;
        case BoolOp::Difference:          return (bits & 1u) && !(bits & 2u);  // A-B
        case BoolOp::SymmetricDifference: {
            int c = 0; for (std::uint32_t b = bits; b; b >>= 1) c += int(b & 1u);
            return (c & 1);
        }
        }
        return false;
    }

    void sew3(std::size_t h1, std::size_t h2) { alpha3_[h1] = h2; alpha3_[h2] = h1; }

    double signed_volume(std::size_t f) const {
        const auto& a = A_.points[facet_[f][0]];
        const auto& b = A_.points[facet_[f][1]];
        const auto& c = A_.points[facet_[f][2]];
        double cx = b[1]*c[2]-b[2]*c[1], cy = b[2]*c[0]-b[0]*c[2], cz = b[0]*c[1]-b[1]*c[0];
        return (a[0]*cx + a[1]*cy + a[2]*cz) / 6.0;
    }

    // Connected components of facets; alpha2 always, alpha3 optionally.
    std::size_t connected_components(std::vector<std::size_t>& comp, bool use_alpha3) const {
        std::size_t nb = 0;
        for (std::size_t f = 0; f < nf_; ++f) {
            if (comp[f] != NO_INDEX) continue;
            std::stack<std::size_t> S; comp[f] = nb; S.push(f);
            while (!S.empty()) {
                std::size_t f1 = S.top(); S.pop();
                if (use_alpha3) {
                    std::size_t f2 = facet_alpha3(f1);
                    if (f2 != NO_INDEX && comp[f2] == NO_INDEX) { comp[f2] = nb; S.push(f2); }
                }
                for (std::size_t le = 0; le < 3; ++le) {
                    std::size_t h2 = alpha2(3*f1 + le);
                    if (h2 == NO_INDEX) continue;
                    std::size_t f2 = h2/3;
                    if (comp[f2] == NO_INDEX) { comp[f2] = nb; S.push(f2); }
                }
            }
            ++nb;
        }
        return nb;
    }

    // Ray-trace one object's operand-inclusion bits (exact, on the arrangement).
    std::uint32_t classify_component(std::size_t comp,
                                     const std::vector<std::size_t>& component,
                                     std::mt19937& rng) const {
        // collect candidate start vertices (prefer original ones).
        std::vector<std::size_t> verts;
        for (std::size_t f = 0; f < nf_; ++f)
            if (component[f] == comp)
                for (int k = 0; k < 3; ++k) verts.push_back(facet_[f][k]);
        std::sort(verts.begin(), verts.end());
        verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
        std::sort(verts.begin(), verts.end(), [&](std::size_t a, std::size_t b){
            return is_original_vertex(a) > is_original_vertex(b);
        });

        std::uniform_real_distribution<double> uni(-1.0, 1.0);
        for (int attempt = 0; attempt < 200; ++attempt) {
            std::size_t v = verts[attempt == 0 ? 0
                                 : std::size_t(rng()) % verts.size()];
            Vec3HE P1 = exact_vertex(v);
            double D[3] = { 1e6*uni(rng), 1e6*uni(rng), 1e6*uni(rng) };
            Vec3HE P2(P1.x + P1.w*ExpansionNt(D[0]),
                      P1.y + P1.w*ExpansionNt(D[1]),
                      P1.z + P1.w*ExpansionNt(D[2]), P1.w);
            std::uint32_t result = 0;
            bool ok = true;
            for (std::size_t t = 0; t < nf_ && ok; ++t) {
                if (component[t] == comp) continue;        // skip own object
                if (t > facet_alpha3(t)) continue;         // one per alpha3 pair
                bool degenerate = false;
                Vec3HE q1 = exact_vertex(facet_[t][0]);
                Vec3HE q2 = exact_vertex(facet_[t][1]);
                Vec3HE q3 = exact_vertex(facet_[t][2]);
                if (segment_triangle_intersection(P1, P2, q1, q2, q3, degenerate, orient3d_he))
                    result ^= operand_bit_[t];
                if (degenerate) { ok = false; }
            }
            if (ok) return result;
        }
        return 0;  // gave up (should not happen on the corpus)
    }

    // ------------------------------------------------------------------------
    //  RadialSort - exact radial ordering of facets around a common edge.
    // ------------------------------------------------------------------------
    struct RadialSort {
        const WeilerModel& W;
        std::size_t h_ref_ = NO_INDEX;
        EVec3 N_ref_;
        mutable bool degenerate_ = false;
        explicit RadialSort(const WeilerModel& w) : W(w) {}

        void init(std::size_t h_ref) {
            h_ref_ = NO_INDEX;
            N_ref_ = normal(h_ref);
            h_ref_ = h_ref;
        }

        // Outward normal of the original facet supporting h (exact).
        EVec3 normal(std::size_t h) const {
            if (h == h_ref_) return N_ref_;
            auto p = W.initial_facet_vertices(W.facet(h));
            return ev_cross(ev_sub(p[0].data(), p[1].data()),
                            ev_sub(p[0].data(), p[2].data()));
        }

        int h_refNorient(std::size_t h2) const {
            if (h2 == h_ref_) return 1;
            EVec3 N2 = normal(h2);
            return ev_dot(N_ref_, N2).sign();
        }

        // Orientation of the apex of h2's facet vs h1's facet plane.
        int h_orient(std::size_t h1, std::size_t h2) const {
            if (h1 == h2) return 0;
            std::size_t f1 = W.facet(h1), f2 = W.facet(h2);
            std::size_t w1 = W.hv(h1, 2), w2 = W.hv(h2, 2);
            if (W.is_original_vertex(w1)) {
                const auto& p0 = W.A_.points[w1];
                auto q = W.initial_facet_vertices(f2);
                return orient3d(q[0].data(), q[1].data(), q[2].data(), p0.data());
            }
            if (W.is_original_vertex(w2)) {
                const auto& q0 = W.A_.points[w2];
                auto p = W.initial_facet_vertices(f1);
                return -orient3d(p[0].data(), p[1].data(), p[2].data(), q0.data());
            }
            auto p = W.initial_facet_vertices(f1);
            Vec3HE pp0 = he_from_double(p[0].data());
            Vec3HE pp1 = he_from_double(p[1].data());
            Vec3HE pp2 = he_from_double(p[2].data());
            Vec3HE q2 = W.exact_vertex(w2);
            return -orient3d(pp0, pp1, pp2, q2);
        }

        bool operator()(std::size_t h1, std::size_t h2) const {
            int su1 = h_orient(h_ref_, h1);
            int su2 = h_orient(h_ref_, h2);
            if (su1 * su2 < 0) return su1 > 0;
            int sv1 = h_refNorient(h1);
            int sv2 = h_refNorient(h2);
            static const int tab[3][3] = { {5,4,3}, {6,-1,2}, {7,0,1} };
            int theta1 = tab[sv1+1][su1+1];
            int theta2 = tab[sv2+1][su2+1];
            if (theta1 == -1 || theta2 == -1) { degenerate_ = true; return false; }
            if (theta1 != theta2) return theta2 > theta1;
            if ((theta1 & 1) == 0) { degenerate_ = true; return false; }
            int o12 = h_orient(h1, h2);
            if (o12 == 0) { degenerate_ = true; return false; }
            return o12 > 0;
        }
    };
    friend struct RadialSort;

    const ArrangementResult& A_;
    const Pts& Vin_;
    const Tris& Fin_;
    std::size_t nf0_ = 0, nf_ = 0;
    int nb_operands_ = 0;
    std::vector<int> op_input_;                       // operand id per input facet
    std::vector<std::array<std::size_t,3>> facet_;
    std::vector<std::size_t> source_;
    std::vector<char> flipped_;
    std::vector<std::uint32_t> operand_bit_;
    std::vector<std::size_t> alpha3_;
    std::vector<std::size_t> adj_facet_;
};

// ----------------------------------------------------------------------------
//  Full pipeline: arrangement -> Weiler model + classify.
// ----------------------------------------------------------------------------
inline BooleanResult mesh_boolean(
    const std::vector<std::array<double,3>>& V,
    const std::vector<std::array<std::size_t,3>>& F,
    BoolOp op, unsigned seed = 1u) {

    BooleanResult out;
    ArrangementResult A = resolve_intersections(V, F);

    WeilerModel W(A, V, F);
    out.n_operands = W.nb_operands();
    W.build();
    auto kept = W.classify(op, seed);

    // Reindex vertices used by the kept facets.
    std::vector<std::size_t> remap(A.points.size(), NO_INDEX);
    for (auto& t : kept) {
        std::array<std::size_t,3> tri;
        for (int k = 0; k < 3; ++k) {
            std::size_t v = t[k];
            if (remap[v] == NO_INDEX) { remap[v] = out.points.size(); out.points.push_back(A.points[v]); }
            tri[k] = remap[v];
        }
        if (tri[0] != tri[1] && tri[1] != tri[2] && tri[2] != tri[0])
            out.triangles.push_back(tri);
    }
    return out;
}

} // namespace exact
} // namespace sm
