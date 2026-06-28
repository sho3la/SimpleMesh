// ============================================================================
//  simplemesh/exact/Expansion.h - arbitrary-precision sums of doubles
// ----------------------------------------------------------------------------
//  An arbitrary-precision number type built from Shewchuk's adaptive-precision
//  routines (two_sum / two_product / fast_expansion_sum_zeroelim /
//  scale_expansion_zeroelim). Each expansion's component array is backed by a
//  std::vector<double> (heap), which keeps the implementation simple; the values
//  produced are exact.
//
//  An "expansion" is a sequence of non-overlapping doubles, ordered by
//  increasing magnitude, whose exact sum is the represented value. Because the
//  components are non-overlapping and sorted, the sign of the value equals the
//  sign of the largest-magnitude (last) component - that is what makes exact
//  predicates possible.
//
//  IMPORTANT: compile any TU that relies on exactness with -ffp-contract=off
//  (GCC/Clang) so the compiler does not fuse a*b+c into an FMA and break the
//  two_product error-free transform. See CMake for the flag.
// ============================================================================
#pragma once

#include <vector>
#include <cstddef>
#include <cmath>

namespace sm {
namespace exact {

// ----------------------------------------------------------------------------
//  Error-free transforms (Shewchuk)
// ----------------------------------------------------------------------------

// x + y = a + b exactly, with x = fl(a+b).
inline void two_sum(double a, double b, double& x, double& y) {
    x = a + b;
    double bvirt = x - a;
    double avirt = x - bvirt;
    double bround = b - bvirt;
    double around = a - avirt;
    y = around + bround;
}

// x + y = a - b exactly, with x = fl(a-b).
inline void two_diff(double a, double b, double& x, double& y) {
    x = a - b;
    double bvirt = a - x;
    double avirt = x + bvirt;
    double bround = bvirt - b;
    double around = a - avirt;
    y = around + bround;
}

// As two_sum but requires |a| >= |b|.
inline void fast_two_sum(double a, double b, double& x, double& y) {
    x = a + b;
    double bvirt = x - a;
    y = b - bvirt;
}

// As two_diff but requires |a| >= |b|.
inline void fast_two_diff(double a, double b, double& x, double& y) {
    x = a - b;
    double bvirt = a - x;
    y = bvirt - b;
}

// Veltkamp split: a = ahi + alo, each with half the mantissa.
inline void split(double a, double& ahi, double& alo) {
    // 2^ceil(p/2)+1 = 134217729 for IEEE double (p = 53).
    const double splitter = 134217729.0;
    double c = splitter * a;
    double abig = c - a;
    ahi = c - abig;
    alo = a - ahi;
}

// x + y = a * b exactly, with x = fl(a*b). Split-based (no FMA) to match
// the non-FP_FAST_FMA path.
inline void two_product(double a, double b, double& x, double& y) {
    x = a * b;
    double ahi, alo, bhi, blo;
    split(a, ahi, alo);
    split(b, bhi, blo);
    double err1 = x - (ahi * bhi);
    double err2 = err1 - (alo * bhi);
    double err3 = err2 - (ahi * blo);
    y = (alo * blo) - err3;
}

inline void two_product_presplit(
    double a, double b, double bhi, double blo, double& x, double& y
) {
    x = a * b;
    double ahi, alo;
    split(a, ahi, alo);
    double err1 = x - (ahi * bhi);
    double err2 = err1 - (alo * bhi);
    double err3 = err2 - (ahi * blo);
    y = (alo * blo) - err3;
}

// x + y = a * a exactly.
inline void square(double a, double& x, double& y) {
    x = a * a;
    double ahi, alo;
    split(a, ahi, alo);
    double err1 = x - (ahi * ahi);
    double err3 = err1 - ((ahi + ahi) * alo);
    y = (alo * alo) - err3;
}

inline int sign_of(double x) { return (x > 0) - (x < 0); }

// ----------------------------------------------------------------------------
//  Expansion: a heap-backed non-overlapping component array.
// ----------------------------------------------------------------------------
class Expansion {
public:
    Expansion() = default;
    explicit Expansion(double a) { assign(a); }

    std::size_t length() const { return x_.size(); }
    bool empty() const { return x_.empty(); }
    double operator[](std::size_t i) const { return x_[i]; }

    // value := a
    Expansion& assign(double a) {
        x_.clear();
        if (a != 0.0) x_.push_back(a);
        return *this;
    }

    // Sign of the represented value: sign of the most significant component.
    int sign() const {
        if (x_.empty()) return 0;
        return sign_of(x_.back());
    }

    // Best double approximation (exact sum is recovered as the rounded total).
    double estimate() const {
        double s = 0.0;
        for (double c : x_) s += c;  // components ascend in magnitude
        return s;
    }

    // negation: value := -value
    Expansion& negate() {
        for (double& c : x_) c = -c;
        return *this;
    }

    // --- factory operations -------------------------------------------------

    static Expansion sum(double a, double b) {
        Expansion h; h.x_.resize(2);
        two_sum(a, b, h.x_[1], h.x_[0]);
        h.zeroelim_pair();
        return h;
    }
    static Expansion diff(double a, double b) {
        Expansion h; h.x_.resize(2);
        two_diff(a, b, h.x_[1], h.x_[0]);
        h.zeroelim_pair();
        return h;
    }
    static Expansion product(double a, double b) {
        Expansion h; h.x_.resize(2);
        two_product(a, b, h.x_[1], h.x_[0]);
        h.zeroelim_pair();
        return h;
    }

    // h := e + f  (fast_expansion_sum_zeroelim)
    static Expansion sum(const Expansion& e, const Expansion& f) {
        Expansion h;
        fast_expansion_sum_zeroelim(e, f, h);
        return h;
    }
    // h := e - f
    static Expansion diff(const Expansion& e, const Expansion& f) {
        Expansion g = f; g.negate();
        return sum(e, g);
    }
    // h := e * b  (scale_expansion_zeroelim)
    static Expansion scale(const Expansion& e, double b) {
        Expansion h;
        scale_expansion_zeroelim(e, b, h);
        return h;
    }
    // h := e * f  (distillation: sum of scaled rows; exact)
    static Expansion product(const Expansion& e, const Expansion& f) {
        Expansion h;  // zero
        if (e.empty() || f.empty()) return h;  // 0 * anything = 0
        for (std::size_t i = 0; i < f.length(); ++i) {
            Expansion row = scale(e, f[i]);
            h = sum(h, row);
        }
        return h;
    }

private:
    std::vector<double> x_;  // non-overlapping, ascending magnitude

    void zeroelim_pair() {
        // Drop a zero low word from a freshly built length-2 expansion.
        if (x_.size() == 2 && x_[0] == 0.0) {
            x_[0] = x_[1];
            x_.pop_back();
        }
        if (x_.size() >= 1 && x_.back() == 0.0 && x_.size() == 1) {
            // value is exactly zero -> empty
            x_.clear();
        }
    }

    // ---- Shewchuk zero-eliminating combinators ------------------------------

    static void scale_expansion_zeroelim(
        const Expansion& e, double b, Expansion& h
    ) {
        double Q, sum, hh, product1, product0, bhi, blo;
        std::size_t eindex, hindex;
        std::size_t elen = e.length();
        if (elen == 0 || b == 0.0) { h.x_.clear(); return; }  // 0 -> empty (zero)
        h.x_.assign(2 * elen, 0.0);

        split(b, bhi, blo);
        two_product_presplit(e[0], b, bhi, blo, Q, hh);
        hindex = 0;
        if (hh != 0) h.x_[hindex++] = hh;
        for (eindex = 1; eindex < elen; eindex++) {
            double enow = e[eindex];
            two_product_presplit(enow, b, bhi, blo, product1, product0);
            two_sum(Q, product0, sum, hh);
            if (hh != 0) h.x_[hindex++] = hh;
            fast_two_sum(product1, sum, Q, hh);
            if (hh != 0) h.x_[hindex++] = hh;
        }
        if ((Q != 0.0) || (hindex == 0)) h.x_[hindex++] = Q;
        h.x_.resize(hindex);
    }

    static void fast_expansion_sum_zeroelim(
        const Expansion& e, const Expansion& f, Expansion& h
    ) {
        double Q, Qnew, hh, enow, fnow;
        std::size_t eindex, findex, hindex;
        std::size_t elen = e.length();
        std::size_t flen = f.length();

        if (elen == 0) { h = f; return; }
        if (flen == 0) { h = e; return; }

        h.x_.assign(elen + flen, 0.0);
        enow = e[0];
        fnow = f[0];
        eindex = findex = 0;
        if ((fnow > enow) == (fnow > -enow)) {
            Q = enow; enow = (++eindex < elen) ? e[eindex] : enow;
        } else {
            Q = fnow; fnow = (++findex < flen) ? f[findex] : fnow;
        }
        hindex = 0;
        if ((eindex < elen) && (findex < flen)) {
            if ((fnow > enow) == (fnow > -enow)) {
                fast_two_sum(enow, Q, Qnew, hh);
                enow = (++eindex < elen) ? e[eindex] : enow;
            } else {
                fast_two_sum(fnow, Q, Qnew, hh);
                fnow = (++findex < flen) ? f[findex] : fnow;
            }
            Q = Qnew;
            if (hh != 0.0) h.x_[hindex++] = hh;
            while ((eindex < elen) && (findex < flen)) {
                if ((fnow > enow) == (fnow > -enow)) {
                    two_sum(Q, enow, Qnew, hh);
                    enow = (++eindex < elen) ? e[eindex] : enow;
                } else {
                    two_sum(Q, fnow, Qnew, hh);
                    fnow = (++findex < flen) ? f[findex] : fnow;
                }
                Q = Qnew;
                if (hh != 0.0) h.x_[hindex++] = hh;
            }
        }
        while (eindex < elen) {
            two_sum(Q, enow, Qnew, hh);
            enow = (++eindex < elen) ? e[eindex] : enow;
            Q = Qnew;
            if (hh != 0.0) h.x_[hindex++] = hh;
        }
        while (findex < flen) {
            two_sum(Q, fnow, Qnew, hh);
            fnow = (++findex < flen) ? f[findex] : fnow;
            Q = Qnew;
            if (hh != 0.0) h.x_[hindex++] = hh;
        }
        if ((Q != 0.0) || (hindex == 0)) h.x_[hindex++] = Q;
        h.x_.resize(hindex);
    }
};

// ----------------------------------------------------------------------------
//  Operator-overloaded arbitrary-precision number type.
//  Minimal Day-1 surface: +, -, *, unary -, and sign-based comparison.
// ----------------------------------------------------------------------------
class ExpansionNt {
public:
    ExpansionNt() : e_(0.0) {}
    ExpansionNt(double a) : e_(a) {}                 // implicit
    explicit ExpansionNt(Expansion e) : e_(std::move(e)) {}

    const Expansion& rep() const { return e_; }
    int sign() const { return e_.sign(); }
    double estimate() const { return e_.estimate(); }

    ExpansionNt operator-() const { Expansion g = e_; g.negate(); return ExpansionNt(g); }

    friend ExpansionNt operator+(const ExpansionNt& a, const ExpansionNt& b) {
        return ExpansionNt(Expansion::sum(a.e_, b.e_));
    }
    friend ExpansionNt operator-(const ExpansionNt& a, const ExpansionNt& b) {
        return ExpansionNt(Expansion::diff(a.e_, b.e_));
    }
    friend ExpansionNt operator*(const ExpansionNt& a, const ExpansionNt& b) {
        return ExpansionNt(Expansion::product(a.e_, b.e_));
    }

    // Comparison via the sign of the exact difference.
    friend int compare(const ExpansionNt& a, const ExpansionNt& b) {
        return Expansion::diff(a.e_, b.e_).sign();
    }
    friend bool operator<(const ExpansionNt& a, const ExpansionNt& b)  { return compare(a, b) < 0; }
    friend bool operator>(const ExpansionNt& a, const ExpansionNt& b)  { return compare(a, b) > 0; }
    friend bool operator==(const ExpansionNt& a, const ExpansionNt& b) { return compare(a, b) == 0; }

private:
    Expansion e_;
};

} // namespace exact
} // namespace sm
