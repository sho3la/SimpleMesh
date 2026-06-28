// ============================================================================
//  simplemesh/exact/Interval.h - round-to-nearest interval arithmetic filter
// ----------------------------------------------------------------------------
//  A "round to nearest" interval after Richard Harris: it operates in the DEFAULT
//  FPU rounding mode and inflates the bounds by a proportional 1+/-0.5eps factor
//  in adjust(). Using this variant (rather than one that switches the FPU
//  rounding mode) means no fragile _MM_SET_ROUNDING_MODE / fesetround control is
//  needed - it is portable across GCC / MSVC without touching the FPU mode.
//
//  Purpose: a FILTER. A predicate evaluates cheaply in interval arithmetic; if
//  the resulting sign is *determined* (the interval does not straddle 0) we trust
//  it; otherwise we fall back to the exact Expansion path. The safety invariant
//  this file must guarantee:
//
//      if sign_is_determined(I.sign())  then  convert_sign(I.sign()) is the
//      TRUE sign of the exact value.   (it may abstain, but never lies)
//
//  Compile with strict FP (-ffp-contract=off / /fp:strict): the interval bound
//  inflation reasons about IEEE round-to-nearest, which FMA contraction breaks.
// ============================================================================
#pragma once

#include "Expansion.h"

#include <limits>
#include <cmath>
#include <utility>

namespace sm {
namespace exact {

// Two-sided sign of an interval (Sign2). Encodes the (sign-of-lower,
// sign-of-upper) pair; the chosen values make the branchless sign() work.
enum Sign2 {
    SIGN2_ERROR = -1,  // NaN / inconsistent
    SIGN2_ZERO  =  0,  // [0,0]
    SIGN2_NP,          // lower<0, upper>0  -> straddles 0 (UNDETERMINED)
    SIGN2_PP,          // both > 0          -> POSITIVE
    SIGN2_ZP,          // lower==0,upper>0  -> undetermined (could be 0)
    SIGN2_NN,          // both < 0          -> NEGATIVE
    SIGN2_NZ           // lower<0,upper==0  -> undetermined
};

inline bool sign_is_determined(Sign2 s) {
    return s == SIGN2_ZERO || s == SIGN2_NN || s == SIGN2_PP;
}
inline bool sign_is_non_zero(Sign2 s) {
    return s == SIGN2_NN || s == SIGN2_PP;
}
// Precondition: sign_is_determined(s). Returns -1/0/+1.
inline int convert_sign(Sign2 s) {
    if (s == SIGN2_NN) return -1;
    if (s == SIGN2_PP) return  1;
    return 0;  // SIGN2_ZERO
}

// ----------------------------------------------------------------------------
//  Interval: round-to-nearest interval arithmetic.
// ----------------------------------------------------------------------------
class Interval {
public:
    Interval() : lb_(0.0), ub_(0.0) {}
    Interval(double x) : lb_(x), ub_(x) {}
    Interval(double l, double u) : lb_(l), ub_(u) {}  // no rounding-mode control
    Interval(const Interval&) = default;
    Interval& operator=(const Interval&) = default;

    // Conversion from an exact Expansion: accumulate components from the one of
    // largest magnitude down, expanding by one ulp when a component is below the
    // current resolution (convert an exact expansion to an interval).
    explicit Interval(const Expansion& e) { assign(e); }

    Interval& operator=(double rhs) { lb_ = ub_ = rhs; return *this; }

    Interval& assign(const Expansion& rhs) {
        std::size_t l = rhs.length();
        if (l == 0) { lb_ = ub_ = 0.0; return *this; }
        lb_ = ub_ = rhs[l - 1];
        for (int comp_idx = int(l) - 2; comp_idx >= 0; --comp_idx) {
            double comp = rhs[std::size_t(comp_idx)];
            if (comp > 0) {
                double nub = ub_ + comp;
                if (nub == ub_) {
                    ub_ = std::nextafter(ub_, std::numeric_limits<double>::infinity());
                    break;
                } else { ub_ = nub; adjust(); }
            } else {
                double nlb = lb_ + comp;
                if (nlb == lb_) {
                    lb_ = std::nextafter(lb_, -std::numeric_limits<double>::infinity());
                    break;
                } else { lb_ = nlb; adjust(); }
            }
        }
        return *this;
    }

    double inf() const { return lb_; }
    double sup() const { return ub_; }
    double estimate() const { return 0.5 * (lb_ + ub_); }
    bool is_nan() const { return !(lb_ == lb_) || !(ub_ == ub_); }

    Sign2 sign() const {
        if (is_nan()) return SIGN2_ERROR;
        int lz = int(lb_ == 0), ln = int(lb_ < 0), lp = int(lb_ > 0);
        int uz = int(ub_ == 0), un = int(ub_ < 0), up = int(ub_ > 0);
        Sign2 result = Sign2(
            ln * up * SIGN2_NP +
            lp * up * SIGN2_PP +
            lz * up * SIGN2_ZP +
            ln * un * SIGN2_NN +
            ln * uz * SIGN2_NZ);
        result = Sign2(int(result) +
                       int(result == SIGN2_ZERO && !(lz && uz)) * SIGN2_ERROR);
        return result;
    }

    Interval& negate() {
        lb_ = -lb_; ub_ = -ub_; std::swap(lb_, ub_);
        return *this;
    }

    Interval& operator+=(const Interval& x) {
        lb_ += x.lb_; ub_ += x.ub_; adjust(); return *this;
    }
    Interval& operator-=(const Interval& x) {
        lb_ -= x.ub_; ub_ -= x.lb_; adjust(); return *this;
    }
    Interval& operator*=(const Interval& x) {
        if (!is_nan() && !x.is_nan()) {
            double ll = lb_ * x.lb_, lu = lb_ * x.ub_;
            double ul = ub_ * x.lb_, uu = ub_ * x.ub_;
            if (!(ll == ll)) ll = 0.0;
            if (!(lu == lu)) lu = 0.0;
            if (!(ul == ul)) ul = 0.0;
            if (!(uu == uu)) uu = 0.0;
            if (lu < ll) std::swap(lu, ll);
            if (ul < ll) std::swap(ul, ll);
            if (uu < ll) std::swap(uu, ll);
            if (lu > uu) uu = lu;
            if (ul > uu) uu = ul;
            lb_ = ll; ub_ = uu;
            adjust();
        } else {
            lb_ = ub_ = std::numeric_limits<double>::quiet_NaN();
        }
        return *this;
    }

private:
    // Inflate [lb_,ub_] outward by the round-to-nearest error bound (Harris).
    void adjust() {
        constexpr double i = std::numeric_limits<double>::infinity();
        constexpr double e = std::numeric_limits<double>::epsilon();
        constexpr double m = std::numeric_limits<double>::min();  // smallest normal
        constexpr double l = 1.0 - e;
        constexpr double u = 1.0 + e;
        constexpr double em = e * m;

        if (lb_ == lb_ && ub_ == ub_ && (lb_ != ub_ || (lb_ != i && lb_ != -i))) {
            if (lb_ > ub_) std::swap(lb_, ub_);
            if (lb_ > m)       lb_ *= l;
            else if (lb_ < -m) lb_ *= u;
            else               lb_ -= em;
            if (ub_ > m)       ub_ *= u;
            else if (ub_ < -m) ub_ *= l;
            else               ub_ += em;
        } else {
            lb_ = ub_ = std::numeric_limits<double>::quiet_NaN();
        }
    }

    double lb_;
    double ub_;
};

inline Interval operator+(const Interval& a, const Interval& b) { Interval r = a; return r += b; }
inline Interval operator-(const Interval& a, const Interval& b) { Interval r = a; return r -= b; }
inline Interval operator*(const Interval& a, const Interval& b) { Interval r = a; return r *= b; }

} // namespace exact
} // namespace sm
