// ============================================================================
//  test_interval.cpp - test for the round-to-nearest interval filter.
// ----------------------------------------------------------------------------
//  The ONE property a filter must never violate:
//
//      if the interval sign is DETERMINED, it equals the TRUE (exact) sign.
//
//  A filter is allowed to abstain (straddle zero) as often as it likes; it is
//  NEVER allowed to return a definite WRONG sign. We verify this against:
//    * an i128 exact oracle for integer expressions a*b - c*d,
//    * the exact Expansion path for Expansion->Interval conversion,
//    * a catastrophic-cancellation family where the filter MUST abstain.
//
//  Because the interval is round-to-nearest, this test needs NO
//  FPU rounding-mode control - which is exactly why that variant was chosen.
// ============================================================================
#include "simplemesh/exact/Interval.h"
#include "simplemesh/exact/Expansion.h"
#include "Int128.h"

#include <cstdint>
#include <iostream>
#include <random>
#include <string>

using namespace sm::exact;

static int failures = 0;
#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } \
        else         { std::cout << "ok:   " << msg << "\n"; }             \
    } while (0)

static int isign(i128 v) { return (v > 0) - (v < 0); }

int main() {
    // --- basics ---------------------------------------------------------------
    CHECK(convert_sign(Interval(3.0).sign()) == 1, "point interval +3 is PP");
    CHECK(convert_sign(Interval(-3.0).sign()) == -1, "point interval -3 is NN");
    CHECK(Interval(0.0).sign() == SIGN2_ZERO, "point interval 0 is ZERO");
    {
        // A computed zero is NOT certified by the filter: adjust() inflates
        // [0,0] to [-eps*min, +eps*min], so the filter ABSTAINS and defers the
        // true-zero decision to the exact path. (A filter proves non-zero; zero
        // always needs exact fallback.) The safety requirement is only that it
        // must not CLAIM a non-zero sign.
        Interval s = Interval(1.0) - Interval(1.0);
        CHECK(!sign_is_non_zero(s.sign()),
              "1-1: filter does not claim a non-zero sign (abstains on true zero)");
    }
    // containment: exact small value lies inside the interval
    {
        Interval r = Interval(7.0) * Interval(11.0) - Interval(70.0);  // = 7
        CHECK(r.inf() <= 7.0 && 7.0 <= r.sup(), "interval contains exact 7");
        CHECK(convert_sign(r.sign()) == 1, "7 is positive (determined)");
    }

    // --- SAFETY INVARIANT on random integer expressions a*b - c*d -------------
    std::mt19937_64 rng(20260628);
    std::uniform_int_distribution<int64_t> D(-(int64_t(1) << 31), (int64_t(1) << 31));

    const int N = 300000;
    int wrong_determined = 0;  // MUST stay 0
    int contain_fail = 0;      // MUST stay 0
    int abstained = 0;         // informational (fast path missed)
    for (int k = 0; k < N; ++k) {
        int64_t a = D(rng), b = D(rng), c = D(rng), d = D(rng);
        i128 exact = (i128)a * b - (i128)c * d;
        int es = isign(exact);

        Interval r = Interval(double(a)) * Interval(double(b))
                   - Interval(double(c)) * Interval(double(d));
        Sign2 s = r.sign();
        if (sign_is_determined(s)) {
            if (convert_sign(s) != es) ++wrong_determined;
        } else {
            ++abstained;
        }
        // containment: the true value (as long double) must lie within bounds.
        long double v = (long double)exact;
        if (!(r.inf() <= v && v <= r.sup())) ++contain_fail;
    }
    CHECK(wrong_determined == 0,
          "no WRONG determined sign in " + std::to_string(N) + " cases (filter safety)");
    CHECK(contain_fail == 0,
          "interval contains the exact value in all " + std::to_string(N) + " cases");
    std::cout << "      filter abstained on " << abstained << "/" << N
              << " random cases (fast path taken otherwise)\n";

    // --- TEETH: catastrophic cancellation. n^2-(n-1)(n+1)=1, products ~2^62.
    //     The filter must abstain (cannot resolve), and must NEVER say <=0.
    int wrong_hard = 0, abstain_hard = 0, trials = 0;
    for (int64_t n = (int64_t(1) << 31); n < (int64_t(1) << 31) + 2000; ++n) {
        ++trials;
        Interval r = Interval(double(n)) * Interval(double(n))
                   - Interval(double(n - 1)) * Interval(double(n + 1));
        Sign2 s = r.sign();
        if (sign_is_determined(s)) {
            if (convert_sign(s) != 1) ++wrong_hard;  // truth is +1
        } else {
            ++abstain_hard;
        }
    }
    CHECK(wrong_hard == 0,
          "filter never returns a wrong determined sign on catastrophic family");
    CHECK(abstain_hard > 0,
          "filter abstains on " + std::to_string(abstain_hard) + "/" +
          std::to_string(trials) + " catastrophic cases (correctly defers to exact)");

    // --- Expansion -> Interval conversion: determined sign must match exact ----
    int conv_wrong = 0, conv_contain = 0;
    for (int k = 0; k < 100000; ++k) {
        int64_t a = D(rng), b = D(rng), c = D(rng), d = D(rng);
        i128 exact = (i128)a * b - (i128)c * d;
        Expansion ex = Expansion::diff(Expansion::product(double(a), double(b)),
                                       Expansion::product(double(c), double(d)));
        Interval r(ex);  // convert exact -> interval
        Sign2 s = r.sign();
        if (sign_is_determined(s) && convert_sign(s) != isign(exact)) ++conv_wrong;
        long double v = (long double)exact;
        if (!(r.inf() <= v && v <= r.sup())) ++conv_contain;
    }
    CHECK(conv_wrong == 0, "Expansion->Interval: no wrong determined sign (100k)");
    CHECK(conv_contain == 0, "Expansion->Interval: exact value always contained (100k)");

    if (failures) {
        std::cerr << "\n" << failures << " interval check(s) FAILED\n";
        return 1;
    }
    std::cout << "\nOK: interval filter is sound (abstains, never lies).\n";
    return 0;
}
