// ============================================================================
//  test_expansion.cpp - unit test for the multi-precision Expansion type.
// ----------------------------------------------------------------------------
//  Oracle: every test value is an integer exactly representable as a double, so
//  the EXACT result of the expressions we build is computable in i128. We
//  compare:
//     * Expansion.sign()      MUST equal the i128 exact sign   (zero misses)
//     * naive double          is ALLOWED to disagree, and we assert it DOES on
//                             some inputs - proving the test exercises the
//                             catastrophic-cancellation path the port exists for.
//
//  Integers are drawn up to 2^31, so products reach ~2^62: double's 53-bit
//  mantissa cannot hold them exactly, forcing rounding. The exact sum of two
//  such products still fits in a signed i128.
// ============================================================================
#include "simplemesh/exact/Expansion.h"
#include "Int128.h"

#include <cstdint>
#include <iostream>
#include <random>
#include <string>

using sm::exact::Expansion;
using sm::exact::ExpansionNt;

static int failures = 0;
#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } \
        else         { std::cout << "ok:   " << msg << "\n"; }             \
    } while (0)

static int isign(i128 v) { return (v > 0) - (v < 0); }

int main() {
    // --- basic identities -----------------------------------------------------
    CHECK(Expansion(0.0).sign() == 0, "zero has sign 0");
    CHECK(Expansion(3.0).sign() == 1, "positive sign");
    CHECK(Expansion(-3.0).sign() == -1, "negative sign");
    CHECK(Expansion::sum(1.0, 2.0).estimate() == 3.0, "1+2 == 3");
    CHECK(Expansion::diff(2.0, 5.0).estimate() == -3.0, "2-5 == -3");
    CHECK(Expansion::product(6.0, 7.0).estimate() == 42.0, "6*7 == 42");

    // A classic exactness case double gets wrong: (2^53 + 1) - 2^53 == 1, but
    // fl(2^53 + 1) == 2^53, so naive double yields 0.
    {
        double big = 9007199254740992.0;  // 2^53
        double naive = (big + 1.0) - big;
        Expansion e = Expansion::diff(Expansion::sum(big, 1.0), Expansion(big));
        CHECK(naive == 0.0, "naive double loses the +1 (sanity: cancellation real)");
        CHECK(e.sign() == 1 && e.estimate() == 1.0, "expansion keeps (2^53+1)-2^53 == 1");
    }

    // --- randomized differential test: sign(a*b - c*d) ------------------------
    std::mt19937_64 rng(20260628);
    std::uniform_int_distribution<int64_t> D(-(int64_t(1) << 31), (int64_t(1) << 31));

    const int N = 200000;
    int exp_mismatch = 0;  // expansion vs exact  (MUST stay 0)

    for (int i = 0; i < N; ++i) {
        int64_t a = D(rng), b = D(rng), c = D(rng), d = D(rng);
        i128 exact = (i128)a * b - (i128)c * d;
        Expansion ab = Expansion::product(double(a), double(b));
        Expansion cd = Expansion::product(double(c), double(d));
        if (Expansion::diff(ab, cd).sign() != isign(exact)) ++exp_mismatch;
    }
    CHECK(exp_mismatch == 0,
          "expansion sign(a*b-c*d) exact on " + std::to_string(N) + " random cases");

    // --- deterministic catastrophic-cancellation family (gives the test teeth) -
    // Identity: n^2 - (n-1)(n+1) = 1, exactly. Near n = 2^31 the products are
    // ~2^62; double's spacing there is ~2^10, so fl(n^2) == fl(n^2-1) and naive
    // double yields 0 (WRONG sign). The expansion must still report +1.
    int exp_fail_hard = 0, double_fail_hard = 0, trials = 0;
    for (int64_t n = (int64_t(1) << 31); n < (int64_t(1) << 31) + 2000; ++n) {
        ++trials;
        i128 exact = (i128)n * n - (i128)(n - 1) * (n + 1);  // == 1
        Expansion nn  = Expansion::product(double(n), double(n));
        Expansion nm  = Expansion::product(double(n - 1), double(n + 1));
        if (Expansion::diff(nn, nm).sign() != isign(exact)) ++exp_fail_hard;
        double naive = double(n) * double(n) - double(n - 1) * double(n + 1);
        if (((naive > 0) - (naive < 0)) != isign(exact)) ++double_fail_hard;
    }
    CHECK(exp_fail_hard == 0,
          "expansion exact on " + std::to_string(trials) + " catastrophic-cancel cases");
    CHECK(double_fail_hard > 0,
          "naive double WRONG on " + std::to_string(double_fail_hard) + "/" +
          std::to_string(trials) + " catastrophic cases (proves the test has teeth)");

    // --- 2x2 determinant form (the shape orient2d uses) ----------------------
    // det = (ax-cx)*(by-cy) - (ay-cy)*(bx-cx), all integer coords.
    int det_mismatch = 0;
    std::uniform_int_distribution<int64_t> C(-(int64_t(1) << 26), (int64_t(1) << 26));
    for (int i = 0; i < 100000; ++i) {
        int64_t ax = C(rng), ay = C(rng), bx = C(rng), by = C(rng), cx = C(rng), cy = C(rng);
        i128 exact = (i128)(ax - cx) * (by - cy) - (i128)(ay - cy) * (bx - cx);

        ExpansionNt EX_ax{double(ax)}, EX_ay{double(ay)};
        ExpansionNt EX_bx{double(bx)}, EX_by{double(by)};
        ExpansionNt EX_cx{double(cx)}, EX_cy{double(cy)};
        ExpansionNt t1 = (EX_ax - EX_cx) * (EX_by - EX_cy);
        ExpansionNt t2 = (EX_ay - EX_cy) * (EX_bx - EX_cx);
        ExpansionNt det = t1 - t2;
        if (det.sign() != isign(exact)) ++det_mismatch;
    }
    CHECK(det_mismatch == 0, "ExpansionNt 2x2 determinant sign exact on 100k cases");

    if (failures) {
        std::cerr << "\n" << failures << " expansion check(s) FAILED\n";
        return 1;
    }
    std::cout << "\nOK: Expansion reproduces exact integer arithmetic.\n";
    return 0;
}
