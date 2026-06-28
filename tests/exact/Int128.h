// ============================================================================
//  Int128.h - portable signed 128-bit integer for the exact-arithmetic oracles
// ----------------------------------------------------------------------------
//  The predicate/expansion tests validate their results against an exact
//  128-bit integer computation. GCC/Clang expose a native __int128, but MSVC
//  does not, so on compilers without it we provide a minimal two's-complement
//  128-bit type supporting exactly the operations the oracles use: construction
//  from a 64-bit integer, +, -, unary -, *, and comparison.
//
//  Multiplication keeps the low 128 bits (sign-agnostic in two's complement),
//  which is the true product whenever it fits in a signed 128-bit value - the
//  range the oracle inputs are chosen to stay within.
// ============================================================================
#pragma once

// Define SM_FORCE_PORTABLE_INT128 to exercise the portable struct even where a
// native __int128 is available (used to test the fallback).
#if defined(__SIZEOF_INT128__) && !defined(SM_FORCE_PORTABLE_INT128)

using i128 = __int128;

#else

#include <cstdint>

struct i128 {
    std::uint64_t lo;
    std::uint64_t hi;

    i128() : lo(0), hi(0) {}
    i128(long long v) : lo(std::uint64_t(v)), hi(v < 0 ? ~std::uint64_t(0) : 0) {}
    i128(int v) : i128(static_cast<long long>(v)) {}
    i128(unsigned v) : lo(v), hi(0) {}
    i128(unsigned long long v) : lo(v), hi(0) {}
    i128(std::uint64_t l, std::uint64_t h) : lo(l), hi(h) {}

    int sign() const {
        if (hi >> 63) return -1;
        return (hi | lo) ? 1 : 0;
    }

    explicit operator long double() const;
    explicit operator double() const { return double(static_cast<long double>(*this)); }
};

inline i128 operator-(const i128& a) {
    i128 r;
    r.lo = ~a.lo + 1;
    r.hi = ~a.hi + (r.lo == 0 ? 1 : 0);
    return r;
}

inline i128::operator long double() const {
    if (sign() < 0) return -static_cast<long double>(-(*this));
    return static_cast<long double>(hi) * 18446744073709551616.0L  // 2^64
         + static_cast<long double>(lo);
}
inline i128 operator+(const i128& a, const i128& b) {
    i128 r;
    r.lo = a.lo + b.lo;
    r.hi = a.hi + b.hi + (r.lo < a.lo ? 1 : 0);
    return r;
}
inline i128 operator-(const i128& a, const i128& b) { return a + (-b); }

inline i128 operator*(const i128& a, const i128& b) {
    // Full 128-bit product of the low limbs via 32-bit halves, plus the low
    // 64 bits of the cross terms (a.lo*b.hi + a.hi*b.lo) shifted into hi.
    const std::uint64_t M = 0xffffffffULL;
    std::uint64_t a0 = a.lo & M, a1 = a.lo >> 32;
    std::uint64_t b0 = b.lo & M, b1 = b.lo >> 32;
    std::uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
    std::uint64_t mid = (p00 >> 32) + (p01 & M) + (p10 & M);
    i128 r;
    r.lo = (p00 & M) | (mid << 32);
    r.hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32) + a.lo * b.hi + a.hi * b.lo;
    return r;
}

inline int cmp(const i128& a, const i128& b) { return (a - b).sign(); }
inline bool operator> (const i128& a, const i128& b) { return cmp(a, b) >  0; }
inline bool operator< (const i128& a, const i128& b) { return cmp(a, b) <  0; }
inline bool operator>=(const i128& a, const i128& b) { return cmp(a, b) >= 0; }
inline bool operator<=(const i128& a, const i128& b) { return cmp(a, b) <= 0; }
inline bool operator==(const i128& a, const i128& b) { return a.hi == b.hi && a.lo == b.lo; }
inline bool operator!=(const i128& a, const i128& b) { return !(a == b); }

#endif
