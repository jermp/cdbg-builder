#pragma once

#include <cassert>
#include <cstdint>
#include <cstring>
#include <functional>
#include <string>

namespace cdgb {

// 2-bit packed k-mer for k <= 63, stored canonical (min of forward / RC).
// Encoding: A=0, C=1, G=2, T=3. Position 0 is the lowest two bits, so the
// first nucleotide of the k-mer ends up in the highest two used bits.
//
// This convention matches SSHash's 2-bit canonical k-mer encoding used by
// Fulgor downstream (forward k-mer is read left-to-right and shifted in
// from the LSB side; RC built by inverting and reversing the 2-bit pairs).

using kmer_int_t = __uint128_t;

constexpr uint32_t MAX_K = 63;

inline uint8_t nuc_to_2bit(char c) {
    switch (c) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return 0xff;
    }
}

inline char twobit_to_nuc(uint8_t x) {
    static constexpr char tab[4] = {'A', 'C', 'G', 'T'};
    return tab[x & 3];
}

inline kmer_int_t kmer_mask(uint32_t k) {
    if (k == 64) return ~kmer_int_t(0);
    return (kmer_int_t(1) << (2 * k)) - 1;
}

// reverse-complement a k-mer. O(k); we use this rarely (we maintain RC
// incrementally during streaming iteration).
inline kmer_int_t reverse_complement(kmer_int_t x, uint32_t k) {
    kmer_int_t y = 0;
    for (uint32_t i = 0; i < k; ++i) {
        uint8_t s = (uint8_t)(x & 3);
        y = (y << 2) | (kmer_int_t)(s ^ 3);
        x >>= 2;
    }
    return y;
}

inline kmer_int_t canonical(kmer_int_t fwd, uint32_t k) {
    kmer_int_t rc = reverse_complement(fwd, k);
    return fwd <= rc ? fwd : rc;
}

inline bool is_canonical(kmer_int_t fwd, uint32_t k) { return fwd <= reverse_complement(fwd, k); }

// shift and append: take a k-mer, drop first symbol, append `nt` (2-bit) at the right.
inline kmer_int_t shift_append(kmer_int_t x, uint8_t nt, uint32_t k) {
    return ((x << 2) | (kmer_int_t)(nt & 3)) & kmer_mask(k);
}

// shift and prepend: drop last symbol, prepend `nt` at the left.
inline kmer_int_t shift_prepend(kmer_int_t x, uint8_t nt, uint32_t k) {
    return ((x >> 2) | ((kmer_int_t)(nt & 3) << (2 * (k - 1)))) & kmer_mask(k);
}

inline std::string kmer_to_string(kmer_int_t x, uint32_t k) {
    std::string s(k, 'N');
    for (uint32_t i = 0; i < k; ++i) {
        s[k - 1 - i] = twobit_to_nuc((uint8_t)(x & 3));
        x >>= 2;
    }
    return s;
}

inline bool string_to_kmer(const char* s, uint32_t k, kmer_int_t& out) {
    kmer_int_t x = 0;
    for (uint32_t i = 0; i < k; ++i) {
        uint8_t v = nuc_to_2bit(s[i]);
        if (v == 0xff) return false;
        x = (x << 2) | v;
    }
    out = x;
    return true;
}

struct kmer_hasher {
    // splitmix64-quality output; tell ankerl::unordered_dense it doesn't
    // need to re-mix.
    using is_avalanching = void;

    size_t operator()(kmer_int_t x) const noexcept {
        // splitmix-style 64x64 mix on the two halves, combined.
        uint64_t lo = (uint64_t)x;
        uint64_t hi = (uint64_t)(x >> 64);
        auto mix = [](uint64_t z) {
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            z = z ^ (z >> 31);
            return z;
        };
        return (size_t)(mix(lo) ^ (mix(hi + 0x9e3779b97f4a7c15ULL) << 1));
    }
};

}  // namespace cdgb
