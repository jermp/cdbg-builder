#pragma once

#include <cassert>
#include <cstdint>
#include <cstring>
#include <functional>
#include <string>

namespace cdbg {

// 2-bit packed k-mer, stored canonical (min of forward / RC).
// Encoding: A=0, C=1, G=2, T=3. Position 0 is the lowest two bits, so the
// first nucleotide of the k-mer ends up in the highest two used bits.
//
// This convention matches SSHash's 2-bit canonical k-mer encoding used by
// Fulgor downstream (forward k-mer is read left-to-right and shifted in
// from the LSB side; RC built by inverting and reversing the 2-bit pairs).
//
// kmer_int_t width is a COMPILE-TIME choice (mirrors SSHash). The default is
// uint64_t, which holds k <= 31 (2*31 = 62 bits) -- the common case -- and uses
// cheaper 64-bit arithmetic plus half the memory in every k-mer-keyed structure
// (notably bucket-process's kmer_info map). Define CDBG_LARGE_K (the CMake
// option of the same name, which adds -DCDBG_LARGE_K) to widen to __uint128_t
// and allow k <= 63. The on-disk formats adapt automatically: each k-mer is
// serialized in ceil(2k/8) bytes or sizeof(kmer_int_t), both derived from the
// active width.
#ifdef CDBG_LARGE_K
using kmer_int_t = __uint128_t;
constexpr uint32_t MAX_K = 63;
#else
using kmer_int_t = uint64_t;
constexpr uint32_t MAX_K = 31;
#endif

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
    // A shift equal to the type width is UB, so special-case the full-width k
    // (k==32 for uint64_t, k==64 for __uint128_t). Valid k <= MAX_K never hits
    // it, but keep it robust.
    if (2 * k >= sizeof(kmer_int_t) * 8) return ~kmer_int_t(0);
    return (kmer_int_t(1) << (2 * k)) - 1;
}

// reverse-complement a k-mer in O(1), using SWAR over the 2-bit packed word
// (the SSHash trick, https://github.com/jermp/sshash/blob/master/include/kmer.hpp).
// `~x` complements every base in parallel (A<->T, C<->G is XOR 3 = XOR all-ones),
// then a constant number of mask-shift steps reverse the order of the 2-bit
// groups across the whole word; the final shift drops the unused high padding so
// the k used groups land in the low 2k bits. This is hot: bucket-process's
// local_ext_mask calls it 1x directly plus 8x through canonical() per k-mer, so
// the old O(k) loop (31 iterations at k=31) dominated the resolve walk.
//
// Templated on the word type so `if constexpr` actually DISCARDS the unused
// width's branch -- otherwise the 128-bit `<< 64` shifts would warn (and be
// ill-formed) when kmer_int_t is uint64_t, and vice-versa.
template <typename T = kmer_int_t>
inline T reverse_complement(T x, uint32_t k) {
    if constexpr (sizeof(T) == sizeof(uint64_t)) {
        x = ~x;
        x = ((x >> 2) & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
        x = ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((x & 0x0F0F0F0F0F0F0F0FULL) << 4);
        x = ((x >> 8) & 0x00FF00FF00FF00FFULL) | ((x & 0x00FF00FF00FF00FFULL) << 8);
        x = ((x >> 16) & 0x0000FFFF0000FFFFULL) | ((x & 0x0000FFFF0000FFFFULL) << 16);
        x = (x >> 32) | (x << 32);
        return x >> (64 - 2 * k);
    } else {
        // 128-bit width (CDBG_LARGE_K): same algorithm with replicated masks.
        auto rep = [](uint64_t b) -> T { return ((T)b << 64) | b; };
        const T M2 = rep(0x3333333333333333ULL);
        const T M4 = rep(0x0F0F0F0F0F0F0F0FULL);
        const T M8 = rep(0x00FF00FF00FF00FFULL);
        const T M16 = rep(0x0000FFFF0000FFFFULL);
        const T M32 = rep(0x00000000FFFFFFFFULL);
        x = ~x;
        x = ((x >> 2) & M2) | ((x & M2) << 2);
        x = ((x >> 4) & M4) | ((x & M4) << 4);
        x = ((x >> 8) & M8) | ((x & M8) << 8);
        x = ((x >> 16) & M16) | ((x & M16) << 16);
        x = ((x >> 32) & M32) | ((x & M32) << 32);
        x = (x >> 64) | (x << 64);
        return x >> (128 - 2 * k);
    }
}

inline kmer_int_t canonical(kmer_int_t fwd, uint32_t k) {
    kmer_int_t rc = reverse_complement(fwd, k);
    return fwd <= rc ? fwd : rc;
}

// shift and append: take a k-mer, drop first symbol, append `nt` (2-bit) at the right.
inline kmer_int_t shift_append(kmer_int_t x, uint8_t nt, uint32_t k) {
    return ((x << 2) | (kmer_int_t)(nt & 3)) & kmer_mask(k);
}

inline std::string kmer_to_string(kmer_int_t x, uint32_t k) {
    std::string s(k, 'N');
    for (uint32_t i = 0; i < k; ++i) {
        s[k - 1 - i] = twobit_to_nuc((uint8_t)(x & 3));
        x >>= 2;
    }
    return s;
}

struct kmer_hasher {
    // splitmix64-quality output; tell ankerl::unordered_dense it doesn't
    // need to re-mix.
    using is_avalanching = void;

    size_t operator()(kmer_int_t x) const noexcept {
        auto mix = [](uint64_t z) {
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            z = z ^ (z >> 31);
            return z;
        };
        uint64_t lo = (uint64_t)x;
        // For the 128-bit width, mix both halves (identical to the prior hash).
        // For uint64_t, x >> 64 is UB, so hash the single word.
        if constexpr (sizeof(kmer_int_t) > sizeof(uint64_t)) {
            uint64_t hi = (uint64_t)(x >> 64);
            return (size_t)(mix(lo) ^ (mix(hi + 0x9e3779b97f4a7c15ULL) << 1));
        } else {
            return (size_t)mix(lo);
        }
    }
};

}  // namespace cdbg
