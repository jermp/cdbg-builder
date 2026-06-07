#pragma once

// Canonical ntHash over m-mers + sliding-window minimum, used for super-k-mer
// minimizer-bucketing ingest (mirrors GGCAT's approach in
// crates/hashes/src/cn_nthash.rs and crates/hashes/src/rolling/batch_minqueue.rs).
//
// We hash an m-mer and its reverse complement simultaneously and take the
// minimum, so the hash is symmetric under RC. That makes the minimizer of a
// k-mer well-defined regardless of orientation: identical canonical k-mers
// always select the same minimizer, hence the same bucket.

#include <cassert>
#include <cstdint>
#include <cstring>
#include <vector>

#include "kmer.hpp"

namespace cdbg {

// Per-base ntHash seeds (the canonical ntHash table by Mohamadi et al.).
inline constexpr uint64_t NT_SEED[4] = {
    0x3c8bfbb395c60474ULL,  // A
    0x3193c18562a02b4cULL,  // C
    0x20323ed082572324ULL,  // G
    0x295549f54be24456ULL   // T
};

inline uint64_t rol64(uint64_t x, unsigned r) {
    r &= 63u;
    return (x << r) | (x >> ((64 - r) & 63));
}
inline uint64_t ror64(uint64_t x, unsigned r) {
    r &= 63u;
    return (x >> r) | (x << ((64 - r) & 63));
}
inline uint64_t comp_seed(uint8_t b) { return NT_SEED[b ^ 3u]; }

// Initialize fwd/rc ntHash for the m-mer at `s[0..m]`. Returns false if any
// base is non-ACGT (in which case fwd/rc are not set). rc-hash of `s[0..m]`
// is defined as the fwd ntHash of revcomp(s[0..m]).
inline bool nthash_init(uint8_t const* s, uint32_t m, uint64_t& fwd, uint64_t& rc) {
    uint64_t f = 0, r = 0;
    for (uint32_t i = 0; i < m; ++i) {
        uint8_t b = s[i];
        if (b > 3) return false;
        f = rol64(f, 1) ^ NT_SEED[b];
    }
    for (uint32_t i = 0; i < m; ++i) {
        uint8_t b = s[m - 1 - i] ^ 3u;
        r = rol64(r, 1) ^ NT_SEED[b];
    }
    fwd = f;
    rc = r;
    return true;
}

// Roll one position right: drop `out_b` (leftmost of old window), append `in_b`.
inline void nthash_roll(uint8_t out_b, uint8_t in_b, uint32_t m, uint64_t& fwd, uint64_t& rc) {
    fwd = rol64(fwd, 1) ^ rol64(NT_SEED[out_b], m) ^ NT_SEED[in_b];
    rc = ror64(rc, 1) ^ ror64(NT_SEED[out_b ^ 3u], 1) ^ rol64(NT_SEED[in_b ^ 3u], m - 1);
}

inline uint64_t canonical_mhash(uint64_t fwd, uint64_t rc) { return fwd <= rc ? fwd : rc; }

// GGCAT's compute_best_m for picking a minimizer length given k.
// (crates/utils/src/lib.rs:29)
inline uint32_t compute_best_m(uint32_t k) {
    if (k <= 13) return (k / 2 > k - 4) ? (k / 2) : (k > 4 ? (k - 4) : 1);
    if (k <= 15) return 9;
    if (k <= 21) return 10;
    if (k <= 30) return 11;
    if (k <= 37) return 12;
    if (k <= 42) return 13;
    if (k <= 64) return 14;
    return (k + 2) / 4;
}

// Sliding-window minimum over a stream of (uint64_t hash, int32_t pos) pairs,
// returning the current minimum's hash. Rescan-on-expiry over a fixed-size ring
// buffer of the last W m-mer hashes plus a cached (cur_min, cur_min_pos):
//
//   - Common path (~all pushes): O(1) — store h at pos % W in the ring,
//     and if h <= cur_min update the cached min.
//   - When the cached min's position falls out of the window, rescan the
//     W ring slots for the new min: O(W).
//
// For random DNA the expected expiration rate is ~1/W per push, so the
// amortised cost stays O(1) with a small constant: no allocations, a dense
// linear rescan, and a branch-predictable common path.
//
// W = k - m + 1, bounded by MAX_K + 1 = 64 for the supported k range.
// The 64-slot static array fits in a single cache line.
struct windowed_min {
    static constexpr int32_t MAX_W = 64;
    uint64_t hashes[MAX_W];
    int32_t window_size = 0;
    uint64_t cur_min = ~uint64_t(0);
    int32_t cur_min_pos = -1;

    void reset(int32_t window) {
        assert(window > 0 and window <= MAX_W);
        window_size = window;
        cur_min = ~uint64_t(0);
        cur_min_pos = -1;
    }

    void push(uint64_t h, int32_t pos) {
        hashes[(uint32_t)pos % (uint32_t)window_size] = h;
        // Did the previously-tracked minimum just slide out of the window?
        if (cur_min_pos + window_size <= pos) {
            int32_t start = pos - window_size + 1;
            if (start < 0) start = 0;
            uint64_t best = hashes[(uint32_t)start % (uint32_t)window_size];
            int32_t best_pos = start;
            for (int32_t p = start + 1; p <= pos; ++p) {
                uint64_t hp = hashes[(uint32_t)p % (uint32_t)window_size];
                if (hp <= best) {  // <= keeps the rightmost (newest) position on ties
                    best = hp;
                    best_pos = p;
                }
            }
            cur_min = best;
            cur_min_pos = best_pos;
        } else if (h <= cur_min) {
            cur_min = h;
            cur_min_pos = pos;
        }
    }

    uint64_t min_hash() const { return cur_min; }
};

}  // namespace cdbg
