#pragma once

// Canonical ntHash over m-mers + sliding-window minimum, used for super-k-mer
// minimizer-bucketing ingest.
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

// Pick "best" minimizer length given k (taken from GGCAT).
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

// Sliding-window minimum of the canonical m-mer hash over a (k-1)-mer's window
// of W = k-m m-mers, returning the current minimum's hash. O(1) SPACE: it keeps
// only the running (min hash, min position) plus the rolling ntHash state -- no
// ring buffer, no modulo. It OWNS the ntHash roll (sshash minimizer_iterator
// style):
//
//   - init():      roll over the first window's W m-mers, cache the min.
//   - advance():   roll one m-mer right. If the cached min is still in the
//                  window (common path): O(1) -- compare the new m-mer and keep
//                  the running min. If the cached min just left the window:
//                  RE-SCAN by re-rolling the W m-mers of the new window from a
//                  local ntHash and taking their min: O(W).
//
// For random DNA the expiry rate is ~1/W per step, so the re-scan amortises to
// O(1) with no per-step division and no buffer. Ties keep the RIGHTMOST (newest)
// position (`<=`), matching the previous ring-buffer implementation exactly, so
// the chosen minimizer -- and thus the bucketing and all downstream output -- is
// unchanged.
struct windowed_min {
    uint8_t const* bases = nullptr;  // 2-bit-encoded run bases
    uint32_t m = 0;
    int32_t W = 0;                   // window size in m-mers (= k - m)
    int32_t pos = 0;                 // absolute index of the current rightmost m-mer
    uint64_t fwd = 0, rc = 0;        // rolling ntHash at `pos`
    uint64_t cur_min = ~uint64_t(0);
    int32_t cur_min_pos = -1;

    // Roll over the first window (m-mers at positions [0, W-1]) and cache its
    // min. Returns false iff a non-ACGT base is hit (matches nthash_init).
    bool init(uint8_t const* b, uint32_t mm, int32_t window) {
        assert(window > 0);
        bases = b;
        m = mm;
        W = window;
        if (!nthash_init(bases, m, fwd, rc)) return false;
        cur_min = canonical_mhash(fwd, rc);
        cur_min_pos = 0;
        for (int32_t i = 1; i < W; ++i) {
            nthash_roll(bases[i - 1], bases[i + (int32_t)m - 1], m, fwd, rc);
            uint64_t h = canonical_mhash(fwd, rc);
            if (h <= cur_min) {  // <= keeps the rightmost (newest) on ties
                cur_min = h;
                cur_min_pos = i;
            }
        }
        pos = W - 1;
        return true;
    }

    // Slide the window right by one m-mer.
    void advance() {
        ++pos;
        nthash_roll(bases[pos - 1], bases[pos + (int32_t)m - 1], m, fwd, rc);
        if (cur_min_pos + W <= pos) {
            // Cached min left the window: re-scan [pos-W+1, pos] with a local
            // ntHash (leaves fwd/rc -- the steady-state roll -- untouched).
            int32_t start = pos - W + 1;
            uint64_t f, r;
            nthash_init(bases + start, m, f, r);
            uint64_t best = canonical_mhash(f, r);
            int32_t best_pos = start;
            for (int32_t p = start + 1; p <= pos; ++p) {
                nthash_roll(bases[p - 1], bases[p + (int32_t)m - 1], m, f, r);
                uint64_t hp = canonical_mhash(f, r);
                if (hp <= best) {
                    best = hp;
                    best_pos = p;
                }
            }
            cur_min = best;
            cur_min_pos = best_pos;
        } else {
            uint64_t h = canonical_mhash(fwd, rc);
            if (h <= cur_min) {
                cur_min = h;
                cur_min_pos = pos;
            }
        }
    }

    uint64_t min_hash() const { return cur_min; }
};

}  // namespace cdbg
