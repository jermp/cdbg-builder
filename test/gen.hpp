#pragma once

// Synthetic data generators for the per-phase tests.
//
// The tests use a "generate-by-construction" oracle: build a known
// correct answer, derive the phase's input from it, run the phase,
// and assert the phase reconstructs the answer. All randomness is
// seeded so failures are reproducible.

#include <algorithm>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "bucket_walker.hpp"  // cdgb::stitchable_unitig, UNITIG_OPEN_*

namespace cdgb_test {

inline char rc_base(char c) {
    switch (c) {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        default:
            return 'N';
    }
}

inline std::string revcomp(std::string const& s) {
    std::string out(s.size(), 'N');
    for (uint64_t i = 0; i < s.size(); ++i) out[i] = rc_base(s[s.size() - 1 - i]);
    return out;
}

// Lexicographic min of a sequence and its reverse complement. Used to
// compare stitched output against truth modulo orientation.
inline std::string canonical(std::string const& s) {
    std::string r = revcomp(s);
    return s <= r ? s : r;
}

inline std::string random_dna(uint64_t len, std::mt19937_64& rng) {
    static const char b[4] = {'A', 'C', 'G', 'T'};
    std::string s(len, 'A');
    for (uint64_t i = 0; i < len; ++i) s[i] = b[rng() & 3u];
    return s;
}

// Split a true unitig `S` into `num_frags` overlapping fragments (each
// adjacent pair overlaps by exactly k-1 bases -- the shared junction
// (k-1)-mer that stitch glues on). Internal split boundaries are
// marked OPEN; the true unitig's two endpoints are CLOSED. Each
// fragment is independently, randomly reverse-complemented (with its
// open-flag sides swapped) so the test exercises both orientations.
//
// `num_frags` is clamped to the number of k-mers in S (you can't have
// more fragments than k-mers). All fragments are >= k bases.
inline std::vector<cdgb::stitchable_unitig> split_unitig(std::string const& S, uint32_t k,
                                                         uint32_t cid, uint64_t num_frags,
                                                         std::mt19937_64& rng) {
    using cdgb::stitchable_unitig;
    const uint64_t L = S.size();
    const uint64_t max_frags = std::max<uint64_t>(1, L >= k ? (L - k + 1) : 1);
    if (num_frags < 1) num_frags = 1;
    if (num_frags > max_frags) num_frags = max_frags;

    // Internal boundaries: num_frags-1 strictly-increasing positions in
    // [k, L-1]. Fragment f spans [start_f, end_f); start_{f+1} =
    // end_f - (k-1) so consecutive fragments share a (k-1)-mer.
    std::vector<uint64_t> bounds;
    if (num_frags > 1) {
        std::vector<uint64_t> cand;
        for (uint64_t b = k; b + 1 <= L; ++b) cand.push_back(b);  // [k, L-1]
        std::shuffle(cand.begin(), cand.end(), rng);
        cand.resize(num_frags - 1);
        std::sort(cand.begin(), cand.end());
        bounds = std::move(cand);
    }

    std::vector<stitchable_unitig> frags;
    frags.reserve(num_frags);
    uint64_t start = 0;
    for (uint64_t f = 0; f < num_frags; ++f) {
        const uint64_t end = (f + 1 < num_frags) ? bounds[f] : L;
        stitchable_unitig u;
        u.seq = S.substr(start, end - start);
        u.cid = cid;
        uint8_t flags = 0;
        if (f > 0) flags |= cdgb::UNITIG_OPEN_LEFT;               // internal left boundary
        if (f + 1 < num_frags) flags |= cdgb::UNITIG_OPEN_RIGHT;  // internal right boundary
        // Randomly flip orientation; swap the open-flag sides to match.
        if (rng() & 1u) {
            u.seq = revcomp(u.seq);
            uint8_t nf = 0;
            if (flags & cdgb::UNITIG_OPEN_LEFT) nf |= cdgb::UNITIG_OPEN_RIGHT;
            if (flags & cdgb::UNITIG_OPEN_RIGHT) nf |= cdgb::UNITIG_OPEN_LEFT;
            flags = nf;
        }
        u.open_flags = flags;
        frags.push_back(std::move(u));
        start = end - (k - 1);  // next fragment overlaps by k-1
    }
    return frags;
}

}  // namespace cdgb_test
