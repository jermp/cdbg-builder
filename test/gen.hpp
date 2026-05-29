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
// adjacent pair overlaps by exactly k bases -- the shared boundary
// k-mer that stitch glues on, GGCAT's k-base super-k-mer overlap).
// Internal split boundaries are marked OPEN; the true unitig's two
// endpoints are CLOSED. Each fragment is independently, randomly
// reverse-complemented (with its open-flag sides swapped) so the test
// exercises both orientations.
//
// Model in k-mer space: S has M = L-k+1 k-mers (indices 0..M-1).
// Fragment f covers an inclusive k-mer range [start_kmer, end_kmer];
// the boundary k-mer splits[f] is the LAST k-mer of fragment f and the
// FIRST k-mer of fragment f+1, so both fragments physically contain all
// k bases of it -> k-base overlap. `num_frags` is clamped so every
// fragment has >= 2 k-mers (the interior split indices live in
// [1, M-2], giving at most M-1 fragments).
inline std::vector<cdgb::stitchable_unitig> split_unitig(std::string const& S, uint32_t k,
                                                         uint32_t cid, uint64_t num_frags,
                                                         std::mt19937_64& rng) {
    using cdgb::stitchable_unitig;
    const uint64_t L = S.size();
    const uint64_t M = (L >= k) ? (L - k + 1) : 1;  // number of k-mers
    const uint64_t max_frags = std::max<uint64_t>(1, (M >= 2) ? (M - 1) : 1);
    if (num_frags < 1) num_frags = 1;
    if (num_frags > max_frags) num_frags = max_frags;

    // num_frags-1 distinct strictly-increasing split k-mer indices in
    // [1, M-2]. Each splits[f] is the boundary k-mer shared (in full) by
    // fragments f and f+1.
    std::vector<uint64_t> splits;
    if (num_frags > 1) {
        std::vector<uint64_t> cand;
        for (uint64_t j = 1; j + 1 <= M - 1; ++j) cand.push_back(j);  // [1, M-2]
        std::shuffle(cand.begin(), cand.end(), rng);
        cand.resize(num_frags - 1);
        std::sort(cand.begin(), cand.end());
        splits = std::move(cand);
    }

    std::vector<stitchable_unitig> frags;
    frags.reserve(num_frags);
    for (uint64_t f = 0; f < num_frags; ++f) {
        // Inclusive k-mer range; sequence spans [start_kmer, end_kmer+k).
        const uint64_t start_kmer = (f == 0) ? 0 : splits[f - 1];
        const uint64_t end_kmer = (f + 1 < num_frags) ? splits[f] : (M - 1);
        stitchable_unitig u;
        u.seq = S.substr(start_kmer, (end_kmer + k) - start_kmer);
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
    }
    return frags;
}

}  // namespace cdgb_test
