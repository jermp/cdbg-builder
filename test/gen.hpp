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

#include "phase2_bucket_process/bucket_walker.hpp"  // cdbg::stitchable_unitig, UNITIG_OPEN_*

namespace cdbg_test {

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
inline std::vector<cdbg::stitchable_unitig> split_unitig(std::string const& S, uint32_t k,
                                                         uint32_t cid, uint64_t num_frags,
                                                         std::mt19937_64& rng) {
    using cdbg::stitchable_unitig;
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
        u.set_mono(cid, k);
        uint8_t flags = 0;
        if (f > 0) flags |= cdbg::UNITIG_OPEN_LEFT;               // internal left boundary
        if (f + 1 < num_frags) flags |= cdbg::UNITIG_OPEN_RIGHT;  // internal right boundary
        // Randomly flip orientation; swap the open-flag sides to match.
        if (rng() & 1u) {
            u.seq = revcomp(u.seq);
            uint8_t nf = 0;
            if (flags & cdbg::UNITIG_OPEN_LEFT) nf |= cdbg::UNITIG_OPEN_RIGHT;
            if (flags & cdbg::UNITIG_OPEN_RIGHT) nf |= cdbg::UNITIG_OPEN_LEFT;
            flags = nf;
        }
        u.open_flags = flags;
        frags.push_back(std::move(u));
    }
    return frags;
}

// Shared-cid branchy oracle.
//
// The distinct-cid oracle above structurally cannot catch same-cid false
// joins: with a unique cid per true unitig the cid gate alone forbids any
// cross-unitig glue. This builds the adversarial case instead -- THREE
// true unitigs that all share ONE cid and meet at a de Bruijn branch:
//
//   stem ........B]              (B = stem's last k-mer; out-degree 2)
//                 \--[Bsuf+x.... ua
//                  \-[Bsuf+y.... ub
//
// At the branch the three meeting ends share their (k-1) junction
// (B[1:] == Bsuf == first k-1 bases of ua/ub's first k-mer) but NOT their
// full boundary k-mer (B != Bsuf+x != Bsuf+y). The meeting ends are
// CLOSED (a branch is a real unitig terminus), so a correct full-k-mer
// stitcher that honors closed ends keeps all three separate. A (k-1)-keyed
// stitcher -- or one that ignores the cid/closed gate -- would falsely
// glue stem onto ua and/or ub. Each unitig is split into k-overlap
// fragments (randomly RC'd) and appended to `frags_out`; the three whole
// unitigs (closed both ends, shared cid) are appended to `truth_out`.
inline void gen_shared_cid_branch(uint32_t k, uint32_t cid, std::mt19937_64& rng,
                                  std::vector<cdbg::stitchable_unitig>& frags_out,
                                  std::vector<cdbg::stitchable_unitig>& truth_out) {
    using cdbg::stitchable_unitig;
    static const char bases[4] = {'A', 'C', 'G', 'T'};

    // Stem long enough to carry a real interior; its last k bases are B.
    const uint64_t stem_len = 2 * (uint64_t)k + (rng() % (3 * (uint64_t)k));
    std::string stem = random_dna(stem_len, rng);
    const std::string Bsuf = stem.substr(stem.size() - (k - 1));  // B[1:], length k-1

    // Two distinct first bases for the two successor unitigs -> out-degree 2.
    char x = bases[rng() & 3u];
    char y = bases[rng() & 3u];
    while (y == x) y = bases[rng() & 3u];

    // Successor unitig: first k-mer is (Bsuf + base), then a random tail.
    auto make_branch = [&](char base) -> std::string {
        const uint64_t tail = (uint64_t)k + (rng() % (2 * (uint64_t)k));
        return Bsuf + base + random_dna(tail, rng);
    };
    std::string ua = make_branch(x);
    std::string ub = make_branch(y);

    std::string* unitigs[3] = {&stem, &ua, &ub};
    for (uint64_t i = 0; i < 3; ++i) {
        std::string const& S = *unitigs[i];
        stitchable_unitig t;
        t.seq = S;
        t.set_mono(cid, k);
        t.open_flags = 0;
        truth_out.push_back(std::move(t));

        const uint64_t M = (S.size() >= k) ? (S.size() - k + 1) : 1;
        const uint64_t nf = 1 + (rng() % std::max<uint64_t>(1, M - 1));
        auto parts = split_unitig(S, k, cid, nf, rng);
        for (auto& p : parts) frags_out.push_back(std::move(p));
    }
}

}  // namespace cdbg_test
