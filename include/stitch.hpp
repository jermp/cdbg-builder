#pragma once

// Cross-bucket unitig stitching.
//
// Per-bucket walking emits StitchableUnitig fragments split at
// minimizer-change boundaries. Each fragment marks each side as OPEN
// (the walk stopped because the next k-mer was in another bucket) or
// closed.
//
// We hash open ends by the canonical (k-1)-mer junction: adjacent k-mers
// K1, K2 satisfy K1[1..] == K2[..k-1], and that (k-1)-mer's canonical
// form is invariant under either fragment's orientation choice. Two
// fragments that should be glued therefore land in the same canonical
// class.
//
// Within a class with exactly two open ends, the pair is glueable iff
// the side+orientation combination is consistent. Concretely, with
// is_canonical_fwd = (junction_in_own_frame == canonical(junction)):
//
//   (R, L) and F1 == F2 : glue (both unflipped)
//   (L, R) and F1 == F2 : glue (both unflipped)
//   (R, R) and F1 != F2 : glue (one flipped, one not)
//   (L, L) and F1 != F2 : glue (one flipped, one not)
//
// At chain-walk time, when a unitig is exposed at its merged-RIGHT through
// own side `s`, its flip state is f = (s == LEFT). Follow the adjacency
// to (nxt, nxt_side); nxt's flip state is f_nxt = (nxt_side == RIGHT).
// Junctions are guaranteed to align by the link rules above.
//
// We also require the two unitigs' color sets to be equal — this is the
// only departure from GGCAT's behaviour (which is uncolored).

#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "bucket_walker.hpp"
#include "kmer.hpp"

namespace cdgb {

namespace detail {

inline kmer_int_t encode_kminus1(const char* s, uint32_t n) {
    kmer_int_t x = 0;
    for (uint32_t i = 0; i < n; ++i) { x = (x << 2) | (kmer_int_t)nuc_to_2bit(s[i]); }
    return x;
}

inline std::string revcomp_string(const std::string& s) {
    std::string out(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[s.size() - 1 - i];
        char rc;
        switch (c) {
            case 'A':
                rc = 'T';
                break;
            case 'C':
                rc = 'G';
                break;
            case 'G':
                rc = 'C';
                break;
            case 'T':
                rc = 'A';
                break;
            default:
                rc = 'N';
                break;
        }
        out[i] = rc;
    }
    return out;
}

constexpr uint8_t SIDE_LEFT = 0;
constexpr uint8_t SIDE_RIGHT = 1;

struct EndRef {
    uint32_t unitig_idx;
    uint8_t side;
    bool is_canonical_fwd;
};

inline kmer_int_t side_junction_canonical(const StitchableUnitig& u, uint32_t k, uint8_t side,
                                          bool& is_canonical_fwd) {
    const char* p = (side == SIDE_LEFT) ? u.seq.data() : u.seq.data() + (u.seq.size() - (k - 1));
    kmer_int_t fwd = encode_kminus1(p, k - 1);
    kmer_int_t rc = reverse_complement(fwd, k - 1);
    if (fwd <= rc) {
        is_canonical_fwd = true;
        return fwd;
    }
    is_canonical_fwd = false;
    return rc;
}

struct Link {
    uint32_t other = UINT32_MAX;
    uint8_t other_side = 0;
};

}  // namespace detail

inline void stitch_unitigs(std::vector<StitchableUnitig>& frag, uint32_t k,
                           std::vector<StitchableUnitig>& out) {
    using detail::SIDE_LEFT;
    using detail::SIDE_RIGHT;
    using detail::EndRef;
    using detail::Link;

    if (k < 2) {
        out = std::move(frag);
        return;
    }

    // 1) Index every open end by its canonical (k-1)-mer junction.
    std::unordered_map<kmer_int_t, std::vector<EndRef>, KmerHasher> by_junction;
    by_junction.reserve(frag.size() * 2);
    for (uint32_t i = 0; i < frag.size(); ++i) {
        const auto& u = frag[i];
        if (u.seq.size() < k) continue;
        if (u.open_flags & UNITIG_OPEN_LEFT) {
            bool is_fwd;
            kmer_int_t key = detail::side_junction_canonical(u, k, SIDE_LEFT, is_fwd);
            by_junction[key].push_back({i, SIDE_LEFT, is_fwd});
        }
        if (u.open_flags & UNITIG_OPEN_RIGHT) {
            bool is_fwd;
            kmer_int_t key = detail::side_junction_canonical(u, k, SIDE_RIGHT, is_fwd);
            by_junction[key].push_back({i, SIDE_RIGHT, is_fwd});
        }
    }

    // 2) Build adjacency. adj[u][side] = (other_unitig, other_side) on the
    //    own-frame `side` of u that the link enters.
    std::vector<std::array<Link, 2>> adj(frag.size());

    auto pair_compatible = [](const EndRef& a, const EndRef& b) {
        // (R,L) or (L,R) with same is_canonical_fwd, OR same side with
        // different is_canonical_fwd.
        if (a.side != b.side) return a.is_canonical_fwd == b.is_canonical_fwd;
        return a.is_canonical_fwd != b.is_canonical_fwd;
    };

    for (auto& kv : by_junction) {
        auto& v = kv.second;
        if (v.size() != 2) continue;
        const EndRef& e1 = v[0];
        const EndRef& e2 = v[1];
        if (frag[e1.unitig_idx].colors != frag[e2.unitig_idx].colors) continue;
        if (!pair_compatible(e1, e2)) continue;
        // Both directions of the link.
        adj[e1.unitig_idx][e1.side] = {e2.unitig_idx, e2.side};
        adj[e2.unitig_idx][e2.side] = {e1.unitig_idx, e1.side};
    }

    // 3) Walk chains.
    std::vector<uint8_t> visited(frag.size(), 0);
    out.reserve(frag.size());

    auto take_seq = [&](uint32_t idx, bool flipped) -> std::string {
        return flipped ? detail::revcomp_string(frag[idx].seq) : std::move(frag[idx].seq);
    };

    auto open_at_side = [&](uint32_t idx, uint8_t side_own) -> bool {
        uint8_t bit = (side_own == SIDE_LEFT) ? UNITIG_OPEN_LEFT : UNITIG_OPEN_RIGHT;
        return (frag[idx].open_flags & bit) != 0;
    };

    auto walk_chain = [&](uint32_t start_idx, bool start_flipped) {
        StitchableUnitig merged;
        merged.colors = frag[start_idx].colors;
        merged.seq = take_seq(start_idx, start_flipped);
        visited[start_idx] = 1;

        // The merged-LEFT side of the chain corresponds to start's own
        // SIDE_LEFT (if !start_flipped) or own SIDE_RIGHT (if start_flipped).
        uint8_t left_side_own = start_flipped ? SIDE_RIGHT : SIDE_LEFT;
        if (open_at_side(start_idx, left_side_own)) merged.open_flags |= UNITIG_OPEN_LEFT;

        uint32_t cur = start_idx;
        bool f_cur = start_flipped;

        for (;;) {
            // The own-frame side of cur exposed at merged-RIGHT.
            uint8_t exit_side_own = f_cur ? SIDE_LEFT : SIDE_RIGHT;
            const Link& link = adj[cur][exit_side_own];
            if (link.other == UINT32_MAX) {
                if (open_at_side(cur, exit_side_own)) merged.open_flags |= UNITIG_OPEN_RIGHT;
                break;
            }
            uint32_t nxt = link.other;
            if (visited[nxt]) {
                // Cycle closure: stop here. The cycle's k-mers have all
                // already been emitted via earlier appends; we must NOT
                // drop bases or we'd lose k-mers from the coverage.
                break;
            }
            // f_nxt is determined by which own side of nxt the link enters.
            // Entering through nxt's RIGHT means nxt is flipped (so its own
            // R becomes its merged-LEFT); entering through nxt's LEFT means
            // nxt is unflipped.
            bool f_nxt = (link.other_side == SIDE_RIGHT);
            std::string add = take_seq(nxt, f_nxt);
            merged.seq.append(add.begin() + (k - 1), add.end());
            visited[nxt] = 1;
            cur = nxt;
            f_cur = f_nxt;
        }

        out.emplace_back(std::move(merged));
    };

    // Pass A: start at unitigs with a free own-LEFT side.
    for (uint32_t i = 0; i < frag.size(); ++i) {
        if (visited[i]) continue;
        if (adj[i][SIDE_LEFT].other == UINT32_MAX) { walk_chain(i, /*start_flipped=*/false); }
    }
    // Pass B: start at unitigs with a free own-RIGHT side (and own-LEFT
    // already linked, otherwise pass A would have caught it).
    for (uint32_t i = 0; i < frag.size(); ++i) {
        if (visited[i]) continue;
        if (adj[i][SIDE_RIGHT].other == UINT32_MAX) { walk_chain(i, /*start_flipped=*/true); }
    }
    // Pass C: pure cross-bucket cycles (both sides linked but the chain
    // closes on itself).
    for (uint32_t i = 0; i < frag.size(); ++i) {
        if (visited[i]) continue;
        walk_chain(i, /*start_flipped=*/false);
    }
}

}  // namespace cdgb
