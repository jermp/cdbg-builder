#pragma once

// Cross-bucket unitig stitching.
//
// Per-bucket walking emits stitchable_unitig fragments split at
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
// We additionally require the two unitigs' color sets to be equal before
// gluing — emitted unitigs are monochromatic. GGCAT does the same; the
// difference is *when* each k-mer's color-class id is assigned. GGCAT
// interns color sets to a small integer id inside its `kmers_merge`
// phase, so its per-bucket walk and its cross-bucket `links_compaction`
// (the analog of our stitch) can both gate on cheap integer
// comparisons. We do the same thing slightly later: process_bucket
// itself interns each emitted fragment's color list into the
// per-bucket compact_color_set_dict and stores the returned local
// cid in stitchable_unitig::cid; process_buckets remaps local cids
// to global cids when merging into the shared streaming dict. From
// that point on, the stitch equality check is a uint32_t compare
// and walk_chain inherits the cid by integer assignment instead of
// copying a thousand-entry vector.

#include <array>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <unordered_dense/unordered_dense.h>

#include "bucket_walker.hpp"
#include "kmer.hpp"

namespace cdgb {

namespace detail {

inline kmer_int_t encode_kminus1(char const* s, uint32_t n) {
    kmer_int_t x = 0;
    for (uint32_t i = 0; i < n; ++i) { x = (x << 2) | (kmer_int_t)nuc_to_2bit(s[i]); }
    return x;
}

inline std::string revcomp_string(std::string const& s) {
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

struct end_ref {
    uint64_t unitig_idx;
    uint8_t side;
    bool is_canonical_fwd;
};

// Fixed-size value for `by_junction`. The walker only acts on junctions
// with exactly two ends, so we keep up to 2 inline and reject the rest
// via `count` going to 3+ (overflow). This avoids one std::vector
// allocation per junction, which dominated the upfront indexing time on
// large inputs (millions of junctions -> millions of heap allocs).
struct junction_ends {
    end_ref a;
    end_ref b;
    uint8_t count = 0;  // 0, 1, 2, or 3 (overflow: more than 2 ends)
};

inline kmer_int_t side_junction_canonical(std::string_view seq, uint32_t k, uint8_t side,
                                          bool& is_canonical_fwd) {
    char const* p = (side == SIDE_LEFT) ? seq.data() : seq.data() + (seq.size() - (k - 1));
    kmer_int_t fwd = encode_kminus1(p, k - 1);
    kmer_int_t rc = reverse_complement(fwd, k - 1);
    if (fwd <= rc) {
        is_canonical_fwd = true;
        return fwd;
    }
    is_canonical_fwd = false;
    return rc;
}

// Adjacency is stored as two 64-bit "link" values per fragment (one
// per side) in an adj_entry. Each link is (other_idx << 1) |
// other_side; the all-ones sentinel LINK_NONE = UINT64_MAX means "no
// link on this side". 16 bytes per fragment.
//
// 64-bit indices are mandatory: fragment counts exceed 2^32 on large
// collections (e.g. the Blackwell 661k pangenome emits ~5.7e9
// fragments). The previous layout packed two 32-bit links into one
// uint64_t (8 B/fragment) and silently overflowed past 4.29e9
// fragments. The 8 B saving is not worth the correctness ceiling;
// the external-memory stitch redesign removes the per-fragment
// in-RAM array entirely anyway.
inline constexpr uint64_t LINK_NONE = UINT64_MAX;

inline uint64_t pack_link(uint64_t other, uint8_t other_side) {
    return (other << 1) | (other_side & 1u);
}
inline uint64_t link_other(uint64_t packed) { return packed >> 1; }
inline uint8_t link_side(uint64_t packed) { return (uint8_t)(packed & 1u); }

struct adj_entry {
    uint64_t left = LINK_NONE;
    uint64_t right = LINK_NONE;
};

inline uint64_t adj_get(adj_entry const& e, uint8_t side) {
    return side == SIDE_LEFT ? e.left : e.right;
}
inline void adj_set(adj_entry& e, uint8_t side, uint64_t packed) {
    if (side == SIDE_LEFT)
        e.left = packed;
    else
        e.right = packed;
}

}  // namespace detail

// Streaming variant: emits each finished unitig via the sink callback
// (`sink(stitchable_unitig&&)`) instead of accumulating into an output
// vector. Use this when the caller spills unitigs to disk so the
// stitched seq strings don't all coexist in RAM. The vector-returning
// stitch_unitigs() is a thin adapter around this.
//
// `Source` is a generic frag accessor with the interface:
//   size_t size() const
//   uint64_t cid(size_t i) const
//   uint8_t  open_flags(size_t i) const
//   std::string_view seq_view(size_t i) const
//
// Production builds pass a frag_unitig_reader (mmap-backed, fragments'
// seq bytes never sit in RAM as std::string's). The
// vector_frag_source adapter below wraps a std::vector<stitchable_unitig>
// so unit tests / dev paths still work without disk spill.
template <typename Source, typename Sink>
inline void stitch_unitigs_streaming(Source& frag, uint32_t k, Sink&& sink,
                                     std::atomic<uint64_t>* done = nullptr) {
    using detail::SIDE_LEFT;
    using detail::SIDE_RIGHT;
    using detail::end_ref;
    using detail::junction_ends;
    using detail::adj_entry;
    using detail::LINK_NONE;
    using detail::adj_get;
    using detail::adj_set;
    using detail::pack_link;
    using detail::link_other;
    using detail::link_side;

    const uint64_t n_frags = (uint64_t)frag.size();

    if (k < 2) {
        for (uint64_t i = 0; i < n_frags; ++i) {
            if (done) done->fetch_add(1, std::memory_order_relaxed);
            stitchable_unitig u;
            u.runs = frag.runs(i);
            u.open_flags = frag.open_flags(i);
            std::string_view sv = frag.seq_view(i);
            u.seq.assign(sv.data(), sv.size());
            sink(std::move(u));
        }
        return;
    }

    // 1) Index every open end by its canonical (k-1)-mer junction.
    //
    // The map uses bucket_type::big so its internal value-index is
    // 64-bit -- the default (standard) bucket caps at 2^32 entries,
    // which overflows on large collections (billions of open-end
    // junctions). NOTE: this whole in-RAM map is the structure the
    // external-memory stitch redesign replaces; it does not fit for
    // multi-billion-fragment inputs regardless of the index width.
    std::cerr << "[stitch] indexing " << n_frags << " fragments...\n";
    using junction_map_t =
        ankerl::unordered_dense::map<kmer_int_t, junction_ends, kmer_hasher,
                                     std::equal_to<kmer_int_t>,
                                     std::allocator<std::pair<kmer_int_t, junction_ends>>,
                                     ankerl::unordered_dense::bucket_type::big>;
    junction_map_t by_junction;
    auto add_end = [&](kmer_int_t key, end_ref ref) {
        auto& je = by_junction[key];
        if (je.count == 0)
            je.a = ref;
        else if (je.count == 1)
            je.b = ref;
        // count >= 2 stays as-is; we just bump the counter so the build
        // step below sees the overflow.
        if (je.count < 3) ++je.count;
    };
    for (uint64_t i = 0; i < n_frags; ++i) {
        std::string_view sv = frag.seq_view(i);
        if (sv.size() < k) continue;
        const uint8_t flags = frag.open_flags(i);
        if (flags & UNITIG_OPEN_LEFT) {
            bool is_fwd;
            kmer_int_t key = detail::side_junction_canonical(sv, k, SIDE_LEFT, is_fwd);
            add_end(key, {i, SIDE_LEFT, is_fwd});
        }
        if (flags & UNITIG_OPEN_RIGHT) {
            bool is_fwd;
            kmer_int_t key = detail::side_junction_canonical(sv, k, SIDE_RIGHT, is_fwd);
            add_end(key, {i, SIDE_RIGHT, is_fwd});
        }
    }

    // 2) Build adjacency. adj[u] holds two 64-bit links (one per
    //    side). Each link encodes (other_unitig << 1) | other_side;
    //    LINK_NONE means "no link on that side." 16 bytes / fragment.
    std::cerr << "[stitch] building adjacency over " << by_junction.size() << " junctions...\n";
    std::vector<adj_entry> adj(n_frags);

    auto pair_compatible = [](end_ref const& a, end_ref const& b) {
        // (R,L) or (L,R) with same is_canonical_fwd, OR same side with
        // different is_canonical_fwd.
        if (a.side != b.side) return a.is_canonical_fwd == b.is_canonical_fwd;
        return a.is_canonical_fwd != b.is_canonical_fwd;
    };

    for (auto& kv : by_junction) {
        auto const& je = kv.second;
        if (je.count != 2) continue;  // skip empty / singleton / overflow
        end_ref const& e1 = je.a;
        end_ref const& e2 = je.b;
        // O(1) color-class comparison (legacy monochromatic gate).
        if (frag.mono_cid(e1.unitig_idx) != frag.mono_cid(e2.unitig_idx)) continue;
        if (!pair_compatible(e1, e2)) continue;
        // Both directions of the link.
        adj_set(adj[e1.unitig_idx], e1.side, pack_link(e2.unitig_idx, e2.side));
        adj_set(adj[e2.unitig_idx], e2.side, pack_link(e1.unitig_idx, e1.side));
    }

    // 3) Walk chains.
    by_junction = {};  // free now; we only need adj from here on
    std::cerr << "[stitch] walking chains...\n";
    std::vector<uint8_t> visited(n_frags, 0);

    // Always copy from source: source seq is either a vector entry's
    // string (move would invalidate it for re-access) or a mmap'd
    // region (can't be moved). Per-chain transient cost; the merged
    // seq is later moved into the sink.
    auto take_seq = [&](uint64_t idx, bool flipped) -> std::string {
        std::string_view sv = frag.seq_view(idx);
        std::string s(sv.data(), sv.size());
        if (flipped) s = detail::revcomp_string(s);
        return s;
    };

    auto open_at_side = [&](uint64_t idx, uint8_t side_own) -> bool {
        uint8_t bit = (side_own == SIDE_LEFT) ? UNITIG_OPEN_LEFT : UNITIG_OPEN_RIGHT;
        return (frag.open_flags(idx) & bit) != 0;
    };

    auto walk_chain = [&](uint64_t start_idx, bool start_flipped) {
        stitchable_unitig merged;
        // Inherit runs from the start fragment. (Legacy in-RAM stitcher:
        // the colorless run model is fully carried only by the external
        // stitch; this test-only path keeps the start fragment's runs.)
        merged.runs = frag.runs(start_idx);
        merged.seq = take_seq(start_idx, start_flipped);
        visited[start_idx] = 1;
        if (done) done->fetch_add(1, std::memory_order_relaxed);

        // The merged-LEFT side of the chain corresponds to start's own
        // SIDE_LEFT (if !start_flipped) or own SIDE_RIGHT (if start_flipped).
        uint8_t left_side_own = start_flipped ? SIDE_RIGHT : SIDE_LEFT;
        if (open_at_side(start_idx, left_side_own)) merged.open_flags |= UNITIG_OPEN_LEFT;

        uint64_t cur = start_idx;
        bool f_cur = start_flipped;

        for (;;) {
            // The own-frame side of cur exposed at merged-RIGHT.
            uint8_t exit_side_own = f_cur ? SIDE_LEFT : SIDE_RIGHT;
            uint64_t lnk = adj_get(adj[cur], exit_side_own);
            if (lnk == LINK_NONE) {
                if (open_at_side(cur, exit_side_own)) merged.open_flags |= UNITIG_OPEN_RIGHT;
                break;
            }
            uint64_t nxt = link_other(lnk);
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
            bool f_nxt = (link_side(lnk) == SIDE_RIGHT);
            std::string add = take_seq(nxt, f_nxt);
            merged.seq.append(add.begin() + (k - 1), add.end());
            visited[nxt] = 1;
            if (done) done->fetch_add(1, std::memory_order_relaxed);
            cur = nxt;
            f_cur = f_nxt;
        }

        sink(std::move(merged));
    };

    // Pass A: start at unitigs with a free own-LEFT side.
    for (uint64_t i = 0; i < n_frags; ++i) {
        if (visited[i]) continue;
        if (adj_get(adj[i], SIDE_LEFT) == LINK_NONE) { walk_chain(i, /*start_flipped=*/false); }
    }
    // Pass B: start at unitigs with a free own-RIGHT side (and own-LEFT
    // already linked, otherwise pass A would have caught it).
    for (uint64_t i = 0; i < n_frags; ++i) {
        if (visited[i]) continue;
        if (adj_get(adj[i], SIDE_RIGHT) == LINK_NONE) { walk_chain(i, /*start_flipped=*/true); }
    }
    // Pass C: pure cross-bucket cycles (both sides linked but the chain
    // closes on itself).
    for (uint64_t i = 0; i < n_frags; ++i) {
        if (visited[i]) continue;
        walk_chain(i, /*start_flipped=*/false);
    }
}

// Source adapter wrapping a std::vector<stitchable_unitig> so the
// templated stitch_unitigs_streaming can accept it. Used by the
// vector-returning stitch_unitigs() compat adapter and any test code
// that already has its frags in memory.
struct vector_frag_source {
    std::vector<stitchable_unitig> const& v;
    explicit vector_frag_source(std::vector<stitchable_unitig> const& vec) : v(vec) {}
    size_t size() const { return v.size(); }
    std::vector<color_run> runs(size_t i) const { return v[i].runs; }
    uint64_t mono_cid(size_t i) const { return v[i].mono_cid(); }
    uint8_t open_flags(size_t i) const { return v[i].open_flags; }
    std::string_view seq_view(size_t i) const {
        return std::string_view(v[i].seq.data(), v[i].seq.size());
    }
};

// Backwards-compatible vector-returning adapter around the streaming
// variant. Useful for unit tests; the production builder pipeline now
// passes a frag_unitig_reader directly to stitch_unitigs_streaming.
inline void stitch_unitigs(std::vector<stitchable_unitig>& frag, uint32_t k,
                           std::vector<stitchable_unitig>& out,
                           std::atomic<uint64_t>* done = nullptr) {
    out.clear();
    out.reserve(frag.size());
    auto sink = [&](stitchable_unitig&& u) { out.emplace_back(std::move(u)); };
    vector_frag_source src(frag);
    stitch_unitigs_streaming(src, k, sink, done);
}

}  // namespace cdgb
