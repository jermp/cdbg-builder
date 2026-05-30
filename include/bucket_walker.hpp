#pragma once

// Per-bucket de Bruijn graph construction and unitig walk, with awareness of
// "phantom" edges that cross into other buckets.
//
// Each canonical k-mer in this bucket is annotated with two bits saying
// whether there is a known dBG edge into a k-mer that lives in a *different*
// bucket. The walker treats these phantom edges as having degree 1 on the
// relevant side: if the only edge on a side is phantom, the walk stops there
// and the unitig is emitted with that end marked OPEN, ready to be stitched
// against another bucket's open end. If the side mixes a local edge and a
// phantom edge, the global degree is >= 2 and the walk stops as a real
// branch (CLOSED) — exactly what GGCAT's `try_extend_function` does.
//
// Closed unitigs go straight to the caller. Open unitigs carry their first
// and last canonical k-mers and side-flags so the cross-bucket stitch phase
// can match them up.

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <unordered_dense/unordered_dense.h>

#include "bucket_io.hpp"
#include "compact_color_set_dict.hpp"
#include "kmer.hpp"
#include "streaming_color_set_dict.hpp"
#include "super_kmer.hpp"
#include "util.hpp"

namespace cdgb {

// Per-k-mer color accumulator. We don't store raw colors here. Instead,
// each bucket maintains a `record_sets` compact_color_set_dict that
// interns the (already sorted+deduped) color list of every super-k-mer
// record read from the bucket file (hybrid-encoded, decoded on demand
// into a scratch). A k-mer entry holds only the *ids* of the records
// that contributed to it.
//
// In practice almost every k-mer is contributed by a single record (the
// super-k-mer it lives inside), so `first_rsid` covers the common case
// and `rest` stays empty -- no heap allocation per k-mer. The exceptions
// are (a) hot super-k-mers spilled multiple times during bucket-write,
// which appear as several records sharing the same bases, and (b) k-mers
// at super-k-mer boundaries within the same bucket. Both add only a few
// extra rsids per affected k-mer.
//
// Memory win: where we used to store `records_per_kmer * colors_per_record
// * 4 B` per k-mer (e.g. 18 KB for a popular k-mer in 4546 colors), we
// now store ~4 B per k-mer plus one shared interned copy of each distinct
// record-color list per bucket.
struct kmer_entry {
    static constexpr uint32_t NO_RSID = UINT32_MAX;
    uint32_t first_rsid = NO_RSID;
    std::vector<uint32_t> rest;  // empty in the common case

    void add(uint32_t rsid) {
        if (first_rsid == NO_RSID)
            first_rsid = rsid;
        else
            rest.push_back(rsid);
    }
};

// Open-end marker, in unitig-sequence orientation.
inline constexpr uint8_t UNITIG_OPEN_LEFT = 1u << 0;
inline constexpr uint8_t UNITIG_OPEN_RIGHT = 1u << 1;

struct stitchable_unitig {
    std::string seq;  // ACGT characters
    // Color-class id. While process_bucket is emitting unitigs, this is
    // a *local* cid into that bucket's local_dict; process_buckets then
    // remaps it to a global cid as it merges each bucket's local_dict
    // into the shared global streaming_color_set_dict. By the time stitch_unitigs
    // and the FASTA emitter run, all cids are global.
    uint64_t cid = UINT64_MAX;
    uint8_t open_flags = 0;  // bits from UNITIG_OPEN_*
};

namespace detail {

// Per-canonical-k-mer flags for phantom edges, in canonical orientation.
inline constexpr uint8_t KMER_PHANTOM_LEFT = 1u << 0;   // phantom predecessor exists
inline constexpr uint8_t KMER_PHANTOM_RIGHT = 1u << 1;  // phantom successor exists

struct bucket_kmer_info {
    kmer_entry colors;
    uint8_t phantom = 0;
    // k-overlap "ignored boundary" model (GGCAT hashmap.rs:382-398).
    // owned == false means this canonical k-mer is present ONLY as the
    // foreign first k-mer of a begin-ignored super-k-mer -- it is the
    // duplicated boundary k-mer whose primary copy (color + walk) lives
    // in the adjacent bucket that ends a super on it. We keep it in the
    // map as a colorless NODE so local degree counting sees it (that is
    // what makes branch detection local, replacing the old global
    // joinable map), but we never color it, never assign it a cid, never
    // start a walk from it, and never emit it as an interior node. The
    // walk stops when a step lands on an ignored node, marking that end
    // OPEN and including the full boundary k-mer there.
    bool owned = false;
};

using bucket_kmer_map = ankerl::unordered_dense::map<kmer_int_t, bucket_kmer_info, kmer_hasher>;

// Build the per-canonical-k-mer info from a bucket's super-k-mer stream.
// Each record's color list is interned into `record_sets` once; the per
// k-mer storage is just the rsid (or short list of rsids), not the colors
// themselves. See the kmer_entry comment for the memory rationale.
inline void load_bucket(std::string const& path, uint32_t k, bucket_kmer_map& out,
                        compact_color_set_dict& record_sets) {
    auto& prof = process_prof();
    bucket_reader reader(path);
    uint8_t flags = 0;
    std::vector<uint32_t> colors;
    std::vector<uint8_t> bases;
    while (reader.next(flags, colors, bases)) {
        prof.n_records.fetch_add(1, std::memory_order_relaxed);
        if (bases.size() < k) continue;
        prof.n_kmers.fetch_add((uint64_t)(bases.size() - (k - 1)), std::memory_order_relaxed);

        // Intern the record's already-sorted-deduped color list. Move
        // into the dict on a miss; on a hit it's just a heterogeneous
        // find, no copy. `colors` is repopulated by reader.next() on
        // the next iteration via resize() + assignment, so a moved-from
        // state is safe.
        uint32_t rsid = record_sets.intern(std::move(colors));

        const kmer_int_t mask = kmer_mask(k);
        const uint32_t k_minus_1_x2 = 2 * (k - 1);

        kmer_int_t fwd = 0, rc = 0;
        for (uint32_t i = 0; i < k - 1; ++i) {
            uint8_t v = bases[i];
            fwd = ((fwd << 2) | v) & mask;
            rc = (rc >> 2) | ((kmer_int_t)(v ^ 3) << k_minus_1_x2);
        }

        // k-overlap "ignored boundary" model (GGCAT hashmap.rs:382-398).
        // With k-base super-k-mer overlap, the super's FIRST k-mer is the
        // boundary k-mer duplicated from the predecessor super in the
        // adjacent bucket whenever this super does not begin an ACGT run
        // (IS_ACGT_BEGIN clear). That boundary k-mer's own minimizer is the
        // OLD minimizer, so its primary (colored, walkable) copy lives in
        // that predecessor bucket; here it is FOREIGN. We add it to the map
        // as a colorless ignored node (owned stays false) so the local
        // degree count sees it -- it is the predecessor that makes the owned
        // start k-mer open-left, and the walk includes it at that open end
        // so the emitted fragment carries the full boundary k-mer for the
        // cross-bucket join. We never color it, never give it a cid, never
        // start a walk from it. The super's last k-mer is always owned (only
        // index 0 can be foreign), and the END boundary is unchanged from
        // the pre-overlap model: the foreign successor is NOT present in
        // this bucket, so it is still signalled by a phantom-right edge.
        const bool begin_ignored = (flags & SK_FLAG_IS_ACGT_BEGIN) == 0;
        const bool end_ignored = (flags & SK_FLAG_IS_ACGT_END) == 0;

        kmer_int_t last_can = 0;
        bool last_is_fwd = true;
        bool have_owned = false;
        uint32_t idx = 0;
        for (uint32_t i = k - 1; i < bases.size(); ++i, ++idx) {
            uint8_t v = bases[i];
            fwd = ((fwd << 2) | v) & mask;
            rc = (rc >> 2) | ((kmer_int_t)(v ^ 3) << k_minus_1_x2);
            bool is_fwd = (fwd <= rc);
            kmer_int_t can = is_fwd ? fwd : rc;
            if (idx == 0 && begin_ignored) {
                // Foreign boundary k-mer: ensure the node exists but do not
                // color it or claim ownership. If the same canonical k-mer
                // is owned via another record in this bucket, that record's
                // else-branch below sets owned=true; touching here never
                // clears it.
                out[can];
            } else {
                bucket_kmer_info& info = out[can];
                info.owned = true;
                info.colors.add(rsid);
                last_can = can;
                last_is_fwd = is_fwd;
                have_owned = true;
            }
        }
        if (!have_owned) continue;  // super was a lone foreign boundary k-mer

        // End boundary: the owned last k-mer has a successor in another
        // bucket iff IS_ACGT_END is not set. That successor is not a node
        // here, so signal it with a phantom-right edge (canonical-oriented).
        if (end_ignored) {
            uint8_t bit = last_is_fwd ? KMER_PHANTOM_RIGHT : KMER_PHANTOM_LEFT;
            out[last_can].phantom |= bit;
        }
        // No begin-side phantom: the foreign predecessor is represented by
        // the ignored node added above, which the walk detects and includes.
    }
}

// True iff `can` is present in the map only as an ignored boundary node
// (the foreign duplicated copy of a k-mer owned by an adjacent bucket).
// Such a node participates in local degree counting (so branch detection
// stays local) but is never colored, walk-started, or emitted as interior.
inline bool node_is_ignored(kmer_int_t can, bucket_kmer_map const& m) {
    auto it = m.find(can);
    return it != m.end() && !it->second.owned;
}

// 8-bit local extension mask, bits 0..3 = forward-side successors of `can`,
// bits 4..7 = back-side predecessors (= forward-side successors of rc(can)).
// Counts EVERY neighbor present in the map (owned or ignored).
inline uint8_t local_ext_mask(kmer_int_t can, uint32_t k, bucket_kmer_map const& m) {
    uint8_t out = 0;
    kmer_int_t fwd = can;
    kmer_int_t rev = reverse_complement(can, k);
    for (uint8_t nt = 0; nt < 4; ++nt) {
        kmer_int_t ext_f = shift_append(fwd, nt, k);
        if (m.find(canonical(ext_f, k)) != m.end()) out |= (uint8_t)(1u << nt);
        kmer_int_t ext_r = shift_append(rev, nt, k);
        if (m.find(canonical(ext_r, k)) != m.end()) out |= (uint8_t)(1u << (4 + nt));
    }
    return out;
}

// Like local_ext_mask but only counts OWNED neighbors. The difference
// (full minus owned) is the set of edges to ignored boundary nodes -- the
// foreign cross-bucket crossings. We separate them because an ignored
// node and a phantom edge on the same side are the SAME physical crossing
// (the foreign successor is both flagged via the super's begin/end-ignored
// phantom AND present as the adjacent super's duplicated boundary k-mer);
// counting both as distinct edges falsely inflated degree and closed
// hundreds of thousands of ends that are really simple-path boundaries.
inline uint8_t owned_ext_mask(kmer_int_t can, uint32_t k, bucket_kmer_map const& m) {
    uint8_t out = 0;
    kmer_int_t fwd = can;
    kmer_int_t rev = reverse_complement(can, k);
    for (uint8_t nt = 0; nt < 4; ++nt) {
        kmer_int_t cf = canonical(shift_append(fwd, nt, k), k);
        auto itf = m.find(cf);
        if (itf != m.end() && itf->second.owned) out |= (uint8_t)(1u << nt);
        kmer_int_t cr = canonical(shift_append(rev, nt, k), k);
        auto itr = m.find(cr);
        if (itr != m.end() && itr->second.owned) out |= (uint8_t)(1u << (4 + nt));
    }
    return out;
}

// In walking orientation (rc=false means we walk through the canonical
// forward), forward-side of `m` is its low nibble; flipped under rc.
inline uint8_t fwd_nibble(uint8_t m, bool rc) { return (uint8_t)((rc ? (m >> 4) : m) & 0xf); }
inline uint8_t back_nibble(uint8_t m, bool rc) { return (uint8_t)((rc ? m : (m >> 4)) & 0xf); }

// Phantom on the forward / back side in walking orientation.
inline bool phantom_on_fwd(uint8_t phantom, bool rc) {
    return rc ? (phantom & KMER_PHANTOM_LEFT) : (phantom & KMER_PHANTOM_RIGHT);
}
inline bool phantom_on_back(uint8_t phantom, bool rc) {
    return rc ? (phantom & KMER_PHANTOM_RIGHT) : (phantom & KMER_PHANTOM_LEFT);
}

// Per-side degree of `can` in walking orientation `rc`, separating real
// interior continuations (owned neighbors) from foreign cross-bucket
// crossings. A crossing can be signalled two redundant ways for the same
// physical edge -- a phantom flag on this k-mer AND the foreign successor
// present as an ignored boundary node -- so we DEDUP: foreign counts the
// ignored-node edges, plus the phantom only if there is no ignored edge on
// that side. This is the fix for the ~670k ends the walker over-closed by
// counting phantom and its ignored-node twin as two distinct edges.
struct side_degree {
    uint8_t owned_nib;  // owned successor nibble (for picking the next base)
    int owned;          // # owned successors
    int foreign;        // # foreign crossings (deduped)
    int total;          // owned + foreign
};
inline side_degree fwd_degree(kmer_int_t can, bool rc, uint32_t k, bucket_kmer_map const& m,
                              uint8_t phantom) {
    uint8_t full = local_ext_mask(can, k, m);
    uint8_t owned = owned_ext_mask(can, k, m);
    uint8_t full_n = fwd_nibble(full, rc);
    uint8_t owned_n = fwd_nibble(owned, rc);
    int ignored = __builtin_popcount((uint8_t)(full_n & ~owned_n));
    int ph = phantom_on_fwd(phantom, rc) ? 1 : 0;
    int foreign = ignored + ((ph && ignored == 0) ? 1 : 0);
    side_degree d;
    d.owned_nib = owned_n;
    d.owned = __builtin_popcount(owned_n);
    d.foreign = foreign;
    d.total = d.owned + foreign;
    return d;
}
inline side_degree back_degree(kmer_int_t can, bool rc, uint32_t k, bucket_kmer_map const& m,
                               uint8_t phantom) {
    uint8_t full = local_ext_mask(can, k, m);
    uint8_t owned = owned_ext_mask(can, k, m);
    uint8_t full_n = back_nibble(full, rc);
    uint8_t owned_n = back_nibble(owned, rc);
    int ignored = __builtin_popcount((uint8_t)(full_n & ~owned_n));
    int ph = phantom_on_back(phantom, rc) ? 1 : 0;
    int foreign = ignored + ((ph && ignored == 0) ? 1 : 0);
    side_degree d;
    d.owned_nib = owned_n;
    d.owned = __builtin_popcount(owned_n);
    d.foreign = foreign;
    d.total = d.owned + foreign;
    return d;
}

struct b_step {
    kmer_int_t next_can;
    bool next_rc;
};
inline b_step bstep(kmer_int_t can, bool rc, uint32_t k, uint8_t nt) {
    kmer_int_t cur = rc ? reverse_complement(can, k) : can;
    kmer_int_t next_fwd = shift_append(cur, nt, k);
    kmer_int_t next_rev = reverse_complement(next_fwd, k);
    b_step r;
    if (next_fwd <= next_rev) {
        r.next_can = next_fwd;
        r.next_rc = false;
    } else {
        r.next_can = next_rev;
        r.next_rc = true;
    }
    return r;
}

// Returns true iff (can, rc) is the start of a maximal monochromatic
// unitig (no straight-line predecessor in walking orientation, considering
// both local and phantom edges + color-set agreement). Also reports whether
// the back-side reason for being a left-end was a phantom-only edge — that
// makes the resulting unitig OPEN_LEFT.
struct left_end_check {
    bool is_left_end;
    bool open_left;        // true => unitig will be OPEN_LEFT
    bool prepend;          // true => prepend `prepend_base` to carry the full boundary k-mer
    uint8_t prepend_base;  // walking-orientation leading base of the boundary k-mer
};
inline left_end_check classify_left_end(
    kmer_int_t can, bool rc, uint32_t cid, uint32_t k, bucket_kmer_map const& m,
    ankerl::unordered_dense::map<kmer_int_t, uint32_t, kmer_hasher> const& cid_of) {
    auto it = m.find(can);
    uint8_t phantom = it->second.phantom;
    side_degree bd = back_degree(can, rc, k, m, phantom);
    if (bd.total != 1) return {true, false, false, 0};  // branch / isolated: closed left
    if (bd.owned == 0) {
        // The single back-side edge is a foreign crossing. If it is present
        // as an ignored boundary node, `can` is open-left and the fragment
        // prepends that node's leading base so its first k bases equal the
        // shared boundary k-mer. For a pure phantom crossing, `can` IS the
        // boundary k-mer already -- open-left, no prepend.
        uint8_t full_n = back_nibble(local_ext_mask(can, k, m), rc);
        uint8_t owned_n = back_nibble(owned_ext_mask(can, k, m), rc);
        uint8_t ign_n = (uint8_t)(full_n & ~owned_n);
        if (!ign_n) return {true, true, false, 0};  // pure phantom: no prepend
        uint8_t nt = (uint8_t)__builtin_ctz(ign_n);
        auto sr = bstep(can, !rc, k, nt);
        kmer_int_t pred_can = sr.next_can;
        bool pred_rc = !sr.next_rc;
        kmer_int_t pred_walk = pred_rc ? reverse_complement(pred_can, k) : pred_can;
        uint8_t base = (uint8_t)((pred_walk >> (2 * (k - 1))) & 3);
        return {true, true, true, base};
    }
    // The single back edge is an owned interior predecessor: continue the
    // straight-line check -- its forward in walking orientation must also be
    // degree 1 and color-match, else `can` is a (closed) left end.
    uint8_t nt = (uint8_t)__builtin_ctz(bd.owned_nib);
    auto sr = bstep(can, !rc, k, nt);
    kmer_int_t pred_can = sr.next_can;
    bool pred_rc = !sr.next_rc;
    auto pit = m.find(pred_can);
    if (pit == m.end()) return {true, false, false, 0};  // shouldn't happen
    side_degree pfd = fwd_degree(pred_can, pred_rc, k, m, pit->second.phantom);
    if (pfd.total != 1) return {true, false, false, 0};
    auto cit = cid_of.find(pred_can);
    if (cit == cid_of.end() || cit->second != cid) return {true, false, false, 0};
    return {false, false, false, 0};
}

inline void process_bucket(std::string const& path, uint32_t k, uint32_t num_colors,
                           std::vector<stitchable_unitig>& out_local,
                           compact_color_set_dict& out_local_dict) {
    auto& prof = process_prof();
    bucket_kmer_map kmer_info;
    ankerl::unordered_dense::map<kmer_int_t, uint32_t, kmer_hasher> cid_of;
    {
        // record_sets and the per-k-mer rsid storage are only needed
        // while we're building cid_of. Scoping them here releases the
        // record-set arena and the inner rsid vectors (via the
        // kmer_entry{} reset below) before the walk phase, so the walk
        // sees only kmer_info's phantom bits + cid_of.
        //
        // record_sets is the compact (hybrid-encoded) variant -- on
        // dense pangenome inputs the live-vector dict grew to hundreds
        // of MB per bucket; with N threads in flight that became the
        // dominant peak contributor. compact stores each color list
        // hybrid-encoded and decodes on access into a scratch buffer.
        compact_color_set_dict record_sets(num_colors);
        auto t_load = bucket_process_prof::clock::now();
        load_bucket(path, k, kmer_info, record_sets);
        prof.ns_load.fetch_add(bucket_process_prof::since(t_load), std::memory_order_relaxed);

        auto t_resolve = bucket_process_prof::clock::now();

        // Common case: a k-mer is contributed by a single record, so its
        // eventual color set is identical to that record's color list.
        // Cache rsid -> local cid so the second, third, ... single-rsid
        // k-mers with the same rsid skip the intern entirely. Boundary or
        // spill-duplicate k-mers (rest non-empty) sort+unique their
        // rsids; if a single distinct rsid remains we use the same cache,
        // otherwise we k-way-merge the referenced color lists into a
        // fresh sorted set and intern that.
        std::vector<uint32_t> rsid_to_cid(record_sets.size(), UINT32_MAX);
        cid_of.reserve(kmer_info.size());

        // Scratch reused across all .at() calls so we don't allocate a
        // fresh decoded-colors buffer per access.
        std::vector<uint32_t> at_scratch;

        auto cid_for_single_rsid = [&](uint32_t rsid) -> uint32_t {
            if (rsid_to_cid[rsid] == UINT32_MAX) {
                record_sets.at(rsid, at_scratch);
                rsid_to_cid[rsid] = out_local_dict.intern(std::move(at_scratch));
            }
            return rsid_to_cid[rsid];
        };

        std::vector<uint32_t> rsids_scratch;
        std::vector<uint32_t> merged_scratch;
        for (auto& kv : kmer_info) {
            // Ignored boundary nodes carry no color (no rsid). Skip them:
            // they get no cid and are never walk-started or emitted. Their
            // primary copy is colored in the bucket that owns them.
            if (!kv.second.owned) continue;
            kmer_entry& e = kv.second.colors;
            uint32_t cid;
            if (e.rest.empty()) {
                cid = cid_for_single_rsid(e.first_rsid);
            } else {
                rsids_scratch.clear();
                rsids_scratch.reserve(1 + e.rest.size());
                rsids_scratch.push_back(e.first_rsid);
                rsids_scratch.insert(rsids_scratch.end(), e.rest.begin(), e.rest.end());
                std::sort(rsids_scratch.begin(), rsids_scratch.end());
                rsids_scratch.erase(std::unique(rsids_scratch.begin(), rsids_scratch.end()),
                                    rsids_scratch.end());
                if (rsids_scratch.size() == 1) {
                    cid = cid_for_single_rsid(rsids_scratch[0]);
                } else {
                    merged_scratch.clear();
                    for (uint32_t r : rsids_scratch) {
                        record_sets.at(r, at_scratch);
                        merged_scratch.insert(merged_scratch.end(), at_scratch.begin(),
                                              at_scratch.end());
                    }
                    std::sort(merged_scratch.begin(), merged_scratch.end());
                    merged_scratch.erase(
                        std::unique(merged_scratch.begin(), merged_scratch.end()),
                        merged_scratch.end());
                    cid = out_local_dict.intern(merged_scratch);
                }
            }
            cid_of.emplace(kv.first, cid);
            // Drop per-k-mer rsid storage immediately; phantom bits stay.
            e = kmer_entry{};
        }
        prof.ns_resolve.fetch_add(bucket_process_prof::since(t_resolve),
                                  std::memory_order_relaxed);
    }

    auto t_walk = bucket_process_prof::clock::now();
    ankerl::unordered_dense::map<kmer_int_t, uint8_t, kmer_hasher> visited;
    visited.reserve(kmer_info.size());

    auto extend_and_emit = [&](kmer_int_t start_can, bool start_rc, uint32_t cid, bool open_left,
                               bool prepend, uint8_t prepend_base) {
        stitchable_unitig u;
        kmer_int_t cur = start_rc ? reverse_complement(start_can, k) : start_can;
        u.seq = kmer_to_string(cur, k);
        // Carry the full boundary k-mer at an open-left end: prepend the
        // ignored boundary k-mer's leading base so the fragment's first k
        // bases equal that shared boundary k-mer (full-k-mer stitch key).
        if (prepend) u.seq.insert(u.seq.begin(), twobit_to_nuc(prepend_base));
        // Local cid; process_buckets remaps to global after merging the
        // local_dict into the shared global dict. The actual color list
        // lives in out_local_dict and is moved into the global dict
        // there -- we never copy a color vector into the unitig.
        u.cid = cid;
        if (open_left) u.open_flags |= UNITIG_OPEN_LEFT;

        kmer_int_t can = start_can;
        bool rc = start_rc;
        visited[can] = 1;

        for (;;) {
            auto it = kmer_info.find(can);
            uint8_t phantom = it->second.phantom;
            side_degree fd = fwd_degree(can, rc, k, kmer_info, phantom);
            if (fd.total != 1) break;  // closed (branch / dead end)
            if (fd.owned == 0) {
                // The single forward continuation is a foreign crossing into
                // the adjacent bucket: open-right. Include the boundary
                // k-mer's trailing base when the foreign successor is present
                // as an ignored node here, so the fragment ends with the full
                // shared boundary k-mer; for a pure phantom crossing `can` is
                // itself the boundary k-mer and no extra base is needed.
                if (fd.foreign == 1) {
                    uint8_t full_n = fwd_nibble(local_ext_mask(can, k, kmer_info), rc);
                    uint8_t owned_n = fwd_nibble(owned_ext_mask(can, k, kmer_info), rc);
                    uint8_t ign_n = (uint8_t)(full_n & ~owned_n);
                    if (ign_n) {
                        uint8_t nt = (uint8_t)__builtin_ctz(ign_n);
                        auto sr = bstep(can, rc, k, nt);
                        kmer_int_t nxt_walk = sr.next_rc ? reverse_complement(sr.next_can, k)
                                                         : sr.next_can;
                        u.seq.push_back(twobit_to_nuc((uint8_t)(nxt_walk & 3)));
                    }
                }
                u.open_flags |= UNITIG_OPEN_RIGHT;
                break;
            }
            uint8_t nt = (uint8_t)__builtin_ctz(fd.owned_nib);
            auto sr = bstep(can, rc, k, nt);
            // Successor must (a) have back-side degree exactly 1 in walking
            // orientation (no incoming branch) and (b) carry the same cid.
            auto nit = kmer_info.find(sr.next_can);
            if (nit == kmer_info.end()) break;
            uint8_t next_phantom = nit->second.phantom;
            side_degree nbd = back_degree(sr.next_can, sr.next_rc, k, kmer_info, next_phantom);
            if (nbd.total != 1) break;
            if (cid_of[sr.next_can] != cid) break;
            if (visited.find(sr.next_can) != visited.end()) break;  // cycle closure
            visited[sr.next_can] = 1;
            u.seq.push_back(twobit_to_nuc(nt));
            can = sr.next_can;
            rc = sr.next_rc;
        }
        out_local.emplace_back(std::move(u));
    };

    // First pass: start from k-mers that look like a left-end of some unitig.
    for (auto& kv : kmer_info) {
        kmer_int_t can = kv.first;
        // Ignored boundary nodes are never walk-started: their primary copy
        // is walked in the owning bucket. (They are also never in cid_of.)
        if (!kv.second.owned) continue;
        if (visited.find(can) != visited.end()) continue;
        uint32_t cid = cid_of[can];
        // Try both orientations.
        auto le_f = classify_left_end(can, false, cid, k, kmer_info, cid_of);
        auto le_r = classify_left_end(can, true, cid, k, kmer_info, cid_of);
        left_end_check le;
        bool start_rc;
        if (le_f.is_left_end) {
            start_rc = false;
            le = le_f;
        } else if (le_r.is_left_end) {
            start_rc = true;
            le = le_r;
        } else
            continue;
        extend_and_emit(can, start_rc, cid, le.open_left, le.prepend, le.prepend_base);
    }

    // Second pass: anything still unvisited is on a pure cycle inside the
    // bucket. Pick an arbitrary k-mer to break it; both ends will be CLOSED
    // (since by definition no left-end exists, but we've torn the cycle).
    for (auto& kv : kmer_info) {
        kmer_int_t can = kv.first;
        if (!kv.second.owned) continue;  // never start from an ignored node
        if (visited.find(can) != visited.end()) continue;
        extend_and_emit(can, false, cid_of[can], /*open_left=*/false,
                        /*prepend=*/false, /*prepend_base=*/0);
    }
    prof.ns_walk.fetch_add(bucket_process_prof::since(t_walk), std::memory_order_relaxed);
    prof.n_unitigs.fetch_add((uint64_t)out_local.size(), std::memory_order_relaxed);
    prof.n_local_classes.fetch_add((uint64_t)out_local_dict.size(),
                                   std::memory_order_relaxed);
}

}  // namespace detail

// Process all buckets in parallel; emit stitchable fragments through a
// caller-supplied sink. Each worker takes the next bucket from a
// shared atomic counter, builds the per-bucket walker state, walks the
// chains, merges its local compact_color_set_dict into the shared
// streaming `global_dict` under `global_mu`, remaps each fragment's
// cid from local to global, then feeds each fragment into
// `sink(stitchable_unitig&&)`.
//
// Per-bucket dedup folds what used to be a separate single-threaded
// "intern color sets" pass over every emitted unitig (~6 M on
// salmonella-4546) into the parallel bucket-process phase: the global
// merge runs once per distinct local class (~few thousand per bucket),
// not per fragment. The global dict encodes each new color set into
// its bit_vector builder during intern() and immediately spills
// complete 64-bit words to the final <basename>.color_sets file --
// peak RAM stays proportional to the per-class metadata count, not to
// the sum of class sizes.
//
// `sink` is required to be safe for concurrent calls from
// num_threads worker threads (use an internal mutex if needed). The
// production pipeline passes a frag_unitig_writer that streams
// fragments to disk so the in-RAM accumulator never reaches its
// multi-GB peak.
template <typename Sink>
inline void process_buckets(bucket_writer const& writer, uint32_t k, uint32_t num_colors,
                            uint32_t num_threads, Sink&& sink,
                            streaming_color_set_dict& global_dict, std::mutex& global_mu,
                            std::atomic<uint64_t>* done = nullptr) {
    if (num_threads == 0) num_threads = 1;
    const uint32_t B = writer.num_buckets();
    std::atomic<uint32_t> next{0};
    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    auto run = [&]() {
        for (;;) {
            uint32_t b = next.fetch_add(1);
            if (b >= B) break;
            std::vector<stitchable_unitig> bucket_unitigs;
            // Hybrid-encoded local dict: per-bucket footprint shrinks
            // ~10-30x vs the live-vector dict, so N threads in flight
            // stay within the per-thread share of -g on dense
            // pangenome inputs.
            compact_color_set_dict local_dict(num_colors);
            try {
                detail::process_bucket(writer.bucket_path(b), k, num_colors, bucket_unitigs,
                                       local_dict);
            } catch (std::exception& e) {
                std::cerr << "error processing bucket " << b << ": " << e.what() << '\n';
            }
            // Merge this bucket's local dict into the shared global
            // dict. The previous version did decode + wyhash + fnv1a
            // + dedup-find + (on miss) encode all under global_mu;
            // with 32 threads on dense pangenome inputs that pinned
            // ~80% of bucket-process wall on lock contention.
            //
            // Split the work into two halves: (a) decode the local
            // class + compute its 128-bit content hash, lock-free, in
            // batches of MERGE_BATCH; then (b) take global_mu and
            // call intern_with_hashes() for each pre-hashed entry --
            // just the hashmap probe + on-miss encode + offset write,
            // which is microseconds vs the tens of microseconds the
            // hash compute takes on a multi-thousand-color list.
            //
            // Batch size keeps the lock-free buffer bounded:
            // MERGE_BATCH * worst-case-color-list. 32 colors at
            // num_colors=50K is ~6 MB peak per thread, well under
            // budget.
            constexpr size_t MERGE_BATCH = 32;
            struct prepared_class {
                std::vector<uint32_t> colors;
                streaming_color_set_dict::precomputed_hash h;
            };
            std::vector<prepared_class> batch(MERGE_BATCH);
            std::vector<uint64_t> local_to_global(local_dict.size());

            for (uint32_t lc_start = 0; lc_start < local_dict.size();
                 lc_start += MERGE_BATCH) {
                uint32_t batch_n = (uint32_t)std::min<size_t>(
                    MERGE_BATCH, (size_t)local_dict.size() - lc_start);

                // Lock-free phase: decode + hash.
                auto t_dec = bucket_process_prof::clock::now();
                for (uint32_t i = 0; i < batch_n; ++i) {
                    local_dict.at(lc_start + i, batch[i].colors);
                }
                process_prof().ns_pre_decode.fetch_add(
                    bucket_process_prof::since(t_dec), std::memory_order_relaxed);

                auto t_hash = bucket_process_prof::clock::now();
                for (uint32_t i = 0; i < batch_n; ++i) {
                    batch[i].h = streaming_color_set_dict::compute_hashes(batch[i].colors);
                }
                process_prof().ns_pre_hash.fetch_add(
                    bucket_process_prof::since(t_hash), std::memory_order_relaxed);

                // Locked phase: intern with pre-computed hashes.
                auto t_wait = bucket_process_prof::clock::now();
                std::unique_lock<std::mutex> lk(global_mu);
                process_prof().ns_merge_lock_wait.fetch_add(
                    bucket_process_prof::since(t_wait), std::memory_order_relaxed);
                auto t_merge = bucket_process_prof::clock::now();
                for (uint32_t i = 0; i < batch_n; ++i) {
                    local_to_global[lc_start + i] = global_dict.intern_with_hashes(
                        std::move(batch[i].colors), batch[i].h);
                }
                process_prof().ns_merge.fetch_add(bucket_process_prof::since(t_merge),
                                                  std::memory_order_relaxed);
            }
            for (auto& u : bucket_unitigs) {
                u.cid = local_to_global[u.cid];
                sink(std::move(u));
            }
            process_prof().n_buckets.fetch_add(1, std::memory_order_relaxed);
            if (done) done->fetch_add(1, std::memory_order_relaxed);
        }
    };

    for (uint32_t t = 0; t < num_threads; ++t) workers.emplace_back(run);
    for (auto& w : workers) w.join();
}

}  // namespace cdgb
