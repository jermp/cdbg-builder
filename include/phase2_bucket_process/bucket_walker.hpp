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
#include <condition_variable>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <unordered_dense/unordered_dense.h>

#include "phase1_bucket_write/bucket_io.hpp"
#include "phase2_bucket_process/compact_color_set_dict.hpp"
#include "kmer.hpp"
#include "phase2_bucket_process/streaming_color_set_dict.hpp"
#include "phase1_bucket_write/super_kmer.hpp"
#include "util.hpp"

namespace cdbg {

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
    // During load/resolve this holds the k-mer's record-set id(s) (first_rsid +
    // any rest). AFTER resolve it is REPURPOSED to hold the k-mer's resolved
    // local color-set id (cid): the rsids are dead once resolve has interned
    // them, so the cid lives here instead of in a separate kmer->cid map. The
    // walk reads first_rsid as the cid for primary k-mers.
    uint32_t first_rsid = NO_RSID;
    std::vector<uint32_t> rest;  // empty in the common case; freed after resolve

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

// A run of consecutive k-mers along a tig that share one color class.
// num_kmers counts k-mers (a length-L tig of L = seq.size()-k+1 k-mers has
// runs summing to L). This is the RLE color sequence GGCAT carries on a
// topological unitig (UnitigColorData); the monochromatic split at emit
// time cuts a tig at its run boundaries.
//
// INVARIANT: sum of num_kmers over a tig's runs == its seq k-mer count. Every
// k-mer carries exactly one run unit, INCLUDING a foreign boundary k-mer at an
// open end (whose color is owned by the adjacent bucket): it gets a run with
// cid == COLOR_RUN_FOREIGN as a placeholder so seq and runs stay aligned. The
// stitch reconciles that placeholder to the real color when the open end joins
// its primary-side partner (ext_concat_runs); a fully closed (sunk) tig has had
// every boundary joined, so it contains no foreign placeholders.
inline constexpr uint64_t COLOR_RUN_FOREIGN = UINT64_MAX;
struct color_run {
    uint64_t cid;        // local during process_bucket, global after remap
    uint32_t num_kmers;  // number of k-mers covered by this run (>= 1)
};

struct stitchable_unitig {
    std::string seq;  // ACGT characters
    // Colorless/topological tig: extension follows graph topology only
    // (GGCAT hashmap.rs:230), never breaking on color. The per-k-mer color
    // classes ride along as an RLE run sequence, joined at each stitch merge
    // and cut into monochromatic unitigs at emit. `runs` covers exactly the
    // tig's k-mers in 5'->3' order; while process_bucket emits, the cids are
    // *local* (into the bucket's local_dict) and process_buckets remaps each
    // to a global cid. For a monochromatic tig `runs` has a single element.
    std::vector<color_run> runs;
    uint8_t open_flags = 0;  // bits from UNITIG_OPEN_*

    // Convenience: a freshly walked or post-split tig is monochromatic.
    uint64_t mono_cid() const { return runs.empty() ? UINT64_MAX : runs.front().cid; }
    // Set a single color run covering all of this tig's k-mers (seq must
    // already be assigned). Used by tests that build monochromatic frags.
    void set_mono(uint64_t cid, uint32_t k) {
        uint32_t nk = (seq.size() >= k) ? (uint32_t)(seq.size() - k + 1) : 1;
        runs.assign(1, color_run{cid, nk});
    }
};

namespace detail {

// Per-canonical-k-mer cross-bucket boundary flags (GGCAT hashmap.rs:382-398),
// accumulated (OR) over every super-k-mer occurrence of the k-mer in this
// bucket. bit0 marks a boundary on the k-mer's canonical-LEFT side, bit1 on its
// canonical-RIGHT side. With k-base super overlap a boundary k-mer is the LAST
// k-mer of one super (end-ignored) and the FIRST of the next (begin-ignored);
// those occurrences set the corresponding side bit.
inline constexpr uint8_t KMER_BOUND_LEFT = 1u << 0;
inline constexpr uint8_t KMER_BOUND_RIGHT = 1u << 1;

// A k-mer is a cross-bucket contig boundary iff EXACTLY one side bit is set
// (flags == 1 or == 2): the unitig ends there and continues into the adjacent
// bucket, so that end is emitted OPEN with the full boundary k-mer included for
// the stitch to match. flags == 3 (both sides, e.g. a length-1 super shared on
// both ends) is NOT a break -- the k-mer is walked THROUGH as interior; flags
// == 0 is a plain interior k-mer. (GGCAT hashmap.rs:291-292.)
inline bool is_contig_break(uint8_t flags) {
    return flags == KMER_BOUND_LEFT or flags == KMER_BOUND_RIGHT;
}

struct bucket_kmer_info {
    kmer_entry colors;
    uint8_t flags = 0;
    // GGCAT primary/foreign distinction. A canonical k-mer is PRIMARY in the
    // bucket of its own minimizer, where every occurrence routes that is not a
    // begin-ignored overlap copy; there it is colored (full union over all its
    // occurrences/genomes) and walk-seeded. In an ADJACENT bucket it appears
    // only as the begin-ignored idx-0 overlap copy of the next super -- FOREIGN:
    // a node (for the walk to find the boundary and carry the full k-mer at an
    // open end) but NOT colored and NOT seeded. At a stitch join the primary
    // side's colored copy is kept and the foreign side's copy is dropped, so the
    // boundary k-mer is colored exactly once with its union color.
    bool primary = false;
};

using bucket_kmer_map = ankerl::unordered_dense::map<kmer_int_t, bucket_kmer_info, kmer_hasher>;

// Build the per-canonical-k-mer info from a bucket's super-k-mer stream.
// Each record's color list is interned into `record_sets` once; the per
// k-mer storage is just the rsid (or short list of rsids), not the colors
// themselves. See the kmer_entry comment for the memory rationale.
inline void load_bucket(std::string const& path, uint32_t k, bucket_kmer_map& out,
                        compact_color_set_dict& record_sets) {
    bucket_reader reader(path);
    uint8_t flags = 0;
    std::vector<uint32_t> colors;
    std::vector<uint8_t> bases;
    // Accumulate the profiling counts locally and publish them to the shared
    // process_prof() atomics ONCE at the end of the bucket -- a per-record
    // fetch_add on a global atomic would be 32-thread cache-line contention
    // across hundreds of millions of records.
    uint64_t loaded_records = 0, loaded_kmers = 0;
    while (reader.next(flags, colors, bases)) {
        ++loaded_records;
        if (bases.size() < k) continue;
        loaded_kmers += bases.size() - (k - 1);

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

        // GGCAT per-k-mer boundary flags (hashmap.rs:382-398). A k-mer is
        // begin-ignored if it is the super's FIRST k-mer (idx 0) and the super
        // does not begin an ACGT run (IS_ACGT_BEGIN clear) -- i.e. it is the
        // overlap k-mer duplicated from the predecessor super in the adjacent
        // bucket. Symmetrically end-ignored is the LAST k-mer of a super whose
        // run does not end here. Each contributes a side bit oriented to the
        // k-mer's canonical frame: a forward (canonical) k-mer's begin is its
        // LEFT side and its end is its RIGHT; a reverse-canonical k-mer's are
        // swapped (GGCAT's `<< (!is_forward)` / `<< is_forward`). Every k-mer --
        // boundary or not -- is inserted and colored uniformly; the flag, not a
        // separate node class, marks the cross-bucket boundary, and the same
        // boundary k-mer is present (and colored) in BOTH adjacent buckets, the
        // stitch reconciling the two copies on the full-k-mer key.
        const bool begin_incl = (flags & SK_FLAG_IS_ACGT_BEGIN) != 0;
        const bool end_incl = (flags & SK_FLAG_IS_ACGT_END) != 0;
        const bool owns_first = (flags & SK_FLAG_OWNS_FIRST) != 0;
        const bool owns_last = (flags & SK_FLAG_OWNS_LAST) != 0;
        const uint32_t n_kmers = (uint32_t)(bases.size() - (k - 1));
        const uint32_t last_idx = n_kmers - 1;

        uint32_t idx = 0;
        for (uint32_t i = k - 1; i < bases.size(); ++i, ++idx) {
            uint8_t v = bases[i];
            fwd = ((fwd << 2) | v) & mask;
            rc = (rc >> 2) | ((kmer_int_t)(v ^ 3) << k_minus_1_x2);
            bool is_fwd = (fwd <= rc);
            kmer_int_t can = is_fwd ? fwd : rc;

            uint8_t contrib = 0;
            if (!begin_incl and idx == 0) contrib |= is_fwd ? KMER_BOUND_LEFT : KMER_BOUND_RIGHT;
            if (!end_incl and idx == last_idx)
                contrib |= is_fwd ? KMER_BOUND_RIGHT : KMER_BOUND_LEFT;

            bucket_kmer_info& info = out[can];
            info.flags |= contrib;
            // Color (mark primary) this occurrence iff THIS bucket owns the
            // k-mer (BCALM2: the boundary k-mer is colored in bucket(min(lmin,
            // rmin)) only). Interior k-mers (neither first nor last) are always
            // owned -- they are not cross-bucket boundaries. The first/last
            // k-mers are owned per the SK_FLAG_OWNS_* bits set at ingest. This
            // makes each k-mer primary in EXACTLY one bucket, so its color is the
            // full union over its occurrences there (no partial-color split).
            const bool is_first = (idx == 0);
            const bool is_last = (idx == last_idx);
            bool owned = true;
            if (is_first and !owns_first) owned = false;
            if (is_last and !owns_last) owned = false;
            if (owned) {
                info.colors.add(rsid);
                info.primary = true;
            }
        }
    }
    auto& prof = process_prof();
    prof.n_records.fetch_add(loaded_records, std::memory_order_relaxed);
    prof.n_kmers.fetch_add(loaded_kmers, std::memory_order_relaxed);
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

// In walking orientation (rc=false means we walk through the canonical
// forward), forward-side of `m` is its low nibble; flipped under rc.
inline uint8_t fwd_nibble(uint8_t m, bool rc) { return (uint8_t)((rc ? (m >> 4) : m) & 0xf); }
inline uint8_t back_nibble(uint8_t m, bool rc) { return (uint8_t)((rc ? m : (m >> 4)) & 0xf); }

// Local degree of `can` (walking orientation `rc`) on the forward / back side.
// Counts neighbours present in THIS bucket; a cross-bucket boundary edge is not
// counted (its partner lives in the adjacent bucket). Used to distinguish a
// genuine cross-bucket simple-path boundary (this-side local degree <= 1, so the
// only continuation crosses buckets -> emit OPEN) from a k-mer that is BOTH a
// cross-bucket boundary AND a local branch (this-side local degree >= 2 -> a real
// branch; must be CLOSED, not stitched, or the stitch joins past the branch).
inline int local_fwd_degree(kmer_int_t can, bool rc, uint32_t k, bucket_kmer_map const& m) {
    return __builtin_popcount(fwd_nibble(local_ext_mask(can, k, m), rc));
}
inline int local_back_degree(kmer_int_t can, bool rc, uint32_t k, bucket_kmer_map const& m) {
    return __builtin_popcount(back_nibble(local_ext_mask(can, k, m), rc));
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

// One clean unitig extension step from (can, rc) in walking orientation.
// `forward` picks the forward side (else the back side). A step is clean iff
// this side has exactly ONE present neighbor (GGCAT count==1) AND that neighbor
// has exactly one present neighbor facing back toward us (ocount==1, no incoming
// branch). Returns false (no extension) otherwise. Boundary flags are NOT
// consulted here -- the caller stops at a contig-break k-mer after stepping.
inline bool walk_step(kmer_int_t can, bool rc, bool forward, uint32_t k, bucket_kmer_map const& m,
                      b_step& out) {
    uint8_t mask = local_ext_mask(can, k, m);
    uint8_t nib = forward ? fwd_nibble(mask, rc) : back_nibble(mask, rc);
    if (__builtin_popcount(nib) != 1) return false;
    uint8_t nt = (uint8_t)__builtin_ctz(nib);
    b_step s;
    if (forward) {
        s = bstep(can, rc, k, nt);
    } else {
        s = bstep(can, !rc, k, nt);
        s.next_rc = !s.next_rc;
    }
    // The neighbor's count on the side facing back toward `can` must be 1.
    uint8_t nmask = local_ext_mask(s.next_can, k, m);
    uint8_t facing = forward ? back_nibble(nmask, s.next_rc) : fwd_nibble(nmask, s.next_rc);
    if (__builtin_popcount(facing) != 1) return false;
    out = s;
    return true;
}

inline void process_bucket(std::string const& path, uint32_t k, uint64_t num_colors,
                           std::vector<stitchable_unitig>& out_local,
                           compact_color_set_dict& out_local_dict) {
    auto& prof = process_prof();
    bucket_kmer_map kmer_info;
    {
        // record_sets and the per-k-mer rsid storage are only needed while we
        // resolve each k-mer's cid. We then stash the cid back INTO kmer_info
        // (reusing kmer_entry.first_rsid, dead after resolve) and free the rest
        // vectors, so the walk needs no separate kmer->cid map -- it reads the
        // cid straight from the kmer_info entry it already looks up.
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

        // Scratch reused across all .at() calls so we don't allocate a
        // fresh decoded-colors buffer per access.
        std::vector<uint32_t> at_scratch;

        auto cid_for_single_rsid = [&](uint32_t rsid) -> uint32_t {
            if (rsid_to_cid[rsid] == UINT32_MAX) {
                // The class IS record_sets[rsid] -- copy its encoded bits and
                // reuse its hashes instead of decode -> re-hash -> re-encode.
                rsid_to_cid[rsid] = out_local_dict.intern_encoded(record_sets, rsid);
            }
            return rsid_to_cid[rsid];
        };

        std::vector<uint32_t> rsids_scratch;
        std::vector<uint32_t> merged_scratch;
        for (auto& kv : kmer_info) {
            // Only PRIMARY k-mers carry color (an rsid); foreign overlap copies
            // have no rsid and get no cid -- their colored copy lives in the
            // adjacent bucket. Skipping them also avoids indexing rsid_to_cid
            // with NO_RSID. The walk reads a cid only for primary k-mers, so a
            // foreign entry's first_rsid stays NO_RSID and is never used as a cid.
            if (!kv.second.primary) continue;
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
                    merged_scratch.erase(std::unique(merged_scratch.begin(), merged_scratch.end()),
                                         merged_scratch.end());
                    cid = out_local_dict.intern(merged_scratch);
                }
            }
            // Stash the resolved cid back into the entry (first_rsid is dead now)
            // and free the rsid overflow -- the walk reads the cid straight from
            // here, so there is no separate kmer->cid map. The phantom/flags bits
            // on bucket_kmer_info stay.
            std::vector<uint32_t>().swap(e.rest);
            e.first_rsid = cid;
        }
        prof.ns_resolve.fetch_add(bucket_process_prof::since(t_resolve), std::memory_order_relaxed);
    }

    auto t_walk = bucket_process_prof::clock::now();
    ankerl::unordered_dense::map<kmer_int_t, uint8_t, kmer_hasher> visited;
    visited.reserve(kmer_info.size());

    // Append one k-mer's cid to the RLE run sequence, merging with the back
    // run if equal. This is GGCAT's extend_forward (colors/managers/multiple.rs).
    auto push_cid = [](std::vector<color_run>& runs, uint64_t cid) {
        if (!runs.empty() and runs.back().cid == cid)
            runs.back().num_kmers += 1;
        else
            runs.push_back({cid, 1});
    };

    // Walk one maximal unitig containing `seed` (treated in its canonical
    // forward orientation), extending backward then forward (GGCAT
    // hashmap.rs:455-545). Degree is uniform over all present neighbors; the
    // boundary k-mer is INCLUDED at an open end and the walk stops there. A
    // side is OPEN iff it ended on a cross-bucket contig break (the seed's own
    // boundary flag, or a contig-break k-mer reached by extension); a side that
    // ended on a branch / dead-end / cycle is CLOSED.
    auto emit_from_seed = [&](kmer_int_t seed) {
        visited[seed] = 1;
        const uint8_t sflags = kmer_info.find(seed)->second.flags;

        // Step lists in walking orientation. `bw` is nearest-first (it is later
        // reversed so the unitig reads left->right); `fw` is in order.
        std::vector<std::pair<kmer_int_t, bool>> bw, fw;
        bool open_left = false, open_right = false;

        // A boundary k-mer should be emitted OPEN (for cross-bucket stitching)
        // ONLY when its local degree on the boundary side is <= 1: then the
        // single continuation truly crosses into the adjacent bucket. If the
        // boundary side also has >= 2 local neighbours, the k-mer is a genuine
        // branch and must be CLOSED -- otherwise the stitch joins one local
        // branch arm past the branch point (the strict-topology internal-branch
        // failures on tandem-repeat k-mers).
        // Backward. Seed shortcut: if the seed itself is a left boundary it is
        // the open-left end already -- do not extend back (subject to the
        // degree guard above).
        if (sflags == KMER_BOUND_LEFT) {
            open_left = (local_back_degree(seed, false, k, kmer_info) <= 1);
        } else {
            kmer_int_t can = seed;
            bool rc = false;
            for (;;) {
                b_step s;
                if (!walk_step(can, rc, /*forward=*/false, k, kmer_info, s)) break;  // closed
                if (visited.find(s.next_can) != visited.end()) break;                // cycle
                visited[s.next_can] = 1;
                bw.push_back({s.next_can, s.next_rc});
                can = s.next_can;
                rc = s.next_rc;
                if (is_contig_break(kmer_info.find(can)->second.flags)) {
                    open_left = (local_back_degree(can, rc, k, kmer_info) <= 1);
                    break;
                }
            }
        }

        // Forward. Symmetric seed shortcut for a right boundary.
        if (sflags == KMER_BOUND_RIGHT) {
            open_right = (local_fwd_degree(seed, false, k, kmer_info) <= 1);
        } else {
            kmer_int_t can = seed;
            bool rc = false;
            for (;;) {
                b_step s;
                if (!walk_step(can, rc, /*forward=*/true, k, kmer_info, s)) break;
                if (visited.find(s.next_can) != visited.end()) break;
                visited[s.next_can] = 1;
                fw.push_back({s.next_can, s.next_rc});
                can = s.next_can;
                rc = s.next_rc;
                if (is_contig_break(kmer_info.find(can)->second.flags)) {
                    open_right = (local_fwd_degree(can, rc, k, kmer_info) <= 1);
                    break;
                }
            }
        }

        // Assemble left->right: reverse(bw), seed, fw. Consecutive entries are
        // forward dBG edges (overlap k-1), so the sequence appends one base per
        // step and the cid RLE pushes one run element per k-mer.
        stitchable_unitig u;
        auto walk_kmer = [&](std::pair<kmer_int_t, bool> const& p) -> kmer_int_t {
            return p.second ? reverse_complement(p.first, k) : p.first;
        };
        bool first = true;
        auto append_kmer = [&](std::pair<kmer_int_t, bool> const& p) {
            if (first) {
                u.seq = kmer_to_string(walk_kmer(p), k);
                first = false;
            } else {
                u.seq.push_back(twobit_to_nuc((uint8_t)(walk_kmer(p) & 3)));
            }
            // Every k-mer gets exactly one run unit so sum(runs)==seq k-mers.
            // A foreign boundary k-mer (open end, colored by the adjacent
            // bucket) gets a COLOR_RUN_FOREIGN placeholder; the stitch replaces
            // it with the real color from the primary-side partner at the join.
            // One lookup: .primary and the resolved cid both live on the entry
            // (cid was stashed into colors.first_rsid by the resolve pass).
            auto const& info = kmer_info.find(p.first)->second;
            push_cid(u.runs, info.primary ? info.colors.first_rsid : COLOR_RUN_FOREIGN);
        };
        for (auto it = bw.rbegin(); it != bw.rend(); ++it) append_kmer(*it);
        append_kmer({seed, false});
        for (auto const& p : fw) append_kmer(p);

        if (open_left) u.open_flags |= UNITIG_OPEN_LEFT;
        if (open_right) u.open_flags |= UNITIG_OPEN_RIGHT;
        out_local.emplace_back(std::move(u));
    };

    // Seed every k-mer present in this bucket once -- NOT just the color-owned
    // (primary) ones. Every super-k-mer record in this bucket lives here (it
    // was written to bucket_of(its own minimizer)), so all of its k-mers must
    // be walked here; ownership gates COLORING only (an unowned/foreign k-mer
    // gets a COLOR_RUN_FOREIGN run unit at emit, reconciled by the stitch).
    //
    // Gating the seed on `primary` strands degenerate supers that own none of
    // their k-mers: a short super [F,X] sandwiched between two smaller-
    // minimizer neighbors has owns_first==owns_last==0 (both k-mers colored in
    // the adjacent buckets), so neither is primary. It is exactly the bridge
    // fragment whose open ends are the stitch partners of F (in F's owner
    // bucket) and X (in X's owner bucket); never walking it leaves both F and X
    // with a lonely, partnerless open end -- the cross-bucket joins then fail
    // and unitigs come out fragmented (salmonella-10: 440k unitigs instead of
    // 87,297). Boundary k-mers are shared between two buckets by construction;
    // each is seeded in both and emitted as the open end of each side's
    // fragment, which is precisely what the full-k-mer stitch joins on. The
    // `visited` set still prevents re-walking within a bucket, so each bucket
    // emits each of its fragments exactly once.
    for (auto& kv : kmer_info) {
        if (visited.find(kv.first) != visited.end()) continue;
        emit_from_seed(kv.first);
    }
    prof.ns_walk.fetch_add(bucket_process_prof::since(t_walk), std::memory_order_relaxed);
    prof.n_unitigs.fetch_add((uint64_t)out_local.size(), std::memory_order_relaxed);
    prof.n_local_classes.fetch_add((uint64_t)out_local_dict.size(), std::memory_order_relaxed);
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
// `mem_budget_bytes` (0 = unbounded) is the TOTAL-RSS target for the phase (a
// fraction of -g). A worker waits to "admit" its bucket until
//   carry + reserved_kmer_info + live_color_dict + this_bucket_cost <= target
// where carry is the RSS at phase start and the color dict is re-read live, so
// the (monotonically growing) dict is reserved for dynamically rather than
// left to overflow on top of a static kmer_info budget. All num_threads
// threads stay alive; under a tight budget FEWER (large) buckets load at once,
// MORE (small) buckets do -- bounding RAM without ever reducing the user's
// thread count. One bucket is always allowed even if it alone exceeds the
// budget (forward progress).
template <typename Sink>
inline void process_buckets(bucket_writer const& writer, uint32_t k, uint64_t num_colors,
                            uint32_t num_threads, Sink&& sink,
                            streaming_color_set_dict& global_dict, std::mutex& global_mu,
                            std::atomic<uint64_t>* done = nullptr, uint64_t mem_budget_bytes = 0) {
    if (num_threads == 0) num_threads = 1;
    // global_dict is now internally sharded/thread-safe, so the caller's
    // global_mu is no longer used to guard the merge. Kept in the signature
    // for API stability (and a possible future use); silence the warning.
    (void)global_mu;
    const uint32_t B = writer.num_buckets();
    std::atomic<uint32_t> next{0};
    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    // Per-bucket resident WORKING SET as a multiple of the bucket's uncompressed
    // on-disk bytes. This must cover everything an in-flight bucket holds at its
    // walk-phase peak, not just kmer_info: kmer_info (cid now stored inline) +
    // visited + bucket_unitigs (emitted fragments) + local_dict. Measured on 100K -g16:
    // peak/reserved was 1.20x at small buckets and 1.37x at large ones with a
    // 16x kmer_info-only estimate -> the full working set is ~24x. Under-
    // reserving here over-admits and busts -g; the live-RSS ceiling below is the
    // backstop for whatever this still misses (fragment volume is graph-shaped).
    constexpr double WORKING_SET_OVERHEAD = 24.0;
    std::mutex mem_mu;
    std::condition_variable mem_cv;
    uint64_t mem_in_use = 0;
    uint32_t mem_n_admitted = 0;  // how many buckets currently admitted
    // mem_budget_bytes is the TOTAL-RSS target for the phase (a fraction of
    // -g), NOT a pre-shrunk kmer_info budget. mem_carry = everything resident
    // at phase start (frag sink, residual heap); the color dict starts ~empty
    // here and is reserved for separately and LIVE in the predicate below, so
    // the admission keeps  carry + reserved_working_set + live_dict <= target.
    const uint64_t mem_carry = (mem_budget_bytes > 0) ? current_rss_bytes() : 0;

    auto run = [&]() {
        for (;;) {
            uint32_t b = next.fetch_add(1);
            if (b >= B) break;
            const uint64_t cost = WORKING_SET_OVERHEAD * writer.bucket_unc_bytes(b);
            // Admission: wait until this bucket fits the budget, or it would be
            // the only one running (so an oversized bucket never deadlocks).
            if (mem_budget_bytes > 0) {
                std::unique_lock<std::mutex> lk(mem_mu);
                mem_cv.wait(lk, [&] {
                    if (mem_n_admitted == 0) return true;  // forward progress
                    // Hard live-RSS ceiling: never admit if the actual resident
                    // set plus this bucket's estimated working set would exceed
                    // the budget. current_rss() already includes carry, the live
                    // dict, every loaded bucket's FULL footprint, and glibc
                    // fragmentation -- so this catches whatever WORKING_SET_
                    // OVERHEAD under-reserves, workload-independently.
                    if (current_rss_bytes() + cost > mem_budget_bytes) return false;
                    // Reservation (lag-safe): bound concurrent admits before
                    // their loads register in RSS. Reserve dynamically for the
                    // live color dict (read lock-free; m_classes.size() is a
                    // monotone gauge, a benign approximate race).
                    const uint64_t dyn = mem_carry + global_dict.resident_bytes();
                    const uint64_t avail = mem_budget_bytes > dyn ? mem_budget_bytes - dyn : 0;
                    return mem_in_use + cost <= avail;
                });
                mem_in_use += cost;
                ++mem_n_admitted;
            }
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
            // Release the admitted memory once kmer_info/local_dict for this
            // bucket are gone (they go out of scope at the end of this loop
            // iteration; release after the merge below, before next bucket).
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

            for (uint32_t lc_start = 0; lc_start < local_dict.size(); lc_start += MERGE_BATCH) {
                uint32_t batch_n =
                    (uint32_t)std::min<size_t>(MERGE_BATCH, (size_t)local_dict.size() - lc_start);

                // Lock-free phase: decode + hash.
                auto t_dec = bucket_process_prof::clock::now();
                for (uint32_t i = 0; i < batch_n; ++i) {
                    local_dict.at(lc_start + i, batch[i].colors);
                }
                process_prof().ns_pre_decode.fetch_add(bucket_process_prof::since(t_dec),
                                                       std::memory_order_relaxed);

                auto t_hash = bucket_process_prof::clock::now();
                for (uint32_t i = 0; i < batch_n; ++i) {
                    batch[i].h = streaming_color_set_dict::compute_hashes(batch[i].colors);
                }
                process_prof().ns_pre_hash.fetch_add(bucket_process_prof::since(t_hash),
                                                     std::memory_order_relaxed);

                // Intern phase: intern_with_hashes is internally sharded and
                // thread-safe -- per-shard dedup locks for the common
                // duplicate case, and for a genuinely new class the expensive
                // hybrid-encode runs lock-free (only the append to the single
                // output stream is briefly serialized). No external global
                // lock, so threads no longer serialize the whole merge on one
                // mutex (this was ~14s of merge_wait / phase on bw20k).
                auto t_merge = bucket_process_prof::clock::now();
                for (uint32_t i = 0; i < batch_n; ++i) {
                    local_to_global[lc_start + i] =
                        global_dict.intern_with_hashes(std::move(batch[i].colors), batch[i].h);
                }
                process_prof().ns_merge.fetch_add(bucket_process_prof::since(t_merge),
                                                  std::memory_order_relaxed);
            }
            for (auto& u : bucket_unitigs)
                for (auto& r : u.runs)
                    if (r.cid != COLOR_RUN_FOREIGN) r.cid = local_to_global[r.cid];
            // One lock acquisition for the whole bucket (see write_batch): the
            // per-fragment sink call was the bucket-process lock-contention /
            // voluntary-context-switch hot spot.
            sink.write_batch(bucket_unitigs);
            // Release this bucket's memory admission: kmer_info + local_dict +
            // bucket_unitigs for bucket b are done (kmer_info/local_dict freed
            // inside process_bucket and at scope end; bucket_unitigs just
            // drained to the sink). Wake any worker waiting to admit.
            if (mem_budget_bytes > 0) {
                std::lock_guard<std::mutex> lk(mem_mu);
                mem_in_use -= cost;
                --mem_n_admitted;
                mem_cv.notify_all();
            }
            process_prof().n_buckets.fetch_add(1, std::memory_order_relaxed);
            if (done) done->fetch_add(1, std::memory_order_relaxed);
        }
    };

    for (uint32_t t = 0; t < num_threads; ++t) workers.emplace_back(run);
    for (auto& w : workers) w.join();
}

}  // namespace cdbg
