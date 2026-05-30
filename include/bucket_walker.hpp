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

// A run of consecutive k-mers along a tig that share one color class.
// num_kmers counts k-mers (a length-L tig of L = seq.size()-k+1 k-mers has
// runs summing to L). This is the RLE color sequence GGCAT carries on a
// topological unitig (UnitigColorData); the monochromatic split at emit
// time cuts a tig at its run boundaries.
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
    return flags == KMER_BOUND_LEFT || flags == KMER_BOUND_RIGHT;
}

struct bucket_kmer_info {
    kmer_entry colors;
    uint8_t flags = 0;
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
            if (!begin_incl && idx == 0) contrib |= is_fwd ? KMER_BOUND_LEFT : KMER_BOUND_RIGHT;
            if (!end_incl && idx == last_idx) contrib |= is_fwd ? KMER_BOUND_RIGHT : KMER_BOUND_LEFT;

            bucket_kmer_info& info = out[can];
            info.flags |= contrib;
            info.colors.add(rsid);
        }
    }
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
inline bool walk_step(kmer_int_t can, bool rc, bool forward, uint32_t k,
                      bucket_kmer_map const& m, b_step& out) {
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
            // Every k-mer (including cross-bucket boundary k-mers) is colored
            // from the super(s) it occurs in; the global dict reconciles the
            // two bucket-local copies of each boundary k-mer on its content.
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

    // Append one k-mer's cid to the RLE run sequence, merging with the back
    // run if equal. This is GGCAT's extend_forward (colors/managers/multiple.rs).
    auto push_cid = [](std::vector<color_run>& runs, uint64_t cid) {
        if (!runs.empty() && runs.back().cid == cid)
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

        // Backward. Seed shortcut: if the seed itself is a left boundary
        // (begin/end-ignored on its canonical-left), it is the open-left end
        // already -- do not extend back.
        if (sflags == KMER_BOUND_LEFT) {
            open_left = true;
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
                    open_left = true;
                    break;
                }
            }
        }

        // Forward. Symmetric seed shortcut for a right boundary.
        if (sflags == KMER_BOUND_RIGHT) {
            open_right = true;
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
                    open_right = true;
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
            push_cid(u.runs, cid_of[p.first]);
        };
        for (auto it = bw.rbegin(); it != bw.rend(); ++it) append_kmer(*it);
        append_kmer({seed, false});
        for (auto const& p : fw) append_kmer(p);

        if (open_left) u.open_flags |= UNITIG_OPEN_LEFT;
        if (open_right) u.open_flags |= UNITIG_OPEN_RIGHT;
        out_local.emplace_back(std::move(u));
    };

    // Seed every k-mer once. Unvisited k-mers on a pure cycle are torn at an
    // arbitrary seed; both their ends come out CLOSED (no boundary reached).
    for (auto& kv : kmer_info) {
        if (visited.find(kv.first) != visited.end()) continue;
        emit_from_seed(kv.first);
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
                for (auto& r : u.runs) r.cid = local_to_global[r.cid];
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
