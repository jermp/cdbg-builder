#pragma once

// Sorting-based unitig extender (GGCAT v2 §3.8, "stabletig" extender), ported to
// cdbg's colored, super-k-mer model.
//
// The legacy phase-2 path (bucket_walker.hpp) builds a per-canonical-k-mer
// hashmap over an entire bucket and walks it node-by-node (up to 8 find()s per
// step). This extender instead processes each MINIMIZER GROUP -- the set of
// super-k-mers sharing one (k-1)-minimizer m-mer, aligned so the m-mer reads
// forward at `min_pos` -- with sorting + centered longest-common-prefix/suffix
// comparisons and O(#super-k-mers) index arrays. No Theta(#distinct k-mers)
// hashmap.
//
// Correctness rationale (why per-group is equivalent to the whole-bucket walk):
// within a cdbg bucket, super-k-mers that must connect into one unitig either
// (a) share the same minimizer m-mer (same group -- a branch's junction (k-1)-
// mer and all its arms carry min(junction), so they co-locate here), or
// (b) meet at a boundary k-mer between two groups, which each group emits as an
// OPEN end that phase-3 stitch rejoins -- exactly how a real cross-bucket
// boundary is handled. So processing each group independently yields the same
// final unitig set; a within-bucket minimizer collision just turns one interior
// k-mer into an open-end + stitch, which is set-equivalent.
//
// This header holds the pure per-group algorithm; the driver that reads a
// bucket, groups by minimizer, aligns records, and routes non-normalizable /
// even-k groups to the hashmap fallback lives in bucket_walker.hpp.

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include "kmer.hpp"
#include "phase2_bucket_process/compact_color_set_dict.hpp"
#include "phase2_bucket_process/stitchable_unitig.hpp"

namespace cdbg {
namespace detail {

// One aligned super-k-mer record within a minimizer group. `bases` live in the
// group's shared 0..3-per-byte buffer, already reverse-complemented so the
// minimizer m-mer reads forward at [min_pos, min_pos+m).
struct se_read {
    uint32_t base_off;  // start of this record's bases in the group buffer
    uint32_t base_len;  // number of bases
    uint32_t min_pos;   // minimizer m-mer offset within [base_off, base_off+base_len)
    uint32_t rsid;      // interned color-set id (into record_sets)
    uint8_t flags;      // SE_* below, in the aligned frame
};

inline constexpr uint8_t SE_INCL_BEGIN = 1u << 0;  // first k-mer is a real run-begin (not a boundary)
inline constexpr uint8_t SE_INCL_END = 1u << 1;    // last k-mer is a real run-end (not a boundary)
inline constexpr uint8_t SE_OWNS_FIRST = 1u << 2;  // this bucket colors the first k-mer
inline constexpr uint8_t SE_OWNS_LAST = 1u << 3;   // this bucket colors the last k-mer

class group_sorting_extender {
public:
    // Process one aligned minimizer group into stitchable fragments appended to
    // `out_local`. `reads` is reordered in place. `record_sets` holds each
    // record's interned color list; stabletig colors are interned into
    // `out_local_dict` (local cids, remapped to global later by process_buckets).
    void process_group(std::vector<se_read>& reads, std::vector<uint8_t> const& gbases, uint32_t k,
                       compact_color_set_dict& record_sets, compact_color_set_dict& out_local_dict,
                       std::vector<stitchable_unitig>& out_local) {
        m_gbases = &gbases;
        m_reads = &reads;
        m_k = k;
        m_record_sets = &record_sets;
        m_local_dict = &out_local_dict;
        m_out = &out_local;
        run();
    }

private:
    // ---- inputs / context ----
    std::vector<uint8_t> const* m_gbases = nullptr;
    std::vector<se_read>* m_reads = nullptr;
    uint32_t m_k = 0;
    compact_color_set_dict* m_record_sets = nullptr;
    compact_color_set_dict* m_local_dict = nullptr;
    std::vector<stitchable_unitig>* m_out = nullptr;

    static constexpr uint32_t NIL = UINT32_MAX;
    static constexpr uint8_t BWD_EXTRA = 1;  // backward (left) open end
    static constexpr uint8_t FWD_EXTRA = 2;  // forward  (right) open end
    static constexpr uint8_t LINKED = 4;

    struct supertig {
        uint32_t base_off;  // stabletig sequence slice in the group buffer
        uint32_t base_len;
        uint64_t cid;
        uint32_t next;
        uint8_t flags;
        uint8_t bwd_owned;  // whether this bucket owns the backward extra-base k-mer
        uint8_t fwd_owned;  // whether this bucket owns the forward extra-base k-mer
        bool backward_extra() const { return flags & BWD_EXTRA; }
        bool forward_extra() const { return flags & FWD_EXTRA; }
        bool is_linked() const { return flags & LINKED; }
    };

    // ---- persistent scratch (reused across groups within a bucket) ----
    std::vector<size_t> forward_skip;
    std::vector<size_t> suffix_sizes;
    std::vector<size_t> elements_mapping;
    std::vector<size_t> prefix_sort_mapping;
    std::vector<size_t> prefix_sort_order;
    std::vector<size_t> lccs_array;
    std::vector<size_t> lccp_array;
    std::vector<supertig> supertigs;
    std::vector<uint32_t> supertigs_mapping;
    std::vector<uint32_t> branching_supertigs;
    std::vector<uint32_t> color_rsids;  // rsids contributing to the current stabletig

    struct sorted_chunk {
        size_t start, end;
        bool is_sorted;
    };
    struct suffix_stack_el {
        size_t index, last_index, sorted_chunks_start;
    };
    struct prefix_stack_el {
        size_t index, last_index, suffix_length;
    };
    struct heap_el {
        size_t index, last;
    };
    std::vector<sorted_chunk> sorted_chunks;
    std::vector<suffix_stack_el> processing_stack;
    std::vector<prefix_stack_el> prefix_stack;
    std::vector<heap_el> heap;
    std::vector<size_t> support_vec;

    // ---- centered comparisons over the aligned group bases ----
    // Longest common centered suffix (bases from min_pos forward); `order` gives
    // the reads' ordering (<0 a before b). On a tie, the shorter suffix is less.
    struct diff {
        size_t lcp;
        int order;
    };
    static int sgn(long v) { return v < 0 ? -1 : (v > 0 ? 1 : 0); }

    diff suffix_diff(se_read const& a, se_read const& b) const {
        uint8_t const* g = m_gbases->data();
        size_t alen = a.base_len - a.min_pos;
        size_t blen = b.base_len - b.min_pos;
        size_t n = std::min(alen, blen);
        uint8_t const* pa = g + a.base_off + a.min_pos;
        uint8_t const* pb = g + b.base_off + b.min_pos;
        for (size_t j = 0; j < n; ++j) {
            if (pa[j] != pb[j]) return {j, pa[j] < pb[j] ? -1 : 1};
        }
        return {n, sgn((long)alen - (long)blen)};
    }
    diff prefix_diff(se_read const& a, se_read const& b) const {
        uint8_t const* g = m_gbases->data();
        size_t alen = a.min_pos;  // bases before the minimizer
        size_t blen = b.min_pos;
        size_t n = std::min(alen, blen);
        uint8_t const* pa = g + a.base_off + a.min_pos;  // read right-to-left from here
        uint8_t const* pb = g + b.base_off + b.min_pos;
        for (size_t j = 0; j < n; ++j) {
            uint8_t ca = pa[-(long)j - 1];
            uint8_t cb = pb[-(long)j - 1];
            if (ca != cb) return {j, ca < cb ? -1 : 1};
        }
        return {n, sgn((long)alen - (long)blen)};
    }

    // ---- color interning ----
    uint64_t cid_for_rsids(std::vector<uint32_t>& rsids) {
        std::sort(rsids.begin(), rsids.end());
        rsids.erase(std::unique(rsids.begin(), rsids.end()), rsids.end());
        if (rsids.size() == 1) return m_local_dict->intern_encoded(*m_record_sets, rsids[0]);
        m_merge_scratch.clear();
        for (uint32_t r : rsids) {
            m_record_sets->at(r, m_at_scratch);
            m_merge_scratch.insert(m_merge_scratch.end(), m_at_scratch.begin(), m_at_scratch.end());
        }
        std::sort(m_merge_scratch.begin(), m_merge_scratch.end());
        m_merge_scratch.erase(std::unique(m_merge_scratch.begin(), m_merge_scratch.end()),
                              m_merge_scratch.end());
        return m_local_dict->intern(m_merge_scratch);
    }
    std::vector<uint32_t> m_at_scratch;
    std::vector<uint32_t> m_merge_scratch;

    // ---- main driver (GGCAT sorting.rs process_reads) ----
    void run() {
        std::vector<se_read>& reads = *m_reads;
        const size_t n = reads.size();
        supertigs.clear();
        if (n == 0) return;

        // 0. Sort by centered suffix (descending) and build the lccs array.
        lccs_array.clear();
        std::sort(reads.begin(), reads.end(), [this](se_read const& a, se_read const& b) {
            return suffix_diff(a, b).order > 0;  // larger suffix first
        });
        for (size_t i = 1; i < n; ++i) lccs_array.push_back(suffix_diff(reads[i - 1], reads[i]).lcp);
        lccs_array.push_back(0);

        // 1. suffix sizes
        suffix_sizes.clear();
        for (auto const& r : reads) suffix_sizes.push_back(r.base_len - r.min_pos);

        // 2. reinit arrays
        elements_mapping.resize(n);
        for (size_t i = 0; i < n; ++i) elements_mapping[i] = i;
        supertigs_mapping.assign(n, NIL);
        forward_skip.assign(n, 1);

        // 3. prefix ordering
        prefix_sort_mapping.resize(n);
        for (size_t i = 0; i < n; ++i) prefix_sort_mapping[i] = i;
        std::sort(prefix_sort_mapping.begin(), prefix_sort_mapping.end(),
                  [this, &reads](size_t a, size_t b) {
                      return prefix_diff(reads[a], reads[b]).order > 0;  // larger prefix first
                  });
        prefix_sort_order.assign(n, 0);
        for (size_t idx = 0; idx < n; ++idx) prefix_sort_order[prefix_sort_mapping[idx]] = idx;

        sorted_chunks.clear();
        processing_stack.clear();
        size_t next_position = 0;

        while (!processing_stack.empty() || next_position < n) {
            if (processing_stack.empty()) {
                processing_stack.push_back({next_position, next_position + 1, sorted_chunks.size()});
                sorted_chunks.push_back({next_position, next_position + 1, false});
            }

            size_t suffix_stack_block_start = processing_stack.back().index;
            size_t sorted_chunks_start = processing_stack.back().sorted_chunks_start;
            const size_t needed_suffix = suffix_sizes[suffix_stack_block_start];

            bool restarted = false;
            while (processing_stack.back().last_index < n &&
                   lccs_array[processing_stack.back().last_index - 1] >= needed_suffix) {
                const size_t last_idx = processing_stack.back().last_index;
                if (suffix_sizes[last_idx] > needed_suffix) {
                    // A deeper (longer-suffix) level: push it and restart.
                    next_position = last_idx;
                    processing_stack.push_back({next_position, next_position + 1, sorted_chunks.size()});
                    sorted_chunks.push_back({next_position, next_position + 1, false});
                    restarted = true;
                    break;
                }

                constexpr size_t MIN_SORTED_SIZE = 64;
                const size_t skip = forward_skip[last_idx];
                const bool next_sorted = skip >= MIN_SORTED_SIZE;
                sorted_chunk& last_chunk = sorted_chunks.back();
                if ((!last_chunk.is_sorted && !next_sorted) ||
                    (last_chunk.end - last_chunk.start + skip) < MIN_SORTED_SIZE) {
                    last_chunk.end = last_idx + skip;
                    last_chunk.is_sorted = false;
                } else {
                    sorted_chunks.push_back({last_idx, last_idx + skip, next_sorted});
                }
                processing_stack.back().last_index = last_idx + skip;
            }
            if (restarted) continue;

            const size_t suffix_stack_block_last = processing_stack.back().last_index;

            // Ensure every chunk is prefix-sorted, then k-way merge them.
            heap.clear();
            for (size_t c = sorted_chunks_start; c < sorted_chunks.size(); ++c) {
                sorted_chunk& chunk = sorted_chunks[c];
                if (!chunk.is_sorted) {
                    std::sort(elements_mapping.begin() + chunk.start,
                              elements_mapping.begin() + chunk.end,
                              [this](size_t a, size_t b) {
                                  return prefix_sort_order[a] < prefix_sort_order[b];
                              });
                    chunk.is_sorted = true;
                }
                heap.push_back({chunk.start, chunk.end});
            }
            if (heap.size() > 1) {
                support_vec.assign(elements_mapping.begin() + suffix_stack_block_start,
                                   elements_mapping.begin() + suffix_stack_block_last);
                auto cursor_key = [this, suffix_stack_block_start](heap_el const& e) {
                    return prefix_sort_order[support_vec[e.index - suffix_stack_block_start]];
                };
                // min-heap on prefix_sort_order
                auto cmp = [&](heap_el const& a, heap_el const& b) {
                    return cursor_key(a) > cursor_key(b);
                };
                std::make_heap(heap.begin(), heap.end(), cmp);
                size_t write_index = suffix_stack_block_start;
                while (!heap.empty()) {
                    std::pop_heap(heap.begin(), heap.end(), cmp);
                    heap_el& top = heap.back();
                    elements_mapping[write_index++] = support_vec[top.index - suffix_stack_block_start];
                    if (top.index + 1 >= top.last) {
                        heap.pop_back();
                    } else {
                        top.index += 1;
                        std::push_heap(heap.begin(), heap.end(), cmp);
                    }
                }
            }

            const size_t max_allowed_suffix_prev =
                processing_stack.size() <= 1
                    ? 0
                    : suffix_sizes[processing_stack[processing_stack.size() - 2].index];
            const size_t max_allowed_suffix_post = lccs_array[suffix_stack_block_last - 1];
            const size_t max_allowed_suffix = std::max(max_allowed_suffix_prev, max_allowed_suffix_post);

            process_reads_block(suffix_stack_block_start, suffix_stack_block_last, needed_suffix,
                                max_allowed_suffix);

            forward_skip[suffix_stack_block_start] = suffix_stack_block_last - suffix_stack_block_start;
            next_position = suffix_stack_block_last;
            if (max_allowed_suffix == max_allowed_suffix_prev) {
                size_t sc_start = processing_stack.back().sorted_chunks_start;
                processing_stack.pop_back();
                sorted_chunks.resize(sc_start);
            }
        }

        emit_chains();
    }

    // GGCAT sorting.rs process_reads_block, stripped of multiplicity / abundance
    // / simplitig / even-k paths (cdbg is colored, cutoff 1, odd-k here).
    void process_reads_block(size_t range_start, size_t range_end, size_t suffix_length,
                             size_t target_suffix_length) {
        std::vector<se_read>& reads = *m_reads;
        const uint32_t k = m_k;

        // lccp over the prefix-sorted elements_mapping[range]
        lccp_array.clear();
        for (size_t idx = range_start + 1; idx < range_end; ++idx) {
            se_read const& a = reads[elements_mapping[idx - 1]];
            se_read const& b = reads[elements_mapping[idx]];
            lccp_array.push_back(prefix_diff(a, b).lcp);
        }
        lccp_array.push_back(0);

        prefix_stack.clear();
        prefix_stack.push_back({range_start, range_end, suffix_length});

        while (!prefix_stack.empty()) {
            prefix_stack_el pe = prefix_stack.back();
            prefix_stack.pop_back();
            size_t element_target_index = pe.index;
            const size_t target_index_end = pe.last_index;
            const size_t shared_suffix = pe.suffix_length;

            size_t supertigs_range_start = supertigs.size();

            while (element_target_index < target_index_end) {
                se_read const& reference = reads[elements_mapping[element_target_index]];
                const size_t reference_index = element_target_index;

                branching_supertigs.push_back(supertigs_mapping[elements_mapping[element_target_index]]);
                element_target_index += 1;

                size_t minimum_prefix_share = reference.min_pos;
                if (minimum_prefix_share + shared_suffix < k) {
                    const size_t next_matching = lccp_array[element_target_index - 1 - range_start];
                    if (next_matching + shared_suffix < (size_t)k - 1 ||
                        element_target_index >= target_index_end) {
                        process_supertigs(supertigs_range_start);
                        branching_supertigs.clear();
                        supertigs_range_start = supertigs.size();
                    }
                    suffix_sizes[reference_index] = 0;
                    continue;
                }

                color_rsids.clear();
                color_rsids.push_back(reference.rsid);

                bool km1mer_break = true;
                while (element_target_index < target_index_end) {
                    se_read const& next_read = reads[elements_mapping[element_target_index]];
                    const size_t left_matching = lccp_array[element_target_index - 1 - range_start];
                    const size_t total_matching = left_matching + shared_suffix;
                    if (total_matching < k) {
                        km1mer_break = total_matching < (size_t)k - 1;
                        break;
                    }
                    branching_supertigs.push_back(
                        supertigs_mapping[elements_mapping[element_target_index]]);
                    minimum_prefix_share = std::min(minimum_prefix_share, left_matching);
                    color_rsids.push_back(next_read.rsid);
                    element_target_index += 1;
                }

                // reference_index..element_target_index all share >= k bases.
                const size_t suffix_limited_suffix = target_suffix_length + 1;
                const size_t prefix_limited_suffix = (size_t)k - minimum_prefix_share;
                const size_t leftmost_allowed_suffix =
                    std::max(prefix_limited_suffix, suffix_limited_suffix);

                bool forward_extra_base =
                    (reference.flags & SE_INCL_END) == 0 &&
                    shared_suffix == (reference.base_len - reference.min_pos);
                bool backward_extra_base = (reference.flags & SE_INCL_BEGIN) == 0 &&
                                           (reference.min_pos + leftmost_allowed_suffix) == k;
                bool fwd_owned = (reference.flags & SE_OWNS_LAST) != 0;
                bool bwd_owned = (reference.flags & SE_OWNS_FIRST) != 0;

                // Double-count guard (GGCAT 502-524): a single-k-mer stabletig
                // with an extra base on both sides is a boundary k-mer shared
                // between two merged buckets; treat it as interior (clear both
                // extra bases) rather than opening it on both ends. cdbg colors
                // it via ownership, so no multiplicity halving is needed.
                if (leftmost_allowed_suffix == shared_suffix &&
                    (forward_extra_base || backward_extra_base)) {
                    for (size_t idx = reference_index; idx < element_target_index; ++idx) {
                        se_read const& rr = reads[elements_mapping[idx]];
                        if ((rr.flags & SE_INCL_END) == 0 &&
                            shared_suffix == (rr.base_len - rr.min_pos)) {
                            forward_extra_base = true;
                            fwd_owned = fwd_owned || (rr.flags & SE_OWNS_LAST) != 0;
                        }
                        if ((rr.flags & SE_INCL_BEGIN) == 0 && rr.min_pos + shared_suffix == k) {
                            backward_extra_base = true;
                            bwd_owned = bwd_owned || (rr.flags & SE_OWNS_FIRST) != 0;
                        }
                    }
                    if (forward_extra_base && backward_extra_base) {
                        forward_extra_base = false;
                        backward_extra_base = false;
                    }
                }

                // Stabletig sequence slice of the reference record.
                const size_t seq_start = reference.min_pos + leftmost_allowed_suffix - k;
                const size_t seq_end = reference.min_pos + shared_suffix;
                if (seq_end - seq_start >= k) {
                    const uint32_t supertig_index = (uint32_t)supertigs.size();
                    supertig st;
                    st.base_off = (uint32_t)(reference.base_off + seq_start);
                    st.base_len = (uint32_t)(seq_end - seq_start);
                    st.cid = cid_for_rsids(color_rsids);
                    st.next = NIL;
                    st.flags = (forward_extra_base ? FWD_EXTRA : 0) |
                               (backward_extra_base ? BWD_EXTRA : 0);
                    st.bwd_owned = bwd_owned ? 1 : 0;
                    st.fwd_owned = fwd_owned ? 1 : 0;
                    supertigs.push_back(st);
                    for (size_t idx = reference_index; idx < element_target_index; ++idx)
                        supertigs_mapping[elements_mapping[idx]] = supertig_index;
                } else {
                    for (size_t idx = reference_index; idx < element_target_index; ++idx)
                        supertigs_mapping[elements_mapping[idx]] = NIL;
                }

                if (km1mer_break) {
                    process_supertigs(supertigs_range_start);
                    branching_supertigs.clear();
                    supertigs_range_start = supertigs.size();
                }

                if (leftmost_allowed_suffix != target_suffix_length + 1) {
                    prefix_stack.push_back(
                        {reference_index, element_target_index, leftmost_allowed_suffix - 1});
                }
            }
        }

        for (size_t idx = range_start; idx < range_end; ++idx) suffix_sizes[idx] = target_suffix_length;
    }

    // GGCAT process_supertigs (unitig mode): link a left and right stabletig only
    // when the junction is unambiguous (1-in / 1-out de Bruijn rule).
    void process_supertigs(size_t supertigs_range_start) {
        std::sort(branching_supertigs.begin(), branching_supertigs.end());
        branching_supertigs.erase(std::unique(branching_supertigs.begin(), branching_supertigs.end()),
                                  branching_supertigs.end());
        if (!branching_supertigs.empty() && branching_supertigs.back() == NIL)
            branching_supertigs.pop_back();

        const size_t left_count = supertigs.size() - supertigs_range_start;
        const bool valid_joining = left_count > 0 && !branching_supertigs.empty() &&
                                   !supertigs[supertigs_range_start].forward_extra() &&
                                   !supertigs[branching_supertigs[0]].backward_extra();
        if (valid_joining) {
            const bool unitig_extendable = branching_supertigs.size() == 1 && left_count == 1;
            if (unitig_extendable) {
                supertigs[branching_supertigs[0]].flags |= LINKED;
                supertigs[supertigs_range_start].next = branching_supertigs[0];
            }
        }
    }

    // GGCAT process_reads join (930-1049): each unlinked chain head assembles its
    // successors (each adds bases_count-(k-1) new bases) into one fragment.
    void emit_chains() {
        const uint32_t k = m_k;
        uint8_t const* g = m_gbases->data();
        for (size_t h = 0; h < supertigs.size(); ++h) {
            if (supertigs[h].is_linked()) continue;

            // Assemble bases and per-stabletig color runs along the chain.
            // Coalesce adjacent same-cid runs (GGCAT push_cid): the downstream
            // monochromatic split emits one unitig PER run, so uncoalesced
            // same-color runs would shatter one tig into many single-run pieces.
            stitchable_unitig u;
            std::vector<color_run>& runs = u.runs;
            auto push_run = [&](uint64_t cid, uint32_t nk) {
                if (!runs.empty() && runs.back().cid == cid)
                    runs.back().num_kmers += nk;
                else
                    runs.push_back({cid, nk});
            };
            uint8_t bwd_owned = supertigs[h].bwd_owned;
            uint8_t fwd_owned = 1;
            uint8_t cum_flags = 0;
            const supertig* st = &supertigs[h];
            bool first = true;
            for (;;) {
                size_t skip = first ? 0 : (k - 1);
                for (size_t i = skip; i < st->base_len; ++i)
                    u.seq.push_back(twobit_to_nuc(g[st->base_off + i]));
                uint32_t nk = st->base_len - k + 1;  // this stabletig's k-mers, all new
                push_run(st->cid, nk);
                cum_flags |= st->flags;
                fwd_owned = st->fwd_owned;
                first = false;
                if (st->next == NIL) break;
                st = &supertigs[st->next];
            }

            const bool open_left = (cum_flags & BWD_EXTRA) != 0;
            const bool open_right = (cum_flags & FWD_EXTRA) != 0;
            if (open_left) u.open_flags |= UNITIG_OPEN_LEFT;
            if (open_right) u.open_flags |= UNITIG_OPEN_RIGHT;

            // Foreign-mark the open-end boundary k-mers this bucket does not own,
            // matching the hashmap walk's primary/foreign split so stitch keeps
            // exactly one colored copy of each boundary k-mer.
            if (open_left && !bwd_owned) mark_first_foreign(runs);
            if (open_right && !fwd_owned) mark_last_foreign(runs);

            m_out->emplace_back(std::move(u));
        }
    }

    static void mark_first_foreign(std::vector<color_run>& runs) {
        if (runs.empty()) return;
        if (runs.front().num_kmers == 1) {
            runs.front().cid = COLOR_RUN_FOREIGN;
        } else {
            runs.front().num_kmers -= 1;
            runs.insert(runs.begin(), color_run{COLOR_RUN_FOREIGN, 1});
        }
    }
    static void mark_last_foreign(std::vector<color_run>& runs) {
        if (runs.empty()) return;
        if (runs.back().num_kmers == 1) {
            runs.back().cid = COLOR_RUN_FOREIGN;
        } else {
            runs.back().num_kmers -= 1;
            runs.push_back(color_run{COLOR_RUN_FOREIGN, 1});
        }
    }
};

}  // namespace detail
}  // namespace cdbg
