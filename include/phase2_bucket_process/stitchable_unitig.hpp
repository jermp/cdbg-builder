#pragma once

// Shared open-ended fragment type produced by phase 2 (bucket-process) and
// consumed by phase 3 (stitch). Both the per-k-mer hashmap walk
// (bucket_walker.hpp) and the sorting-based extender (sorting_extender.hpp)
// emit these, so the type lives in its own header to avoid a walker<->extender
// include cycle.

#include <cstdint>
#include <string>
#include <vector>

namespace cdbg {

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

}  // namespace cdbg
