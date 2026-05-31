#pragma once

// Shared cross-bucket stitch helpers.
//
// The cross-bucket stitch itself lives in stitch_extmem.hpp (the
// external-memory, bounded-RAM, GGCAT/BCALM2-faithful implementation used by
// the production pipeline). This header now holds only the small pieces that
// both that stitch and the test harness share:
//
//   - SIDE_LEFT / SIDE_RIGHT      side tags
//   - encode_kminus1              pack a (k-1)-mer from chars
//   - revcomp_string              reverse-complement an ACGT string
//   - side_junction_canonical     canonical (k-1)-mer junction at a side
//   - vector_frag_source          in-memory frag accessor for tests/dev
//
// History: an older in-RAM stitch (stitch_unitigs / stitch_unitigs_streaming)
// lived here. It indexed every open end in one global junction map + a
// per-fragment adjacency array -- O(num_fragments) resident, hundreds of GB on
// the 661k pangenome -- and keyed joins on the (k-1)-mer junction, the
// convention the k-overlap redesign replaced. It was retired once the
// external-memory stitch (full-k-mer keyed, bounded RAM) became the only path
// the builder uses; the test harness now cross-checks ext_mem against ext_file.

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "bucket_walker.hpp"
#include "kmer.hpp"

namespace cdbg {

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

// Canonical (k-1)-mer junction at one side of a fragment. SIDE_LEFT takes the
// first (k-1) bases, SIDE_RIGHT the last (k-1). `is_canonical_fwd` reports
// whether the junction-in-own-frame equals its canonical form. (The
// external-memory stitch keys joins on the full boundary k-mer; this (k-1)
// helper is retained for orientation bookkeeping that both paths reuse.)
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

}  // namespace detail

// Source adapter wrapping a std::vector<stitchable_unitig> so the templated
// stitch entry points (stitch_unitigs_extmem / _file) can accept in-memory
// frags. Used by unit tests and any dev path that already has frags in RAM;
// the production builder passes a frag_unitig_reader instead.
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

}  // namespace cdbg
