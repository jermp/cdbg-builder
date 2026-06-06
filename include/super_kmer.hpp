#pragma once

// Binary record format for super-k-mers stored in a bucket file:
//
//   [color_count : varint]
//   [color_delta_0 : varint] ... [color_delta_{n-1} : varint]
//   [length_with_flags : varint]
//   [2-bit packed bases]
//
// Records are *compacted*: a super-k-mer that occurs in multiple input files
// is stored once with the union of its colors. Colors are emitted in
// ascending order as deltas from the previous color (delta_0 is absolute).
//
// length_with_flags is `(length_in_bases << 2) | flags`, where flags is:
//   bit 0 (IS_ACGT_BEGIN) = this super-k-mer's first k-mer (in input order) is
//                           the actual first k-mer of an ACGT-only run, i.e. it
//                           has no predecessor in the input.
//   bit 1 (IS_ACGT_END)   = analogously, this super-k-mer's last k-mer is the
//                           last k-mer of its ACGT run.
//
// When the same super-k-mer is contributed by multiple files, the compactor
// AND-merges flag bits across contributors. So IS_ACGT_BEGIN is set in the
// stored record iff *every* contributing input had this super-k-mer at the
// start of an ACGT-only run; otherwise at least one contributor had a
// predecessor in another bucket and the bucket walker correctly inserts a
// phantom edge on that side. Same for IS_ACGT_END.
//
// Bases are packed 4 per byte, lowest 2 bits = first base of the run. Records
// are concatenated in the file with no padding between them.

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

namespace cdbg {

// ---- 2-bit packing helpers ---------------------------------------------------

// Append `n` 2-bit-encoded bases from `src` (each entry in [0,3]) to `dst`.
// Bases are packed into bytes with base i occupying bits (2*(i%4)..2*(i%4)+1)
// of byte (i/4).
inline void pack_2bit(uint8_t const* src, size_t n, std::vector<uint8_t>& dst) {
    size_t bytes_needed = (n + 3) / 4;
    size_t off = dst.size();
    dst.resize(off + bytes_needed, 0);
    for (size_t i = 0; i < n; ++i) {
        dst[off + (i >> 2)] |= (uint8_t)((src[i] & 3u) << (2 * (i & 3)));
    }
}

// Unpack `n` bases from `src` (starting at base index 0 in src's first byte)
// into `dst` (cleared and resized to n).
inline void unpack_2bit(uint8_t const* src, size_t n, std::vector<uint8_t>& dst) {
    dst.assign(n, 0);
    for (size_t i = 0; i < n; ++i) { dst[i] = (src[i >> 2] >> (2 * (i & 3))) & 3u; }
}

// ---- LEB128 varint -----------------------------------------------------------

inline void varint_write(uint64_t v, std::vector<uint8_t>& out) {
    while (v >= 0x80) {
        out.push_back((uint8_t)((v & 0x7f) | 0x80));
        v >>= 7;
    }
    out.push_back((uint8_t)v);
}

// Read a varint starting at out[pos], advancing pos. Returns 0 on overflow.
inline uint64_t varint_read(uint8_t const* buf, size_t buf_len, size_t& pos) {
    uint64_t v = 0;
    unsigned shift = 0;
    while (pos < buf_len) {
        uint8_t b = buf[pos++];
        v |= (uint64_t)(b & 0x7f) << shift;
        if (!(b & 0x80)) return v;
        shift += 7;
        if (shift >= 64) break;
    }
    return v;
}

// ---- Super-k-mer record ------------------------------------------------------

inline constexpr uint8_t SK_FLAG_IS_ACGT_BEGIN = 1u << 0;
inline constexpr uint8_t SK_FLAG_IS_ACGT_END = 1u << 1;
// BCALM2 boundary-k-mer ownership. A boundary k-mer X (shared by k-overlap
// between two adjacent super-k-mers) is colored in exactly ONE bucket: the
// bucket of min(lmin(X), rmin(X)) where lmin/rmin are the minimizers of X's
// prefix/suffix (k-1)-mers (intrinsic to X, so all of X's occurrences agree on
// it). At a split between an ending super (run minimizer cur_min) and a starting
// super (new_min), the shared k-mer is the ending super's LAST k-mer and the
// starting super's FIRST k-mer; the ending super owns it iff cur_min < new_min.
// These two bits record, per super, whether it owns its first / last boundary
// k-mer; the walker colors a boundary k-mer only in the bucket that owns it.
inline constexpr uint8_t SK_FLAG_OWNS_FIRST = 1u << 2;
inline constexpr uint8_t SK_FLAG_OWNS_LAST = 1u << 3;

// Serialize a compacted super-k-mer record (one or more colors) into `out`.
// `colors` must be sorted ascending and contain no duplicates.
inline void write_super_kmer(uint8_t flags, uint32_t const* colors, uint32_t num_colors,
                             uint8_t const* bases, uint32_t len, std::vector<uint8_t>& out) {
    varint_write(num_colors, out);
    uint32_t prev = 0;
    for (uint32_t i = 0; i < num_colors; ++i) {
        uint32_t c = colors[i];
        varint_write((uint64_t)(c - prev), out);
        prev = c;
    }
    // 4 low bits carry flags (BEGIN, END, OWNS_FIRST, OWNS_LAST); len is shifted
    // up by 4. (Was 2 bits / 0x3 -- the OWNS ownership bits were being dropped.)
    varint_write(((uint64_t)len << 4) | (flags & 0xfu), out);
    pack_2bit(bases, len, out);
}

// Returns the byte length consumed; sets out_flags, out_colors and copies
// bases into out_bases. Returns 0 on malformed/EOF.
inline size_t read_super_kmer(uint8_t const* buf, size_t buf_len, uint8_t& out_flags,
                              std::vector<uint32_t>& out_colors, std::vector<uint8_t>& out_bases) {
    if (buf_len < 1) return 0;
    size_t p = 0;
    uint64_t num_colors = varint_read(buf, buf_len, p);
    if (p > buf_len) return 0;
    out_colors.resize(num_colors);
    uint32_t prev = 0;
    for (uint64_t i = 0; i < num_colors; ++i) {
        uint64_t delta = varint_read(buf, buf_len, p);
        if (p > buf_len) return 0;
        prev += (uint32_t)delta;
        out_colors[i] = prev;
    }
    if (p >= buf_len) return 0;
    uint64_t lenflags = varint_read(buf, buf_len, p);
    uint64_t len = lenflags >> 4;
    out_flags = (uint8_t)(lenflags & 0xfu);
    size_t base_bytes = (size_t)((len + 3) / 4);
    if (p + base_bytes > buf_len) return 0;
    unpack_2bit(buf + p, (size_t)len, out_bases);
    return p + base_bytes;
}

}  // namespace cdbg
