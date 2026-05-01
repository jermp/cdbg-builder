#pragma once

// Binary record format for super-k-mers stored in a bucket file:
//
//   [flags : u8] [color_id : varint] [length_in_bases : varint] [2-bit packed bases]
//
// Flags:
//   bit 0 (IS_ACGT_BEGIN) = this super-k-mer's first k-mer (in input order) is
//                           the actual first k-mer of an ACGT-only run, i.e. it
//                           has no predecessor in the input.
//   bit 1 (IS_ACGT_END)   = analogously, this super-k-mer's last k-mer is the
//                           last k-mer of its ACGT run.
//
// When IS_ACGT_BEGIN is unset, this super-k-mer's first k-mer has a predecessor
// k-mer in the input that lives in a different bucket; the cross-bucket
// stitcher uses that fact later. Same for IS_ACGT_END on the right side.
//
// Bases are packed 4 per byte, lowest 2 bits = first base of the run. Records
// are concatenated in the file with no padding between them.

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

namespace cdgb {

// ---- 2-bit packing helpers ---------------------------------------------------

// Append `n` 2-bit-encoded bases from `src` (each entry in [0,3]) to `dst`.
// Bases are packed into bytes with base i occupying bits (2*(i%4)..2*(i%4)+1)
// of byte (i/4).
inline void pack_2bit(const uint8_t* src, size_t n, std::vector<uint8_t>& dst) {
    size_t bytes_needed = (n + 3) / 4;
    size_t off = dst.size();
    dst.resize(off + bytes_needed, 0);
    for (size_t i = 0; i < n; ++i) {
        dst[off + (i >> 2)] |= (uint8_t)((src[i] & 3u) << (2 * (i & 3)));
    }
}

// Unpack `n` bases from `src` (starting at base index 0 in src's first byte)
// into `dst` (cleared and resized to n).
inline void unpack_2bit(const uint8_t* src, size_t n, std::vector<uint8_t>& dst) {
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
inline uint64_t varint_read(const uint8_t* buf, size_t buf_len, size_t& pos) {
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
inline constexpr uint8_t SK_FLAG_IS_ACGT_END   = 1u << 1;

inline void write_super_kmer(uint8_t flags, uint32_t color, const uint8_t* bases, uint32_t len,
                             std::vector<uint8_t>& out) {
    out.push_back(flags);
    varint_write(color, out);
    varint_write(len, out);
    pack_2bit(bases, len, out);
}

// Returns the byte length consumed; sets out_flags, out_color and copies bases
// into out_bases. Returns 0 on malformed/EOF.
inline size_t read_super_kmer(const uint8_t* buf, size_t buf_len, uint8_t& out_flags,
                              uint32_t& out_color, std::vector<uint8_t>& out_bases) {
    if (buf_len < 1) return 0;
    size_t p = 0;
    out_flags = buf[p++];
    uint64_t color = varint_read(buf, buf_len, p);
    uint64_t len = varint_read(buf, buf_len, p);
    size_t base_bytes = (size_t)((len + 3) / 4);
    if (p + base_bytes > buf_len) return 0;
    out_color = (uint32_t)color;
    unpack_2bit(buf + p, (size_t)len, out_bases);
    return p + base_bytes;
}

}  // namespace cdgb
