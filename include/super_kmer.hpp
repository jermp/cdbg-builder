#pragma once

// Binary record format for super-k-mers stored in a bucket file:
//
//   [color_id : varint] [length_in_bases : varint] [2-bit packed bases]
//
// Bases are packed 4 per byte, lowest 2 bits = first base of the run. Records
// are concatenated in the file with no padding between them (the byte
// containing the last base may have unused high bits, which is fine because
// the next record starts on the next byte).
//
// We deliberately omit GGCAT's INCL_BEGIN/INCL_END flag bits and minimizer
// position: this first cut does not stitch unitigs across buckets, so we
// don't need that information.

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
    for (size_t i = 0; i < n; ++i) {
        dst[i] = (src[i >> 2] >> (2 * (i & 3))) & 3u;
    }
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

inline void write_super_kmer(uint32_t color, const uint8_t* bases, uint32_t len,
                             std::vector<uint8_t>& out) {
    varint_write(color, out);
    varint_write(len, out);
    pack_2bit(bases, len, out);
}

// Returns the byte length consumed; sets out_color, out_len, and copies bases
// into out_bases. Returns 0 on malformed/EOF.
inline size_t read_super_kmer(const uint8_t* buf, size_t buf_len,
                              uint32_t& out_color, std::vector<uint8_t>& out_bases) {
    size_t p = 0;
    uint64_t color = varint_read(buf, buf_len, p);
    uint64_t len = varint_read(buf, buf_len, p);
    size_t base_bytes = (size_t)((len + 3) / 4);
    if (p + base_bytes > buf_len) return 0;
    out_color = (uint32_t)color;
    unpack_2bit(buf + p, (size_t)len, out_bases);
    return p + base_bytes;
}

}  // namespace cdgb
