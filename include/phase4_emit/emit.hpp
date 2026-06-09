#pragma once

// Phase 4: emit. Turns the stitch's cid-range unitig buckets into the two
// final on-disk outputs -- <out>.fa (cid-ascending FASTA) + <out>.u2c
// (unitig->color-set bit_vector) -- and finalizes <out>.color_sets (whose
// bit_vector was already streamed out during bucket-process). These free
// functions were extracted from builder so phase 4 has its own home; the
// builder orchestrator calls them after stitch.

#include <charconv>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <bit_vector.hpp>
#include <essentials.hpp>

#include "util.hpp"  // timer, format_bytes
#include "phase2_bucket_process/streaming_color_set_dict.hpp"
#include "phase2_bucket_process/unitig_spill.hpp"  // unitig_bucket_writer

namespace cdbg {

// Streams a bits::bit_vector straight to disk, set bit-by-bit in ASCENDING
// position order, holding only the current 64-bit partial word in RAM instead
// of the whole num_bits/8-byte bitmap. The on-disk layout is byte-identical to
// essentials::save(bits::bit_vector): [u64 num_bits][u64 num_words][num_words *
// u64] (bit_vector::visit emits m_num_bits then the m_data word vector, which
// essentials serializes as an 8-byte length + the raw words). Both totals are
// known up front (num_bits is the unitig count), so the header is written
// immediately -- no seek-back. Used for the u2c run-end bitmap, whose set bits
// are produced strictly in .fa emission order, so it never needs to sit in RAM.
struct streaming_bit_vector_writer {
    streaming_bit_vector_writer(std::string const& path, uint64_t num_bits)
        : m_num_words((num_bits + 63) / 64) {
        m_file = std::fopen(path.c_str(), "wb");
        if (!m_file) throw std::runtime_error("cannot open " + path);
        write_pod_(num_bits);               // bit_vector::m_num_bits
        write_pod_((uint64_t)m_num_words);  // m_data vector length (8 B on LP64)
    }
    ~streaming_bit_vector_writer() {
        if (m_file) std::fclose(m_file);
    }
    streaming_bit_vector_writer(streaming_bit_vector_writer const&) = delete;
    streaming_bit_vector_writer& operator=(streaming_bit_vector_writer const&) = delete;

    // Set bit `pos`. Positions must arrive non-decreasing (ascending across the
    // whole stream); a position in a later word flushes the words before it.
    void set(uint64_t pos) {
        uint64_t word_idx = pos >> 6;
        while (m_flushed_words < word_idx) flush_cur_();
        m_cur_word |= uint64_t(1) << (pos & 63);
    }

    // Flush the trailing word + any remaining (zero) words and close.
    void finish() {
        while (m_flushed_words < m_num_words) flush_cur_();
        std::fflush(m_file);
        std::fclose(m_file);
        m_file = nullptr;
    }

private:
    void flush_cur_() {
        if (std::fwrite(&m_cur_word, sizeof(m_cur_word), 1, m_file) != 1)
            throw std::runtime_error("short write of u2c word");
        m_cur_word = 0;
        ++m_flushed_words;
    }
    template <typename T>
    void write_pod_(T v) {
        if (std::fwrite(&v, sizeof(T), 1, m_file) != 1)
            throw std::runtime_error("short header write of u2c");
    }
    std::FILE* m_file = nullptr;
    uint64_t m_num_words;
    uint64_t m_flushed_words = 0;
    uint64_t m_cur_word = 0;
};

// Walk the K cid-range unitig buckets in ascending order; within each, sort by
// cid (so the .fa is strictly cid-ascending) and write `>cid\n<seq>\n` with a
// hand-rolled 1 MiB buffer + std::to_chars (far faster than std::ofstream on
// millions of small records). In the same pass build the u2c bit_vector: bit i
// is set iff unitig i (in .fa emission order) is the last unitig of a color-set
// run -- length = num_unitigs, popcount = num_color_classes -- which Fulgor
// reads via rank1(unitig_id). The per-bucket cid-sort is RAM-capped
// (read_bucket_sorted external merge-sorts a bucket that exceeds the cap), so
// emit peak is bounded regardless of cid skew.
inline void emit_fasta(unitig_bucket_writer& uwriter, std::string const& out_basename,
                       uint64_t num_unitigs, double max_ram_gb) {
    timer _("emit fasta");

    FILE* fa = std::fopen((out_basename + ".fa").c_str(), "wb");
    if (!fa) throw std::runtime_error("cannot open " + out_basename + ".fa");
    constexpr size_t BUF_BYTES = 1 << 20;
    std::vector<char> buf(BUF_BYTES);
    size_t pos = 0;
    auto flush_buf = [&] {
        if (pos == 0) return;
        if (std::fwrite(buf.data(), 1, pos, fa) != pos) {
            std::fclose(fa);
            throw std::runtime_error("short write to " + out_basename + ".fa");
        }
        pos = 0;
    };
    auto reserve = [&](size_t n) {
        if (pos + n > BUF_BYTES) flush_buf();
    };

    // u2c is streamed straight to <out>.u2c word-by-word: its set bits (run
    // ends) are produced in ascending unitig order, so it never needs to sit in
    // RAM as a num_unitigs-bit bitmap (a non-spillable -g violator at scale).
    streaming_bit_vector_writer u2c_writer(out_basename + ".u2c", num_unitigs);
    size_t emitted = 0;
    uint64_t prev_cid = 0;

    // Per-bucket RAM cap for the cid-sort: a modest share of -g (fixed fallback
    // with no -g). read_bucket_sorted sorts in RAM when a bucket fits, else
    // external merge-sorts, so the emit peak is bounded REGARDLESS of cid skew.
    const uint64_t emit_mem_cap = max_ram_gb > 0
                                      ? (uint64_t)(max_ram_gb * 1024.0 * 1024.0 * 1024.0 * 0.10)
                                      : (uint64_t)(1ull << 30);  // 1 GiB default

    auto emit_record = [&](uint64_t cid, std::string_view seq) {
        // cid is globally non-decreasing (buckets ascending by cid range, sorted
        // within), so cid != prev_cid is exactly a color-set group boundary.
        if (emitted > 0 and cid != prev_cid) u2c_writer.set(emitted - 1);
        prev_cid = cid;
        ++emitted;

        reserve(22);  // '>' + up to 20 digits (uint64_t) + '\n'
        buf[pos++] = '>';
        auto rr = std::to_chars(buf.data() + pos, buf.data() + pos + 20, cid);
        pos = (size_t)(rr.ptr - buf.data());
        buf[pos++] = '\n';
        size_t s_pos = 0;
        while (s_pos < seq.size()) {
            if (pos == BUF_BYTES) flush_buf();
            size_t take = std::min(BUF_BYTES - pos, seq.size() - s_pos);
            std::memcpy(buf.data() + pos, seq.data() + s_pos, take);
            pos += take;
            s_pos += take;
        }
        reserve(1);
        buf[pos++] = '\n';
    };

    for (uint32_t b = 0; b < uwriter.num_buckets(); ++b) {
        uwriter.read_bucket_sorted(b, emit_mem_cap, emit_record);
    }
    flush_buf();
    std::fclose(fa);

    std::cout << "  [emit-fasta] " << uwriter.num_buckets() << " buckets ("
              << uwriter.spilled_buckets() << " spilled to disk, rest read from RAM), cid-sort RAM cap "
              << format_bytes(emit_mem_cap) << "\n";

    // Close out the very last run, then flush the trailing/zero words.
    if (emitted > 0) u2c_writer.set(emitted - 1);
    u2c_writer.finish();

    uwriter.close_and_unlink();
}

// Finalize <out>.color_sets. The bit_vector was streamed out word-by-word
// during bucket-process (each global_dict.intern flushed complete 64-bit
// words); finalize() flushes the trailing partial word, builds + appends the
// Elias-Fano over per-class bit-offsets, then fseeks back to write the
// now-known header totals.
inline void emit_colors(streaming_color_set_dict& global_dict) {
    timer _("emit color_sets");
    global_dict.finalize();
}

}  // namespace cdbg
