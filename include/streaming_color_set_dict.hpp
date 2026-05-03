#pragma once

// Streaming color-set dict that writes its output incrementally to the
// final on-disk artifact, never holding the whole compressed bit_vector
// or the EF offsets in memory at once.
//
// On intern() we (a) compute a 128-bit content hash of the candidate
// color list, (b) look it up in m_index, and either return the existing
// id or (c) encode the candidate's bits into an in-memory
// bits::bit_vector::builder (m_bvb) using the hybrid sparse / dense /
// complementary-dense rules, then immediately flush every COMPLETE 64-
// bit word from m_bvb into the output file. We never keep more than
// the trailing partial word + the bits being actively encoded for the
// current candidate in memory.
//
// When finalize() is called, we
//   - flush the trailing partial word (zero-padded to a full 64-bit word),
//   - build a bits::elias_fano over the per-class bit_offsets we tracked
//     during interning,
//   - serialize the EF onto the file's tail,
//   - fseek back to file start and overwrite the placeholder header
//     with the now-known totals (bit_vector_num_bits / num_words /
//     num_color_sets / etc.).
//
// On-disk layout of the resulting file (matches what the user's
// downstream consumer expects: header -> color_sets -> EF offsets):
//
//   [u32 num_colors]
//   [u32 sparse_threshold]
//   [u32 dense_threshold]
//   [u64 num_color_sets]
//   [u64 bit_vector_num_bits]
//   [u64 bit_vector_num_words]   == ceil(bit_vector_num_bits / 64)
//   [bit_vector_num_words * u64] color-set bit_vector words (LE host order)
//   [serialized bits::elias_fano<false,false>]   per-class offsets
//
// Memory per class: 32 bytes of metadata (bit_offset, bit_length,
// primary_hash, secondary_hash). Memory across classes: just
// num_classes * 32 bytes + the dedup hashtable + the trailing partial
// word in m_bvb. The compressed color-set bits themselves are NOT held
// in memory beyond the flush threshold.
//
// Dedup uses 128 bits of hash:
//   primary   = wyhash(bytes)   (from ankerl::unordered_dense::detail)
//   secondary = FNV-1a(bytes)
// Two independent hash families -> birthday collision over the run is
// ~2^-64 per pair. We rely on this to skip the byte-level equality
// check (we don't have the original colors any more after encoding).

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include <bit_vector.hpp>
#include <elias_fano.hpp>
#include <essentials.hpp>
#include <integer_codes.hpp>
#include <util.hpp>
#include <unordered_dense/unordered_dense.h>

#include "hybrid_color_sets.hpp"
#include "util.hpp"

namespace cdgb {

struct streaming_color_set_dict {
    streaming_color_set_dict(uint32_t num_colors, std::string output_path)
        : m_num_colors(num_colors)
        , m_sparse_threshold((uint32_t)(0.25 * num_colors))
        , m_dense_threshold((uint32_t)(0.75 * num_colors))
        , m_index(0, hasher{&m_classes}, key_eq{&m_classes})
        , m_output_path(std::move(output_path)) {
        m_file = std::fopen(m_output_path.c_str(), "wb+");
        if (!m_file)
            throw std::runtime_error("cannot open color-set output: " + m_output_path);
        // Reserve space for the fixed-size header. We don't know the
        // final values (num_color_sets, bit_vector_num_bits, num_words)
        // until finalize(), so write zeros for now and overwrite at the
        // end via fseek.
        char hdr[HEADER_BYTES] = {};
        if (std::fwrite(hdr, 1, HEADER_BYTES, m_file) != HEADER_BYTES)
            throw std::runtime_error("short write of header to " + m_output_path);
    }

    ~streaming_color_set_dict() {
        if (m_file) std::fclose(m_file);
    }

    streaming_color_set_dict(streaming_color_set_dict const&) = delete;
    streaming_color_set_dict& operator=(streaming_color_set_dict const&) = delete;
    streaming_color_set_dict(streaming_color_set_dict&&) = delete;
    streaming_color_set_dict& operator=(streaming_color_set_dict&&) = delete;

    uint32_t intern(std::vector<uint32_t>&& candidate) {
        uint64_t primary = wyhash_(candidate);
        uint64_t secondary = fnv1a_(candidate);

        auto it = m_index.find(hash_pair{primary, secondary});
        if (it != m_index.end()) return *it;

        uint64_t bit_offset = m_flushed_words * 64 + m_bvb.num_bits();
        hybrid_builder::encode_one(m_bvb, candidate.data(), candidate.size(), m_num_colors,
                                   m_sparse_threshold, m_dense_threshold);
        uint64_t bit_length = (m_flushed_words * 64 + m_bvb.num_bits()) - bit_offset;

        uint32_t id = (uint32_t)m_classes.size();
        m_classes.push_back({bit_offset, bit_length, primary, secondary});
        m_index.insert(id);
        m_total_integers += candidate.size();

        // Flush every COMPLETE 64-bit word to the output file. The
        // trailing partial word stays in m_bvb so the next encode can
        // continue at the same bit position. With this flushed every
        // intern() call, in-memory bit storage is bounded by one word.
        spill_complete_words_();

        return id;
    }

    uint32_t size() const { return (uint32_t)m_classes.size(); }
    uint64_t total_integers() const { return m_total_integers; }
    uint64_t total_bits() const { return m_flushed_words * 64 + m_bvb.num_bits(); }

    // Finalize the on-disk file: flush the trailing partial word,
    // build & serialize the EF over per-class bit-offsets, then
    // fseek back to overwrite the placeholder header.
    void finalize() {
        if (m_finalized)
            throw std::runtime_error("streaming_color_set_dict::finalize called twice");
        m_finalized = true;

        const uint64_t total_bit_count = total_bits();
        const uint64_t total_word_count = (total_bit_count + 63) / 64;

        // 1) Flush any trailing partial word (zero-padded to 64 bits)
        //    so the file's bit_vector words section has exactly
        //    total_word_count u64s.
        if (m_bvb.num_bits() > 0) {
            auto const& words = m_bvb.data();
            // bit_vector::builder always keeps at least
            // ceil(num_bits/64) words; we want exactly that many.
            size_t want = (m_bvb.num_bits() + 63) / 64;
            if (want > words.size()) want = words.size();
            if (want > 0) {
                size_t bytes = want * sizeof(uint64_t);
                if (std::fwrite(words.data(), 1, bytes, m_file) != bytes)
                    throw std::runtime_error("short write of trailing words to " + m_output_path);
            }
            m_bvb.clear();
        }

        // 2) Build per-class offsets array (last entry = sentinel
        //    total_bit_count) and encode it as elias_fano. The EF
        //    itself is small relative to the bit_vector (typically a
        //    few bytes per class), so we hold it in memory and
        //    serialize via essentials onto the file tail.
        std::vector<uint64_t> offsets;
        offsets.reserve(m_classes.size() + 1);
        for (auto const& e : m_classes) offsets.push_back(e.bit_offset);
        offsets.push_back(total_bit_count);  // sentinel

        bits::elias_fano<false, false> ef;
        ef.encode(offsets.begin(), offsets.size(), offsets.back());

        // 3) Serialize EF to bytes via essentials (it walks ef.visit
        //    and writes pods + vec sizes), then fwrite the bytes onto
        //    the file's tail. essentials::generic_saver wants a
        //    std::ostream so we round-trip through a stringstream.
        std::ostringstream oss(std::ios::binary);
        {
            essentials::generic_saver gs(oss);
            gs.visit(ef);
        }
        std::string ef_bytes = oss.str();
        if (!ef_bytes.empty()) {
            if (std::fwrite(ef_bytes.data(), 1, ef_bytes.size(), m_file) != ef_bytes.size())
                throw std::runtime_error("short write of EF to " + m_output_path);
        }

        // 4) Update the placeholder header with the now-known totals.
        std::fflush(m_file);
        if (std::fseek(m_file, 0, SEEK_SET) != 0)
            throw std::runtime_error("fseek to header failed on " + m_output_path);
        write_pod_(m_num_colors);
        write_pod_(m_sparse_threshold);
        write_pod_(m_dense_threshold);
        write_pod_((uint64_t)m_classes.size());
        write_pod_(total_bit_count);
        write_pod_(total_word_count);

        std::fflush(m_file);
        std::fclose(m_file);
        m_file = nullptr;

        std::cout << "  num_color_sets = " << m_classes.size() << "\n";
        std::cout << "  num_total_integers = " << m_total_integers << "\n";
        std::cout << "  total bits for ints  = " << total_bit_count << "\n";
        std::cout << "  total bits for offs  = " << 8 * ef_bytes.size() << "\n";
    }

private:
    // Header layout: see top-of-file comment.
    //   u32 num_colors
    //   u32 sparse_threshold
    //   u32 dense_threshold
    //   u64 num_color_sets
    //   u64 bit_vector_num_bits
    //   u64 bit_vector_num_words
    static constexpr size_t HEADER_BYTES = 4 + 4 + 4 + 8 + 8 + 8;

    struct class_entry {
        uint64_t bit_offset;
        uint64_t bit_length;
        uint64_t primary_hash;
        uint64_t secondary_hash;
    };

    struct hash_pair {
        uint64_t primary;
        uint64_t secondary;
    };

    struct hasher {
        std::vector<class_entry> const* classes;
        using is_transparent = void;
        using is_avalanching = void;
        size_t operator()(uint32_t id) const noexcept { return (*classes)[id].primary_hash; }
        size_t operator()(hash_pair const& h) const noexcept { return h.primary; }
    };

    struct key_eq {
        std::vector<class_entry> const* classes;
        using is_transparent = void;
        bool operator()(uint32_t a, uint32_t b) const noexcept {
            auto const& ea = (*classes)[a];
            auto const& eb = (*classes)[b];
            return ea.primary_hash == eb.primary_hash && ea.secondary_hash == eb.secondary_hash;
        }
        bool operator()(uint32_t a, hash_pair const& h) const noexcept {
            auto const& ea = (*classes)[a];
            return ea.primary_hash == h.primary && ea.secondary_hash == h.secondary;
        }
        bool operator()(hash_pair const& h, uint32_t a) const noexcept { return (*this)(a, h); }
    };

    static uint64_t wyhash_(std::vector<uint32_t> const& v) noexcept {
        return ankerl::unordered_dense::detail::wyhash::hash(
            v.data(), v.size() * sizeof(uint32_t));
    }
    static uint64_t fnv1a_(std::vector<uint32_t> const& v) noexcept {
        uint64_t h = 1469598103934665603ULL;
        unsigned char const* p = (unsigned char const*)v.data();
        size_t n = v.size() * sizeof(uint32_t);
        for (size_t i = 0; i < n; ++i) {
            h ^= p[i];
            h *= 1099511628211ULL;
        }
        return h;
    }

    template <typename T>
    void write_pod_(T const& v) {
        static_assert(std::is_trivially_copyable<T>::value, "POD only");
        if (std::fwrite(&v, sizeof(T), 1, m_file) != 1)
            throw std::runtime_error("short header write to " + m_output_path);
    }

    // Flush every COMPLETE 64-bit word to the output file; keep the
    // trailing partial word (if any) in m_bvb so the next encode
    // continues exactly where this one left off.
    void spill_complete_words_() {
        uint64_t bits = m_bvb.num_bits();
        uint64_t complete_words = bits / 64;
        if (complete_words == 0) return;

        auto& words = m_bvb.data();
        size_t bytes = complete_words * sizeof(uint64_t);
        if (std::fwrite(words.data(), 1, bytes, m_file) != bytes)
            throw std::runtime_error("short write of color-set words to " + m_output_path);
        m_flushed_words += complete_words;

        uint64_t trailing_bits = bits - complete_words * 64;
        uint64_t partial_word = trailing_bits ? words[complete_words] : 0;
        m_bvb.clear();
        if (trailing_bits) m_bvb.append_bits(partial_word, trailing_bits);
    }

    uint32_t m_num_colors;
    uint32_t m_sparse_threshold;
    uint32_t m_dense_threshold;

    bits::bit_vector::builder m_bvb;
    uint64_t m_flushed_words = 0;
    std::vector<class_entry> m_classes;
    ankerl::unordered_dense::set<uint32_t, hasher, key_eq> m_index;
    uint64_t m_total_integers = 0;

    std::FILE* m_file = nullptr;
    std::string m_output_path;
    bool m_finalized = false;
};

}  // namespace cdgb
