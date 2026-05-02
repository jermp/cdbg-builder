#pragma once

// Streaming color-set dict.
//
// On intern() we (a) compute a 128-bit content hash of the candidate
// color list, (b) look it up in m_index, and either return the existing
// id or (c) encode the candidate's bits *immediately* into an
// in-memory bits::bit_vector::builder using the same sparse / dense /
// complementary-dense rules as hybrid_builder. We never keep the
// uncompressed std::vector<uint32_t> alive past the call.
//
// When `spill_threshold_bytes` is set (> 0) and the in-memory bvb
// crosses that threshold, the dict flushes its complete 64-bit words
// to a sidecar file and keeps the trailing partial word in memory.
// Bit-offset tracking is `m_flushed_words * 64 + m_bvb.num_bits()`,
// which stays consistent across spills.
//
// Memory per class: 32 bytes of metadata (bit_offset, bit_length,
// primary_hash, secondary_hash) plus the per-class bits inside the
// bvb (capped by spill_threshold_bytes once spilling is enabled).
//
// Dedup uses 128 bits of hash:
//   primary   = wyhash(bytes)   (from ankerl::unordered_dense::detail)
//   secondary = FNV-1a(bytes)
// Two independent hash families -> birthday collision over the run
// is ~2^-64 per pair. We rely on this to skip the byte-level equality
// check (we don't have the original colors any more).
//
// finalize() builds the full hybrid struct in `out`: bit_vector is
// loaded from sidecar (if any) + the remaining bvb; offsets are
// turned into an elias_fano from the per-class bit_offsets we
// tracked along the way.

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include <bit_vector.hpp>
#include <elias_fano.hpp>
#include <integer_codes.hpp>
#include <util.hpp>
#include <unordered_dense/unordered_dense.h>

#include "hybrid_color_sets.hpp"
#include "util.hpp"

namespace cdgb {

struct streaming_color_set_dict {
    streaming_color_set_dict(uint32_t num_colors)
        : m_num_colors(num_colors)
        , m_sparse_threshold((uint32_t)(0.25 * num_colors))
        , m_dense_threshold((uint32_t)(0.75 * num_colors))
        , m_index(0, hasher{&m_classes}, key_eq{&m_classes}) {}

    ~streaming_color_set_dict() {
        if (m_sidecar) std::fclose(m_sidecar);
    }

    streaming_color_set_dict(streaming_color_set_dict const&) = delete;
    streaming_color_set_dict& operator=(streaming_color_set_dict const&) = delete;
    streaming_color_set_dict(streaming_color_set_dict&&) = delete;
    streaming_color_set_dict& operator=(streaming_color_set_dict&&) = delete;

    // Optional: enable disk spilling. If `spill_threshold_bytes == 0`,
    // never spill (keep everything in memory). The sidecar file at
    // `path` is created on the first spill and deleted by finalize()
    // after we read it back.
    void enable_spill(std::string path, uint64_t spill_threshold_bytes) {
        m_sidecar_path = std::move(path);
        m_spill_threshold_bytes = spill_threshold_bytes;
    }

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

        // Spill complete words if over threshold. We always retain the
        // trailing partial word in memory so the next encode can
        // continue at the same bit position.
        if (m_spill_threshold_bytes != 0 &&
            m_bvb.data().size() * sizeof(uint64_t) >= m_spill_threshold_bytes) {
            spill_();
        }

        return id;
    }

    uint32_t size() const { return (uint32_t)m_classes.size(); }
    uint64_t total_integers() const { return m_total_integers; }
    uint64_t total_bits() const { return m_flushed_words * 64 + m_bvb.num_bits(); }
    bool spilled() const { return m_flushed_words > 0; }

    // Build the full hybrid struct from the encoded bits + per-class
    // offsets. If we spilled to a sidecar during intern(), that file
    // is read back here, combined with any in-memory bits, and the
    // sidecar is removed.
    void finalize(hybrid& out) {
        out.m_num_colors = m_num_colors;
        out.m_sparse_set_threshold_size = m_sparse_threshold;
        out.m_very_dense_set_threshold_size = m_dense_threshold;

        uint64_t total_bit_count = total_bits();

        if (m_flushed_words == 0) {
            // Never spilled -- bvb already has everything.
            m_bvb.build(out.m_color_sets);
        } else {
            // Append remaining bvb words to sidecar, then read the
            // whole file back into a fresh bit_vector::builder of the
            // correct bit count.
            flush_remaining_to_sidecar_();

            uint64_t total_words = (total_bit_count + 63) / 64;
            std::fflush(m_sidecar);
            std::rewind(m_sidecar);

            bits::bit_vector::builder b;
            b.resize(total_bit_count);
            auto& dst = b.data();
            if (total_words > 0) {
                size_t got = std::fread(dst.data(), sizeof(uint64_t), total_words, m_sidecar);
                if (got != total_words)
                    throw std::runtime_error("sidecar short read at finalize");
            }
            b.build(out.m_color_sets);

            std::fclose(m_sidecar);
            m_sidecar = nullptr;
            std::error_code ec;
            std::filesystem::remove(m_sidecar_path, ec);
        }

        std::vector<uint64_t> offsets;
        offsets.reserve(m_classes.size() + 1);
        for (auto const& e : m_classes) offsets.push_back(e.bit_offset);
        offsets.push_back(out.m_color_sets.num_bits());  // sentinel
        out.m_offsets.encode(offsets.begin(), offsets.size(), offsets.back());

        std::cout << "  num_color_sets = " << m_classes.size() << "\n";
        std::cout << "  num_total_integers = " << m_total_integers << "\n";
        std::cout << "  total bits for ints  = " << 8 * out.m_color_sets.num_bytes() << "\n";
        std::cout << "  total bits for offs  = " << 8 * out.m_offsets.num_bytes() << "\n";
        if (m_flushed_words > 0) {
            std::cout << "  color-set bits spilled to disk: yes ("
                      << format_bytes(m_flushed_words * sizeof(uint64_t)) << ")\n";
        }
    }

private:
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

    void open_sidecar_if_needed_() {
        if (m_sidecar) return;
        if (m_sidecar_path.empty())
            throw std::runtime_error(
                "streaming_color_set_dict: spill triggered but no sidecar path set");
        m_sidecar = std::fopen(m_sidecar_path.c_str(), "wb+");
        if (!m_sidecar)
            throw std::runtime_error("cannot open sidecar: " + m_sidecar_path);
    }

    // Flush every COMPLETE 64-bit word to the sidecar; keep the
    // trailing partial word (if any) in memory so the next encode
    // continues exactly where this one left off.
    void spill_() {
        uint64_t bits = m_bvb.num_bits();
        uint64_t complete_words = bits / 64;
        if (complete_words == 0) return;

        open_sidecar_if_needed_();
        auto& words = m_bvb.data();
        size_t bytes = complete_words * sizeof(uint64_t);
        if (std::fwrite(words.data(), 1, bytes, m_sidecar) != bytes)
            throw std::runtime_error("sidecar short write");
        m_flushed_words += complete_words;

        uint64_t trailing_bits = bits - complete_words * 64;
        uint64_t partial_word = trailing_bits ? words[complete_words] : 0;
        m_bvb.clear();
        if (trailing_bits) m_bvb.append_bits(partial_word, trailing_bits);
    }

    // At finalize time, write the remainder of the bvb (whole words,
    // including the trailing partial word as a final word with the
    // unused upper bits zero) to the sidecar.
    void flush_remaining_to_sidecar_() {
        if (m_bvb.num_bits() == 0) return;
        open_sidecar_if_needed_();
        auto const& words = m_bvb.data();
        size_t n_words = words.size();
        if (n_words == 0) return;
        size_t bytes = n_words * sizeof(uint64_t);
        if (std::fwrite(words.data(), 1, bytes, m_sidecar) != bytes)
            throw std::runtime_error("sidecar short write at finalize");
        m_flushed_words += n_words;  // (note: we've now over-counted by trailing zeros, but
                                     // the bit count is tracked separately via total_bit_count
                                     // computed *before* this call.)
        m_bvb.clear();
    }

    uint32_t m_num_colors;
    uint32_t m_sparse_threshold;
    uint32_t m_dense_threshold;

    bits::bit_vector::builder m_bvb;
    uint64_t m_flushed_words = 0;  // 64-bit words already on the sidecar
    std::vector<class_entry> m_classes;
    ankerl::unordered_dense::set<uint32_t, hasher, key_eq> m_index;
    uint64_t m_total_integers = 0;

    std::string m_sidecar_path;
    FILE* m_sidecar = nullptr;
    uint64_t m_spill_threshold_bytes = 0;  // 0 = never spill
};

}  // namespace cdgb
