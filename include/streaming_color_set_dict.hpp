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
// Memory per class: 32 bytes of metadata (bit_offset, bit_length,
// primary_hash, secondary_hash) plus the per-class bits inside the
// bvb -- which is the same content that ends up in the .colors file
// at emit time. So peak RAM for the dict is ~ |.colors| + 32 B/class
// instead of "sum of class sizes × 4 bytes" (which on the
// 4546-genome run was ~8.6 GiB vs ~191 MiB for the bvb).
//
// Dedup uses 128 bits of hash:
//   primary   = wyhash(bytes)   (from ankerl::unordered_dense::detail)
//   secondary = FNV-1a(bytes)
// Two independent hash families -> birthday collision over the run
// is ~2^-64 per pair. We rely on this to skip the byte-level equality
// check (we don't have the original colors any more).
//
// finalize() builds the full hybrid struct in `out`: bit_vector is
// moved out of the builder; offsets are turned into an elias_fano
// from the per-class bit_offsets we tracked along the way.

#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

#include <bit_vector.hpp>
#include <elias_fano.hpp>
#include <integer_codes.hpp>
#include <util.hpp>
#include <unordered_dense/unordered_dense.h>

#include "hybrid_color_sets.hpp"

namespace cdgb {

struct streaming_color_set_dict {
    streaming_color_set_dict(uint32_t num_colors)
        : m_num_colors(num_colors)
        , m_sparse_threshold((uint32_t)(0.25 * num_colors))
        , m_dense_threshold((uint32_t)(0.75 * num_colors))
        , m_index(0, hasher{&m_classes}, key_eq{&m_classes}) {}

    streaming_color_set_dict(streaming_color_set_dict const&) = delete;
    streaming_color_set_dict& operator=(streaming_color_set_dict const&) = delete;
    streaming_color_set_dict(streaming_color_set_dict&&) = delete;
    streaming_color_set_dict& operator=(streaming_color_set_dict&&) = delete;

    uint32_t intern(std::vector<uint32_t>&& candidate) {
        uint64_t primary = wyhash_(candidate);
        uint64_t secondary = fnv1a_(candidate);

        // Heterogeneous lookup: hashes/compares against the (primary,
        // secondary) pair without going through m_classes for the
        // candidate side.
        auto it = m_index.find(hash_pair{primary, secondary});
        if (it != m_index.end()) return *it;

        uint64_t bit_offset = m_bvb.num_bits();
        hybrid_builder::encode_one(m_bvb, candidate.data(), candidate.size(), m_num_colors,
                                   m_sparse_threshold, m_dense_threshold);
        uint64_t bit_length = m_bvb.num_bits() - bit_offset;

        uint32_t id = (uint32_t)m_classes.size();
        m_classes.push_back({bit_offset, bit_length, primary, secondary});
        m_index.insert(id);
        m_total_integers += candidate.size();
        return id;
    }

    uint32_t size() const { return (uint32_t)m_classes.size(); }
    uint64_t total_integers() const { return m_total_integers; }
    uint64_t total_bits() const { return m_bvb.num_bits(); }

    // Build the full hybrid struct from the encoded bvb + the offsets
    // tracked during build. Consumes the internal bit_vector::builder
    // (after this call the dict is empty).
    void finalize(hybrid& out) {
        out.m_num_colors = m_num_colors;
        out.m_sparse_set_threshold_size = m_sparse_threshold;
        out.m_very_dense_set_threshold_size = m_dense_threshold;

        m_bvb.build(out.m_color_sets);

        std::vector<uint64_t> offsets;
        offsets.reserve(m_classes.size() + 1);
        for (auto const& e : m_classes) offsets.push_back(e.bit_offset);
        offsets.push_back(out.m_color_sets.num_bits());  // sentinel
        out.m_offsets.encode(offsets.begin(), offsets.size(), offsets.back());

        std::cout << "  num_color_sets = " << m_classes.size() << "\n";
        std::cout << "  num_total_integers = " << m_total_integers << "\n";
        std::cout << "  total bits for ints  = " << 8 * out.m_color_sets.num_bytes() << "\n";
        std::cout << "  total bits for offs  = " << 8 * out.m_offsets.num_bytes() << "\n";
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

    uint32_t m_num_colors;
    uint32_t m_sparse_threshold;
    uint32_t m_dense_threshold;

    bits::bit_vector::builder m_bvb;
    std::vector<class_entry> m_classes;
    ankerl::unordered_dense::set<uint32_t, hasher, key_eq> m_index;
    uint64_t m_total_integers = 0;
};

}  // namespace cdgb
