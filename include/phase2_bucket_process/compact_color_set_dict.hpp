#pragma once

// In-memory hybrid-encoded color-set dict for the bucket walker.
//
// Same dedup-by-128-bit-hash story as streaming_color_set_dict, but
// the encoded bits stay in memory (one per dict instance) rather than
// being flushed to a file. Used inside process_bucket for both
// record_sets (per-record interned color lists) and local_dict (the
// per-bucket final color classes assigned to each k-mer). With 16
// threads in flight, holding the source color lists as live
// std::vector<uint32_t>'s grew to ~10 GB on the 25K-genome workload;
// hybrid-encoding them here cuts that to ~hundreds of MB.
//
// API:
//   - intern(colors)          -> id, identical hashing rules as the
//                                 streaming dict
//   - at(id, scratch)         -> fills the caller-supplied scratch
//                                 vector with the decoded color list.
//                                 Returns no reference because the
//                                 bits would otherwise be re-decoded
//                                 on every access.
//   - size()                  -> number of distinct classes in the dict
//
// No mutable_at / classes() iterator: those would force decoding all
// classes at once, defeating the memory win.
//
// Dedup is 128-bit hash only; we don't keep the original colors after
// encoding, so byte-level confirmation isn't possible. Collision rate
// per pair is ~2^-64; with ~5000 entries per bucket and ~2048 buckets
// that's ~3e-9 chance of any collision over a run.
//
// Per-class metadata: 24 B (primary hash, secondary hash, bit_offset).
// Plus the hybrid-encoded bits in m_bvb. For 5000 classes per bucket
// at ~14000 colors each, expect ~10-30 MB per bucket.

#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <bit_vector.hpp>
#include <integer_codes.hpp>
#include <unordered_dense/unordered_dense.h>

#include "phase2_bucket_process/hybrid_color_sets.hpp"

namespace cdbg {

struct compact_color_set_dict {
    explicit compact_color_set_dict(uint32_t num_colors)
        : m_num_colors(num_colors)
        , m_sparse_threshold((uint32_t)(0.25 * num_colors))
        , m_dense_threshold((uint32_t)(0.75 * num_colors))
        , m_index(0, hasher{&m_classes}, key_eq{&m_classes}) {}

    compact_color_set_dict(compact_color_set_dict const&) = delete;
    compact_color_set_dict& operator=(compact_color_set_dict const&) = delete;
    compact_color_set_dict(compact_color_set_dict&&) = delete;
    compact_color_set_dict& operator=(compact_color_set_dict&&) = delete;

    uint32_t intern(std::vector<uint32_t> const& colors) { return intern_impl_(colors); }
    uint32_t intern(std::vector<uint32_t>&& colors) { return intern_impl_(colors); }

    uint32_t size() const { return (uint32_t)m_classes.size(); }
    uint32_t num_colors() const { return m_num_colors; }

    // Decode class `id` into `out`. `out` is cleared then filled in
    // ascending-color order, matching the original input.
    void at(uint32_t id, std::vector<uint32_t>& out) const {
        if (id >= m_classes.size())
            throw std::out_of_range("compact_color_set_dict::at id out of range");
        out.clear();
        bits::bit_vector::iterator it(
            m_bvb.data().empty() ? nullptr : m_bvb.data().data(),
            m_bvb.data().size(),
            m_classes[id].bit_offset);
        decode_one_(it, out);
    }

private:
    struct class_meta {
        uint64_t primary;
        uint64_t secondary;
        uint64_t bit_offset;
    };

    struct hash_pair {
        uint64_t primary;
        uint64_t secondary;
    };

    struct hasher {
        std::vector<class_meta> const* classes;
        using is_transparent = void;
        using is_avalanching = void;
        size_t operator()(uint32_t id) const noexcept { return (*classes)[id].primary; }
        size_t operator()(hash_pair const& h) const noexcept { return h.primary; }
    };

    struct key_eq {
        std::vector<class_meta> const* classes;
        using is_transparent = void;
        bool operator()(uint32_t a, uint32_t b) const noexcept {
            auto const& ea = (*classes)[a];
            auto const& eb = (*classes)[b];
            return ea.primary == eb.primary && ea.secondary == eb.secondary;
        }
        bool operator()(uint32_t a, hash_pair const& h) const noexcept {
            auto const& ea = (*classes)[a];
            return ea.primary == h.primary && ea.secondary == h.secondary;
        }
        bool operator()(hash_pair const& h, uint32_t a) const noexcept { return (*this)(a, h); }
    };

    uint32_t intern_impl_(std::vector<uint32_t> const& colors) {
        const uint64_t primary = wyhash_(colors);
        const uint64_t secondary = fnv1a_(colors);
        auto it = m_index.find(hash_pair{primary, secondary});
        if (it != m_index.end()) return *it;

        const uint64_t bit_offset = m_bvb.num_bits();
        hybrid_builder::encode_one(m_bvb, colors.data(), colors.size(), m_num_colors,
                                   m_sparse_threshold, m_dense_threshold);

        const uint32_t id = (uint32_t)m_classes.size();
        m_classes.push_back({primary, secondary, bit_offset});
        m_index.insert(id);
        return id;
    }

    // Inverse of hybrid_builder::encode_one. Reads exactly the bits
    // written by encode_one starting at the iterator's current position.
    void decode_one_(bits::bit_vector::iterator& it, std::vector<uint32_t>& out) const {
        const uint64_t N = bits::util::read_delta(it);
        if (N == 0) return;
        out.reserve(N);
        if (N < m_sparse_threshold) {
            // Sparse: write_delta(set[0]); write_delta(set[i]-set[i-1]-1)
            uint32_t prev = (uint32_t)bits::util::read_delta(it);
            out.push_back(prev);
            for (uint64_t i = 1; i < N; ++i) {
                const uint32_t gap = (uint32_t)bits::util::read_delta(it);
                const uint32_t v = prev + gap + 1;
                out.push_back(v);
                prev = v;
            }
        } else if (N < m_dense_threshold) {
            // Plain bitmap of m_num_colors bits.
            for (uint64_t i = 0; i < m_num_colors; ++i) {
                if (it.take(1)) out.push_back((uint32_t)i);
            }
        } else {
            // Complementary delta-gaps over (m_num_colors - N) absent values.
            const uint64_t comp = m_num_colors - N;
            // Reuse one allocation for the absent set; we throw it away
            // before returning so reserve is enough.
            std::vector<uint32_t> absent;
            absent.reserve(comp);
            if (comp > 0) {
                uint32_t prev = (uint32_t)bits::util::read_delta(it);
                absent.push_back(prev);
                for (uint64_t i = 1; i < comp; ++i) {
                    const uint32_t gap = (uint32_t)bits::util::read_delta(it);
                    const uint32_t v = prev + gap + 1;
                    absent.push_back(v);
                    prev = v;
                }
            }
            // Emit positions in [0, m_num_colors) NOT in absent.
            size_t a = 0;
            for (uint32_t v = 0; v < m_num_colors; ++v) {
                if (a < absent.size() && absent[a] == v) {
                    ++a;
                } else {
                    out.push_back(v);
                }
            }
        }
    }

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
    std::vector<class_meta> m_classes;
    ankerl::unordered_dense::set<uint32_t, hasher, key_eq> m_index;
};

}  // namespace cdbg
