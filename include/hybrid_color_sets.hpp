#pragma once

// Encoder for color sets in Fulgor's "hybrid" on-disk format.
//
// Byte-for-byte compatible with `fulgor::hybrid` as serialized by
// essentials::save(...). Member layout (visited in this exact order):
//
//   uint32_t m_num_colors
//   uint32_t m_sparse_set_threshold_size       = floor(0.25 * num_colors)
//   uint32_t m_very_dense_set_threshold_size   = floor(0.75 * num_colors)
//   bits::elias_fano<false, false> m_offsets   // size = num_color_sets + 1
//   bits::bit_vector               m_color_sets
//
// Encoding rules per color set (sorted, deduped uint32_t list of size N):
//   (1) write_delta(N)
//   (2) if N < sparse_threshold:
//          delta-gaps:
//            write_delta(set[0])
//            for i in 1..N: write_delta(set[i] - set[i-1] - 1)
//   (3) elif N < very_dense_threshold:
//          plain bitmap of `num_colors` bits with those positions set
//   (4) else:
//          delta-gaps over the *complementary* set (the (num_colors - N)
//          absent values), encoded identically to (2).
//
// Offsets array stores the bit offset into m_color_sets at which each color
// set begins, plus a sentinel = total bit count.

#include <bit_vector.hpp>
#include <elias_fano.hpp>
#include <integer_codes.hpp>
#include <util.hpp>
#include <essentials.hpp>

#include <cassert>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>

#include "color_set_dict.hpp"

namespace cdgb {

struct hybrid {
    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }
    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    uint32_t m_num_colors = 0;
    uint32_t m_sparse_set_threshold_size = 0;
    uint32_t m_very_dense_set_threshold_size = 0;
    bits::elias_fano<false, false> m_offsets;
    bits::bit_vector m_color_sets;

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_num_colors);
        visitor.visit(t.m_sparse_set_threshold_size);
        visitor.visit(t.m_very_dense_set_threshold_size);
        visitor.visit(t.m_offsets);
        visitor.visit(t.m_color_sets);
    }
};

struct hybrid_builder {
    hybrid_builder() = default;
    explicit hybrid_builder(uint64_t num_colors) { init(num_colors); }

    void init(uint64_t num_colors) {
        m_num_colors = num_colors;
        m_sparse_set_threshold_size = 0.25 * num_colors;
        m_very_dense_set_threshold_size = 0.75 * num_colors;
        m_offsets.clear();
        m_offsets.push_back(0);
        m_bvb.clear();
        m_num_color_sets = 0;
        m_num_total_integers = 0;
    }

    void encode_color_set(uint32_t const* color_set, const uint64_t size) {
        encode_one(m_bvb, color_set, size, m_num_colors, m_sparse_set_threshold_size,
                   m_very_dense_set_threshold_size);
        m_offsets.push_back(m_bvb.num_bits());
        m_num_total_integers += size;
        ++m_num_color_sets;
    }

    // Encode every color set in `dict` into this builder using up to
    // `num_threads` worker threads. Each worker takes a contiguous range
    // of class ids, encodes its range into a private bit_vector::builder
    // + relative-offsets array, then a single-thread sequential merge
    // appends the partials into this builder's m_bvb / m_offsets. The
    // encode is the heavy work (write_delta over the union of all color
    // sets); the merge is just one bit_vector::append per partial plus a
    // base-offset translation, and runs in negligible time.
    //
    // Caller can still alternate this with sequential encode_color_set
    // calls if needed -- merge_part is the same primitive used here.
    void encode_parallel(color_set_dict const& dict, uint32_t num_threads) {
        if (num_threads == 0) num_threads = 1;
        const size_t num_classes = dict.size();
        if (num_classes == 0) return;

        struct part {
            bits::bit_vector::builder bvb;
            std::vector<uint64_t> offsets;
            uint64_t total_ints = 0;
        };
        std::vector<part> parts(num_threads);

        const uint64_t num_colors = m_num_colors;
        const uint64_t sparse_thresh = m_sparse_set_threshold_size;
        const uint64_t dense_thresh = m_very_dense_set_threshold_size;

        std::vector<std::thread> workers;
        workers.reserve(num_threads);
        for (uint32_t t = 0; t < num_threads; ++t) {
            workers.emplace_back([&, t] {
                auto& p = parts[t];
                p.offsets.push_back(0);
                size_t lo = (size_t)t * num_classes / num_threads;
                size_t hi = (size_t)(t + 1) * num_classes / num_threads;
                for (size_t cid = lo; cid < hi; ++cid) {
                    auto const& cs = dict.at((uint32_t)cid);
                    encode_one(p.bvb, cs.data(), cs.size(), num_colors, sparse_thresh,
                               dense_thresh);
                    p.offsets.push_back(p.bvb.num_bits());
                    p.total_ints += cs.size();
                }
            });
        }
        for (auto& w : workers) w.join();
        for (auto& p : parts) merge_part(p.bvb, p.offsets, p.total_ints);
    }

    // Pure encoder: appends one color set's bits to `out_bvb` using the
    // same sparse / dense / complementary-dense rules as encode_color_set.
    // No shared state; safe to call from many threads each with its own
    // bit_vector::builder. The thresholds are passed in (vs. read from a
    // member) so this can be called without a builder instance, e.g. from
    // the parallel encoders inside encode_parallel.
    static void encode_one(bits::bit_vector::builder& out_bvb, uint32_t const* color_set,
                           uint64_t size, uint64_t num_colors, uint64_t sparse_threshold,
                           uint64_t dense_threshold) {
        bits::util::write_delta(out_bvb, size);
        if (size < sparse_threshold) {
            uint32_t prev = color_set[0];
            bits::util::write_delta(out_bvb, prev);
            for (uint64_t i = 1; i < size; ++i) {
                uint32_t v = color_set[i];
                assert(v >= prev + 1);
                bits::util::write_delta(out_bvb, v - (prev + 1));
                prev = v;
            }
        } else if (size < dense_threshold) {
            bits::bit_vector::builder bvb;
            bvb.resize(num_colors);
            for (uint64_t i = 0; i < size; ++i) bvb.set(color_set[i]);
            out_bvb.append(bvb);
        } else {
            bool first = true;
            uint32_t val = 0;
            uint32_t prev = uint32_t(-1);
            uint32_t written = 0;
            for (uint64_t i = 0; i < size; ++i) {
                uint32_t x = color_set[i];
                while (val < x) {
                    if (first) {
                        bits::util::write_delta(out_bvb, val);
                        first = false;
                    } else {
                        bits::util::write_delta(out_bvb, val - (prev + 1));
                    }
                    prev = val;
                    ++val;
                    ++written;
                }
                ++val;
            }
            while (val < num_colors) {
                bits::util::write_delta(out_bvb, val - (prev + 1));
                prev = val;
                ++val;
                ++written;
            }
            assert(written == num_colors - size);
            (void)written;
        }
    }

    // Append one thread's worth of pre-encoded color sets into this
    // builder. `part_offsets` must be relative bit positions: front = 0,
    // back = part_bvb.num_bits(); one entry per encoded set + the
    // sentinel. Used by the parallel emit-colors path to merge per-thread
    // partials into a single hybrid_builder before build().
    void merge_part(bits::bit_vector::builder& part_bvb, std::vector<uint64_t> const& part_offsets,
                    uint64_t part_num_integers) {
        if (part_offsets.empty()) return;
        uint64_t base = m_bvb.num_bits();
        m_bvb.append(part_bvb);
        // Skip the leading 0; m_offsets already has its base entry.
        for (size_t i = 1; i < part_offsets.size(); ++i) {
            m_offsets.push_back(base + part_offsets[i]);
        }
        m_num_total_integers += part_num_integers;
        m_num_color_sets += part_offsets.size() - 1;
    }

    void build(hybrid& h) {
        h.m_num_colors = m_num_colors;
        h.m_sparse_set_threshold_size = m_sparse_set_threshold_size;
        h.m_very_dense_set_threshold_size = m_very_dense_set_threshold_size;
        h.m_offsets.encode(m_offsets.begin(), m_offsets.size(), m_offsets.back());
        m_bvb.build(h.m_color_sets);

        std::cout << "  num_color_sets = " << m_num_color_sets << "\n";
        std::cout << "  num_total_integers = " << m_num_total_integers << "\n";
        std::cout << "  total bits for ints  = " << 8 * h.m_color_sets.num_bytes() << "\n";
        std::cout << "  total bits for offs  = " << 8 * h.m_offsets.num_bytes() << "\n";
    }

    uint64_t num_color_sets() const { return m_num_color_sets; }

private:
    uint64_t m_num_colors = 0;
    uint64_t m_sparse_set_threshold_size = 0;
    uint64_t m_very_dense_set_threshold_size = 0;
    uint64_t m_num_color_sets = 0;
    uint64_t m_num_total_integers = 0;
    bits::bit_vector::builder m_bvb;
    std::vector<uint64_t> m_offsets;
};

}  // namespace cdgb
