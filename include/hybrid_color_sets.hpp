#pragma once

// Hybrid color-set encoder primitive.
//
// Encodes a single sorted-deduped color list into a bits::bit_vector
// using the same sparse / dense / very-dense (complementary delta-gap)
// rules as Fulgor's `hybrid` color-set format.
//
// Per color set (sorted, deduped uint32_t list of size N, alphabet
// [0, num_colors), with thresholds sparse_threshold = 0.25 * num_colors
// and very_dense_threshold = 0.75 * num_colors):
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
// This file used to also define the full `hybrid` struct + builder
// (with elias_fano offsets, a `hybrid` instance built via the
// essentials visitor pattern, and a parallel encode helper). That
// code path is gone: the streaming color-set dict writes the on-
// disk format incrementally to <basename>.color_sets, and the
// per-bucket compact dict holds its hybrid-encoded bits in memory.
// Both call this static `encode_one` directly. No instances of
// `hybrid_builder` are constructed anywhere; only the static method
// matters.

#include <bit_vector.hpp>
#include <integer_codes.hpp>

#include <cassert>
#include <cstdint>

namespace cdgb {

struct hybrid_builder {
    // Append one color set's bits to `out_bvb` using the sparse /
    // dense / complementary-dense rules described at the top of this
    // file. No shared state; safe to call from many threads each
    // with its own bit_vector::builder. The thresholds are passed in
    // explicitly (vs. read from a member) so callers don't need to
    // construct a hybrid_builder instance.
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
};

}  // namespace cdgb
