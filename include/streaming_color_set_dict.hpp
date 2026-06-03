#pragma once

// Streaming color-set dict that writes its output incrementally to the
// final on-disk artifact, never holding the whole compressed bit_vector
// or the EF offsets array in memory at once.
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
// When a new class is interned, its bit_offset is appended (as a
// uint64_t) to a sidecar offsets file at <output_path>.tmp_offsets.
// At finalize() we
//   - flush the trailing partial word (zero-padded to a full 64-bit word),
//   - append the sentinel total_bit_count to the offsets sidecar,
//   - rewind the sidecar and feed it through bits::elias_fano::encode
//     via a single-pass file-backed iterator (so the offsets array
//     never sits in RAM as a vector),
//   - serialize the EF onto the file's tail,
//   - fseek back to file start and overwrite the placeholder header
//     with the now-known totals,
//   - delete the offsets sidecar.
//
// On-disk layout of the resulting file:
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
// Memory per class: 16 bytes of metadata (primary + secondary hash for
// dedup). Memory across classes: num_classes * 16 bytes + the dedup
// hashtable + the trailing partial word in m_bvb. The compressed
// color-set bits and the EF offsets array are NOT held in RAM. The EF
// internal structures (compact_vector for low bits + bit_vector for
// high bits) are still built in memory at finalize, sized
// approximately num_classes * ceil(log2(universe / num_classes)) bits
// for the low part and ~num_classes + universe/2^l bits for the high
// part -- both compact, and dominated by the universe / num_classes
// ratio rather than by num_classes itself.
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
#include <iterator>
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

namespace cdbg {

struct streaming_color_set_dict {
    streaming_color_set_dict(uint32_t num_colors, std::string output_path)
        : m_num_colors(num_colors)
        , m_sparse_threshold((uint32_t)(0.25 * num_colors))
        , m_dense_threshold((uint32_t)(0.75 * num_colors))
        , m_index(0, hasher{&m_classes}, key_eq{&m_classes})
        , m_output_path(std::move(output_path))
        , m_offsets_path(m_output_path + ".tmp_offsets") {
        m_file = std::fopen(m_output_path.c_str(), "wb+");
        if (!m_file)
            throw std::runtime_error("cannot open color-set output: " + m_output_path);
        char hdr[HEADER_BYTES] = {};
        if (std::fwrite(hdr, 1, HEADER_BYTES, m_file) != HEADER_BYTES)
            throw std::runtime_error("short write of header to " + m_output_path);

        m_offsets_file = std::fopen(m_offsets_path.c_str(), "wb+");
        if (!m_offsets_file)
            throw std::runtime_error("cannot open color-set offsets sidecar: " +
                                     m_offsets_path);
    }

    ~streaming_color_set_dict() {
        if (m_file) std::fclose(m_file);
        if (m_offsets_file) std::fclose(m_offsets_file);
        if (!m_offsets_path.empty()) {
            std::error_code ec;
            std::filesystem::remove(m_offsets_path, ec);
        }
    }

    streaming_color_set_dict(streaming_color_set_dict const&) = delete;
    streaming_color_set_dict& operator=(streaming_color_set_dict const&) = delete;
    streaming_color_set_dict(streaming_color_set_dict&&) = delete;
    streaming_color_set_dict& operator=(streaming_color_set_dict&&) = delete;

    // 128-bit content hash for a candidate color list. Caller can
    // compute this lock-free and pass it to intern_with_hashes(),
    // saving the per-call wyhash + fnv1a work (the dominant cost
    // of intern() on multi-thousand-color lists).
    struct precomputed_hash {
        uint64_t primary;
        uint64_t secondary;
    };
    static precomputed_hash compute_hashes(std::vector<uint32_t> const& v) noexcept {
        return {wyhash_(v), fnv1a_(v)};
    }

    uint64_t intern(std::vector<uint32_t>&& candidate) {
        return intern_with_hashes(std::move(candidate), compute_hashes(candidate));
    }

    // Variant of intern() that uses caller-supplied hashes. The two
    // hashes MUST be wyhash + fnv1a of `candidate` (i.e. produced by
    // compute_hashes()) -- internal dedup relies on this. Lets the
    // bucket-process merge phase compute hashes outside the global
    // lock so threads only serialise on the actual hashmap-find +
    // encode work.
    uint64_t intern_with_hashes(std::vector<uint32_t>&& candidate, precomputed_hash h) {
        assert(!m_released && "intern() after release_index()");
        auto it = m_index.find(hash_pair{h.primary, h.secondary});
        if (it != m_index.end()) return *it;

        uint64_t bit_offset = m_flushed_words * 64 + m_bvb.num_bits();
        hybrid_builder::encode_one(m_bvb, candidate.data(), candidate.size(), m_num_colors,
                                   m_sparse_threshold, m_dense_threshold);

        // Spill the per-class bit_offset to the sidecar (8 B per
        // intern). Read back at finalize via a single-pass iterator;
        // never lives in RAM as a vector.
        if (std::fwrite(&bit_offset, sizeof(bit_offset), 1, m_offsets_file) != 1)
            throw std::runtime_error("short write of class offset to " + m_offsets_path);

        uint64_t id = m_classes.size();
        m_classes.push_back({h.primary, h.secondary});
        m_index.insert(id);
        m_total_integers += candidate.size();

        spill_complete_words_();

        return id;
    }

    uint64_t size() const { return m_released ? m_released_count : m_classes.size(); }
    uint64_t total_integers() const { return m_total_integers; }
    uint64_t total_bits() const { return m_flushed_words * 64 + m_bvb.num_bits(); }

    // Free the dedup index + per-class hashes once interning is DONE (after
    // bucket-process). Neither m_index nor m_classes is read by finalize() --
    // it only needs the class COUNT (stashed here) and the on-disk bits +
    // sidecar offsets. This dict otherwise stays fully resident through stitch
    // AND emit (it's freed at finalize today), and at high class counts that's
    // the dominant cross-phase CARRY-IN that pushed stitch over budget (100K:
    // ~15 GiB carried into stitch; 661k would be far worse). Releasing here
    // hands that RAM back before stitch starts. After release, intern() must
    // NOT be called again (asserted); size()/total_*()/finalize() still work.
    void release_index() {
        if (m_released) return;
        m_released_count = m_classes.size();
        decltype(m_index){}.swap(m_index);     // free the dedup hashset
        std::vector<hash_pair>{}.swap(m_classes);  // free the per-class hashes
        m_released = true;
    }

    // Approx RAM held by the dedup structures: m_classes (16 B/class) + the
    // m_index hashset (ankerl flat backing ~ 1.6 * 9 B/entry). This is what
    // release_index() frees. Reported at phase boundaries to attribute budget.
    uint64_t resident_bytes() const {
        if (m_released) return 0;
        const uint64_t n = m_classes.size();
        const uint64_t classes_bytes = n * sizeof(hash_pair);             // 16 B/class
        const uint64_t index_bytes = (uint64_t)((double)n * 1.6 * 9.0);   // hashset est
        return classes_bytes + index_bytes;
    }

    // Finalize the on-disk file: flush trailing partial word, build &
    // serialize the EF over per-class bit-offsets streamed back from
    // the sidecar, then fseek back to overwrite the placeholder header.
    void finalize() {
        if (m_finalized)
            throw std::runtime_error("streaming_color_set_dict::finalize called twice");
        m_finalized = true;

        const uint64_t total_bit_count = total_bits();
        const uint64_t total_word_count = (total_bit_count + 63) / 64;

        // 1) Flush the trailing partial word (zero-padded to 64 bits).
        if (m_bvb.num_bits() > 0) {
            auto const& words = m_bvb.data();
            size_t want = (m_bvb.num_bits() + 63) / 64;
            if (want > words.size()) want = words.size();
            if (want > 0) {
                size_t bytes = want * sizeof(uint64_t);
                if (std::fwrite(words.data(), 1, bytes, m_file) != bytes)
                    throw std::runtime_error("short write of trailing words to " + m_output_path);
            }
            m_bvb.clear();
        }

        // 2) Append the sentinel total_bit_count to the offsets
        //    sidecar so EF::encode reads num_classes + 1 values.
        if (std::fwrite(&total_bit_count, sizeof(total_bit_count), 1, m_offsets_file) != 1)
            throw std::runtime_error("short write of EF sentinel to " + m_offsets_path);
        std::fflush(m_offsets_file);
        std::rewind(m_offsets_file);

        // 3) Build EF via a file-backed input iterator. EF::encode
        //    walks the sequence exactly once when universe is given,
        //    so a single-pass iterator suffices.
        const uint64_t n = (uint64_t)size() + 1;
        bits::elias_fano<false, false> ef;
        offset_file_iterator begin(m_offsets_file);
        ++begin;  // load the first value into operator*
        ef.encode(begin, n, total_bit_count);

        std::fclose(m_offsets_file);
        m_offsets_file = nullptr;
        std::error_code ec;
        std::filesystem::remove(m_offsets_path, ec);
        m_offsets_path.clear();

        // 4) Serialize EF to the file's tail via essentials.
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

        // 5) Update the placeholder header with the now-known totals.
        std::fflush(m_file);
        if (std::fseek(m_file, 0, SEEK_SET) != 0)
            throw std::runtime_error("fseek to header failed on " + m_output_path);
        write_pod_(m_num_colors);
        write_pod_(m_sparse_threshold);
        write_pod_(m_dense_threshold);
        write_pod_((uint64_t)size());
        write_pod_(total_bit_count);
        write_pod_(total_word_count);

        std::fflush(m_file);
        std::fclose(m_file);
        m_file = nullptr;

        std::cout << "  num_color_sets = " << size() << "\n";
        std::cout << "  num_total_integers = " << m_total_integers << "\n";
        std::cout << "  total bits for ints  = " << total_bit_count << "\n";
        std::cout << "  total bits for offs  = " << 8 * ef_bytes.size() << "\n";
    }

private:
    static constexpr size_t HEADER_BYTES = 4 + 4 + 4 + 8 + 8 + 8;

    struct hash_pair {
        uint64_t primary;
        uint64_t secondary;
    };

    struct hasher {
        std::vector<hash_pair> const* classes;
        using is_transparent = void;
        using is_avalanching = void;
        size_t operator()(uint64_t id) const noexcept { return (*classes)[id].primary; }
        size_t operator()(hash_pair const& h) const noexcept { return h.primary; }
    };

    struct key_eq {
        std::vector<hash_pair> const* classes;
        using is_transparent = void;
        bool operator()(uint64_t a, uint64_t b) const noexcept {
            auto const& ea = (*classes)[a];
            auto const& eb = (*classes)[b];
            return ea.primary == eb.primary && ea.secondary == eb.secondary;
        }
        bool operator()(uint64_t a, hash_pair const& h) const noexcept {
            auto const& ea = (*classes)[a];
            return ea.primary == h.primary && ea.secondary == h.secondary;
        }
        bool operator()(hash_pair const& h, uint64_t a) const noexcept { return (*this)(a, h); }
    };

    // Forward-input iterator over a sequence of u64s on disk. Used at
    // finalize to feed bits::elias_fano::encode without ever
    // materialising the offsets array in RAM. Single-pass: ++ reads
    // the next u64 from the file; * returns the most recently read
    // value. No equality / sentinel needed -- EF::encode iterates a
    // known number of items.
    struct offset_file_iterator {
        using iterator_category = std::input_iterator_tag;
        using value_type = uint64_t;
        using difference_type = std::ptrdiff_t;
        using pointer = uint64_t const*;
        using reference = uint64_t const&;

        std::FILE* f;
        uint64_t cur = 0;

        explicit offset_file_iterator(std::FILE* file) : f(file) {}

        uint64_t operator*() const { return cur; }
        offset_file_iterator& operator++() {
            if (std::fread(&cur, sizeof(cur), 1, f) != 1) cur = 0;
            return *this;
        }
        offset_file_iterator operator++(int) {
            offset_file_iterator prev = *this;
            ++(*this);
            return prev;
        }
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
    bool m_released = false;          // release_index() called (interning done)
    uint64_t m_released_count = 0;    // class count stashed before release
    std::vector<hash_pair> m_classes;  // 16 B per class (dedup hashes only)
    // bucket_type::big -> 64-bit internal value-index; the default
    // standard bucket caps at 2^32 entries, which overflows past
    // 4.29e9 color classes.
    ankerl::unordered_dense::set<uint64_t, hasher, key_eq, std::allocator<uint64_t>,
                                 ankerl::unordered_dense::bucket_type::big>
        m_index;
    uint64_t m_total_integers = 0;

    std::FILE* m_file = nullptr;
    std::string m_output_path;
    std::FILE* m_offsets_file = nullptr;
    std::string m_offsets_path;
    bool m_finalized = false;
};

}  // namespace cdbg
