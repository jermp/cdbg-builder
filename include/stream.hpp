#pragma once

#include <fstream>
#include <vector>
#include <algorithm>
#include "bit_vector.hpp"
#include "elias_fano.hpp"
#include "integer_codes.hpp"
#include "util.hpp"

namespace cdbg {
class streaming_bit_vector {
public:
    // chunk_size_bytes defaults to 64MB for high-throughput block reads
    explicit streaming_bit_vector(std::ifstream&& stream, const std::streamoff skip = 0,
                                  size_t chunk_size_bytes = 64 * 1024 * 1024)
        : m_stream(std::move(stream)), m_chunk_size_words(chunk_size_bytes / sizeof(uint64_t)) {
        m_stream.seekg(skip);
        init();
    }

    /**
     * Call this before starting to decode a new color set.
     * If the remaining bits in the current bit_vector chunk fall below the safe margin,
     * it shifts unconsumed bits to the front and appends new words from the stream.
     */
    void ensure_bits_available(uint64_t safe_margin_bits = 16 * 1024 * 1024) {
        const uint64_t current_pos = m_it.position();
        const uint64_t total_bits = m_bv.num_bits();

        if (total_bits - current_pos < safe_margin_bits && !m_stream.eof()) {
            uint64_t remaining_bits = total_bits - current_pos;
            bits::bit_vector::builder builder;

            // DRY: Reuse the iterator's native take() to seamlessly transfer
            // remaining unconsumed bits to the front of the next chunk.
            while (remaining_bits > 0) {
                const uint64_t grab = std::min<uint64_t>(64, remaining_bits);
                const uint64_t val = m_it.take(grab);
                builder.append_bits(val, grab);
                remaining_bits -= grab;
            }

            // Calculate how many free words are available in the chunk to fill up
            const size_t words_to_read = m_chunk_size_words - builder.data().size();
            if (words_to_read > 0) {
                std::vector<uint64_t> read_buf(words_to_read, 0);
                m_stream.read(reinterpret_cast<char*>(read_buf.data()),
                              words_to_read * sizeof(uint64_t));
                const size_t words_read = m_stream.gcount() / sizeof(uint64_t);

                for (size_t i = 0; i < words_read; ++i) { builder.append_bits(read_buf[i], 64); }
            }

            // Build the new chunk and point our iterator back to the beginning (0)
            builder.build(m_bv);
            m_it = m_bv.begin();
        }
    }

    uint64_t read_delta() {
        ensure_bits_available(64);
        m_it.skip_to(m_it.position());
        return bits::util::read_delta(m_it);
    }

    uint64_t take(const uint64_t l) {
        ensure_bits_available(l);
        return m_it.take(l);
    }

    uint64_t skip_zeros() { return m_it.skip_zeros(); }

    uint64_t next() { return m_it.next(); }
    uint64_t position() const { return m_it.position(); }
    void skip_to(const uint64_t pos) { m_it.skip_to(pos); }

    bits::bit_vector::iterator& iterator() { return m_it; }
    void span(bits::bit_vector& bv, uint64_t len) {
        bits::bit_vector::builder builder;
        while (len >= 64) {
            builder.append_bits(take(64), 64);
            len -= 64;
        }
        builder.append_bits(take(len), len);

        builder.build(bv);
    }

    bool eof() const { return m_stream.eof() && (m_it.position() >= m_bv.num_bits()); }

private:
    void init() {
        bits::bit_vector::builder builder;
        builder.data().resize(m_chunk_size_words);
        m_stream.seekg(16, std::ios::cur);  // skip bv metadata (num_bits, num_words)

        // Highly optimized word-aligned block read from file
        m_stream.read(reinterpret_cast<char*>(builder.data().data()),
                      m_chunk_size_words * sizeof(uint64_t));
        const size_t words_read = m_stream.gcount() / sizeof(uint64_t);

        builder.data().resize(words_read);
        builder.resize(words_read * 64);
        builder.build(m_bv);

        m_it = m_bv.begin();
    }

    std::ifstream m_stream;
    size_t m_chunk_size_words;
    bits::bit_vector m_bv;
    bits::bit_vector::iterator m_it;
};

struct unitigs_color_set {
    uint64_t unitig_start{};
    uint64_t num_unitigs{};
    uint64_t color_set_id{};
    std::vector<uint32_t> color_set;
};

class unitigs_color_set_stream {
    streaming_bit_vector m_u2c_reader;
    streaming_bit_vector m_color_reader;
    const metadata m_metadata;

    uint32_t num_colors{};
    uint32_t sparse_threshold{};
    uint32_t very_dense_threshold{};

    std::atomic<uint64_t> m_parsed_sets = 0;

    bits::elias_fano<false, false> m_offsets;
    std::mutex m_mutex;

public:
    explicit unitigs_color_set_stream(const std::string& base_filename,
                                      const size_t chunk_size_bytes = 64 * 1024 * 1024)
        : m_u2c_reader(std::ifstream(u2c_filename(base_filename), std::ios::binary), 0,
                       chunk_size_bytes)
        , m_color_reader(std::ifstream(cs_filename(base_filename), std::ios::binary), 12,
                         chunk_size_bytes)
        , m_metadata(metadata_filename(base_filename)) {
        auto cs_file = std::ifstream(cs_filename(base_filename), std::ios::binary);
        cs_file.read(reinterpret_cast<char*>(&num_colors), sizeof(num_colors));
        cs_file.read(reinterpret_cast<char*>(&sparse_threshold), sizeof(sparse_threshold));
        cs_file.read(reinterpret_cast<char*>(&very_dense_threshold), sizeof(very_dense_threshold));
        assert(num_colors == m_metadata.num_colors);
        uint64_t bit_vector_num_bits;
        uint64_t bit_vector_num_words;
        cs_file.read(reinterpret_cast<char*>(&bit_vector_num_bits), sizeof(bit_vector_num_bits));
        cs_file.read(reinterpret_cast<char*>(&bit_vector_num_words), sizeof(bit_vector_num_words));

        cs_file.seekg(28 + bit_vector_num_words * 8, std::ios::beg);
        essentials::generic_loader loader(cs_file);
        loader.visit(m_offsets);
    }

    std::optional<unitigs_color_set> get() {
        bits::bit_vector encoded_cs;
        uint64_t cs_id = 0;
        uint64_t u2c_begin = 0;
        uint64_t u2c_num_unitigs = 0;
        {
            std::lock_guard guard(m_mutex);
            if (m_parsed_sets >= m_metadata.num_color_sets) { return std::nullopt; }

            cs_id = m_parsed_sets++;
            u2c_begin = m_u2c_reader.position();
            u2c_num_unitigs = m_u2c_reader.skip_zeros() + 1;
            assert(m_u2c_reader.position() - u2c_begin == u2c_num_unitigs);

            const uint64_t num_bits = m_offsets.access(cs_id + 1) - m_offsets.access(cs_id);
            m_color_reader.ensure_bits_available(num_bits);
            m_color_reader.span(encoded_cs, num_bits);
        }

        unitigs_color_set rec;
        rec.unitig_start = u2c_begin;
        rec.num_unitigs = u2c_num_unitigs;
        rec.color_set_id = cs_id;

        rec.color_set = decode_color_set(std::move(encoded_cs));

        return rec;
    }

private:
    std::vector<uint32_t> decode_color_set(bits::bit_vector&& encoded_bv) const {
        std::vector<uint32_t> color_set;
        auto it = encoded_bv.begin();

        const uint64_t cs_size = bits::util::read_delta(it);
        color_set.reserve(cs_size);

        if (cs_size < sparse_threshold) {
            uint32_t prev_val = bits::util::read_delta(it);
            color_set.push_back(static_cast<uint32_t>(prev_val));

            for (uint64_t i = 1; i < cs_size; ++i) {
                auto val = static_cast<uint32_t>(bits::util::read_delta(it) + prev_val + 1);
                color_set.push_back(val);
                prev_val = val;
            }
        } else if (cs_size < very_dense_threshold) {
            const uint32_t init_pos = it.position();
            for (uint32_t i = 0; i < cs_size; ++i) { color_set.push_back(it.next() - init_pos); }
            it.skip_to(init_pos + num_colors);
        } else {
            const uint64_t size = num_colors - cs_size;
            uint32_t curr = 0;
            uint32_t next = 0;
            for (uint64_t i = 0; i < size; ++i) {
                next += bits::util::read_delta(it);
                while (curr < next) {
                    color_set.push_back(curr);
                    curr++;
                }
                curr = next + 1;
                ++next;
            }
            while (curr < num_colors) {
                color_set.push_back(curr);
                curr++;
            }
        }

        assert(color_set.size() == cs_size);
        return color_set;
    }
};

}  // namespace cdbg