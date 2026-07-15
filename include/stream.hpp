#pragma once

#include <fstream>
#include <vector>
#include <algorithm>
#include "bit_vector.hpp"
#include "integer_codes.hpp"
#include "util.hpp"

namespace cdbg {
class StreamBitVectorReader {
public:
    // chunk_size_bytes defaults to 64MB for high-throughput block reads
    explicit StreamBitVectorReader(std::ifstream&& stream, const std::streamoff skip = 0,
                                   size_t chunk_size_bytes = 64 * 1024 * 1024)
        : m_stream(std::move(stream)), m_chunk_size_words(chunk_size_bytes / sizeof(uint64_t)) {
        m_stream.seekg(skip);
        refill_initial();
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

    // Exposes a reference to the native iterator if external querying is needed
    bits::bit_vector::iterator& iterator() { return m_it; }

    bool eof() const { return m_stream.eof() && (m_it.position() >= m_bv.num_bits()); }

private:
    void refill_initial() {
        bits::bit_vector::builder builder;
        builder.data().resize(m_chunk_size_words);
        m_stream.seekg(16, std::ios::cur);

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

struct ColorRecord {
    uint64_t unitig_start;
    uint64_t unitig_end;
    std::vector<uint32_t> color_set;
};

class JointColorStreamer {
    StreamBitVectorReader m_u2c_reader;
    StreamBitVectorReader m_color_reader;
    uint32_t num_colors;
    uint32_t sparse_threshold;
    uint32_t very_dense_threshold;
    uint64_t num_color_sets;

    uint64_t m_current_bit_index = 0;
    uint64_t m_start_unitig = 0;
    uint64_t m_safe_margin_bits = 0;
    uint64_t m_parsed_sets = 0;

public:
    JointColorStreamer(const std::string& base_filename,
                       const size_t chunk_size_bytes = 64 * 1024 * 1024)
        : m_u2c_reader(std::ifstream(u2c_filename(base_filename), std::ios::binary), 0,
                       chunk_size_bytes)
        , m_color_reader(std::ifstream(cs_filename(base_filename), std::ios::binary), 20,
                         chunk_size_bytes) {
        auto cs_file = std::ifstream(cs_filename(base_filename), std::ios::binary);
        cs_file.read(reinterpret_cast<char*>(&num_colors), sizeof(num_colors));
        cs_file.read(reinterpret_cast<char*>(&sparse_threshold), sizeof(sparse_threshold));
        cs_file.read(reinterpret_cast<char*>(&very_dense_threshold), sizeof(very_dense_threshold));
        cs_file.read(reinterpret_cast<char*>(&num_color_sets), sizeof(num_color_sets));

        // Define a safe margin bound for the color stream buffer.
        // The dense/very dense formats can consume up to m_num_colors bits per set.
        // We ensure at least 16MB or enough space for a massive set exists.
        m_safe_margin_bits =
            std::max<uint64_t>(16 * 1024 * 1024, static_cast<uint64_t>(num_colors) * 64);
    }

    // Fills the user-provided buffer container up to max_records.
    // Returns false when the out.u2c stream is completely exhausted.
    bool get_next_batch(std::vector<ColorRecord>& record_buffer, const size_t max_records) {
        record_buffer.clear();

        while (record_buffer.size() < max_records && m_parsed_sets < num_color_sets) {
            // Standard small check to ensure the bitvector stream has data available
            m_u2c_reader.ensure_bits_available(64 * 1024);

            if (m_u2c_reader.take(1)) {
                m_color_reader.ensure_bits_available(m_safe_margin_bits);

                ColorRecord rec;
                rec.unitig_start = m_start_unitig;
                rec.unitig_end = m_current_bit_index;
                m_start_unitig = m_current_bit_index + 1;

                rec.color_set = decode_next_color_set();
                record_buffer.push_back(std::move(rec));
                ++m_parsed_sets;
            }
            ++m_current_bit_index;
        }

        return !record_buffer.empty();
    }

private:
    std::vector<uint32_t> decode_next_color_set() {
        std::vector<uint32_t> color_set;
        const uint64_t cs_size = m_color_reader.read_delta();
        color_set.reserve(cs_size);

        if (cs_size < sparse_threshold) {
            uint32_t prev_val = m_color_reader.read_delta();
            color_set.push_back(static_cast<uint32_t>(prev_val));

            for (uint64_t i = 1; i < cs_size; ++i) {
                auto val = static_cast<uint32_t>(m_color_reader.read_delta() + prev_val + 1);
                color_set.push_back(val);
                prev_val = val;
            }
        } else if (cs_size < very_dense_threshold) {
            uint32_t init_pos = m_color_reader.position();
            for (uint32_t i = 0; i < cs_size; ++i) {
                color_set.push_back(m_color_reader.next() - init_pos);
            }
            m_color_reader.skip_to(init_pos + num_colors);
        } else {
            const uint64_t size = num_colors - cs_size;
            uint32_t curr = 0;
            for (uint64_t i = 0; i < size; ++i) {
                const uint32_t next = m_color_reader.read_delta();
                while (curr < next) {
                    color_set.push_back(curr);
                    curr++;
                }
                curr = next + 1;
            }
            while (curr < num_colors) {
                color_set.push_back(curr);
                curr++;
            }
        }

        return color_set;
    }
};

}  // namespace cdbg