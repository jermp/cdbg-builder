#pragma once

#include <fstream>
#include <vector>
#include <algorithm>
#include "bit_vector.hpp"
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

    // Exposes a reference to the native iterator if external querying is needed
    bits::bit_vector::iterator& iterator() { return m_it; }

    bool eof() const { return m_stream.eof() && (m_it.position() >= m_bv.num_bits()); }

private:
    void init() {
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

struct unitigs_color_set {
    uint64_t unitig_start{};
    uint64_t num_unitigs{};
    uint64_t color_set_id{};
    std::vector<uint32_t> color_set;
};

class unitigs_color_set_stream {
    streaming_bit_vector m_u2c_reader;
    streaming_bit_vector m_color_reader;
    uint32_t num_colors{};
    uint32_t sparse_threshold{};
    uint32_t very_dense_threshold{};
    uint64_t num_color_sets{};
    uint64_t m_safe_margin_bits{};

    uint64_t m_start_unitig = 0;
    uint64_t m_parsed_sets = 0;

    uint64_t max_queue_size_;
    std::queue<unitigs_color_set> queue_;

    std::mutex mutex_;
    std::condition_variable cv_not_full_;
    std::condition_variable cv_not_empty_;

    std::atomic<bool> stop_requested_{false};
    std::thread worker_thread_;

    void producer_loop() {
        while (!stop_requested_.load() && m_parsed_sets < num_color_sets) {
            auto ucs = parse_next();
            std::unique_lock lock(mutex_);

            cv_not_full_.wait(
                lock, [this] { return queue_.size() < max_queue_size_ || stop_requested_.load(); });

            queue_.push(std::move(ucs));
            cv_not_empty_.notify_one();
        }

        cv_not_empty_.notify_all();
    }

    unitigs_color_set parse_next() {
        const uint64_t begin = m_u2c_reader.position();
        const uint64_t num_unitigs = m_u2c_reader.skip_zeros() + 1;
        assert(m_u2c_reader.position() - begin == num_unitigs);

        unitigs_color_set rec;
        rec.unitig_start = begin;
        rec.num_unitigs = num_unitigs;
        rec.color_set_id = m_parsed_sets;
        m_start_unitig = rec.unitig_start;

        rec.color_set = decode_next_color_set();
        ++m_parsed_sets;

        return rec;
    }

public:
    explicit unitigs_color_set_stream(const std::string& base_filename,
                                      const uint64_t max_queue_size,
                                      const size_t chunk_size_bytes = 64 * 1024 * 1024)
        : m_u2c_reader(std::ifstream(u2c_filename(base_filename), std::ios::binary), 0,
                       chunk_size_bytes)
        , m_color_reader(std::ifstream(cs_filename(base_filename), std::ios::binary), 20,
                         chunk_size_bytes)
        , max_queue_size_(max_queue_size) {
        auto cs_file = std::ifstream(cs_filename(base_filename), std::ios::binary);
        cs_file.read(reinterpret_cast<char*>(&num_colors), sizeof(num_colors));
        cs_file.read(reinterpret_cast<char*>(&sparse_threshold), sizeof(sparse_threshold));
        cs_file.read(reinterpret_cast<char*>(&very_dense_threshold), sizeof(very_dense_threshold));
        cs_file.read(reinterpret_cast<char*>(&num_color_sets), sizeof(num_color_sets));

        // Each color_set should take at most num_colors bit. Times 2 to be sure.
        m_safe_margin_bits = num_colors * 2;
    }

    void start() {
        m_start_unitig = 0;
        m_parsed_sets = 0;
        worker_thread_ = std::thread([this] { producer_loop(); });
    }

    std::optional<unitigs_color_set> get() {
        std::unique_lock lock(mutex_);

        cv_not_empty_.wait(lock, [this] {
            return !queue_.empty() || m_parsed_sets == num_color_sets || stop_requested_.load();
        });

        if (stop_requested_.load() || (queue_.empty() && m_parsed_sets == num_color_sets)) {
            return std::nullopt;
        }
        auto ucs = std::move(queue_.front());
        queue_.pop();
        cv_not_full_.notify_one();

        return ucs;
    }

    void request_stop_and_join() {
        stop_requested_.store(true);

        cv_not_full_.notify_all();
        cv_not_empty_.notify_all();

        if (worker_thread_.joinable()) { worker_thread_.join(); }
    }

    ~unitigs_color_set_stream() { request_stop_and_join(); }

private:
    std::vector<uint32_t> decode_next_color_set() {
        m_color_reader.ensure_bits_available(m_safe_margin_bits);

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
            const uint32_t init_pos = m_color_reader.position();
            for (uint32_t i = 0; i < cs_size; ++i) {
                color_set.push_back(m_color_reader.next() - init_pos);
            }
            m_color_reader.skip_to(init_pos + num_colors);
        } else {
            const uint64_t size = num_colors - cs_size;
            uint32_t curr = 0;
            uint32_t next = 0;
            for (uint64_t i = 0; i < size; ++i) {
                next += m_color_reader.read_delta();
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