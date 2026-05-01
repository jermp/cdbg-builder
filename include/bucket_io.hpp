#pragma once

// Disk-backed bucket I/O for the minimizer-bucketed ingest path.
//
// The writer side is sharded by minimizer-derived bucket id. Each ingest
// thread owns a per-bucket buffer; when a buffer fills (or the thread
// finishes a file), the thread takes the bucket's mutex and appends its
// buffer to the bucket's gzip stream. Bucket files are streaming gzip
// (level 1) so per-bucket on-disk size matches the GGCAT-style compressed
// format; super-k-mers in the same bucket share their minimizer m-mer,
// which compresses very well.
//
// The reader side decompresses the whole bucket file into memory and yields
// records via super_kmer.hpp. Buckets are bounded in size by 1/B of the
// dataset, so for our target inputs (a few-thousand bacterial genomes)
// per-bucket size stays comfortably in RAM.

#include <atomic>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <mutex>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>
#include <zlib.h>

#include "super_kmer.hpp"

namespace cdgb {

class BucketWriter {
public:
    BucketWriter(const std::string& dir, uint32_t num_buckets, size_t flush_bytes = 64 * 1024)
        : m_dir(dir)
        , m_num_buckets(num_buckets)
        , m_flush_bytes(flush_bytes)
        , m_mus(num_buckets)
        , m_files(num_buckets, nullptr) {
        std::filesystem::create_directories(m_dir);
        for (uint32_t b = 0; b < num_buckets; ++b) {
            std::string p = bucket_path(b);
            // Level 1 (Z_BEST_SPEED): fast streaming gzip; on dense
            // 2-bit-packed DNA records the cost is small but the on-disk
            // savings are large (super-k-mers in the same bucket share
            // their minimizer, so they compress strongly together).
            gzFile f = gzopen(p.c_str(), "wb1");
            if (!f)
                throw std::runtime_error("cannot open bucket file: " + p + ": " +
                                         std::strerror(errno));
            // 256 KiB internal buffer reduces syscall overhead per chunk.
            gzbuffer(f, 256 * 1024);
            m_files[b] = f;
        }
    }

    ~BucketWriter() { close(); }

    BucketWriter(const BucketWriter&) = delete;
    BucketWriter& operator=(const BucketWriter&) = delete;

    uint32_t num_buckets() const { return m_num_buckets; }
    size_t flush_bytes() const { return m_flush_bytes; }
    std::string bucket_path(uint32_t b) const {
        return m_dir + "/bucket_" + std::to_string(b) + ".bin";
    }

    // Append `buf` to bucket `b`'s file under the bucket's mutex. The caller
    // typically passes its per-thread buffer here once it has reached the
    // flush threshold.
    void flush(uint32_t b, std::vector<uint8_t>& buf) {
        if (buf.empty()) return;
        std::lock_guard<std::mutex> lk(m_mus[b]);
        int n = gzwrite(m_files[b], buf.data(), (unsigned)buf.size());
        if (n <= 0 || (size_t)n != buf.size()) {
            throw std::runtime_error("short write to " + bucket_path(b));
        }
        m_total_uncompressed.fetch_add(buf.size(), std::memory_order_relaxed);
        buf.clear();
    }

    void close() {
        for (uint32_t b = 0; b < m_num_buckets; ++b) {
            if (m_files[b]) {
                gzclose(m_files[b]);
                m_files[b] = nullptr;
            }
        }
        if (m_total_compressed.load(std::memory_order_relaxed) == 0) {
            // First close call: stat each file to measure on-disk bytes.
            uint64_t total = 0;
            for (uint32_t b = 0; b < m_num_buckets; ++b) {
                struct stat st;
                if (::stat(bucket_path(b).c_str(), &st) == 0) total += (uint64_t)st.st_size;
            }
            m_total_compressed.store(total, std::memory_order_relaxed);
        }
    }

    // Bytes actually written to disk (compressed).
    uint64_t total_bytes() const { return m_total_compressed.load(std::memory_order_relaxed); }
    // Logical record bytes fed into the compressor (pre-compression).
    uint64_t total_uncompressed_bytes() const {
        return m_total_uncompressed.load(std::memory_order_relaxed);
    }

private:
    std::string m_dir;
    uint32_t m_num_buckets;
    size_t m_flush_bytes;
    std::vector<std::mutex> m_mus;
    std::vector<gzFile> m_files;
    std::atomic<uint64_t> m_total_uncompressed{0};
    std::atomic<uint64_t> m_total_compressed{0};
};

// Per-thread sidecar that batches records destined for each bucket.
struct PerThreadBucketBuffers {
    std::vector<std::vector<uint8_t>> bufs;
    BucketWriter* sink = nullptr;

    explicit PerThreadBucketBuffers(BucketWriter& w) : bufs(w.num_buckets()), sink(&w) {}

    // Append a record to bucket `b`; flush if the buffer crossed the threshold.
    void append(uint32_t b, const uint8_t* rec, size_t rec_len) {
        auto& buf = bufs[b];
        size_t off = buf.size();
        buf.resize(off + rec_len);
        std::memcpy(buf.data() + off, rec, rec_len);
        if (buf.size() >= sink->flush_bytes()) sink->flush(b, buf);
    }

    void flush_all() {
        for (uint32_t b = 0; b < (uint32_t)bufs.size(); ++b) {
            if (!bufs[b].empty()) sink->flush(b, bufs[b]);
        }
    }
};

class BucketReader {
public:
    explicit BucketReader(const std::string& path) {
        gzFile f = gzopen(path.c_str(), "rb");
        if (!f) throw std::runtime_error("cannot open " + path + ": " + std::strerror(errno));
        gzbuffer(f, 256 * 1024);
        constexpr size_t CHUNK = 64 * 1024;
        size_t off = 0;
        for (;;) {
            if (m_buf.size() < off + CHUNK) m_buf.resize(off + CHUNK);
            int n = gzread(f, m_buf.data() + off, (unsigned)CHUNK);
            if (n < 0) {
                gzclose(f);
                throw std::runtime_error("gzread failed on " + path);
            }
            if (n == 0) break;
            off += (size_t)n;
        }
        m_buf.resize(off);
        gzclose(f);
    }

    // Iterate records in order. Returns false when no more records are available.
    bool next(uint8_t& flags, uint32_t& color, std::vector<uint8_t>& bases) {
        if (m_pos >= m_buf.size()) return false;
        size_t consumed =
            read_super_kmer(m_buf.data() + m_pos, m_buf.size() - m_pos, flags, color, bases);
        if (consumed == 0) return false;
        m_pos += consumed;
        return true;
    }

    size_t size_bytes() const { return m_buf.size(); }

private:
    std::vector<uint8_t> m_buf;
    size_t m_pos = 0;
};

}  // namespace cdgb
