#pragma once

// Disk-backed bucket I/O for the minimizer-bucketed ingest path, with
// inline compaction.
//
// Each bucket on disk holds a stream of *compacted* super-k-mer records
// (see super_kmer.hpp): a super-k-mer that occurs in many input files is
// stored once with the union of its colors. Compaction happens online in
// per-bucket hashmaps; ingest threads flush their per-thread buffers into
// the bucket's hashmap under the bucket mutex, and the hashmap spills its
// contents to disk (gzip level 1) whenever it crosses a memory budget. The
// bucket walker correctly merges color sets across multiple spilled
// records for the same super-k-mer, so spilling is purely a memory bound.
//
// On dense pangenome inputs (e.g. many closely-related bacterial genomes)
// most super-k-mers are shared across colors, so compaction is the main
// disk-write reduction over a naive one-record-per-occurrence layout.

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sys/stat.h>
#include <vector>
#include <zlib.h>

#include <unordered_dense/unordered_dense.h>

#include "super_kmer.hpp"

namespace cdgb {

// ---- Per-bucket in-memory compactor + gzip-streaming output -----------------

// Default per-bucket hashmap memory budget (estimated). With 1024 buckets
// at this budget the global ingest peak from the compactors is ~256 MiB.
inline constexpr size_t DEFAULT_COMPACTOR_SPILL_BYTES = 256 * 1024;

class BucketCompactor {
public:
    BucketCompactor(std::string path, size_t spill_bytes)
        : m_path(std::move(path)), m_spill_bytes(spill_bytes) {
        m_file = gzopen(m_path.c_str(), "wb1");
        if (!m_file)
            throw std::runtime_error("cannot open bucket file: " + m_path + ": " +
                                     std::strerror(errno));
        gzbuffer(m_file, 256 * 1024);
    }

    ~BucketCompactor() { close(); }

    BucketCompactor(const BucketCompactor&) = delete;
    BucketCompactor& operator=(const BucketCompactor&) = delete;

    const std::string& path() const { return m_path; }
    uint64_t total_uncompressed_bytes() const {
        return m_total_uncompressed.load(std::memory_order_relaxed);
    }

    // Insert a batch of records (parsed). `bases_storage` holds the 2-bit
    // values (one per byte) for every record; each record points into it
    // via [bases_off, bases_off + bases_len). Caller may clear/reset its
    // buffers after this returns.
    struct PendingRecord {
        uint32_t color;
        uint32_t bases_off;
        uint32_t bases_len;
        uint8_t flags;
    };

    void insert_batch(const std::vector<PendingRecord>& recs,
                      const std::vector<uint8_t>& bases_storage) {
        if (recs.empty()) return;
        std::lock_guard<std::mutex> lk(m_mu);
        for (const auto& r : recs) {
            std::string_view key((const char*)bases_storage.data() + r.bases_off, r.bases_len);
            // Transparent find avoids allocating a std::string on the hot
            // (already-seen) path. The std::string is constructed only on
            // a miss, when we actually have to insert into the map.
            auto it = m_dedup.find(key);
            if (it == m_dedup.end()) {
                Entry e;
                // Within a single batch from one input thread, a fresh
                // entry sees only one set of flags; subsequent merges
                // across batches/colors AND-shrink the begin/end bits.
                e.flags = r.flags;
                e.colors.push_back(r.color);
                m_bytes += key.size() + sizeof(uint32_t);
                m_dedup.emplace(std::string(key), std::move(e));
            } else {
                Entry& e = it->second;
                e.flags &= r.flags;  // AND across contributors
                // Append unsorted; sort+unique runs once at spill time.
                // For high-redundancy inputs (a popular super-k-mer hit by
                // every input file) this turns the per-color cost from
                // O(prev_color_count) into O(1).
                e.colors.push_back(r.color);
                m_bytes += sizeof(uint32_t);
            }
        }
        if (m_bytes >= m_spill_bytes) spill_locked();
    }

    void close() {
        std::lock_guard<std::mutex> lk(m_mu);
        if (m_file) {
            spill_locked();
            gzclose(m_file);
            m_file = nullptr;
        }
    }

private:
    struct Entry {
        std::vector<uint32_t> colors;  // ascending, deduped
        uint8_t flags = 0;
    };

    void spill_locked() {
        if (m_dedup.empty()) return;
        std::vector<uint8_t> rec_buf;
        rec_buf.reserve(256);
        uint64_t spilled = 0;
        for (auto& kv : m_dedup) {
            const std::string& key = kv.first;
            Entry& e = kv.second;
            // Sort+unique once per spill; the insert-time path is just
            // push_back, so duplicates are common when the same color
            // recurs across batches.
            std::sort(e.colors.begin(), e.colors.end());
            e.colors.erase(std::unique(e.colors.begin(), e.colors.end()), e.colors.end());
            rec_buf.clear();
            write_super_kmer(e.flags, e.colors.data(), (uint32_t)e.colors.size(),
                             (const uint8_t*)key.data(), (uint32_t)key.size(), rec_buf);
            int n = gzwrite(m_file, rec_buf.data(), (unsigned)rec_buf.size());
            if (n <= 0 || (size_t)n != rec_buf.size())
                throw std::runtime_error("short write to " + m_path);
            spilled += rec_buf.size();
        }
        m_total_uncompressed.fetch_add(spilled, std::memory_order_relaxed);
        m_dedup.clear();
        m_bytes = 0;
    }

    // Transparent hash/eq so find() can take a string_view directly without
    // building a std::string. Saves a heap allocation + copy on every
    // already-seen super-k-mer, which is the common case once the first
    // input file has populated each bucket's hashmap.
    struct StringHash {
        using is_transparent = void;
        using is_avalanching = void;
        size_t operator()(std::string_view sv) const noexcept {
            return ankerl::unordered_dense::hash<std::string_view>{}(sv);
        }
    };
    struct StringEq {
        using is_transparent = void;
        bool operator()(std::string_view a, std::string_view b) const noexcept { return a == b; }
    };

    std::string m_path;
    size_t m_spill_bytes;
    std::mutex m_mu;
    gzFile m_file = nullptr;
    ankerl::unordered_dense::map<std::string, Entry, StringHash, StringEq> m_dedup;
    size_t m_bytes = 0;
    std::atomic<uint64_t> m_total_uncompressed{0};
};

// ---- BucketWriter: owns one BucketCompactor per bucket ----------------------

class BucketWriter {
public:
    BucketWriter(const std::string& dir, uint32_t num_buckets, size_t flush_bases = 64 * 1024,
                 size_t spill_bytes = DEFAULT_COMPACTOR_SPILL_BYTES)
        : m_dir(dir), m_num_buckets(num_buckets), m_flush_bases(flush_bases) {
        std::filesystem::create_directories(m_dir);
        m_compactors.reserve(num_buckets);
        for (uint32_t b = 0; b < num_buckets; ++b) {
            m_compactors.emplace_back(
                std::make_unique<BucketCompactor>(bucket_path(b), spill_bytes));
        }
    }

    ~BucketWriter() { close(); }

    BucketWriter(const BucketWriter&) = delete;
    BucketWriter& operator=(const BucketWriter&) = delete;

    uint32_t num_buckets() const { return m_num_buckets; }
    // Per-thread, per-bucket buffer threshold in *bases* (2-bit values).
    size_t flush_bases() const { return m_flush_bases; }

    std::string bucket_path(uint32_t b) const {
        return m_dir + "/bucket_" + std::to_string(b) + ".bin";
    }

    void flush(uint32_t b, std::vector<BucketCompactor::PendingRecord>& recs,
               std::vector<uint8_t>& bases_buf) {
        if (recs.empty()) return;
        m_compactors[b]->insert_batch(recs, bases_buf);
        recs.clear();
        bases_buf.clear();
    }

    void close() {
        if (m_closed) return;
        for (auto& c : m_compactors) c->close();
        // Stat each file once for the on-disk byte total.
        uint64_t total_compressed = 0;
        uint64_t total_uncompressed = 0;
        for (uint32_t b = 0; b < m_num_buckets; ++b) {
            struct stat st;
            if (::stat(bucket_path(b).c_str(), &st) == 0) total_compressed += (uint64_t)st.st_size;
            total_uncompressed += m_compactors[b]->total_uncompressed_bytes();
        }
        m_total_compressed.store(total_compressed, std::memory_order_relaxed);
        m_total_uncompressed.store(total_uncompressed, std::memory_order_relaxed);
        m_closed = true;
    }

    uint64_t total_bytes() const { return m_total_compressed.load(std::memory_order_relaxed); }
    uint64_t total_uncompressed_bytes() const {
        return m_total_uncompressed.load(std::memory_order_relaxed);
    }

private:
    std::string m_dir;
    uint32_t m_num_buckets;
    size_t m_flush_bases;
    std::vector<std::unique_ptr<BucketCompactor>> m_compactors;
    std::atomic<uint64_t> m_total_compressed{0};
    std::atomic<uint64_t> m_total_uncompressed{0};
    bool m_closed = false;
};

// ---- Per-thread batching sidecar -------------------------------------------
//
// Each ingest worker owns one of these. For each bucket it accumulates a
// flat list of pending records plus a parallel 2-bit-value buffer for the
// bases. When a bucket's bases-buffer crosses `flush_bases`, the thread
// flushes that bucket's batch into the writer (which delegates to the
// per-bucket compactor under that bucket's mutex).

struct PerThreadBucketBuffers {
    std::vector<std::vector<BucketCompactor::PendingRecord>> recs;
    std::vector<std::vector<uint8_t>> bases;
    BucketWriter* sink = nullptr;

    explicit PerThreadBucketBuffers(BucketWriter& w)
        : recs(w.num_buckets()), bases(w.num_buckets()), sink(&w) {}

    void append(uint32_t b, uint8_t flags, uint32_t color, const uint8_t* sk_bases, uint32_t len) {
        auto& bbuf = bases[b];
        uint32_t off = (uint32_t)bbuf.size();
        bbuf.insert(bbuf.end(), sk_bases, sk_bases + len);
        recs[b].push_back({color, off, len, (uint8_t)(flags & 0x3u)});
        if (bbuf.size() >= sink->flush_bases()) sink->flush(b, recs[b], bbuf);
    }

    void flush_all() {
        for (uint32_t b = 0; b < (uint32_t)recs.size(); ++b) {
            if (!recs[b].empty()) sink->flush(b, recs[b], bases[b]);
        }
    }
};

// ---- BucketReader: streams compacted records from disk ----------------------

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

    // Iterate records in order. Returns false when no more records remain.
    bool next(uint8_t& flags, std::vector<uint32_t>& colors, std::vector<uint8_t>& bases) {
        if (m_pos >= m_buf.size()) return false;
        size_t consumed =
            read_super_kmer(m_buf.data() + m_pos, m_buf.size() - m_pos, flags, colors, bases);
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
