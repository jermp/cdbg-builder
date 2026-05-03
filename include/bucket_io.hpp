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
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <string_view>
#include <sys/stat.h>
#include <vector>
#include <zlib.h>

#include <unordered_dense/unordered_dense.h>

#include "super_kmer.hpp"
#include "util.hpp"

namespace cdgb {

// ---- Per-bucket in-memory compactor + gzip-streaming output -----------------

// Default per-bucket hashmap memory budget (estimated). With 1024 buckets
// at this budget the global ingest peak from the compactors is ~256 MiB.
inline constexpr size_t DEFAULT_COMPACTOR_SPILL_BYTES = 256 * 1024;
// Default zlib output buffer per bucket. Each gzFile also carries an
// internal deflate state of similar size, so the per-bucket zlib RAM
// footprint is roughly 2 * gzbuffer.
inline constexpr size_t DEFAULT_COMPACTOR_GZBUFFER_BYTES = 256 * 1024;

class bucket_compactor {
public:
    // `under_pressure` (optional) is a writer-owned atomic flag set by the
    // RSS watcher when the process approaches the --max-ram budget. When
    // set, insert_batch spills at the end of every batch (ignoring the
    // m_spill_bytes threshold), and try_spill() lets the watcher itself
    // force-flush quiet buckets.
    bucket_compactor(std::string path, size_t spill_bytes, size_t gzbuffer_bytes,
                     std::atomic<bool> const* under_pressure = nullptr)
        : m_path(std::move(path)), m_spill_bytes(spill_bytes), m_under_pressure(under_pressure) {
        m_file = gzopen(m_path.c_str(), "wb1");
        if (!m_file)
            throw std::runtime_error("cannot open bucket file: " + m_path + ": " +
                                     std::strerror(errno));
        gzbuffer(m_file, (unsigned)gzbuffer_bytes);
    }

    ~bucket_compactor() { close(); }

    bucket_compactor(bucket_compactor const&) = delete;
    bucket_compactor& operator=(bucket_compactor const&) = delete;

    std::string const& path() const { return m_path; }
    uint64_t total_uncompressed_bytes() const {
        return m_total_uncompressed.load(std::memory_order_relaxed);
    }

    // Insert a batch of records (parsed). `bases_storage` holds the 2-bit
    // values (one per byte) for every record; each record points into it
    // via [bases_off, bases_off + bases_len). Caller may clear/reset its
    // buffers after this returns.
    struct pending_record {
        uint32_t color;
        uint32_t bases_off;
        uint32_t bases_len;
        uint8_t flags;
    };

    void insert_batch(std::vector<pending_record> const& recs,
                      std::vector<uint8_t> const& bases_storage) {
        if (recs.empty()) return;
        auto& prof = bucket_prof();
        auto t_lock = bucket_write_prof::clock::now();
        std::lock_guard<std::mutex> lk(m_mu);
        prof.ns_lock_wait.fetch_add(bucket_write_prof::since(t_lock), std::memory_order_relaxed);
        auto t_map = bucket_write_prof::clock::now();
        for (auto const& r : recs) {
            std::string_view key((char const*)bases_storage.data() + r.bases_off, r.bases_len);
            // Transparent find avoids allocating a std::string on the hot
            // (already-seen) path. The std::string is constructed only on
            // a miss, when we actually have to insert into the map.
            auto it = m_dedup.find(key);
            if (it == m_dedup.end()) {
                entry e;
                // Within a single batch from one input thread, a fresh
                // entry sees only one set of flags; subsequent merges
                // across batches/colors AND-shrink the begin/end bits.
                e.flags = r.flags;
                e.colors.push_back(r.color);
                m_bytes += key.size() + sizeof(uint32_t);
                m_dedup.emplace(std::string(key), std::move(e));
                prof.n_inserts.fetch_add(1, std::memory_order_relaxed);
            } else {
                entry& e = it->second;
                e.flags &= r.flags;  // AND across contributors
                // Append unsorted; sort+unique runs once at spill time.
                // For high-redundancy inputs (a popular super-k-mer hit by
                // every input file) this turns the per-color cost from
                // O(prev_color_count) into O(1).
                e.colors.push_back(r.color);
                m_bytes += sizeof(uint32_t);
            }
        }
        prof.ns_hashmap.fetch_add(bucket_write_prof::since(t_map), std::memory_order_relaxed);
        prof.n_records.fetch_add(recs.size(), std::memory_order_relaxed);
        // Spill at the configured threshold only -- don't shred dedup
        // state by spilling on every batch under pressure. The
        // bucket_writer's RSS watcher provides backpressure
        // independently by sweeping all compactors whenever live RSS
        // crosses the budget threshold; per-batch backpressure on top
        // of that just produces tiny spills (worse dedup, more disk
        // records, more rsids per k-mer in bucket-process). On the
        // 4546-genome workload this kept spill count at the natural
        // ~500K-700K range instead of the 4M+ explosion seen with
        // per-batch pressure spilling.
        if (m_bytes >= m_spill_bytes) {
            auto t_sp = bucket_write_prof::clock::now();
            spill_locked();
            prof.ns_spill.fetch_add(bucket_write_prof::since(t_sp), std::memory_order_relaxed);
            prof.n_spills.fetch_add(1, std::memory_order_relaxed);
        }
    }

    // Force-spill any pending entries, no-op if already empty. Called by
    // the writer's RSS watcher to drain quiet buckets that aren't seeing
    // ingest traffic but still hold accumulated state.
    void try_spill() {
        std::lock_guard<std::mutex> lk(m_mu);
        if (m_dedup.empty()) return;
        auto& prof = bucket_prof();
        auto t_sp = bucket_write_prof::clock::now();
        spill_locked();
        prof.ns_spill.fetch_add(bucket_write_prof::since(t_sp), std::memory_order_relaxed);
        prof.n_spills.fetch_add(1, std::memory_order_relaxed);
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
    struct entry {
        std::vector<uint32_t> colors;  // ascending, deduped
        uint8_t flags = 0;
    };

    void spill_locked() {
        if (m_dedup.empty()) return;
        std::vector<uint8_t> rec_buf;
        rec_buf.reserve(256);
        uint64_t spilled = 0;
        for (auto& kv : m_dedup) {
            std::string const& key = kv.first;
            entry& e = kv.second;
            // Sort+unique once per spill; the insert-time path is just
            // push_back, so duplicates are common when the same color
            // recurs across batches.
            std::sort(e.colors.begin(), e.colors.end());
            e.colors.erase(std::unique(e.colors.begin(), e.colors.end()), e.colors.end());
            rec_buf.clear();
            write_super_kmer(e.flags, e.colors.data(), (uint32_t)e.colors.size(),
                             (uint8_t const*)key.data(), (uint32_t)key.size(), rec_buf);
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
    struct string_hash {
        using is_transparent = void;
        using is_avalanching = void;
        size_t operator()(std::string_view sv) const noexcept {
            return ankerl::unordered_dense::hash<std::string_view>{}(sv);
        }
    };
    struct string_eq {
        using is_transparent = void;
        bool operator()(std::string_view a, std::string_view b) const noexcept { return a == b; }
    };

    std::string m_path;
    size_t m_spill_bytes;
    std::atomic<bool> const* m_under_pressure = nullptr;
    std::mutex m_mu;
    gzFile m_file = nullptr;
    ankerl::unordered_dense::map<std::string, entry, string_hash, string_eq> m_dedup;
    size_t m_bytes = 0;
    std::atomic<uint64_t> m_total_uncompressed{0};
};

// ---- bucket_writer: owns one bucket_compactor per bucket ----------------------

class bucket_writer {
public:
    bucket_writer(std::string const& dir, uint32_t num_buckets, size_t flush_bases = 64 * 1024,
                  size_t spill_bytes = DEFAULT_COMPACTOR_SPILL_BYTES,
                  size_t gzbuffer_bytes = DEFAULT_COMPACTOR_GZBUFFER_BYTES)
        : m_dir(dir), m_num_buckets(num_buckets), m_flush_bases(flush_bases) {
        std::filesystem::create_directories(m_dir);
        m_compactors.reserve(num_buckets);
        for (uint32_t b = 0; b < num_buckets; ++b) {
            m_compactors.emplace_back(std::make_unique<bucket_compactor>(
                bucket_path(b), spill_bytes, gzbuffer_bytes, &m_under_pressure));
        }
    }

    ~bucket_writer() {
        stop_rss_watcher();
        close();
    }

    bucket_writer(bucket_writer const&) = delete;
    bucket_writer& operator=(bucket_writer const&) = delete;

    uint32_t num_buckets() const { return m_num_buckets; }
    // Per-thread, per-bucket buffer threshold in *bases* (2-bit values).
    size_t flush_bases() const { return m_flush_bases; }

    std::string bucket_path(uint32_t b) const {
        return m_dir + "/bucket_" + std::to_string(b) + ".bin";
    }

    void flush(uint32_t b, std::vector<bucket_compactor::pending_record>& recs,
               std::vector<uint8_t>& bases_buf) {
        if (recs.empty()) return;
        auto& prof = bucket_prof();
        auto t = bucket_write_prof::clock::now();
        m_compactors[b]->insert_batch(recs, bases_buf);
        prof.ns_flush.fetch_add(bucket_write_prof::since(t), std::memory_order_relaxed);
        prof.n_flushes.fetch_add(1, std::memory_order_relaxed);
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

    // ---- RSS pressure watcher --------------------------------------------
    //
    // Background thread that polls live RSS (via /proc/self/status:VmRSS
    // when available, else getrusage's monotonic peak as fallback). On
    // OFF -> ON transition (RSS crosses high_threshold_bytes), the
    // watcher sets m_under_pressure and sweeps all compactors once,
    // force-spilling any with pending state. On ON -> OFF transition
    // (RSS drops below low_threshold_bytes), the watcher clears the
    // flag and ingest threads return to their normal m_bytes-based
    // spill cadence.
    //
    // Hysteresis prevents the death-spiral seen with sticky pressure:
    // without a release condition, every ingest call would spill for
    // the rest of the phase, exploding spill counts and bucket file
    // sizes (3.4M spills / 6.8 GB written observed on the 4546-genome
    // workload). With hysteresis the steady state is "spill enough to
    // stay under the cap, then resume normal cadence."
    //
    // Fallback path: if current_rss_bytes() returns 0 (no /proc), the
    // watcher uses process_peak_rss_bytes(). That is monotonic so
    // there's no release; pressure stays sticky once tripped. Less
    // efficient but still correct.
    //
    // Per-thread bucket buffers and zlib state are NOT spillable --
    // they form a structural floor of approximately
    //   T*B*flush_bases*overhead + B*(gzbuffer + ~384 KiB)
    // bytes that no amount of pressure response can reduce. The auto-
    // tune in builder picks flush_bases / gzbuffer / num_buckets to
    // keep that floor under the budget share for bucket-write.
    void start_rss_watcher(uint64_t high_threshold_bytes, uint64_t low_threshold_bytes,
                           std::chrono::milliseconds interval = std::chrono::milliseconds(100)) {
        if (m_watcher_running.exchange(true)) return;  // already started
        m_high_threshold_bytes = high_threshold_bytes;
        m_low_threshold_bytes = low_threshold_bytes;
        m_watcher_thread = std::thread([this, interval] { watcher_run(interval); });
    }

    void stop_rss_watcher() {
        if (!m_watcher_running.exchange(false)) return;
        if (m_watcher_thread.joinable()) m_watcher_thread.join();
    }

    bool under_pressure() const {
        return m_under_pressure.load(std::memory_order_relaxed);
    }

    // Highest live RSS observed by the watcher during bucket-write.
    // Useful for verifying that the cap actually held (vs the lifetime
    // peak from getrusage, which can include earlier spikes).
    uint64_t observed_rss_high() const {
        return m_observed_rss_high.load(std::memory_order_relaxed);
    }

    bool pressure_was_engaged() const {
        return m_pressure_was_engaged.load(std::memory_order_relaxed);
    }

private:
    void watcher_run(std::chrono::milliseconds interval) {
        while (m_watcher_running.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(interval);
            if (!m_watcher_running.load(std::memory_order_relaxed)) break;
            uint64_t rss = current_rss_bytes();
            if (rss == 0) {
                // /proc unavailable: fall back to lifetime peak. Sweep
                // semantics still work with peak (it's monotonic, so
                // once we're over we keep sweeping every poll), just
                // less precise.
                rss = process_peak_rss_bytes();
                if (rss == 0) continue;  // can't enforce on this platform
            }
            uint64_t prev_high = m_observed_rss_high.load(std::memory_order_relaxed);
            while (rss > prev_high &&
                   !m_observed_rss_high.compare_exchange_weak(prev_high, rss,
                                                              std::memory_order_relaxed)) {}
            // Hysteresis on the pressure flag (informational; useful for
            // the post-phase log line). Sweep semantics are independent:
            // we sweep whenever RSS is at/above HIGH, regardless of the
            // flag's current state, because if pressure is already ON
            // and RSS climbs *back* to HIGH between sweeps we still
            // need to spill. Otherwise pressure could stay ON, RSS
            // climbs unchecked past HIGH, and the cap is lost.
            if (rss >= m_high_threshold_bytes) {
                m_under_pressure.store(true, std::memory_order_relaxed);
                m_pressure_was_engaged.store(true, std::memory_order_relaxed);
                for (auto& c : m_compactors) {
                    if (!m_watcher_running.load(std::memory_order_relaxed)) break;
                    c->try_spill();
                }
            } else if (rss <= m_low_threshold_bytes) {
                m_under_pressure.store(false, std::memory_order_relaxed);
            }
        }
    }

    std::string m_dir;
    uint32_t m_num_buckets;
    size_t m_flush_bases;
    std::vector<std::unique_ptr<bucket_compactor>> m_compactors;
    std::atomic<uint64_t> m_total_compressed{0};
    std::atomic<uint64_t> m_total_uncompressed{0};
    bool m_closed = false;
    // Pressure flag set by the RSS watcher; consulted by every
    // bucket_compactor::insert_batch via the pointer we hand them at
    // construction.
    std::atomic<bool> m_under_pressure{false};
    std::atomic<bool> m_pressure_was_engaged{false};
    std::atomic<uint64_t> m_observed_rss_high{0};
    std::atomic<bool> m_watcher_running{false};
    std::thread m_watcher_thread;
    uint64_t m_high_threshold_bytes = 0;
    uint64_t m_low_threshold_bytes = 0;
};

// ---- Per-thread batching sidecar -------------------------------------------
//
// Each ingest worker owns one of these. For each bucket it accumulates a
// flat list of pending records plus a parallel 2-bit-value buffer for the
// bases. When a bucket's bases-buffer crosses `flush_bases`, the thread
// flushes that bucket's batch into the writer (which delegates to the
// per-bucket compactor under that bucket's mutex).

struct per_thread_bucket_buffers {
    std::vector<std::vector<bucket_compactor::pending_record>> recs;
    std::vector<std::vector<uint8_t>> bases;
    bucket_writer* sink = nullptr;

    explicit per_thread_bucket_buffers(bucket_writer& w)
        : recs(w.num_buckets()), bases(w.num_buckets()), sink(&w) {}

    void append(uint32_t b, uint8_t flags, uint32_t color, uint8_t const* sk_bases, uint32_t len) {
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

// ---- bucket_reader: streams compacted records from disk ----------------------

class bucket_reader {
public:
    explicit bucket_reader(std::string const& path) {
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
