#pragma once

// Disk-backed bucket I/O for the minimizer-bucketed ingest path, with
// inline compaction.
//
// Each bucket on disk holds a stream of *compacted* super-k-mer records
// (see super_kmer.hpp): a super-k-mer that occurs in many input files is
// stored once with the union of its colors. Compaction happens online in
// per-bucket hashmaps; ingest threads flush their per-thread buffers into
// the bucket's hashmap under the bucket mutex, and the hashmap spills its
// contents to disk whenever it crosses a memory budget. The bucket walker
// correctly merges color sets across multiple spilled records for the
// same super-k-mer, so spilling is purely a memory bound.
//
// On dense pangenome inputs (e.g. many closely-related bacterial genomes)
// most super-k-mers are shared across colors, so compaction is the main
// disk-write reduction over a naive one-record-per-occurrence layout.
//
// Bucket file format (LZ4 block API + custom per-spill framing):
//
//   repeat:
//     [u32 uncompressed_size]   little-endian; 0 marks end-of-stream
//     [u32 compressed_size]     bytes of LZ4-compressed data that follow
//     [compressed bytes]        LZ4_compress_default output
//
// Each spill produces ONE frame containing all of that spill's
// super-k-mer records concatenated. Batching per spill (instead of
// per-record) is essential: LZ4 has a few-byte per-block framing
// overhead, and per-record compression of ~100-byte inputs ended up
// producing files *larger* than the input (overhead exceeded
// compression savings). Per-spill batching gives proper LZ4
// compression ratio on the joint super-k-mer record stream.
//
// We use the LZ4 block API (LZ4_compress_default / LZ4_decompress_safe)
// rather than the frame API because there's no per-bucket compression
// context: LZ4_compress_default is stateless from the caller's view
// (its internal hash table lives on the call's stack). With thousands
// of buckets that means zero persistent compressor state per bucket --
// the per-bucket footprint is just the batch + output scratch buffers,
// already accounted for in spill_bytes.

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
#include <lz4.h>

#include <unordered_dense/unordered_dense.h>

#include "super_kmer.hpp"
#include "util.hpp"

namespace cdbg {

// ---- Per-bucket in-memory compactor + LZ4-framed disk output ---------------

// Default per-bucket hashmap memory budget (estimated). With 1024 buckets
// at this budget the global ingest peak from the compactors is ~256 MiB.
inline constexpr size_t DEFAULT_COMPACTOR_SPILL_BYTES = 256 * 1024;

class bucket_compactor {
public:
    // try_spill() is the writer's RSS-watcher entry point: when the
    // process is over the -g budget the watcher walks every
    // compactor and force-spills any with pending data. Per-batch
    // ingest itself is unaware of pressure; spill cadence stays at
    // m_bytes >= m_spill_bytes so dedup state per spill remains
    // worthwhile.
    //
    // The compactor doesn't hold any persistent LZ4 state -- all
    // compression happens via stateless LZ4_compress_default calls
    // inside spill_locked(). The per-bucket buffers we *do* keep
    // (m_batch_buf, m_out_buf) grow up to one spill's worth of bytes
    // and are reused across spills. They're already counted as part
    // of the compactor's structural cost in the auto-tune model.
    bucket_compactor(std::string path, size_t spill_bytes)
        : m_path(std::move(path)), m_spill_bytes(spill_bytes) {
        m_file = std::fopen(m_path.c_str(), "wb");
        if (!m_file)
            throw std::runtime_error("cannot open bucket file: " + m_path + ": " +
                                     std::strerror(errno));
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
    // via [bases_off, bases_off + bases_len). The `hash` field is a
    // wyhash of the bases, computed by the producer at append time
    // while the bases are still hot in L1; we feed it directly to the
    // dedup map's find() instead of rehashing the bases bytes here
    // (where they are typically cold). Caller may clear/reset its
    // buffers after this returns.
    struct pending_record {
        uint64_t hash;
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
            // Heterogeneous lookup with a precomputed hash. Our
            // string_hash::operator()(hashed_view) just returns h.hash,
            // so find() doesn't re-hash the bases. The bases bytes are
            // only re-touched on a hash collision (then string_eq does
            // the byte compare).
            hashed_view hv{key, r.hash};
            auto it = m_dedup.find(hv);
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
            // End-of-stream marker: a [u32 0] uncompressed_size with
            // no following bytes. bucket_reader stops on this.
            uint32_t eof = 0;
            std::fwrite(&eof, sizeof(eof), 1, m_file);
            std::fclose(m_file);
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
        m_batch_buf.clear();
        for (auto& kv : m_dedup) {
            std::string const& key = kv.first;
            entry& e = kv.second;
            // Sort+unique once per spill; the insert-time path is just
            // push_back, so duplicates are common when the same color
            // recurs across batches.
            std::sort(e.colors.begin(), e.colors.end());
            e.colors.erase(std::unique(e.colors.begin(), e.colors.end()), e.colors.end());
            // Append the serialized record directly into m_batch_buf.
            // We compress the whole batch with one LZ4 call below;
            // amortising the per-block framing overhead across all of
            // the spill's records is what makes compression actually
            // pay off on inputs this small (~50-200 B per record).
            write_super_kmer(e.flags, e.colors.data(), (uint32_t)e.colors.size(),
                             (uint8_t const*)key.data(), (uint32_t)key.size(), m_batch_buf);
        }
        if (m_batch_buf.empty()) {
            m_dedup.clear();
            m_bytes = 0;
            return;
        }
        // LZ4 block API: stateless, fast mode (acceleration=1). The
        // function expects int parameters, so we cap each frame at
        // ~1 GiB; in practice spills are kilobytes to a few MiB so
        // this never trips.
        if (m_batch_buf.size() > (size_t)LZ4_MAX_INPUT_SIZE)
            throw std::runtime_error("spill batch exceeds LZ4_MAX_INPUT_SIZE on " + m_path);
        int src_size = (int)m_batch_buf.size();
        int bound = LZ4_compressBound(src_size);
        if (bound <= 0)
            throw std::runtime_error("LZ4_compressBound failed on " + m_path);
        if (m_out_buf.size() < (size_t)bound) m_out_buf.resize((size_t)bound);
        int compressed = LZ4_compress_default((char const*)m_batch_buf.data(),
                                              (char*)m_out_buf.data(), src_size,
                                              (int)m_out_buf.size());
        if (compressed <= 0)
            throw std::runtime_error("LZ4_compress_default failed on " + m_path);
        // Write per-spill frame: [u32 uncompressed][u32 compressed][bytes].
        uint32_t u = (uint32_t)src_size;
        uint32_t c = (uint32_t)compressed;
        if (std::fwrite(&u, sizeof(u), 1, m_file) != 1 ||
            std::fwrite(&c, sizeof(c), 1, m_file) != 1 ||
            std::fwrite(m_out_buf.data(), 1, (size_t)compressed, m_file) != (size_t)compressed) {
            throw std::runtime_error("short write to " + m_path);
        }
        m_total_uncompressed.fetch_add((uint64_t)src_size, std::memory_order_relaxed);
        m_dedup.clear();
        m_bytes = 0;
    }

    // Heterogeneous lookup key bundling bases + a precomputed wyhash.
    // The producer (per_thread_bucket_buffers::append) computes the
    // hash while the bases are still hot in L1; the compactor's
    // find(hashed_view) returns that hash without touching the bases
    // bytes again. Bytes are only re-read on a hash collision, when
    // string_eq runs the byte compare against the stored std::string.
    struct hashed_view {
        std::string_view sv;
        uint64_t hash;
    };

public:
    // wyhash function exposed so producers (per_thread_bucket_buffers
    // ::append) can compute the same hash the compactor's hashmap
    // would have computed internally on find(). Same call ankerl uses
    // for its default std::string_view hash, so the values match.
    static uint64_t hash_bases(std::string_view sv) noexcept {
        return ankerl::unordered_dense::hash<std::string_view>{}(sv);
    }
    static uint64_t hash_bases(uint8_t const* p, uint32_t n) noexcept {
        return hash_bases(std::string_view((char const*)p, n));
    }

private:
    struct string_hash {
        using is_transparent = void;
        using is_avalanching = void;
        size_t operator()(std::string_view sv) const noexcept {
            return ankerl::unordered_dense::hash<std::string_view>{}(sv);
        }
        size_t operator()(hashed_view const& h) const noexcept { return (size_t)h.hash; }
    };
    struct string_eq {
        using is_transparent = void;
        bool operator()(std::string_view a, std::string_view b) const noexcept { return a == b; }
        bool operator()(std::string_view a, hashed_view const& b) const noexcept {
            return a == b.sv;
        }
        bool operator()(hashed_view const& a, std::string_view b) const noexcept {
            return a.sv == b;
        }
        bool operator()(hashed_view const& a, hashed_view const& b) const noexcept {
            return a.sv == b.sv;
        }
    };

    std::string m_path;
    size_t m_spill_bytes;
    std::mutex m_mu;
    std::FILE* m_file = nullptr;
    // m_batch_buf accumulates the serialized super-k-mer records of one
    // spill so we can compress them as a single LZ4 block. m_out_buf
    // holds the LZ4 output. Both grow lazily to one spill's high-water
    // and stay there for the rest of the phase (cleared but not
    // shrunk). For our typical spill_bytes ~ 16-256 KiB the combined
    // cost is a few hundred KiB per bucket.
    std::vector<uint8_t> m_batch_buf;
    std::vector<uint8_t> m_out_buf;
    ankerl::unordered_dense::map<std::string, entry, string_hash, string_eq> m_dedup;
    size_t m_bytes = 0;
    std::atomic<uint64_t> m_total_uncompressed{0};
};

// ---- bucket_writer: owns one bucket_compactor per bucket ----------------------

class bucket_writer {
public:
    bucket_writer(std::string const& dir, uint32_t num_buckets, size_t flush_bases = 64 * 1024,
                  size_t spill_bytes = DEFAULT_COMPACTOR_SPILL_BYTES)
        : m_dir(dir), m_num_buckets(num_buckets), m_flush_bases(flush_bases) {
        std::filesystem::create_directories(m_dir);
        m_compactors.reserve(num_buckets);
        for (uint32_t b = 0; b < num_buckets; ++b) {
            m_compactors.emplace_back(
                std::make_unique<bucket_compactor>(bucket_path(b), spill_bytes));
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
        // Stat each file once for the on-disk byte total. Also collect the
        // per-bucket uncompressed sizes so we can report the bucket-size
        // distribution: bucket-process holds ONE bucket's kmer_info resident
        // per in-flight thread, so the peak working set is driven by the
        // LARGEST buckets, not the average. The distribution (and the sum of
        // the top-`report_threads` buckets) tells us how the bucket-process
        // admission gate will behave -- how many buckets fit in -g at once, and
        // thus the effective concurrency. Cheap: one stat()/accessor per bucket.
        uint64_t total_compressed = 0;
        uint64_t total_uncompressed = 0;
        std::vector<uint64_t> unc;
        unc.reserve(m_num_buckets);
        for (uint32_t b = 0; b < m_num_buckets; ++b) {
            struct stat st;
            if (::stat(bucket_path(b).c_str(), &st) == 0) total_compressed += (uint64_t)st.st_size;
            uint64_t u = m_compactors[b]->total_uncompressed_bytes();
            total_uncompressed += u;
            unc.push_back(u);
        }
        m_total_compressed.store(total_compressed, std::memory_order_relaxed);
        m_total_uncompressed.store(total_uncompressed, std::memory_order_relaxed);
        m_bucket_unc_sizes = std::move(unc);
        m_closed = true;
    }

    // Free the per-bucket compactor objects (dedup hashmaps + LZ4 scratch
    // buffers) after close(). The on-disk bucket files and the cached per-bucket
    // sizes (m_bucket_unc_sizes) remain, so bucket_path()/bucket_unc_bytes()/
    // num_buckets() still work and process_buckets can read every bucket. Call
    // this before bucket-process: the compactors hold hundreds of MB - GBs of
    // buffers and pin glibc's churned arenas, so leaving them alive inflates the
    // RSS that bucket-process's live-RSS admission gate reads as its baseline
    // carry -- which starves admission (it would admit ~one bucket at a time).
    // Requires close() first; flush()/close() must not be called afterwards.
    void release_compactors() {
        m_compactors.clear();
        m_compactors.shrink_to_fit();
    }

    // Print the per-bucket (uncompressed) size distribution and the estimated
    // bucket-process working set. `concurrency` is the number of buckets
    // processed at once (num_threads). Call after close().
    void report_bucket_size_distribution(uint32_t concurrency) const {
        if (m_bucket_unc_sizes.empty()) return;
        std::vector<uint64_t> v = m_bucket_unc_sizes;  // copy; sort ascending
        std::sort(v.begin(), v.end());
        const size_t n = v.size();
        auto pct = [&](double p) -> uint64_t {
            size_t i = (size_t)(p * (double)(n - 1));
            return v[i];
        };
        uint64_t total = 0;
        for (uint64_t x : v) total += x;
        // Sum of the largest `concurrency` buckets = worst-case set of buckets
        // resident together (each thread grabs the next bucket; the slowest /
        // biggest can co-reside). This is the bucket-process peak driver, in
        // *on-disk uncompressed* bytes -- the in-RAM kmer_info is a further
        // ~constant multiple of this (record bytes -> hashmap entries).
        uint64_t top_sum = 0;
        for (size_t i = (n > concurrency ? n - concurrency : 0); i < n; ++i) top_sum += v[i];
        uint32_t nonempty = 0;
        for (uint64_t x : v) if (x) ++nonempty;
        std::fprintf(stderr,
            "[bucket size dist] buckets=%zu nonempty=%u  uncompressed bytes:\n"
            "  mean=%.2f MiB  p50=%.2f  p90=%.2f  p99=%.2f  p999=%.2f  max=%.2f MiB\n"
            "  largest-%u-sum=%.2f GiB (worst-case co-resident bucket bytes; "
            "kmer_info RAM is a ~constant multiple of this)\n"
            "  max-bucket=%llu bytes  total=%.2f GiB\n",
            n, nonempty,
            (double)total / n / 1048576.0,
            (double)pct(0.50) / 1048576.0, (double)pct(0.90) / 1048576.0,
            (double)pct(0.99) / 1048576.0, (double)pct(0.999) / 1048576.0,
            (double)v[n - 1] / 1048576.0,
            concurrency, (double)top_sum / 1073741824.0,
            (unsigned long long)v[n - 1], (double)total / 1073741824.0);
    }

    uint64_t total_bytes() const { return m_total_compressed.load(std::memory_order_relaxed); }
    uint64_t total_uncompressed_bytes() const {
        return m_total_uncompressed.load(std::memory_order_relaxed);
    }

    // Sum of the `n` largest buckets' uncompressed bytes = worst-case set of
    // buckets that can be co-resident when `n` threads each process one bucket
    // (an atomic counter hands the biggest ones out; in the worst case the n
    // largest are in flight together). The builder multiplies this by an
    // Uncompressed on-disk bytes of bucket b (0 if unknown). The bucket-process
    // memory-admission gate uses this to charge each bucket's kmer_info RAM
    // (~a fixed multiple of these bytes) against the budget before loading it.
    uint64_t bucket_unc_bytes(uint32_t b) const {
        return b < m_bucket_unc_sizes.size() ? m_bucket_unc_sizes[b] : 0;
    }
    // Largest single bucket's uncompressed bytes (one bucket's worst case).
    uint64_t max_bucket_unc_bytes() const {
        uint64_t m = 0;
        for (uint64_t s : m_bucket_unc_sizes) m = std::max(m, s);
        return m;
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
    // Per-thread bucket buffers and per-compactor scratch (m_batch_buf,
    // m_out_buf) are NOT spillable -- they form a structural floor of
    // approximately
    //   T * B * flush_bases * BUFFER_OVERHEAD
    //     + B * spill_bytes * COMPACTOR_OVERHEAD
    // bytes that no amount of pressure response can reduce. The auto-
    // tune in builder picks flush_bases / spill_bytes / num_buckets so
    // that floor stays under the bucket-write share of -g.
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
    std::vector<uint64_t> m_bucket_unc_sizes;  // per-bucket uncompressed bytes (set in close())
    bool m_closed = false;
    // Pressure flag set by the RSS watcher; informational only
    // (surfaced via under_pressure() / pressure_was_engaged() for the
    // post-phase log line). Backpressure is implemented by the watcher
    // sweeping all compactors with try_spill(), not by ingest threads
    // consulting this flag.
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
        // Hash the bases now while they're still hot in L1 from the
        // emit_super_kmers buffer; the compactor reuses this hash on
        // its dedup-map find() instead of re-hashing the (typically
        // cold) bytes from the writer's bases_storage.
        uint64_t h = bucket_compactor::hash_bases(sk_bases, len);
        recs[b].push_back({h, color, off, len, (uint8_t)(flags & 0xfu)});
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
    // Slurps the whole bucket file (decompressed) into m_buf so that
    // bucket_walker can iterate records via a simple byte cursor.
    //
    // File format (matches bucket_compactor on the write side):
    //   repeat:
    //     [u32 uncompressed_size]   little-endian; 0 marks end-of-stream
    //     [u32 compressed_size]
    //     [compressed bytes]
    //
    // Each frame is one spill. We read uncompressed_size + compressed_size,
    // realloc m_buf if needed, and call LZ4_decompress_safe directly into
    // m_buf at the current write offset. Output buffer grows geometrically.
    explicit bucket_reader(std::string const& path) {
        std::FILE* f = std::fopen(path.c_str(), "rb");
        if (!f) throw std::runtime_error("cannot open " + path + ": " + std::strerror(errno));

        std::vector<uint8_t> comp_buf;
        size_t out_off = 0;

        for (;;) {
            uint32_t u = 0;
            size_t got = std::fread(&u, sizeof(u), 1, f);
            if (got != 1) {
                // Bucket file written without an EOF marker (shouldn't
                // happen with our writer, but handle it anyway): EOF
                // here is fine as long as the stream ends cleanly.
                if (std::feof(f)) break;
                std::fclose(f);
                throw std::runtime_error("short read of frame header on " + path);
            }
            if (u == 0) break;  // explicit end-of-stream marker
            uint32_t c = 0;
            if (std::fread(&c, sizeof(c), 1, f) != 1) {
                std::fclose(f);
                throw std::runtime_error("short read of compressed_size on " + path);
            }
            if (comp_buf.size() < c) comp_buf.resize(c);
            if (std::fread(comp_buf.data(), 1, c, f) != c) {
                std::fclose(f);
                throw std::runtime_error("short read of compressed payload on " + path);
            }
            if (m_buf.size() < out_off + u) m_buf.resize(out_off + u);
            int decoded = LZ4_decompress_safe((char const*)comp_buf.data(),
                                              (char*)m_buf.data() + out_off, (int)c, (int)u);
            if (decoded < 0 || (uint32_t)decoded != u) {
                std::fclose(f);
                throw std::runtime_error("LZ4_decompress_safe failed on " + path);
            }
            out_off += u;
        }
        m_buf.resize(out_off);
        std::fclose(f);
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

}  // namespace cdbg
