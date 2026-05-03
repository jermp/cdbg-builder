#pragma once

// Disk-backed bucket I/O for the minimizer-bucketed ingest path, with
// inline compaction.
//
// Each bucket on disk holds a stream of *compacted* super-k-mer records
// (see super_kmer.hpp): a super-k-mer that occurs in many input files is
// stored once with the union of its colors. Compaction happens online in
// per-bucket hashmaps; ingest threads flush their per-thread buffers into
// the bucket's hashmap under the bucket mutex, and the hashmap spills its
// contents to disk (LZ4 frame, fast mode) whenever it crosses a memory
// budget. The bucket walker correctly merges color sets across multiple
// spilled records for the same super-k-mer, so spilling is purely a
// memory bound.
//
// On dense pangenome inputs (e.g. many closely-related bacterial genomes)
// most super-k-mers are shared across colors, so compaction is the main
// disk-write reduction over a naive one-record-per-occurrence layout.
//
// LZ4 (frame API) replaced zlib level-1 here for two reasons:
//   1. Per-bucket compression-context state is much smaller (~24 KiB
//      with autoFlush=1 and max64KB blocks vs ~256 KiB for zlib's
//      deflate state + hash chains + sliding window). On runs with
//      thousands of buckets this slashes the unspillable structural
//      floor of bucket-write.
//   2. LZ4 fast mode is 3-5x faster than zlib level-1 at similar (or
//      slightly better) compression ratios on our small-record stream,
//      so spill time drops too.
//
// File format: an LZ4 frame written via LZ4F_compressBegin /
// LZ4F_compressUpdate / LZ4F_compressEnd. Bucket files are ~1.5-2x
// larger than they were under zlib level-1 in practice; that's the
// tradeoff for the much smaller in-RAM state.

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
#include <lz4frame.h>

#include <unordered_dense/unordered_dense.h>

#include "super_kmer.hpp"
#include "util.hpp"

namespace cdgb {

// ---- Per-bucket in-memory compactor + gzip-streaming output -----------------

// Default per-bucket hashmap memory budget (estimated). With 1024 buckets
// at this budget the global ingest peak from the compactors is ~256 MiB.
inline constexpr size_t DEFAULT_COMPACTOR_SPILL_BYTES = 256 * 1024;
// Per-bucket LZ4 output scratch buffer initial size. Records are
// typically <1 KiB and autoFlush=1 produces output proportional to
// input, so 16 KiB covers the common case; compress_chunk_to_file_locked
// resizes (rarely) for atypical large records.
inline constexpr size_t DEFAULT_COMPACTOR_LZ4_OUT_BUF_BYTES = 16 * 1024;

class bucket_compactor {
public:
    // try_spill() is the writer's RSS-watcher entry point: when the
    // process is over the --max-ram budget the watcher walks every
    // compactor and force-spills any with pending data. Per-batch
    // ingest itself is unaware of pressure; spill cadence stays at
    // m_bytes >= m_spill_bytes so dedup state per spill remains
    // worthwhile.
    //
    // The spill stream is an LZ4 frame: blockSize=max64KB,
    // blockMode=independent, contentChecksum=off, blockChecksum=off,
    // compressionLevel=0 (fast). autoFlush=1 keeps the cctx's internal
    // tmpBuff small at the cost of slightly worse compression on
    // sub-block writes -- a good trade for thousands of compactors
    // alive at once.
    bucket_compactor(std::string path, size_t spill_bytes,
                     size_t lz4_out_buf_bytes = DEFAULT_COMPACTOR_LZ4_OUT_BUF_BYTES)
        : m_path(std::move(path)), m_spill_bytes(spill_bytes),
          m_out_buf(lz4_out_buf_bytes) {
        m_file = std::fopen(m_path.c_str(), "wb");
        if (!m_file)
            throw std::runtime_error("cannot open bucket file: " + m_path + ": " +
                                     std::strerror(errno));
        LZ4F_errorCode_t e = LZ4F_createCompressionContext(&m_cctx, LZ4F_VERSION);
        if (LZ4F_isError(e)) {
            std::fclose(m_file);
            m_file = nullptr;
            throw std::runtime_error(std::string("LZ4F_createCompressionContext failed: ") +
                                     LZ4F_getErrorName(e));
        }
        LZ4F_preferences_t prefs;
        std::memset(&prefs, 0, sizeof(prefs));
        prefs.frameInfo.blockSizeID = LZ4F_max64KB;
        prefs.frameInfo.blockMode = LZ4F_blockIndependent;
        prefs.frameInfo.contentChecksumFlag = LZ4F_noContentChecksum;
        prefs.frameInfo.blockChecksumFlag = LZ4F_noBlockChecksum;
        prefs.compressionLevel = 0;
        prefs.autoFlush = 1;
        m_prefs = prefs;
        size_t hdr = LZ4F_compressBegin(m_cctx, m_out_buf.data(), m_out_buf.size(), &prefs);
        if (LZ4F_isError(hdr))
            throw std::runtime_error(std::string("LZ4F_compressBegin failed on ") + m_path +
                                     ": " + LZ4F_getErrorName(hdr));
        if (hdr > 0 && std::fwrite(m_out_buf.data(), 1, hdr, m_file) != hdr)
            throw std::runtime_error("short write to " + m_path);
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
            // Finalise the LZ4 frame (writes the end-of-frame marker)
            // before closing the file. Errors here are fatal: the
            // bucket file would be unreadable otherwise.
            size_t end = LZ4F_compressEnd(m_cctx, m_out_buf.data(), m_out_buf.size(), nullptr);
            if (LZ4F_isError(end))
                throw std::runtime_error(std::string("LZ4F_compressEnd failed on ") + m_path +
                                         ": " + LZ4F_getErrorName(end));
            if (end > 0 && std::fwrite(m_out_buf.data(), 1, end, m_file) != end)
                throw std::runtime_error("short write to " + m_path);
            std::fclose(m_file);
            m_file = nullptr;
        }
        if (m_cctx) {
            LZ4F_freeCompressionContext(m_cctx);
            m_cctx = nullptr;
        }
    }

private:
    struct entry {
        std::vector<uint32_t> colors;  // ascending, deduped
        uint8_t flags = 0;
    };

    // Compress `src_size` bytes from `src` into the frame and write to
    // disk. Grows m_out_buf if needed -- worst-case for autoFlush=1 is
    // bounded by LZ4F_compressBound(src_size, &prefs), which for our
    // typical record sizes (<1 KiB) fits comfortably in the default
    // 80 KiB scratch.
    void compress_chunk_to_file_locked(uint8_t const* src, size_t src_size) {
        size_t bound = LZ4F_compressBound(src_size, &m_prefs);
        if (m_out_buf.size() < bound) m_out_buf.resize(bound);
        size_t out = LZ4F_compressUpdate(m_cctx, m_out_buf.data(), m_out_buf.size(), src,
                                         src_size, nullptr);
        if (LZ4F_isError(out))
            throw std::runtime_error(std::string("LZ4F_compressUpdate failed on ") + m_path +
                                     ": " + LZ4F_getErrorName(out));
        if (out > 0 && std::fwrite(m_out_buf.data(), 1, out, m_file) != out)
            throw std::runtime_error("short write to " + m_path);
    }

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
            compress_chunk_to_file_locked(rec_buf.data(), rec_buf.size());
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
    std::mutex m_mu;
    std::FILE* m_file = nullptr;
    LZ4F_cctx* m_cctx = nullptr;
    LZ4F_preferences_t m_prefs{};
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
    // Per-thread bucket buffers and LZ4 cctx state are NOT spillable --
    // they form a structural floor of approximately
    //   T*B*flush_bases*overhead + B*(lz4_cctx_state + lz4_out_buf)
    // bytes that no amount of pressure response can reduce. The auto-
    // tune in builder picks flush_bases / num_buckets to keep that
    // floor under the budget share for bucket-write.
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
    // Slurps the whole bucket file (decompressed) into m_buf so that
    // bucket_walker can iterate records via a simple byte cursor. The
    // input file is an LZ4 frame produced by bucket_compactor; we
    // decompress it via LZ4F_decompress in chunks. Output buffer grows
    // geometrically.
    explicit bucket_reader(std::string const& path) {
        std::FILE* f = std::fopen(path.c_str(), "rb");
        if (!f) throw std::runtime_error("cannot open " + path + ": " + std::strerror(errno));

        LZ4F_dctx* dctx = nullptr;
        LZ4F_errorCode_t e = LZ4F_createDecompressionContext(&dctx, LZ4F_VERSION);
        if (LZ4F_isError(e)) {
            std::fclose(f);
            throw std::runtime_error(std::string("LZ4F_createDecompressionContext failed: ") +
                                     LZ4F_getErrorName(e));
        }

        constexpr size_t IN_CHUNK = 64 * 1024;
        constexpr size_t OUT_CHUNK = 256 * 1024;
        std::vector<uint8_t> in_buf(IN_CHUNK);
        size_t out_off = 0;

        for (;;) {
            size_t in_size = std::fread(in_buf.data(), 1, IN_CHUNK, f);
            if (in_size == 0) {
                if (std::ferror(f)) {
                    LZ4F_freeDecompressionContext(dctx);
                    std::fclose(f);
                    throw std::runtime_error("read failed on " + path);
                }
                break;  // EOF
            }
            size_t in_pos = 0;
            while (in_pos < in_size) {
                if (m_buf.size() < out_off + OUT_CHUNK) m_buf.resize(out_off + OUT_CHUNK);
                size_t out_capacity = m_buf.size() - out_off;
                size_t in_remaining = in_size - in_pos;
                size_t hint = LZ4F_decompress(dctx, m_buf.data() + out_off, &out_capacity,
                                              in_buf.data() + in_pos, &in_remaining, nullptr);
                if (LZ4F_isError(hint)) {
                    LZ4F_freeDecompressionContext(dctx);
                    std::fclose(f);
                    throw std::runtime_error(std::string("LZ4F_decompress failed on ") + path +
                                             ": " + LZ4F_getErrorName(hint));
                }
                in_pos += in_remaining;
                out_off += out_capacity;
                if (hint == 0) {
                    // End of frame. Drain any trailing input and stop.
                    in_pos = in_size;
                    break;
                }
                if (out_capacity == 0 && in_remaining == 0) {
                    // Nothing consumed and nothing produced; LZ4 wants
                    // more input -- break out of inner loop to refill.
                    break;
                }
            }
        }
        m_buf.resize(out_off);
        LZ4F_freeDecompressionContext(dctx);
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

}  // namespace cdgb
