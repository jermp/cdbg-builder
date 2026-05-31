#pragma once

// Small project-wide utilities: build configuration, progress reporter,
// and bucket-write profiling counters. Bundled together because each is
// a few-line struct with no internal dependencies.

#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/resource.h>
#include <thread>
#include <unistd.h>

namespace cdbg {

// ---- build configuration ----------------------------------------------------

struct build_config {
    std::string filenames_list;  // text file: one input path per line
    std::string out_basename;    // produces <basename>.fa and <basename>.colors
    uint32_t k = 31;
    uint32_t num_threads = 1;
    uint32_t m = 0;             // minimizer length, 0 = auto (compute_best_m(k))
    uint32_t bucket_log2 = 0;   // 0 = auto (derived from max_ram_gb if set, else 10)
    std::string tmp_dir;        // scratch dir; empty -> mkdtemp under $TMPDIR
    // Soft RAM budget in GiB. 0 = no budget. When set, the builder
    // auto-picks bucket_log2 (more buckets -> smaller per-bucket data
    // structures) and streams the encoded color bit_vector to a
    // sidecar file when it would exceed a fraction of the budget.
    // The peak RSS is reported at end-of-build; if it exceeded the
    // budget, the report says so (the build is not killed -- the
    // budget is a soft target, not a hard cap).
    double max_ram_gb = 0.0;
    bool verbose = false;
};

// ---- timer ------------------------------------------------------------------

// RAII phase timer: prints "[label] X.XX s" on destruction. Used to bracket
// individual stages inside cdbg_builder::build().
class timer {
public:
    timer(char const* label) : m_label(label), m_t0(std::chrono::steady_clock::now()) {}
    ~timer() {
        auto t1 = std::chrono::steady_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - m_t0).count();
        std::cout << "[" << m_label << "] " << (ms / 1000.0) << " s\n";
    }

    timer(timer const&) = delete;
    timer& operator=(timer const&) = delete;

private:
    char const* m_label;
    std::chrono::steady_clock::time_point m_t0;
};

// ---- process-memory query ---------------------------------------------------
//
// Peak resident set size in bytes via getrusage(RUSAGE_SELF) — POSIX,
// portable across Linux / macOS / BSD. The kernel maintains the peak
// continuously, so a single read at end of build() gives a true
// high-water mark with zero hot-path overhead.
//
// Unit gotcha: Linux's man page says ru_maxrss is in KiB; Darwin
// returns bytes. The #ifdef below normalises to bytes on both.

inline uint64_t process_peak_rss_bytes() {
    struct rusage ru;
    if (::getrusage(RUSAGE_SELF, &ru) != 0) return 0;
#ifdef __APPLE__
    return (uint64_t)ru.ru_maxrss;
#else
    return (uint64_t)ru.ru_maxrss * 1024ULL;
#endif
}

// Pretty-printer: 3.42 GiB / 728 MiB / 12 KiB, picking the largest
// unit at which the number is >= 1.
inline std::string format_bytes(uint64_t b) {
    static char const* units[] = {"B", "KiB", "MiB", "GiB", "TiB"};
    double v = (double)b;
    int u = 0;
    while (v >= 1024.0 && u + 1 < (int)(sizeof(units) / sizeof(units[0]))) {
        v /= 1024.0;
        ++u;
    }
    char buf[64];
    if (u == 0)
        std::snprintf(buf, sizeof(buf), "%llu B", (unsigned long long)b);
    else
        std::snprintf(buf, sizeof(buf), "%.2f %s", v, units[u]);
    return buf;
}

// Current resident set size in bytes (live RSS, not the lifetime peak).
// Linux: parses VmRSS from /proc/self/status using a single read()
// syscall to be robust against environments where stdio fopen() is
// flaky. Returns 0 if unavailable -- callers should fall back to
// process_peak_rss_bytes() (lifetime monotonic peak) in that case.
//
// Used by the bucket-write RSS watcher for hysteresis: trip pressure
// at a high threshold, release at a low threshold once the spills
// have actually brought live RSS back down. process_peak_rss_bytes()
// is monotonic so it can't observe the release.
inline uint64_t current_rss_bytes() {
#if defined(__linux__)
    int fd = ::open("/proc/self/status", O_RDONLY | O_CLOEXEC);
    if (fd < 0) return 0;
    char buf[4096];
    ssize_t total = 0;
    for (;;) {
        ssize_t n = ::read(fd, buf + total, sizeof(buf) - 1 - (size_t)total);
        if (n <= 0) break;
        total += n;
        if ((size_t)total >= sizeof(buf) - 1) break;
    }
    ::close(fd);
    if (total <= 0) return 0;
    buf[total] = 0;
    char const* p = std::strstr(buf, "VmRSS:");
    if (!p) return 0;
    p += 6;
    while (*p == ' ' || *p == '\t') ++p;
    char* end = nullptr;
    unsigned long kb = std::strtoul(p, &end, 10);
    if (end == p) return 0;
    return (uint64_t)kb * 1024ULL;
#else
    return 0;
#endif
}

// ---- per-phase RSS marker --------------------------------------------------
//
// process_peak_rss_bytes() above is the lifetime high-water mark via
// getrusage; it is monotonically non-decreasing across the run. We
// piggy-back on it: at the end of each phase, snapshot it and print
//
//     [<phase> peak RSS] X.XX GiB (+Y.YY GiB)
//
// where the parenthesised delta is the increase since the previous
// snapshot. A non-zero delta means this phase pushed the high-water
// mark; a zero delta means everything this phase touched fit under
// the previous peak. The signal is coarser than an in-phase sampler
// (a phase that allocates and frees within itself without crossing
// the prior watermark shows +0.00) but it works wherever getrusage
// works -- no /proc parsing, no sampler thread.

class phase_rss_marker {
public:
    explicit phase_rss_marker(std::string label) : m_label(std::move(label)) {
        m_baseline = process_peak_rss_bytes();
    }
    ~phase_rss_marker() { stop(); }

    phase_rss_marker(phase_rss_marker const&) = delete;
    phase_rss_marker& operator=(phase_rss_marker const&) = delete;

    void stop() {
        if (m_stopped) return;
        m_stopped = true;
        uint64_t now = process_peak_rss_bytes();
        if (now == 0) {
            std::cout << "  [" << m_label << " peak RSS] unavailable\n";
            return;
        }
        std::cout << "  [" << m_label << " peak RSS] " << format_bytes(now);
        if (now > m_baseline) {
            std::cout << " (+" << format_bytes(now - m_baseline) << ")";
        } else {
            std::cout << " (+0)";
        }
        std::cout << "\n";
    }

private:
    std::string m_label;
    uint64_t m_baseline = 0;
    bool m_stopped = false;
};

// ---- progress reporter ------------------------------------------------------

// Lightweight progress printer for long-running phases.
//
// Spawns a background thread that periodically samples a worker-owned atomic
// counter and prints a status line to stderr. When stderr is a TTY it updates
// a single line in place via '\r'; otherwise it emits one line per tick so
// progress is still visible in piped logs.
//
// Workers don't talk to progress directly — they bump `counter`, and progress
// only reads it. This keeps the hot path lock-free and progress optional.
class progress {
public:
    progress(std::string label, std::atomic<uint64_t>& counter, uint64_t total,
             std::ostream& os = std::cerr,
             std::chrono::milliseconds interval = std::chrono::milliseconds(500))
        : m_label(std::move(label))
        , m_counter(counter)
        , m_total(total)
        , m_os(os)
        , m_interval(interval)
        , m_is_tty(::isatty(2) != 0)
        , m_start(std::chrono::steady_clock::now()) {
        m_thread = std::thread([this] { run(); });
    }

    ~progress() { stop(); }

    progress(progress const&) = delete;
    progress& operator=(progress const&) = delete;

    void stop() {
        bool was = m_stopped.exchange(true);
        if (was) return;
        if (m_thread.joinable()) m_thread.join();
        print_line(true);
    }

private:
    void run() {
        // Print one line up front so the user sees the phase started.
        print_line(false);
        while (!m_stopped.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(m_interval);
            if (m_stopped.load(std::memory_order_relaxed)) break;
            print_line(false);
        }
    }

    void print_line(bool final) {
        uint64_t done = m_counter.load(std::memory_order_relaxed);
        if (final) done = m_total ? m_total : done;

        auto now = std::chrono::steady_clock::now();
        double sec =
            std::chrono::duration_cast<std::chrono::milliseconds>(now - m_start).count() / 1000.0;

        char buf[256];
        if (m_total) {
            double pct = 100.0 * (double)done / (double)m_total;
            std::snprintf(buf, sizeof(buf), "[%s] %llu/%llu (%.1f%%) %.1fs", m_label.c_str(),
                          (unsigned long long)done, (unsigned long long)m_total, pct, sec);
        } else {
            std::snprintf(buf, sizeof(buf), "[%s] %llu %.1fs", m_label.c_str(),
                          (unsigned long long)done, sec);
        }

        if (m_is_tty) {
            m_os << '\r' << buf << "      ";
            if (final) m_os << '\n';
        } else {
            // Skip mid-flight prints if nothing changed since last line.
            if (!final && done == m_last_done) return;
            m_last_done = done;
            m_os << buf << '\n';
        }
        m_os.flush();
    }

    std::string m_label;
    std::atomic<uint64_t>& m_counter;
    uint64_t m_total;
    std::ostream& m_os;
    std::chrono::milliseconds m_interval;
    bool m_is_tty;
    std::chrono::steady_clock::time_point m_start;
    std::atomic<bool> m_stopped{false};
    std::thread m_thread;
    uint64_t m_last_done = (uint64_t)-1;
};

// ---- bucket-write profiling counters ----------------------------------------
//
// Cumulative bucket-write timings, summed across ingest threads. Reported
// once at the end of the bucket-write phase. Each value is a total of all
// per-thread time spent in that stage; dividing by num_threads gives a
// wall-clock-equivalent under perfect parallelism.

struct bucket_write_prof {
    using clock = std::chrono::steady_clock;

    // ---- ingest-loop stages (covers the ingest_file_bucketed body) ----
    // seq_reader::next: kseq parsing (libdeflate gzip decompression
    // happens upfront in seq_reader's ctor and is amortised across the
    // file's records, not counted here).
    std::atomic<uint64_t> ns_seq_read{0};
    // Whole loop body (ACGT scan, 2-bit conversion, emit_super_kmers,
    // including its calls into the per-thread buffer and any flushes).
    // Subtract ns_flush below to isolate pure compute (parse + minimizer).
    std::atomic<uint64_t> ns_loop_body{0};

    // ---- bucket-writer / compactor stages (subsets of ns_loop_body) ----
    // Total time inside bucket_writer::flush() (= insert_batch + spills).
    std::atomic<uint64_t> ns_flush{0};
    // Time waiting on the per-bucket compactor mutex.
    std::atomic<uint64_t> ns_lock_wait{0};
    // Hashmap find/insert/update inside insert_batch (excludes spill).
    std::atomic<uint64_t> ns_hashmap{0};
    // sort+unique colors + write_super_kmer + gzwrite for each entry.
    std::atomic<uint64_t> ns_spill{0};

    // ---- counters ----
    std::atomic<uint64_t> n_files{0};
    std::atomic<uint64_t> n_records{0};
    std::atomic<uint64_t> n_flushes{0};
    std::atomic<uint64_t> n_spills{0};
    std::atomic<uint64_t> n_inserts{0};

    // Helper: nanoseconds since `t0`.
    static inline uint64_t since(clock::time_point t0) {
        return (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now() - t0)
            .count();
    }

    void print(uint32_t num_threads) const {
        auto load = [](std::atomic<uint64_t> const& a) {
            return a.load(std::memory_order_relaxed);
        };
        double per_thread = num_threads > 0 ? (double)num_threads : 1.0;
        auto s = [&](std::atomic<uint64_t> const& a) { return (double)load(a) / 1e9 / per_thread; };
        uint64_t loop_body = load(ns_loop_body);
        uint64_t flush = load(ns_flush);
        double s_compute =
            (loop_body > flush ? (double)(loop_body - flush) : 0.0) / 1e9 / per_thread;
        std::fprintf(
            stderr,
            "[bucket-write profile] (per-thread time, ns/threads -> wall-equiv):\n"
            "  seq_read     %6.2fs   (libdeflate gzip + kseq parsing)\n"
            "  compute      %6.2fs   (ACGT scan + 2-bit + ntHash + minimizer + per-thread append)\n"
            "  flush        %6.2fs   (writer.flush total = lock + hashmap + spill)\n"
            "    lock_wait  %6.2fs\n"
            "    hashmap    %6.2fs\n"
            "    spill      %6.2fs   (sort+unique + write_super_kmer + gzwrite)\n"
            "  counts: files=%llu records=%llu flushes=%llu spills=%llu inserts=%llu\n",
            s(ns_seq_read), s_compute, s(ns_flush), s(ns_lock_wait), s(ns_hashmap), s(ns_spill),
            (unsigned long long)load(n_files), (unsigned long long)load(n_records),
            (unsigned long long)load(n_flushes), (unsigned long long)load(n_spills),
            (unsigned long long)load(n_inserts));
    }
};

inline bucket_write_prof& bucket_prof() {
    static bucket_write_prof p;
    return p;
}

// ---- bucket-process profiling counters --------------------------------------
//
// Mirrors bucket_write_prof for the bucket-process phase. Phase-level
// timers only (per-record / per-kmer timers would dwarf the work
// they're trying to measure: 24.5B records for the 50K workload at
// ~20ns/clock_now() == ~10 minutes of overhead). Phase counts let us
// see whether the 1141s is mostly load_bucket vs rsid->cid resolve
// vs walk vs the global-dict merge under mutex.
struct bucket_process_prof {
    using clock = std::chrono::steady_clock;

    // ---- per-bucket phases inside process_bucket ----
    // load_bucket: read every super-k-mer record, intern its color list
    // into record_sets, roll k-mers into kmer_info. Most likely the
    // dominant line because it touches every record + every k-mer.
    std::atomic<uint64_t> ns_load{0};
    // rsid->cid resolve: for each k-mer, decode its (single or merged)
    // record-color list and intern into local_dict. ~1 intern per
    // distinct local class per bucket.
    std::atomic<uint64_t> ns_resolve{0};
    // walk: classify_left_end + extend_and_emit passes. Repeated
    // hashmap lookups, no color-list work.
    std::atomic<uint64_t> ns_walk{0};

    // ---- merge work split into pre-* (lock-free) and locked phases ----
    // Lock-free pre-merge: decode each local_dict class into a vector
    // (compact_color_set_dict::at).
    std::atomic<uint64_t> ns_pre_decode{0};
    // Lock-free pre-merge: wyhash + fnv1a of the decoded color list.
    // Heavy on multi-thousand-color lists; previously ran inside the
    // critical section.
    std::atomic<uint64_t> ns_pre_hash{0};

    // ---- global merge (in process_buckets driver, under global_mu) ----
    // Time waiting to acquire global_mu before merging this bucket's
    // local_dict into the streaming global_dict.
    std::atomic<uint64_t> ns_merge_lock_wait{0};
    // Time held under global_mu: dedup-find + on-miss encode + offset
    // sidecar write. With pre-decode + pre-hash done lock-free this is
    // just the hashmap probe and (rarely) the encoder work.
    std::atomic<uint64_t> ns_merge{0};

    // ---- counters ----
    std::atomic<uint64_t> n_buckets{0};
    std::atomic<uint64_t> n_records{0};       // total super-k-mer records read
    std::atomic<uint64_t> n_kmers{0};         // total k-mers rolled (sum of bases.size()-k+1)
    std::atomic<uint64_t> n_local_classes{0}; // sum of local_dict.size() across buckets
    std::atomic<uint64_t> n_unitigs{0};       // total bucket_unitigs emitted

    static inline uint64_t since(clock::time_point t0) {
        return (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now() - t0)
            .count();
    }

    void print(uint32_t num_threads) const {
        auto load_a = [](std::atomic<uint64_t> const& a) {
            return a.load(std::memory_order_relaxed);
        };
        double per_thread = num_threads > 0 ? (double)num_threads : 1.0;
        auto s = [&](std::atomic<uint64_t> const& a) {
            return (double)load_a(a) / 1e9 / per_thread;
        };
        std::fprintf(
            stderr,
            "[bucket-process profile] (per-thread time, ns/threads -> wall-equiv):\n"
            "  load        %7.2fs   (bucket_reader + LZ4 + record intern + kmer hashmap roll)\n"
            "  resolve     %7.2fs   (rsid -> local cid; record_sets.at + local_dict.intern)\n"
            "  walk        %7.2fs   (classify_left_end + extend_and_emit)\n"
            "  pre_decode  %7.2fs   (local_dict.at, lock-free)\n"
            "  pre_hash    %7.2fs   (wyhash + fnv1a on decoded class, lock-free)\n"
            "  merge_wait  %7.2fs   (waiting on global_mu)\n"
            "  merge       %7.2fs   (global_dict.intern_with_hashes under global_mu)\n"
            "  counts: buckets=%llu records=%llu kmers=%llu local_classes=%llu unitigs=%llu\n",
            s(ns_load), s(ns_resolve), s(ns_walk), s(ns_pre_decode), s(ns_pre_hash),
            s(ns_merge_lock_wait), s(ns_merge),
            (unsigned long long)load_a(n_buckets), (unsigned long long)load_a(n_records),
            (unsigned long long)load_a(n_kmers), (unsigned long long)load_a(n_local_classes),
            (unsigned long long)load_a(n_unitigs));
    }
};

inline bucket_process_prof& process_prof() {
    static bucket_process_prof p;
    return p;
}

}  // namespace cdbg
