#pragma once

// Small project-wide utilities: build configuration, progress reporter,
// and bucket-write profiling counters. Bundled together because each is
// a few-line struct with no internal dependencies.

#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <string>
#include <thread>
#include <unistd.h>

namespace cdgb {

// ---- build configuration ----------------------------------------------------

struct build_config {
    std::string filenames_list;  // text file: one input path per line
    std::string out_basename;    // produces <basename>.fa and <basename>.colors
    uint32_t k = 31;
    uint32_t num_threads = 1;
    uint32_t m = 0;             // minimizer length, 0 = auto (compute_best_m(k))
    uint32_t bucket_log2 = 10;  // 2^10 = 1024 minimizer buckets
    std::string tmp_dir;        // scratch dir; empty -> mkdtemp under $TMPDIR
    bool verbose = false;
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
    // seq_reader::next: gzip decode + kseq parsing.
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
            "  seq_read     %6.2fs   (gzip + kseq parsing)\n"
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

}  // namespace cdgb
