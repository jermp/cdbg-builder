#pragma once

// Cumulative bucket-write timings, summed across ingest threads.
// Reported once at the end of the bucket-write phase. Each value is a
// total of all per-thread time spent in that stage; dividing by
// num_threads gives a wall-clock-equivalent under perfect parallelism.

#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>

namespace cdgb {

struct bucket_write_prof {
    using clock = std::chrono::steady_clock;

    // ---- ingest-loop stages (covers the ingest_file_bucketed body) ----
    // SeqReader::next: gzip decode + kseq parsing.
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
        return (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(
                   clock::now() - t0)
            .count();
    }

    void print(uint32_t num_threads) const {
        auto load = [](const std::atomic<uint64_t>& a) { return a.load(std::memory_order_relaxed); };
        double per_thread = num_threads > 0 ? (double)num_threads : 1.0;
        auto s = [&](const std::atomic<uint64_t>& a) {
            return (double)load(a) / 1e9 / per_thread;
        };
        uint64_t loop_body = load(ns_loop_body);
        uint64_t flush = load(ns_flush);
        double s_compute = (loop_body > flush ? (double)(loop_body - flush) : 0.0) / 1e9 / per_thread;
        std::fprintf(stderr,
                     "[bucket-write profile] (per-thread time, ns/threads -> wall-equiv):\n"
                     "  seq_read     %6.2fs   (gzip + kseq parsing)\n"
                     "  compute      %6.2fs   (ACGT scan + 2-bit + ntHash + minimizer + per-thread append)\n"
                     "  flush        %6.2fs   (writer.flush total = lock + hashmap + spill)\n"
                     "    lock_wait  %6.2fs\n"
                     "    hashmap    %6.2fs\n"
                     "    spill      %6.2fs   (sort+unique + write_super_kmer + gzwrite)\n"
                     "  counts: files=%llu records=%llu flushes=%llu spills=%llu inserts=%llu\n",
                     s(ns_seq_read), s_compute, s(ns_flush), s(ns_lock_wait), s(ns_hashmap),
                     s(ns_spill), (unsigned long long)load(n_files),
                     (unsigned long long)load(n_records),
                     (unsigned long long)load(n_flushes),
                     (unsigned long long)load(n_spills),
                     (unsigned long long)load(n_inserts));
    }
};

inline bucket_write_prof& bucket_prof() {
    static bucket_write_prof p;
    return p;
}

}  // namespace cdgb
