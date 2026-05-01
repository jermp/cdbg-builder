#pragma once

// Minimizer-based bucketing ingest. Each input file is one color.
//
// For every ACGT-only run in every sequence, we slide a window of (k-m+1)
// m-mers and find the canonical-ntHash minimum. Whenever the window's minimum
// changes, the current super-k-mer ends and a new one begins. Each super-k-mer
// is appended to the bucket file selected by `(min_hash >> 1) & (B-1)`,
// matching GGCAT's bucket-id derivation in cn_nthash.rs.
//
// The output is a directory of B bucket files, each holding a stream of
// super-k-mer records (see super_kmer.hpp). Thread safety: each ingest worker
// holds a private PerThreadBucketBuffers; concurrency on the shared bucket
// files is one mutex per bucket.

#include <atomic>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

#include "bucket_io.hpp"
#include "kmer.hpp"
#include "minimizer.hpp"
#include "seq_reader.hpp"
#include "super_kmer.hpp"

namespace cdgb {

namespace detail {

// Emit super-k-mers for a single ACGT-only run of length L (>= k).
// `bases` holds 2-bit-encoded bases. Records are appended via `sink` to the
// appropriate bucket files.
inline void emit_super_kmers(const uint8_t* bases, uint32_t L, uint32_t k, uint32_t m,
                             uint32_t color, uint32_t bucket_log2,
                             PerThreadBucketBuffers& sink) {
    if (L < k) return;
    const uint32_t K = L - k + 1;  // number of k-mers
    const uint32_t W = k - m + 1;  // window size, in m-mer indices
    const uint64_t bucket_mask = (uint64_t(1) << bucket_log2) - 1;

    auto bucket_of = [&](uint64_t h) -> uint32_t {
        // Skip the bottom bit to match GGCAT's "uniqueness-flag" reservation.
        return (uint32_t)((h >> 1) & bucket_mask);
    };

    uint64_t fwd = 0, rc = 0;
    if (!nthash_init(bases, m, fwd, rc)) return;  // run is ACGT-only, shouldn't happen

    MinQueue mq;
    mq.reset((int32_t)W);
    mq.push(canonical_mhash(fwd, rc), 0);

    // Fill the first k-mer's window: m-mers at positions 0..k-m.
    for (uint32_t i = 1; i <= k - m; ++i) {
        uint8_t out_b = bases[i - 1];
        uint8_t in_b = bases[i + m - 1];
        nthash_roll(out_b, in_b, m, fwd, rc);
        mq.push(canonical_mhash(fwd, rc), (int32_t)i);
    }

    uint32_t super_start = 0;  // first k-mer index in the running super-k-mer
    bool is_first_super = true;
    uint64_t cur_min = mq.min_hash();
    uint32_t cur_bucket = bucket_of(cur_min);

    auto emit_super = [&](uint32_t first_kmer, uint32_t last_kmer, uint32_t bucket,
                          bool is_run_begin, bool is_run_end) {
        uint32_t base_start = first_kmer;
        uint32_t base_len = (last_kmer - first_kmer) + k;
        uint8_t flags = 0;
        if (is_run_begin) flags |= SK_FLAG_IS_ACGT_BEGIN;
        if (is_run_end) flags |= SK_FLAG_IS_ACGT_END;
        sink.append(bucket, flags, color, bases + base_start, base_len);
    };

    for (uint32_t i = 1; i < K; ++i) {
        uint32_t mpos = i + k - m;
        uint8_t out_b = bases[mpos - 1];
        uint8_t in_b = bases[mpos + m - 1];
        nthash_roll(out_b, in_b, m, fwd, rc);
        mq.push(canonical_mhash(fwd, rc), (int32_t)mpos);

        uint64_t new_min = mq.min_hash();
        if (new_min != cur_min) {
            emit_super(super_start, i - 1, cur_bucket, is_first_super, /*is_run_end=*/false);
            super_start = i;
            is_first_super = false;
            cur_min = new_min;
            cur_bucket = bucket_of(new_min);
        }
    }
    emit_super(super_start, K - 1, cur_bucket, is_first_super, /*is_run_end=*/true);
}

inline void ingest_file_bucketed(const std::string& path, uint32_t k, uint32_t m,
                                 uint32_t bucket_log2, uint32_t color,
                                 PerThreadBucketBuffers& sink) {
    SeqReader r(path);
    const char* s = nullptr;
    size_t l = 0;
    std::vector<uint8_t> bases_buf;
    while (r.next(s, l)) {
        size_t pos = 0;
        while (pos < l) {
            // Find an ACGT-only run starting at pos.
            size_t end = pos;
            while (end < l && nuc_to_2bit(s[end]) != 0xff) ++end;
            size_t run_len = end - pos;
            if (run_len >= k) {
                bases_buf.resize(run_len);
                for (size_t i = 0; i < run_len; ++i) { bases_buf[i] = nuc_to_2bit(s[pos + i]); }
                emit_super_kmers(bases_buf.data(), (uint32_t)run_len, k, m, color, bucket_log2,
                                 sink);
            }
            pos = end;
            while (pos < l && nuc_to_2bit(s[pos]) == 0xff) ++pos;
        }
    }
}

}  // namespace detail

// Parallel driver. Spawns `num_threads` workers, each pulling files from a
// shared queue. `done` (if non-null) is incremented after each file finishes.
inline void ingest_bucketed(const std::vector<std::string>& files, uint32_t k, uint32_t m,
                            uint32_t bucket_log2, BucketWriter& writer, uint32_t num_threads,
                            std::atomic<uint64_t>* done = nullptr) {
    if (num_threads == 0) num_threads = 1;
    std::atomic<size_t> next{0};
    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    auto run = [&]() {
        PerThreadBucketBuffers bufs(writer);
        for (;;) {
            size_t i = next.fetch_add(1);
            if (i >= files.size()) break;
            try {
                detail::ingest_file_bucketed(files[i], k, m, bucket_log2, (uint32_t)i, bufs);
            } catch (std::exception& e) {
                std::cerr << "error ingesting " << files[i] << ": " << e.what() << '\n';
            }
            if (done) done->fetch_add(1, std::memory_order_relaxed);
        }
        bufs.flush_all();
    };

    for (uint32_t t = 0; t < num_threads; ++t) workers.emplace_back(run);
    for (auto& w : workers) w.join();
}

}  // namespace cdgb
