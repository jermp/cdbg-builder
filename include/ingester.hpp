#pragma once

#include <atomic>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

#include "concurrent_kmer_map.hpp"
#include "kmer.hpp"
#include "seq_reader.hpp"

namespace cdgb {

// Per-thread, per-shard outgoing buffer. We batch (kmer, color) pairs and
// flush in bulk under one lock acquisition per flush.
//
// FLUSH_THRESHOLD is the per-shard buffer size before a flush is triggered.
struct ShardBuffer {
    std::vector<std::pair<kmer_int_t, uint32_t>> pairs;
};

// Slide canonical k-mers across `seq` and route them to shard buffers.
// Maintains forward and reverse-complement of the current window in parallel
// to avoid per-step O(k) recomputation.
inline void emit_kmers_to_buffers(const char* seq, size_t len, uint32_t k, uint32_t color,
                                  uint64_t shard_mask, std::vector<ShardBuffer>& bufs,
                                  ConcurrentKmerMap& dst, size_t flush_threshold) {
    if (len < k) return;

    auto flush_one = [&](uint64_t s) {
        auto& b = bufs[s];
        if (b.pairs.empty()) return;
        auto& shard = dst.shard(s);
        std::lock_guard<std::mutex> lk(shard.mu);
        for (auto& p : b.pairs) shard.map[p.first].add_color(p.second);
        b.pairs.clear();
    };

    kmer_int_t fwd = 0, rc = 0;
    uint32_t valid = 0;  // how many of the last `valid` characters were valid nucleotides
    const kmer_int_t mask = kmer_mask(k);
    const uint32_t k_minus_1_x2 = 2 * (k - 1);

    for (size_t i = 0; i < len; ++i) {
        uint8_t v = nuc_to_2bit(seq[i]);
        if (v == 0xff) {
            valid = 0;
            fwd = rc = 0;
            continue;
        }
        fwd = ((fwd << 2) | v) & mask;
        rc  = (rc >> 2) | ((kmer_int_t)(v ^ 3) << k_minus_1_x2);
        if (valid + 1 < k) {
            ++valid;
            continue;
        }
        valid = k;
        kmer_int_t can = fwd <= rc ? fwd : rc;
        uint64_t s = KmerHasher{}(can) & shard_mask;
        bufs[s].pairs.emplace_back(can, color);
        if (bufs[s].pairs.size() >= flush_threshold) flush_one(s);
    }
}

// Process one file (one color) into the shared map.
inline void ingest_file(const std::string& path, uint32_t color, uint32_t k,
                        ConcurrentKmerMap& dst, size_t flush_threshold = 4096) {
    const uint64_t shard_mask = dst.num_shards() - 1;
    std::vector<ShardBuffer> bufs(dst.num_shards());

    SeqReader r(path);
    const char* s = nullptr;
    size_t l = 0;
    while (r.next(s, l)) {
        emit_kmers_to_buffers(s, l, k, color, shard_mask, bufs, dst, flush_threshold);
    }
    // final flush
    for (uint64_t i = 0; i < dst.num_shards(); ++i) {
        auto& b = bufs[i];
        if (b.pairs.empty()) continue;
        auto& shard = dst.shard(i);
        std::lock_guard<std::mutex> lk(shard.mu);
        for (auto& p : b.pairs) shard.map[p.first].add_color(p.second);
        b.pairs.clear();
    }
}

// Parallel driver: distribute files across threads. Each thread fully processes
// one file at a time. Files are independent inputs (one color each).
//
// If `done` is non-null, it is incremented by 1 after each file completes; a
// caller-side Progress can then report file-level ingest progress.
inline void ingest_parallel(const std::vector<std::string>& files, uint32_t k,
                            ConcurrentKmerMap& dst, uint32_t num_threads,
                            std::atomic<uint64_t>* done = nullptr) {
    if (num_threads == 0) num_threads = 1;
    std::atomic<size_t> next{0};
    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    auto run = [&]() {
        for (;;) {
            size_t i = next.fetch_add(1);
            if (i >= files.size()) return;
            try {
                ingest_file(files[i], (uint32_t)i, k, dst);
            } catch (std::exception& e) {
                std::cerr << "error ingesting " << files[i] << ": " << e.what() << '\n';
            }
            if (done) done->fetch_add(1, std::memory_order_relaxed);
        }
    };

    for (uint32_t t = 0; t < num_threads; ++t) workers.emplace_back(run);
    for (auto& w : workers) w.join();
}

}  // namespace cdgb
