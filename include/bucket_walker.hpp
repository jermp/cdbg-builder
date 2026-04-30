#pragma once

// Per-bucket de Bruijn graph construction and unitig walk.
//
// One bucket file holds all super-k-mers whose minimizer hashes into this
// bucket id. By construction every canonical k-mer lives in exactly one
// bucket, so a bucket-local walker produces a correct (but possibly
// fragmented at minimizer-change boundaries) compacted colored dBG.
//
// Workflow per bucket:
//   1. Read super-k-mer records from the bucket file.
//   2. Slide canonical k-mers, accumulating per-k-mer color lists.
//   3. Build a local ColorSetDict and a 1-shard FinalKmerMap.
//   4. Run UnitigWalker on it.
//   5. Re-attach each unitig's color set inline (BucketUnitig) so the
//      caller can globally intern at the end.

#include <atomic>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "bucket_io.hpp"
#include "color_set_dict.hpp"
#include "concurrent_kmer_map.hpp"
#include "kmer.hpp"
#include "super_kmer.hpp"
#include "unitig_walker.hpp"

namespace cdgb {

struct BucketUnitig {
    std::string seq;
    std::vector<uint32_t> colors;  // inline (sorted, deduped) color set
};

namespace detail {

inline void process_bucket(const std::string& path, uint32_t k,
                           std::vector<BucketUnitig>& out_local) {
    BucketReader reader(path);

    // Per-k-mer color accumulator. KmerEntry already does sorted-dedupe insertion.
    std::unordered_map<kmer_int_t, KmerEntry, KmerHasher> kmer_colors;

    uint32_t color = 0;
    std::vector<uint8_t> bases;
    while (reader.next(color, bases)) {
        if (bases.size() < k) continue;
        // Slide canonical k-mers across the super-k-mer.
        kmer_int_t fwd = 0, rc = 0;
        const kmer_int_t mask = kmer_mask(k);
        const uint32_t k_minus_1_x2 = 2 * (k - 1);
        for (uint32_t i = 0; i < k - 1; ++i) {
            uint8_t v = bases[i];
            fwd = ((fwd << 2) | v) & mask;
            rc  = (rc >> 2) | ((kmer_int_t)(v ^ 3) << k_minus_1_x2);
        }
        for (uint32_t i = k - 1; i < bases.size(); ++i) {
            uint8_t v = bases[i];
            fwd = ((fwd << 2) | v) & mask;
            rc  = (rc >> 2) | ((kmer_int_t)(v ^ 3) << k_minus_1_x2);
            kmer_int_t can = fwd <= rc ? fwd : rc;
            kmer_colors[can].add_color(color);
        }
    }

    // Local color-set dict + 1-shard FinalKmerMap.
    ColorSetDict local_dict;
    FinalKmerMap fkm(0);  // 0 -> 1 shard
    auto& shard = fkm.shard(0);
    shard.map.reserve(kmer_colors.size());
    for (auto& kv : kmer_colors) {
        std::vector<uint32_t> sorted = kv.second.to_sorted();
        uint32_t cid = local_dict.intern(std::move(sorted));
        shard.map.emplace(kv.first, cid);
    }
    // Free the accumulator now.
    kmer_colors = {};

    // Walk unitigs.
    UnitigWalker w(fkm, k);
    std::vector<Unitig> unitigs;
    std::mutex mu;
    w.walk_all(1, unitigs, mu);

    // Re-attach color sets inline so the caller can globally intern later.
    out_local.reserve(out_local.size() + unitigs.size());
    for (auto& u : unitigs) {
        BucketUnitig b;
        b.seq = std::move(u.seq);
        b.colors = local_dict.at(u.color_set_id);
        out_local.push_back(std::move(b));
    }
}

}  // namespace detail

// Parallel driver. Each worker processes one bucket at a time and appends its
// unitigs to the shared output vector under `out_mu`.
inline void process_buckets(const BucketWriter& writer, uint32_t k,
                            uint32_t num_threads,
                            std::vector<BucketUnitig>& out, std::mutex& out_mu,
                            std::atomic<uint64_t>* done = nullptr) {
    if (num_threads == 0) num_threads = 1;
    const uint32_t B = writer.num_buckets();
    std::atomic<uint32_t> next{0};
    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    auto run = [&]() {
        std::vector<BucketUnitig> local;
        for (;;) {
            uint32_t b = next.fetch_add(1);
            if (b >= B) break;
            try {
                detail::process_bucket(writer.bucket_path(b), k, local);
            } catch (std::exception& e) {
                std::cerr << "error processing bucket " << b << ": " << e.what() << '\n';
            }
            if (done) done->fetch_add(1, std::memory_order_relaxed);
        }
        std::lock_guard<std::mutex> lk(out_mu);
        for (auto& u : local) out.emplace_back(std::move(u));
    };

    for (uint32_t t = 0; t < num_threads; ++t) workers.emplace_back(run);
    for (auto& w : workers) w.join();
}

}  // namespace cdgb
