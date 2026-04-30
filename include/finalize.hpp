#pragma once

#include <atomic>
#include <thread>
#include <vector>

#include "color_set_dict.hpp"
#include "concurrent_kmer_map.hpp"

namespace cdgb {

// Fold a per-k-mer set of colors into a global color-set dictionary, then
// drop the colors from the per-k-mer payload (replacing them with a single
// `color_set_id`).
//
// Strategy:
//   1. Each worker thread takes a range of shards, builds a *local*
//      ColorSetDict and a *local* FinalKmerMap shard, recording per-shard
//      local ids.
//   2. Sequentially merge each local dict into the global ColorSetDict,
//      recording a local_id -> global_id translation per thread.
//   3. Each worker rewrites its assigned FinalKmerMap shards to use global
//      ids.
inline void finalize(ConcurrentKmerMap& src, FinalKmerMap& dst, ColorSetDict& global_dict,
                     uint32_t num_threads) {
    if (num_threads == 0) num_threads = 1;
    const uint64_t S = src.num_shards();

    struct ThreadState {
        ColorSetDict local_dict;
        // For each shard owned by this thread: temp map kmer -> local_id.
        std::vector<std::vector<std::pair<kmer_int_t, uint32_t>>> shard_records;
        std::vector<uint64_t> shard_indices;  // which shards this thread owned
    };
    std::vector<ThreadState> states(num_threads);

    // 1) parallel per-shard local interning.
    std::atomic<uint64_t> next_shard{0};
    auto worker_local = [&](uint32_t tid) {
        auto& st = states[tid];
        for (;;) {
            uint64_t i = next_shard.fetch_add(1);
            if (i >= S) return;
            st.shard_indices.push_back(i);
            st.shard_records.emplace_back();
            auto& records = st.shard_records.back();
            auto& shard = src.shard(i);
            records.reserve(shard.map.size());
            for (auto& kv : shard.map) {
                std::vector<uint32_t> sorted = kv.second.to_sorted();
                uint32_t lid = st.local_dict.intern(std::move(sorted));
                records.emplace_back(kv.first, lid);
            }
            // Free the source shard's memory now that we've extracted what we need.
            std::unordered_map<kmer_int_t, KmerEntry, KmerHasher> tmp;
            shard.map.swap(tmp);
        }
    };
    {
        std::vector<std::thread> ts;
        for (uint32_t t = 0; t < num_threads; ++t) ts.emplace_back(worker_local, t);
        for (auto& t : ts) t.join();
    }

    // 2) sequential merge of local dicts -> global, build local-to-global translation.
    std::vector<std::vector<uint32_t>> local_to_global(num_threads);
    for (uint32_t t = 0; t < num_threads; ++t) {
        const auto& classes = states[t].local_dict.classes();
        local_to_global[t].resize(classes.size());
        for (uint32_t lid = 0; lid < classes.size(); ++lid) {
            std::vector<uint32_t> copy = classes[lid];
            uint32_t gid = global_dict.intern(std::move(copy));
            local_to_global[t][lid] = gid;
        }
    }

    // 3) parallel rewrite into FinalKmerMap.
    std::atomic<uint32_t> next_tid{0};
    auto worker_emit = [&]() {
        for (;;) {
            uint32_t t = next_tid.fetch_add(1);
            if (t >= num_threads) return;
            auto& st = states[t];
            const auto& xlate = local_to_global[t];
            for (size_t k = 0; k < st.shard_records.size(); ++k) {
                uint64_t shard_idx = st.shard_indices[k];
                auto& records = st.shard_records[k];
                auto& dst_shard = dst.shard(shard_idx);
                dst_shard.map.reserve(records.size());
                for (auto& p : records) dst_shard.map.emplace(p.first, xlate[p.second]);
                std::vector<std::pair<kmer_int_t, uint32_t>>().swap(records);
            }
        }
    };
    {
        std::vector<std::thread> ts;
        for (uint32_t t = 0; t < num_threads; ++t) ts.emplace_back(worker_emit);
        for (auto& t : ts) t.join();
    }
}

}  // namespace cdgb
