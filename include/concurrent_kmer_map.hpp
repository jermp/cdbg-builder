#pragma once

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <vector>

#include "kmer.hpp"

namespace cdgb {

// Sharded concurrent hashmap from canonical k-mer to a per-k-mer payload.
//
// The payload accumulates the colors that contain this k-mer. Colors are
// inserted as they are discovered; they end up sorted+deduped in `finalize`.
// We use a small-vector representation that stays compact for the common case
// of a k-mer appearing in only a handful of colors, and falls back to a heap
// vector once the count grows.
struct KmerEntry {
    // First two colors are stored inline to avoid allocation for the common
    // case (a singleton or doubleton color set). After that, we spill into
    // `extra`. A 0xffffffff sentinel marks "unused" inline slots.
    uint32_t inline_a = INVALID;
    uint32_t inline_b = INVALID;
    std::vector<uint32_t> extra;  // empty until 3rd distinct color arrives.

    static constexpr uint32_t INVALID = 0xffffffffu;

    // Insert a color; idempotent.
    void add_color(uint32_t c) {
        if (inline_a == INVALID) {
            inline_a = c;
            return;
        }
        if (inline_a == c) return;
        if (inline_b == INVALID) {
            inline_b = c;
            return;
        }
        if (inline_b == c) return;
        for (uint32_t x : extra)
            if (x == c) return;
        extra.push_back(c);
    }

    // Materialize the (sorted, deduped) list of colors.
    std::vector<uint32_t> to_sorted() const {
        std::vector<uint32_t> v;
        v.reserve(2 + extra.size());
        if (inline_a != INVALID) v.push_back(inline_a);
        if (inline_b != INVALID) v.push_back(inline_b);
        for (uint32_t x : extra) v.push_back(x);
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
        return v;
    }
};

struct ConcurrentKmerMap {
    // We pick a power-of-two number of shards. Lock contention is low because
    // each insertion touches one shard and we fan out by hash bits.
    struct Shard {
        std::mutex mu;
        std::unordered_map<kmer_int_t, KmerEntry, KmerHasher> map;
    };

    explicit ConcurrentKmerMap(uint32_t num_shards_log2 = 8)
        : m_mask((uint64_t(1) << num_shards_log2) - 1) {
        m_shards.resize(uint64_t(1) << num_shards_log2);
        for (auto& s : m_shards) s = std::make_unique<Shard>();
    }

    uint64_t num_shards() const { return m_shards.size(); }

    // Insert (k-mer, color). Thread-safe.
    void insert(kmer_int_t key, uint32_t color) {
        size_t h = KmerHasher{}(key);
        uint64_t shard = h & m_mask;
        Shard& s = *m_shards[shard];
        std::lock_guard<std::mutex> lk(s.mu);
        auto& e = s.map[key];
        e.add_color(color);
    }

    Shard& shard(uint64_t i) { return *m_shards[i]; }
    const Shard& shard(uint64_t i) const { return *m_shards[i]; }

    uint64_t shard_of(kmer_int_t key) const { return KmerHasher{}(key)&m_mask; }

    uint64_t num_kmers() const {
        uint64_t n = 0;
        for (auto& s : m_shards) n += s->map.size();
        return n;
    }

private:
    uint64_t m_mask;
    std::vector<std::unique_ptr<Shard>> m_shards;
};

}  // namespace cdgb
