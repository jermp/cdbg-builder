#pragma once

#include <cstdint>
#include <cstring>
#include <memory>
#include <unordered_map>
#include <vector>

#include "kmer.hpp"

namespace cdgb {

// Maps a sorted, deduped color list to a small integer id. Building this
// dictionary collapses the per-k-mer "set of colors" into "color_set_id".
//
// Hashing/equality work on the contiguous byte view of the vector to keep the
// hot path cheap.
struct ColorSetDict {
    struct Hash {
        size_t operator()(const std::vector<uint32_t>& v) const noexcept {
            // FNV-1a over the bytes, fast and decent for small vectors.
            uint64_t h = 1469598103934665603ULL;
            const unsigned char* p = (const unsigned char*)v.data();
            size_t n = v.size() * sizeof(uint32_t);
            for (size_t i = 0; i < n; ++i) {
                h ^= p[i];
                h *= 1099511628211ULL;
            }
            return (size_t)h;
        }
    };

    uint32_t intern(std::vector<uint32_t>&& key) {
        auto it = m_index.find(key);
        if (it != m_index.end()) return it->second;
        uint32_t id = (uint32_t)m_classes.size();
        m_classes.push_back(key);
        m_index.emplace(std::move(key), id);
        return id;
    }

    uint32_t size() const { return (uint32_t)m_classes.size(); }
    const std::vector<uint32_t>& at(uint32_t id) const { return m_classes[id]; }
    const std::vector<std::vector<uint32_t>>& classes() const { return m_classes; }

private:
    std::vector<std::vector<uint32_t>> m_classes;
    std::unordered_map<std::vector<uint32_t>, uint32_t, Hash> m_index;
};

// Final per-k-mer map: canonical k-mer -> color_set_id. Sharded for cache
// friendliness during the parallel unitig walk.
struct FinalKmerMap {
    struct Shard {
        std::unordered_map<kmer_int_t, uint32_t, KmerHasher> map;
    };

    explicit FinalKmerMap(uint32_t num_shards_log2)
        : m_log2(num_shards_log2)
        , m_mask((uint64_t(1) << num_shards_log2) - 1) {
        m_shards.resize(uint64_t(1) << num_shards_log2);
        for (auto& s : m_shards) s = std::make_unique<Shard>();
    }

    uint64_t num_shards() const { return m_shards.size(); }
    Shard& shard(uint64_t i) { return *m_shards[i]; }
    const Shard& shard(uint64_t i) const { return *m_shards[i]; }
    uint64_t shard_of(kmer_int_t k) const { return KmerHasher{}(k) & m_mask; }

    // Lookup; returns UINT32_MAX if not present.
    uint32_t lookup(kmer_int_t k) const {
        const auto& s = *m_shards[KmerHasher{}(k) & m_mask];
        auto it = s.map.find(k);
        return it == s.map.end() ? 0xffffffffu : it->second;
    }

    bool contains(kmer_int_t k) const {
        const auto& s = *m_shards[KmerHasher{}(k) & m_mask];
        return s.map.find(k) != s.map.end();
    }

    uint64_t size() const {
        uint64_t n = 0;
        for (auto& s : m_shards) n += s->map.size();
        return n;
    }

private:
    uint32_t m_log2;
    uint64_t m_mask;
    std::vector<std::unique_ptr<Shard>> m_shards;
};

}  // namespace cdgb
