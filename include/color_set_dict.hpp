#pragma once

#include <cstdint>
#include <cstring>
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

}  // namespace cdgb
