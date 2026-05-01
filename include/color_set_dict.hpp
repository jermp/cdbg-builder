#pragma once

#include <cstdint>
#include <cstring>
#include <vector>

#include <unordered_dense/unordered_dense.h>

#include "kmer.hpp"

namespace cdgb {

// Maps a sorted, deduped color list to a small integer id. Building this
// dictionary collapses the per-k-mer "set of colors" into "color_set_id".
//
// Hashing uses wyhash over the byte view of the vector — ~4-5x faster
// than the FNV-1a we used previously, which mattered once color sets
// reached the thousands-of-entries range (e.g. 4546-genome pangenomes).
// `intern` does a single try_emplace so we hash the key once per call,
// not twice (the previous find + emplace did two hashes per insert).
struct color_set_dict {
    struct hash {
        // wyhash output is already avalanched; tell ankerl not to re-mix.
        using is_avalanching = void;
        size_t operator()(const std::vector<uint32_t>& v) const noexcept {
            return (size_t)ankerl::unordered_dense::detail::wyhash::hash(
                v.data(), v.size() * sizeof(uint32_t));
        }
    };

    uint32_t intern(std::vector<uint32_t>&& key) {
        // try_emplace: if the key is already present, no allocation, no
        // move; if not, the rvalue is moved into the map. Single hash.
        // The placeholder id is corrected after we know it's a fresh
        // insert (we can't read m_classes.size() in the value position
        // and have it match the post-insert id — that's fine, we patch
        // it_->second below).
        uint32_t pending_id = (uint32_t)m_classes.size();
        auto [it, inserted] = m_index.try_emplace(std::move(key), pending_id);
        if (inserted) {
            // The rvalue was moved INTO m_index as the key; copy it
            // into m_classes so we can later look up by id.
            m_classes.push_back(it->first);
        }
        return it->second;
    }

    uint32_t size() const { return (uint32_t)m_classes.size(); }
    const std::vector<uint32_t>& at(uint32_t id) const { return m_classes[id]; }
    const std::vector<std::vector<uint32_t>>& classes() const { return m_classes; }

private:
    std::vector<std::vector<uint32_t>> m_classes;
    ankerl::unordered_dense::map<std::vector<uint32_t>, uint32_t, hash> m_index;
};

}  // namespace cdgb
