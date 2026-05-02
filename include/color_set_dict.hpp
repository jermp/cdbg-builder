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
// Single-store layout: each color list is stored exactly once, in
// `m_classes`. The dedup index `m_index` is a hash *set* whose key type
// is the integer id; its hash and equality dispatch through `m_classes`
// to access the actual vector data, so we never duplicate the bytes.
//
// Lookup of a not-yet-interned vector goes through the transparent
// (heterogeneous) overloads, so `m_index.find(some_vector)` works
// without copying or moving the candidate into the index first.
//
// Pointer / reference stability of m_classes[i].data() is guaranteed:
// std::vector's move is O(1) (it transfers the data buffer pointer),
// so even if the outer vector reallocates during push_back, the data
// pointer of any inner vector stays put.
struct color_set_dict {
    color_set_dict()
        : m_index(0, hasher{&m_classes}, key_eq{&m_classes}) {}

    color_set_dict(color_set_dict const&) = delete;
    color_set_dict& operator=(color_set_dict const&) = delete;
    color_set_dict(color_set_dict&&) = delete;
    color_set_dict& operator=(color_set_dict&&) = delete;

    uint32_t intern(std::vector<uint32_t>&& key) {
        // Heterogeneous find: hashes/compares directly against `key` via
        // the hasher/key_eq overloads that take std::vector const&.
        auto it = m_index.find(key);
        if (it != m_index.end()) return *it;
        uint32_t id = (uint32_t)m_classes.size();
        m_classes.push_back(std::move(key));  // single-store: one move, no copy
        m_index.insert(id);
        return id;
    }

    uint32_t size() const { return (uint32_t)m_classes.size(); }
    std::vector<uint32_t> const& at(uint32_t id) const { return m_classes[id]; }
    std::vector<std::vector<uint32_t>> const& classes() const { return m_classes; }

    // Move-out of the i-th color list; the dict entry is left empty.
    // Used by the per-bucket merge in process_buckets to avoid copying
    // each bucket's classes into the global dict.
    std::vector<uint32_t>& mutable_at(uint32_t id) { return m_classes[id]; }

private:
    // FNV-1a over the bytes of a sorted-deduped color list. Used by both
    // the id-based and vector-based hasher overloads below.
    static size_t hash_vec(std::vector<uint32_t> const& v) noexcept {
        uint64_t h = 1469598103934665603ULL;
        unsigned char const* p = (unsigned char const*)v.data();
        size_t n = v.size() * sizeof(uint32_t);
        for (size_t i = 0; i < n; ++i) {
            h ^= p[i];
            h *= 1099511628211ULL;
        }
        return (size_t)h;
    }

    // Hasher + key_eq carry a back-pointer to m_classes so they can
    // resolve an integer id to its color list. is_transparent enables
    // heterogeneous lookup with std::vector<uint32_t>; is_avalanching
    // tells unordered_dense to skip its own mix step (FNV-1a is good
    // enough on its own).
    struct hasher {
        std::vector<std::vector<uint32_t>> const* classes;
        using is_transparent = void;
        using is_avalanching = void;
        size_t operator()(uint32_t id) const noexcept { return hash_vec((*classes)[id]); }
        size_t operator()(std::vector<uint32_t> const& v) const noexcept { return hash_vec(v); }
    };

    struct key_eq {
        std::vector<std::vector<uint32_t>> const* classes;
        using is_transparent = void;
        bool operator()(uint32_t a, uint32_t b) const noexcept {
            return (*classes)[a] == (*classes)[b];
        }
        bool operator()(uint32_t a, std::vector<uint32_t> const& v) const noexcept {
            return (*classes)[a] == v;
        }
        bool operator()(std::vector<uint32_t> const& v, uint32_t a) const noexcept {
            return (*classes)[a] == v;
        }
    };

    std::vector<std::vector<uint32_t>> m_classes;
    ankerl::unordered_dense::set<uint32_t, hasher, key_eq> m_index;
};

}  // namespace cdgb
