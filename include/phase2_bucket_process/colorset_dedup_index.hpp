#pragma once

// colorset_dedup_index: the in-RAM exact dedup index for color sets.
//
// Maps a color set's 128-bit content hash (primary wyhash + secondary FNV-1a)
// to its dense cid. This is the "color-sets-dedup-map" of algorithm.md §9.1 --
// the dominant, fastest-growing, non-spillable structure in bucket-process
// (~0.5 GiB @ 20k genomes -> ~14 GiB @ 661k). It is encapsulated here as the
// clean seam for a planned RAM-first / spill-overflow externalization:
// today it is purely in-RAM and exact; the overflow path would plug in
// behind this interface.
//
// Sharded NUM_SHARDS ways (selected by the secondary hash) so the bucket-process
// local->global merge dedups across shards in parallel instead of serialising on
// one lock. Each shard maps the 128-bit hash directly to a cid -- no shared
// per-class array, so shards are fully independent.
//
// THREAD-SAFE usage: the caller locks one shard (via shard_for()) across the
// find -> (on miss) encode + cid-assign -> emplace region, so distinct shards
// proceed in parallel while a single shard's dedup stays consistent. The lock is
// the caller's because the on-miss cid assignment needs the dict's output lock
// nested inside it (see streaming_color_set_dict::intern_with_hashes).

#include <array>
#include <cstdint>
#include <mutex>

#include <unordered_dense/unordered_dense.h>

namespace cdbg {

struct colorset_dedup_index {
    // 128-bit content hash key: primary = wyhash, secondary = FNV-1a. Two
    // independent hash families -> birthday collision over the run is ~2^-64 per
    // pair, which is why we skip the byte-level equality check.
    struct key {
        uint64_t primary;
        uint64_t secondary;
    };

    struct hasher {
        using is_avalanching = void;  // primary is already avalanched
        size_t operator()(key const& k) const noexcept { return k.primary; }
    };
    struct eq {
        bool operator()(key const& a, key const& b) const noexcept {
            return a.primary == b.primary and a.secondary == b.secondary;
        }
    };

    static constexpr size_t NUM_SHARDS = 256;  // power of two

    // alignas(64): keep each shard's mutex on its own cache line so locking one
    // shard doesn't false-share with its neighbors.
    struct alignas(64) shard_t {
        std::mutex mu;
        ankerl::unordered_dense::map<key, uint64_t, hasher, eq> map;
    };

    // The shard owning a key, selected by its secondary hash. Caller locks
    // shard.mu across find/emplace.
    shard_t& shard_for(uint64_t secondary) noexcept {
        return m_shards[secondary & (NUM_SHARDS - 1)];
    }

    // Approx RAM held: NUM_SHARDS maps, each entry being the 128-bit key (16 B) +
    // the cid (8 B); ankerl's flat backing adds ~0.6x slack at the default load
    // factor. `n` is the live class count (the dict tracks it atomically).
    static uint64_t resident_bytes(uint64_t n) noexcept { return (uint64_t)(n * 24.0 * 1.6); }

    // Free all shard maps (interning is done). Idempotent.
    void release() {
        for (auto& s : m_shards) {
            decltype(s.map){}.swap(s.map);
        }
    }

private:
    std::array<shard_t, NUM_SHARDS> m_shards;
};

}  // namespace cdbg
