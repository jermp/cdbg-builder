#pragma once

// Per-bucket de Bruijn graph construction and unitig walk, with awareness of
// "phantom" edges that cross into other buckets.
//
// Each canonical k-mer in this bucket is annotated with two bits saying
// whether there is a known dBG edge into a k-mer that lives in a *different*
// bucket. The walker treats these phantom edges as having degree 1 on the
// relevant side: if the only edge on a side is phantom, the walk stops there
// and the unitig is emitted with that end marked OPEN, ready to be stitched
// against another bucket's open end. If the side mixes a local edge and a
// phantom edge, the global degree is >= 2 and the walk stops as a real
// branch (CLOSED) — exactly what GGCAT's `try_extend_function` does.
//
// Closed unitigs go straight to the caller. Open unitigs carry their first
// and last canonical k-mers and side-flags so the cross-bucket stitch phase
// can match them up.

#include <algorithm>
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
#include "kmer.hpp"
#include "super_kmer.hpp"

namespace cdgb {

// Per-k-mer color accumulator. Stores the first two colors inline to avoid
// allocation for the common singleton/doubleton case; a 0xffffffff sentinel
// marks "unused" inline slots.
struct KmerEntry {
    static constexpr uint32_t INVALID = 0xffffffffu;
    uint32_t inline_a = INVALID;
    uint32_t inline_b = INVALID;
    std::vector<uint32_t> extra;

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

// Open-end marker, in unitig-sequence orientation.
inline constexpr uint8_t UNITIG_OPEN_LEFT = 1u << 0;
inline constexpr uint8_t UNITIG_OPEN_RIGHT = 1u << 1;

struct StitchableUnitig {
    std::string seq;               // ACGT characters
    std::vector<uint32_t> colors;  // sorted, deduped color set
    uint8_t open_flags = 0;        // bits from UNITIG_OPEN_*
};

namespace detail {

// Per-canonical-k-mer flags for phantom edges, in canonical orientation.
inline constexpr uint8_t KMER_PHANTOM_LEFT = 1u << 0;   // phantom predecessor exists
inline constexpr uint8_t KMER_PHANTOM_RIGHT = 1u << 1;  // phantom successor exists

struct BucketKmerInfo {
    KmerEntry colors;
    uint8_t phantom = 0;
};

using BucketKmerMap = std::unordered_map<kmer_int_t, BucketKmerInfo, KmerHasher>;

// Build the per-canonical-k-mer info from a bucket's super-k-mer stream.
inline void load_bucket(const std::string& path, uint32_t k, BucketKmerMap& out) {
    BucketReader reader(path);
    uint8_t flags = 0;
    uint32_t color = 0;
    std::vector<uint8_t> bases;
    while (reader.next(flags, color, bases)) {
        if (bases.size() < k) continue;

        const kmer_int_t mask = kmer_mask(k);
        const uint32_t k_minus_1_x2 = 2 * (k - 1);

        kmer_int_t fwd = 0, rc = 0;
        for (uint32_t i = 0; i < k - 1; ++i) {
            uint8_t v = bases[i];
            fwd = ((fwd << 2) | v) & mask;
            rc = (rc >> 2) | ((kmer_int_t)(v ^ 3) << k_minus_1_x2);
        }

        kmer_int_t first_can = 0, last_can = 0;
        bool first_is_fwd = true, last_is_fwd = true;
        bool have_first = false;

        for (uint32_t i = k - 1; i < bases.size(); ++i) {
            uint8_t v = bases[i];
            fwd = ((fwd << 2) | v) & mask;
            rc = (rc >> 2) | ((kmer_int_t)(v ^ 3) << k_minus_1_x2);
            bool is_fwd = (fwd <= rc);
            kmer_int_t can = is_fwd ? fwd : rc;
            out[can].colors.add_color(color);
            if (!have_first) {
                first_can = can;
                first_is_fwd = is_fwd;
                have_first = true;
            }
            last_can = can;
            last_is_fwd = is_fwd;
        }
        if (!have_first) continue;

        // Boundary phantom edges. The super-k-mer's first k-mer (in input
        // order) has a predecessor in another bucket iff IS_ACGT_BEGIN is
        // not set. Likewise the last k-mer for IS_ACGT_END.
        if ((flags & SK_FLAG_IS_ACGT_BEGIN) == 0) {
            // The phantom edge enters first_kmer's input "back" side.
            // In canonical orientation: input-back maps to canonical-left if
            // input == canonical, else canonical-right.
            uint8_t bit = first_is_fwd ? KMER_PHANTOM_LEFT : KMER_PHANTOM_RIGHT;
            out[first_can].phantom |= bit;
        }
        if ((flags & SK_FLAG_IS_ACGT_END) == 0) {
            // Phantom edge exits last_kmer on its input "fwd" side.
            uint8_t bit = last_is_fwd ? KMER_PHANTOM_RIGHT : KMER_PHANTOM_LEFT;
            out[last_can].phantom |= bit;
        }
    }
}

// 8-bit local extension mask, bits 0..3 = forward-side successors of `can`,
// bits 4..7 = back-side predecessors (= forward-side successors of rc(can)).
inline uint8_t local_ext_mask(kmer_int_t can, uint32_t k, const BucketKmerMap& m) {
    uint8_t out = 0;
    kmer_int_t fwd = can;
    kmer_int_t rev = reverse_complement(can, k);
    for (uint8_t nt = 0; nt < 4; ++nt) {
        kmer_int_t ext_f = shift_append(fwd, nt, k);
        if (m.find(canonical(ext_f, k)) != m.end()) out |= (uint8_t)(1u << nt);
        kmer_int_t ext_r = shift_append(rev, nt, k);
        if (m.find(canonical(ext_r, k)) != m.end()) out |= (uint8_t)(1u << (4 + nt));
    }
    return out;
}

// In walking orientation (rc=false means we walk through the canonical
// forward), forward-side of `m` is its low nibble; flipped under rc.
inline uint8_t fwd_nibble(uint8_t m, bool rc) { return (uint8_t)((rc ? (m >> 4) : m) & 0xf); }
inline uint8_t back_nibble(uint8_t m, bool rc) { return (uint8_t)((rc ? m : (m >> 4)) & 0xf); }

// Phantom on the forward / back side in walking orientation.
inline bool phantom_on_fwd(uint8_t phantom, bool rc) {
    return rc ? (phantom & KMER_PHANTOM_LEFT) : (phantom & KMER_PHANTOM_RIGHT);
}
inline bool phantom_on_back(uint8_t phantom, bool rc) {
    return rc ? (phantom & KMER_PHANTOM_RIGHT) : (phantom & KMER_PHANTOM_LEFT);
}

struct BStep {
    kmer_int_t next_can;
    bool next_rc;
};
inline BStep bstep(kmer_int_t can, bool rc, uint32_t k, uint8_t nt) {
    kmer_int_t cur = rc ? reverse_complement(can, k) : can;
    kmer_int_t next_fwd = shift_append(cur, nt, k);
    kmer_int_t next_rev = reverse_complement(next_fwd, k);
    BStep r;
    if (next_fwd <= next_rev) {
        r.next_can = next_fwd;
        r.next_rc = false;
    } else {
        r.next_can = next_rev;
        r.next_rc = true;
    }
    return r;
}

// Returns true iff (can, rc) is the start of a maximal monochromatic
// unitig (no straight-line predecessor in walking orientation, considering
// both local and phantom edges + color-set agreement). Also reports whether
// the back-side reason for being a left-end was a phantom-only edge — that
// makes the resulting unitig OPEN_LEFT.
struct LeftEndCheck {
    bool is_left_end;
    bool back_is_phantom_only;  // true => unitig will be open-left
};
inline LeftEndCheck classify_left_end(
    kmer_int_t can, bool rc, uint32_t cid, uint32_t k, const BucketKmerMap& m,
    const std::unordered_map<kmer_int_t, uint32_t, KmerHasher>& cid_of) {
    auto it = m.find(can);
    uint8_t phantom = it->second.phantom;
    uint8_t mask = local_ext_mask(can, k, m);
    uint8_t back = back_nibble(mask, rc);
    int local_back_n = __builtin_popcount(back);
    int phantom_back_n = phantom_on_back(phantom, rc) ? 1 : 0;
    int total_back = local_back_n + phantom_back_n;
    if (total_back != 1) return {true, false};
    if (local_back_n == 0) {
        // Only edge is phantom: this is a left-end and the unitig is open.
        return {true, true};
    }
    // total_back == 1 and edge is local. Walk to predecessor and check that
    // *its* forward in walking orientation is also degree 1 and color matches.
    uint8_t nt = (uint8_t)__builtin_ctz(back);
    auto sr = bstep(can, !rc, k, nt);
    kmer_int_t pred_can = sr.next_can;
    bool pred_rc = !sr.next_rc;
    auto pit = m.find(pred_can);
    if (pit == m.end()) return {true, false};  // shouldn't happen: local_back_n said it does
    uint8_t pred_phantom = pit->second.phantom;
    uint8_t pred_mask = local_ext_mask(pred_can, k, m);
    int pred_local_fwd = __builtin_popcount(fwd_nibble(pred_mask, pred_rc));
    int pred_phantom_fwd = phantom_on_fwd(pred_phantom, pred_rc) ? 1 : 0;
    if (pred_local_fwd + pred_phantom_fwd != 1) return {true, false};
    auto cit = cid_of.find(pred_can);
    if (cit == cid_of.end() || cit->second != cid) return {true, false};
    return {false, false};
}

inline void process_bucket(const std::string& path, uint32_t k,
                           std::vector<StitchableUnitig>& out_local) {
    BucketKmerMap kmer_info;
    load_bucket(path, k, kmer_info);

    // Build local color-set dict and a parallel map cid_of[can] for quick lookups.
    ColorSetDict local_dict;
    std::unordered_map<kmer_int_t, uint32_t, KmerHasher> cid_of;
    cid_of.reserve(kmer_info.size());
    for (auto& kv : kmer_info) {
        std::vector<uint32_t> sorted = kv.second.colors.to_sorted();
        uint32_t cid = local_dict.intern(std::move(sorted));
        cid_of.emplace(kv.first, cid);
    }
    // Free per-k-mer color storage; we keep phantom flags in kmer_info.
    for (auto& kv : kmer_info) kv.second.colors = KmerEntry{};

    std::unordered_map<kmer_int_t, uint8_t, KmerHasher> visited;
    visited.reserve(kmer_info.size());

    auto extend_and_emit = [&](kmer_int_t start_can, bool start_rc, uint32_t cid, bool open_left) {
        StitchableUnitig u;
        kmer_int_t cur = start_rc ? reverse_complement(start_can, k) : start_can;
        u.seq = kmer_to_string(cur, k);
        u.colors = local_dict.at(cid);
        if (open_left) u.open_flags |= UNITIG_OPEN_LEFT;

        kmer_int_t can = start_can;
        bool rc = start_rc;
        visited[can] = 1;

        for (;;) {
            auto it = kmer_info.find(can);
            uint8_t phantom = it->second.phantom;
            uint8_t mask = local_ext_mask(can, k, kmer_info);
            uint8_t fwd_n = fwd_nibble(mask, rc);
            int local_fwd = __builtin_popcount(fwd_n);
            int phantom_fwd = phantom_on_fwd(phantom, rc) ? 1 : 0;
            if (local_fwd + phantom_fwd != 1) break;  // closed (branch / dead end)
            if (local_fwd == 0) {                     // open right
                u.open_flags |= UNITIG_OPEN_RIGHT;
                break;
            }
            uint8_t nt = (uint8_t)__builtin_ctz(fwd_n);
            auto sr = bstep(can, rc, k, nt);
            // Successor must (a) have back-side degree exactly 1 in walking
            // orientation (no incoming branch) and (b) carry the same cid.
            auto nit = kmer_info.find(sr.next_can);
            if (nit == kmer_info.end()) break;
            uint8_t next_phantom = nit->second.phantom;
            uint8_t next_mask = local_ext_mask(sr.next_can, k, kmer_info);
            int next_local_back = __builtin_popcount(back_nibble(next_mask, sr.next_rc));
            int next_phantom_back = phantom_on_back(next_phantom, sr.next_rc) ? 1 : 0;
            if (next_local_back + next_phantom_back != 1) break;
            if (cid_of[sr.next_can] != cid) break;
            if (visited.find(sr.next_can) != visited.end()) break;  // cycle closure
            visited[sr.next_can] = 1;
            u.seq.push_back(twobit_to_nuc(nt));
            can = sr.next_can;
            rc = sr.next_rc;
        }
        out_local.emplace_back(std::move(u));
    };

    // First pass: start from k-mers that look like a left-end of some unitig.
    for (auto& kv : kmer_info) {
        kmer_int_t can = kv.first;
        if (visited.find(can) != visited.end()) continue;
        uint32_t cid = cid_of[can];
        // Try both orientations.
        auto le_f = classify_left_end(can, false, cid, k, kmer_info, cid_of);
        auto le_r = classify_left_end(can, true, cid, k, kmer_info, cid_of);
        bool start_rc;
        bool open_left;
        if (le_f.is_left_end) {
            start_rc = false;
            open_left = le_f.back_is_phantom_only;
        } else if (le_r.is_left_end) {
            start_rc = true;
            open_left = le_r.back_is_phantom_only;
        } else
            continue;
        extend_and_emit(can, start_rc, cid, open_left);
    }

    // Second pass: anything still unvisited is on a pure cycle inside the
    // bucket. Pick an arbitrary k-mer to break it; both ends will be CLOSED
    // (since by definition no left-end exists, but we've torn the cycle).
    for (auto& kv : kmer_info) {
        kmer_int_t can = kv.first;
        if (visited.find(can) != visited.end()) continue;
        extend_and_emit(can, false, cid_of[can], /*open_left=*/false);
    }
}

}  // namespace detail

// Parallel driver. Each worker processes one bucket at a time and appends its
// stitchable unitigs to the shared output vector under `out_mu`.
inline void process_buckets(const BucketWriter& writer, uint32_t k, uint32_t num_threads,
                            std::vector<StitchableUnitig>& out, std::mutex& out_mu,
                            std::atomic<uint64_t>* done = nullptr) {
    if (num_threads == 0) num_threads = 1;
    const uint32_t B = writer.num_buckets();
    std::atomic<uint32_t> next{0};
    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    auto run = [&]() {
        std::vector<StitchableUnitig> local;
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
