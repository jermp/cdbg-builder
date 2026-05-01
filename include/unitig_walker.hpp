#pragma once

#include <atomic>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "color_set_dict.hpp"
#include "kmer.hpp"

namespace cdgb {

// A monochromatic colored unitig: a maximal path in the canonical de Bruijn
// graph along which every internal node has unique forward/back extension
// AND every node shares the same color_set_id.
struct Unitig {
    std::string seq;
    uint32_t color_set_id;
};

namespace detail {

// 8-bit "extension presence" mask centered on a canonical k-mer K:
//   bits 0..3: nucleotides A,C,G,T appended to fwd(K) -> next k-mer exists?
//   bits 4..7: nucleotides A,C,G,T appended to rc(K)  -> next k-mer exists?
inline uint8_t compute_ext_mask(kmer_int_t can, uint32_t k, const FinalKmerMap& fkm) {
    uint8_t m = 0;
    kmer_int_t fwd = can;
    kmer_int_t rev = reverse_complement(can, k);
    for (uint8_t nt = 0; nt < 4; ++nt) {
        kmer_int_t ext_f = shift_append(fwd, nt, k);
        if (fkm.contains(canonical(ext_f, k))) m |= (uint8_t)(1u << nt);
        kmer_int_t ext_r = shift_append(rev, nt, k);
        if (fkm.contains(canonical(ext_r, k))) m |= (uint8_t)(1u << (4 + nt));
    }
    return m;
}

struct StepResult {
    kmer_int_t next_can;
    bool next_rc;  // orientation of next k-mer relative to its canonical
};

inline StepResult step(kmer_int_t can, bool rc, uint32_t k, uint8_t nt) {
    kmer_int_t cur = rc ? reverse_complement(can, k) : can;
    kmer_int_t next_fwd = shift_append(cur, nt, k);
    kmer_int_t next_rev = reverse_complement(next_fwd, k);
    StepResult r;
    if (next_fwd <= next_rev) {
        r.next_can = next_fwd;
        r.next_rc = false;
    } else {
        r.next_can = next_rev;
        r.next_rc = true;
    }
    return r;
}

}  // namespace detail

struct UnitigWalker {
    UnitigWalker(const FinalKmerMap& fkm, uint32_t k) : m_fkm(fkm), m_k(k) {}

    void walk_all(uint32_t num_threads, std::vector<Unitig>& out, std::mutex& out_mu,
                  std::atomic<uint64_t>* done = nullptr) {
        if (num_threads == 0) num_threads = 1;
        const uint64_t S = m_fkm.num_shards();

        m_shard_keys.resize(S);
        m_visited.resize(S);
        m_pos.resize(S);

        // Build shard_keys, visited flags, and a kmer->index lookup per shard,
        // all in parallel — shards are disjoint so there is no contention.
        std::atomic<uint64_t> next_init{0};
        auto init_worker = [&]() {
            for (;;) {
                uint64_t i = next_init.fetch_add(1);
                if (i >= S) return;
                const auto& shard = m_fkm.shard(i);
                auto& ks = m_shard_keys[i];
                ks.reserve(shard.map.size());
                m_pos[i].reserve(shard.map.size() * 2);
                size_t idx = 0;
                for (auto& kv : shard.map) {
                    ks.push_back(kv.first);
                    m_pos[i].emplace(kv.first, idx);
                    ++idx;
                }
                m_visited[i] =
                    std::unique_ptr<std::atomic<uint8_t>[]>(new std::atomic<uint8_t>[ks.size()]);
                for (size_t j = 0; j < ks.size(); ++j) m_visited[i][j].store(0);
            }
        };
        {
            std::vector<std::thread> ts;
            for (uint32_t t = 0; t < num_threads; ++t) ts.emplace_back(init_worker);
            for (auto& t : ts) t.join();
        }

        // Walk shards in parallel.
        std::atomic<uint64_t> next_shard{0};
        auto worker = [&]() {
            std::vector<Unitig> local;
            for (;;) {
                uint64_t s = next_shard.fetch_add(1);
                if (s >= S) break;
                walk_shard(s, local);
                if (done) done->fetch_add(1, std::memory_order_relaxed);
            }
            std::lock_guard<std::mutex> lk(out_mu);
            for (auto& u : local) out.emplace_back(std::move(u));
        };
        std::vector<std::thread> ts;
        for (uint32_t t = 0; t < num_threads; ++t) ts.emplace_back(worker);
        for (auto& t : ts) t.join();
    }

private:
    const FinalKmerMap& m_fkm;
    uint32_t m_k;
    std::vector<std::vector<kmer_int_t>> m_shard_keys;
    std::vector<std::unique_ptr<std::atomic<uint8_t>[]>> m_visited;
    std::vector<std::unordered_map<kmer_int_t, size_t, KmerHasher>> m_pos;

    bool claim(uint64_t s, size_t idx) {
        uint8_t expected = 0;
        return m_visited[s][idx].compare_exchange_strong(expected, 1, std::memory_order_acq_rel);
    }

    bool locate(kmer_int_t can, uint64_t& s_out, size_t& idx_out) const {
        uint64_t s = m_fkm.shard_of(can);
        const auto& pos = m_pos[s];
        auto it = pos.find(can);
        if (it == pos.end()) return false;
        s_out = s;
        idx_out = it->second;
        return true;
    }

    uint32_t color_id(kmer_int_t can) const { return m_fkm.lookup(can); }

    uint8_t fwd_side(uint8_t mask, bool rc) const {
        return (uint8_t)((rc ? (mask >> 4) : mask) & 0xf);
    }
    uint8_t back_side(uint8_t mask, bool rc) const {
        // back-extensions in orientation rc = fwd-extensions in orientation !rc
        return (uint8_t)((rc ? mask : (mask >> 4)) & 0xf);
    }

    // Filter `bits` (a 4-bit fwd or back nibble computed from compute_ext_mask
    // in walking orientation) to only count neighbors whose color_set_id
    // equals `cid`. The neighbor reached by base `nt` from `(can, rc)` along
    // direction `is_back` is looked up via step.
    int colored_count(uint8_t bits, kmer_int_t can, bool rc, uint32_t cid,
                      bool is_back) const {
        int count = 0;
        for (uint8_t nt = 0; nt < 4; ++nt) {
            if (!(bits & (1u << nt))) continue;
            auto sr = detail::step(can, is_back ? !rc : rc, m_k, nt);
            if (color_id(sr.next_can) == cid) ++count;
        }
        return count;
    }

    bool is_left_end(kmer_int_t can, bool rc, uint32_t cid) const {
        uint8_t m = detail::compute_ext_mask(can, m_k, m_fkm);
        uint8_t back = back_side(m, rc);
        // Branches in the colored dBG only matter when the branching
        // neighbors share our color set; predecessors with different colors
        // are part of a different unitig and don't fragment ours.
        int back_colored = colored_count(back, can, rc, cid, /*is_back=*/true);
        if (back_colored != 1) return true;
        uint8_t nt = 0;
        for (uint8_t b = 0; b < 4; ++b) {
            if (!(back & (1u << b))) continue;
            auto sr = detail::step(can, !rc, m_k, b);
            if (color_id(sr.next_can) == cid) { nt = b; break; }
        }
        auto sr = detail::step(can, !rc, m_k, nt);
        kmer_int_t pred = sr.next_can;
        bool pred_rc = !sr.next_rc;
        uint8_t mp = detail::compute_ext_mask(pred, m_k, m_fkm);
        uint8_t pf = fwd_side(mp, pred_rc);
        int pf_colored = colored_count(pf, pred, pred_rc, cid, /*is_back=*/false);
        if (pf_colored != 1) return true;
        return false;
    }

    void walk_shard(uint64_t s, std::vector<Unitig>& out) {
        const auto& ks = m_shard_keys[s];
        for (size_t i = 0; i < ks.size(); ++i) {
            kmer_int_t can = ks[i];
            uint32_t cid = color_id(can);
            bool start_rc;
            if (is_left_end(can, false, cid))
                start_rc = false;
            else if (is_left_end(can, true, cid))
                start_rc = true;
            else
                continue;
            if (!claim(s, i)) continue;
            extend_and_emit(can, start_rc, cid, out);
        }
        // Pure monochromatic cycles: no left-end exists. Pick any unvisited
        // k-mer and break the cycle there.
        for (size_t i = 0; i < ks.size(); ++i) {
            if (m_visited[s][i].load(std::memory_order_relaxed)) continue;
            if (!claim(s, i)) continue;
            extend_and_emit(ks[i], false, color_id(ks[i]), out);
        }
    }

    void extend_and_emit(kmer_int_t start_can, bool start_rc, uint32_t cid,
                         std::vector<Unitig>& out) {
        Unitig u;
        u.color_set_id = cid;
        kmer_int_t cur = start_rc ? reverse_complement(start_can, m_k) : start_can;
        u.seq = kmer_to_string(cur, m_k);

        kmer_int_t can = start_can;
        bool rc = start_rc;
        for (;;) {
            uint8_t m = detail::compute_ext_mask(can, m_k, m_fkm);
            uint8_t fs = fwd_side(m, rc);
            int fs_colored = colored_count(fs, can, rc, cid, /*is_back=*/false);
            if (fs_colored != 1) break;
            uint8_t nt = 0;
            for (uint8_t b = 0; b < 4; ++b) {
                if (!(fs & (1u << b))) continue;
                auto sr = detail::step(can, rc, m_k, b);
                if (color_id(sr.next_can) == cid) { nt = b; break; }
            }
            auto sr = detail::step(can, rc, m_k, nt);
            uint32_t next_cid = color_id(sr.next_can);
            if (next_cid != cid) break;
            uint8_t mn = detail::compute_ext_mask(sr.next_can, m_k, m_fkm);
            uint8_t bn = back_side(mn, sr.next_rc);
            int bn_colored = colored_count(bn, sr.next_can, sr.next_rc, cid, /*is_back=*/true);
            if (bn_colored != 1) break;
            uint64_t s2;
            size_t idx2;
            if (!locate(sr.next_can, s2, idx2)) break;
            if (!claim(s2, idx2)) break;  // cycle closure
            u.seq.push_back(twobit_to_nuc(nt));
            can = sr.next_can;
            rc = sr.next_rc;
        }
        out.emplace_back(std::move(u));
    }
};

}  // namespace cdgb
