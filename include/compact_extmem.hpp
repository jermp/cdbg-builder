#pragma once

// GGCAT-style id-only cross-bucket compaction + single base/color assembly.
// Companion to stitch_extmem.hpp; see parallel-stitch-plan.md.
//
// Same doubling algorithm as stitch_extmem (boundary-k-mer keying, colorless
// join, structural branch separation, O(log L) rounds), but the doubling loop
// carries only fragment-ID chains -- NO bases and NO colors. A final pass then
// assembles each closed chain's bases+colors exactly ONCE (like GGCAT's
// build_unitigs), which is just a fold of the SAME ext_concat_runs + seq-append
// the rounds use today. This removes the per-round base/color shuffle (LZ4 of
// growing sequences every round) that dominates stitch_extmem.
//
// STEP 1 (this file): in-RAM round store + in-RAM assembly, so it is directly
// comparable to stitch_unitigs_extmem in test_stitch (same Source/Sink). The
// file-backed store and the scalable (sort/bucket) assembly come next; the
// id-only doubling and the assembly fold validated here are reused unchanged.

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include <lz4.h>
#include <unordered_dense/unordered_dense.h>

#include "kmer.hpp"           // kmer_int_t, reverse_complement, kmer_hasher, nuc_to_2bit
#include "stitch.hpp"         // SIDE_*, UNITIG_OPEN_*, revcomp_string
#include "stitch_extmem.hpp"  // ext_mix_bit/ext_rng_next/ext_pair_compatible/ext_end,
                              // ext_concat_runs/ext_reverse_runs/ext_split_monochromatic

namespace cdbg {
namespace detail {

// --- chain entries -----------------------------------------------------------
// An entry packs a fragment id with its orientation: (frag_id << 1) | rc.
// rc = 1 means the fragment is reverse-complemented when laid into the chain's
// LEFT->RIGHT frame.
inline uint64_t id_entry(uint64_t frag_id, bool rc) { return (frag_id << 1) | (rc ? 1ull : 0ull); }
inline uint64_t id_entry_frag(uint64_t e) { return e >> 1; }
inline bool id_entry_rc(uint64_t e) { return (e & 1ull) != 0; }

// Reverse a chain: flip element order AND each fragment's orientation bit.
inline void id_reverse_entries(std::vector<uint64_t>& e) {
    std::reverse(e.begin(), e.end());
    for (auto& x : e) x ^= 1ull;
}

// --- the id-only tig ---------------------------------------------------------
// A partial chain with up to two open ends. kl/kr are the boundary k-mers read
// FORWARD in the current chain frame at the LEFT/RIGHT end (valid iff that end
// is open). We store the forward-frame value (not the canonical) and derive the
// canonical + is_canonical_fwd on demand, so palindromic boundary k-mers (only
// possible for even k) stay correct under reversal -- unlike a stored-canonical
// flag, which would mis-flip on a palindrome.
struct id_tig {
    kmer_int_t kl = 0, kr = 0;
    uint8_t open_flags = 0;
    uint64_t rng = 0;
    std::vector<uint64_t> entries;  // LEFT->RIGHT in the current frame
};

// Deterministic per-fragment rng seed. Unlike stitch_extmem (which seeds from
// sequence content because it has the bases), the compaction loop has only ids,
// so we seed from frag_id. Only requirement: distinct both-open tigs get
// distinct streams so they don't flip their presented end in lockstep -- this is
// a convergence (round-count) property, not a correctness one, so a per-id mix
// is at least as good as per-content. (splitmix64 finalizer, never zero.)
inline uint64_t id_rng_seed(uint64_t frag_id) {
    uint64_t h = frag_id + 0x9e3779b97f4a7c15ull;
    h ^= h >> 30;
    h *= 0xbf58476d1ce4e5b9ull;
    h ^= h >> 27;
    h *= 0x94d049bb133111ebull;
    h ^= h >> 31;
    return h ? h : 0x9e3779b97f4a7c15ull;
}

// Boundary k-mer read FORWARD in `seq`'s frame at the given side (NOT canonical).
inline kmer_int_t id_fwd_kmer(std::string_view seq, uint32_t k, uint8_t side) {
    char const* p = (side == SIDE_LEFT) ? seq.data() : seq.data() + (seq.size() - k);
    kmer_int_t v = 0;
    for (uint32_t i = 0; i < k; ++i) v = (v << 2) | (kmer_int_t)nuc_to_2bit(p[i]);
    return v;
}

// canonical(kf) + whether the forward-frame value IS the canonical one.
inline kmer_int_t id_canon(kmer_int_t kf, uint32_t k, bool& fwd) {
    kmer_int_t rc = reverse_complement(kf, k);
    if (kf <= rc) {
        fwd = true;
        return kf;
    }
    fwd = false;
    return rc;
}

// Choose the presented open side from open_flags + rng (mirrors ext_choose_side,
// which is itself base-agnostic but takes an ext_tig).
inline bool id_choose_side(uint8_t open_flags, uint64_t rng, uint8_t& side) {
    bool l = (open_flags & UNITIG_OPEN_LEFT) != 0;
    bool r = (open_flags & UNITIG_OPEN_RIGHT) != 0;
    if (l && r) {
        side = ext_mix_bit(rng) ? SIDE_RIGHT : SIDE_LEFT;
        return true;
    }
    if (l) {
        side = SIDE_LEFT;
        return true;
    }
    if (r) {
        side = SIDE_RIGHT;
        return true;
    }
    return false;
}

// Build the ext_end (canonical junction / side / is_canonical_fwd) for t's
// chosen side. Returns false if t is fully closed. Reuses ext_end +
// ext_pair_compatible from stitch_extmem unchanged (both are base-free).
inline bool id_end(id_tig const& t, uint32_t k, ext_end& out) {
    uint8_t side;
    if (!id_choose_side(t.open_flags, t.rng, side)) return false;
    kmer_int_t kf = (side == SIDE_LEFT) ? t.kl : t.kr;
    bool fwd;
    kmer_int_t j = id_canon(kf, k, fwd);
    out = ext_end{0, j, side, fwd};
    return true;
}

// Join two id-tigs at a shared, compatible boundary k-mer. `a` is oriented with
// its keyed end at the RIGHT and `b` with its keyed end at the LEFT, then the
// chains concatenate. Mirrors ext_join's end/orientation bookkeeping exactly
// (see its comments) but on entries + end-kmers instead of seq + runs. NO entry
// is dropped: both fragments keep their shared boundary k-mer; the k-base
// overlap is removed once, at assembly.
inline id_tig id_join(id_tig const& a, uint8_t a_side, id_tig const& b, uint8_t b_side,
                      uint32_t k) {
    id_tig m;
    // a oriented keyed-end-RIGHT: merged LEFT = a's other (non-keyed) end.
    std::vector<uint64_t> ae = a.entries;
    bool left_open;
    if (a_side == SIDE_RIGHT) {  // frame unchanged
        m.kl = a.kl;
        left_open = (a.open_flags & UNITIG_OPEN_LEFT) != 0;
    } else {  // a reversed
        id_reverse_entries(ae);
        m.kl = reverse_complement(a.kr, k);
        left_open = (a.open_flags & UNITIG_OPEN_RIGHT) != 0;
    }
    // b oriented keyed-end-LEFT: merged RIGHT = b's other end.
    std::vector<uint64_t> be = b.entries;
    bool right_open;
    if (b_side == SIDE_LEFT) {  // frame unchanged
        m.kr = b.kr;
        right_open = (b.open_flags & UNITIG_OPEN_RIGHT) != 0;
    } else {  // b reversed
        id_reverse_entries(be);
        m.kr = reverse_complement(b.kl, k);
        right_open = (b.open_flags & UNITIG_OPEN_LEFT) != 0;
    }
    m.entries = std::move(ae);
    m.entries.insert(m.entries.end(), be.begin(), be.end());
    m.open_flags = (uint8_t)((left_open ? UNITIG_OPEN_LEFT : 0) |
                             (right_open ? UNITIG_OPEN_RIGHT : 0));
    m.rng = a.rng ^ (b.rng * 0x9e3779b97f4a7c15ull);
    if (m.rng == 0) m.rng = 0x9e3779b97f4a7c15ull;
    return m;
}

// --- in-RAM round store (Step 1) --------------------------------------------
// Same interface as round_store_mem (emit / take_input_bucket / advance), holding
// id_tigs instead of ext_tigs.
struct id_round_store_mem {
    explicit id_round_store_mem(uint32_t num_buckets)
        : m_in(num_buckets), m_out(num_buckets), m_locks(num_buckets),
          m_num_buckets(num_buckets) {}

    uint32_t num_buckets() const { return m_num_buckets; }

    void emit(uint32_t b, id_tig&& t) {
        std::lock_guard<std::mutex> lk(m_locks[b]);
        m_out[b].push_back(std::move(t));
    }

    std::vector<id_tig> take_input_bucket(uint32_t b) {
        std::vector<id_tig> v = std::move(m_in[b]);
        m_in[b].clear();
        return v;
    }

    void advance() {
        m_in.swap(m_out);
        for (auto& v : m_out) v.clear();
    }

private:
    std::vector<std::vector<id_tig>> m_in, m_out;
    std::vector<std::mutex> m_locks;
    uint32_t m_num_buckets;
};

// --- id-tig (de)serialization for the file store ----------------------------
// Layout: [u8 open_flags][u64 rng][kmer kl][kmer kr][u32 n_entries][u64 entries...].
// Mirrors ext_tig_serialize but carries ids + end-kmers instead of seq + runs.
inline void id_tig_serialize(id_tig const& t, std::vector<uint8_t>& out) {
    auto put = [&](void const* p, size_t n) {
        uint8_t const* b = (uint8_t const*)p;
        out.insert(out.end(), b, b + n);
    };
    put(&t.open_flags, 1);
    put(&t.rng, sizeof(t.rng));
    put(&t.kl, sizeof(t.kl));
    put(&t.kr, sizeof(t.kr));
    uint32_t n = (uint32_t)t.entries.size();
    put(&n, sizeof(n));
    if (n) put(t.entries.data(), (size_t)n * sizeof(uint64_t));
}

// Returns bytes consumed, 0 on malformed/EOF.
inline size_t id_tig_deserialize(uint8_t const* buf, size_t buf_len, id_tig& t) {
    constexpr size_t HDR = 1 + sizeof(uint64_t) + 2 * sizeof(kmer_int_t) + sizeof(uint32_t);
    if (buf_len < HDR) return 0;
    size_t p = 0;
    t.open_flags = buf[p];
    p += 1;
    std::memcpy(&t.rng, buf + p, sizeof(t.rng));
    p += sizeof(t.rng);
    std::memcpy(&t.kl, buf + p, sizeof(t.kl));
    p += sizeof(t.kl);
    std::memcpy(&t.kr, buf + p, sizeof(t.kr));
    p += sizeof(t.kr);
    uint32_t n;
    std::memcpy(&n, buf + p, sizeof(n));
    p += sizeof(n);
    if (buf_len < p + (size_t)n * sizeof(uint64_t)) return 0;
    t.entries.resize(n);
    if (n) std::memcpy(t.entries.data(), buf + p, (size_t)n * sizeof(uint64_t));
    p += (size_t)n * sizeof(uint64_t);
    return p;
}

// File-backed id-tig round store: one LZ4-framed file per bucket per round, one
// bucket resident at a time. Clone of round_store_file (stitch_extmem.hpp) with
// the id_tig codec; see there for the frame format and concurrency notes.
class id_round_store_file {
public:
    id_round_store_file(std::string dir, uint32_t num_buckets)
        : m_dir(std::move(dir)), m_num_buckets(num_buckets), m_batch(num_buckets),
          m_files(num_buckets, nullptr), m_locks(num_buckets) {}

    ~id_round_store_file() {
        for (auto* f : m_files)
            if (f) std::fclose(f);
    }

    id_round_store_file(id_round_store_file const&) = delete;
    id_round_store_file& operator=(id_round_store_file const&) = delete;

    uint32_t num_buckets() const { return m_num_buckets; }

    void emit(uint32_t b, id_tig&& t) {
        std::lock_guard<std::mutex> lk(m_locks[b]);
        id_tig_serialize(t, m_batch[b]);
        if (m_batch[b].size() >= FRAME_BUDGET) flush_frame(b);
    }

    std::vector<id_tig> take_input_bucket(uint32_t b) {
        std::vector<id_tig> out;
        std::string path = bucket_path(m_in_round, b);
        std::FILE* f = std::fopen(path.c_str(), "rb");
        if (!f) return out;  // empty bucket: no file was created
        std::vector<uint8_t> comp, raw;
        for (;;) {
            uint32_t u = 0;
            if (std::fread(&u, sizeof(u), 1, f) != 1) break;
            if (u == 0) break;
            uint32_t c = 0;
            if (std::fread(&c, sizeof(c), 1, f) != 1) break;
            if (comp.size() < c) comp.resize(c);
            if (std::fread(comp.data(), 1, c, f) != c) break;
            size_t base = raw.size();
            raw.resize(base + u);
            int decoded = LZ4_decompress_safe((char const*)comp.data(), (char*)raw.data() + base,
                                              (int)c, (int)u);
            if (decoded < 0 || (uint32_t)decoded != u)
                throw std::runtime_error("id-stitch round bucket LZ4 decode failed: " + path);
        }
        std::fclose(f);
        std::remove(path.c_str());
        size_t pos = 0;
        while (pos < raw.size()) {
            id_tig t;
            size_t got = id_tig_deserialize(raw.data() + pos, raw.size() - pos, t);
            if (got == 0) break;
            pos += got;
            out.push_back(std::move(t));
        }
        return out;
    }

    void advance() {
        for (uint32_t b = 0; b < m_num_buckets; ++b) {
            if (!m_batch[b].empty()) flush_frame(b);
            if (m_files[b]) {
                uint32_t eof = 0;
                std::fwrite(&eof, sizeof(eof), 1, m_files[b]);
                std::fclose(m_files[b]);
                m_files[b] = nullptr;
            }
        }
        m_in_round = m_out_round;
        m_out_round = m_in_round + 1;
    }

private:
    static constexpr size_t FRAME_BUDGET = 4u * 1024 * 1024;

    std::string bucket_path(uint32_t round, uint32_t b) const {
        return m_dir + "/idstitch_r" + std::to_string(round) + "_b" + std::to_string(b) + ".bin";
    }

    void flush_frame(uint32_t b) {  // caller holds m_locks[b], or single-threaded
        auto& batch = m_batch[b];
        if (batch.empty()) return;
        if (!m_files[b]) {
            std::string path = bucket_path(m_out_round, b);
            m_files[b] = std::fopen(path.c_str(), "wb");
            if (!m_files[b])
                throw std::runtime_error("cannot open id-stitch round file: " + path + ": " +
                                         std::strerror(errno));
        }
        int src = (int)batch.size();
        int bound = LZ4_compressBound(src);
        std::vector<uint8_t> scratch((size_t)bound);
        int comp = LZ4_compress_default((char const*)batch.data(), (char*)scratch.data(), src,
                                        (int)scratch.size());
        if (comp <= 0) throw std::runtime_error("id-stitch LZ4 compress failed");
        uint32_t u = (uint32_t)src, c = (uint32_t)comp;
        std::FILE* f = m_files[b];
        if (std::fwrite(&u, sizeof(u), 1, f) != 1 || std::fwrite(&c, sizeof(c), 1, f) != 1 ||
            std::fwrite(scratch.data(), 1, (size_t)comp, f) != (size_t)comp)
            throw std::runtime_error("short write to id-stitch round file");
        batch.clear();
    }

    std::string m_dir;
    uint32_t m_num_buckets;
    std::vector<std::vector<uint8_t>> m_batch;
    std::vector<std::FILE*> m_files;
    std::vector<std::mutex> m_locks;
    uint32_t m_in_round = 0;
    uint32_t m_out_round = 0;
};

// Route an id-tig to a round-output bucket by its presented end's canonical
// junction. Returns false (does nothing) if fully closed -- caller sinks it.
template <typename Store>
inline bool id_route(id_tig&& t, uint32_t k, Store& store) {
    ext_end e;
    if (!id_end(t, k, e)) return false;
    uint32_t b = (uint32_t)((kmer_hasher{}(e.junction) >> 1) % store.num_buckets());
    store.emit(b, std::move(t));
    return true;
}

// A closed (or terminal-flushed) chain handed to assembly: the ordered entry
// list plus the residual open_flags (0 for a genuinely closed chain; the
// unmatched end's flag for a terminal-flushed one, so assembly drops the right
// extremal foreign run -- exactly as the rounds pass open_flags to
// ext_split_monochromatic today).
struct id_chain {
    std::vector<uint64_t> entries;
    uint8_t open_flags;
};

// Seed round 0: build one id-tig per fragment, route open ones, sink closed ones
// as length-1 chains. Mirrors ext_seed_round0 but reads only open_flags (+ the
// boundary k-mers of open ends) -- not the whole sequence into the loop.
// Fragments are independent, so the pass fans across num_threads: store.emit
// (per-bucket locked) and chain_sink (caller-guarded) are both thread-safe.
template <typename Source, typename Store, typename ChainSink>
inline void id_seed(Source& frag, uint32_t k, Store& store, ChainSink&& chain_sink,
                    uint32_t num_threads = 1) {
    const uint64_t n = (uint64_t)frag.size();
    auto seed_one = [&](uint64_t i) {
        uint8_t of = frag.open_flags(i);
        if (k < 2 || of == 0) {
            chain_sink(id_chain{std::vector<uint64_t>{id_entry(i, false)}, 0});
            return;
        }
        std::string_view sv = frag.seq_view(i);
        id_tig t;
        t.open_flags = of;
        t.rng = id_rng_seed(i);
        t.entries.push_back(id_entry(i, false));
        if (of & UNITIG_OPEN_LEFT) t.kl = id_fwd_kmer(sv, k, SIDE_LEFT);
        if (of & UNITIG_OPEN_RIGHT) t.kr = id_fwd_kmer(sv, k, SIDE_RIGHT);
        id_route(std::move(t), k, store);
    };
    if (num_threads <= 1) {
        for (uint64_t i = 0; i < n; ++i) seed_one(i);
    } else {
        std::atomic<uint64_t> next{0};
        constexpr uint64_t CHUNK = 4096;  // batch ids to amortize the atomic
        std::vector<std::thread> ws;
        ws.reserve(num_threads);
        for (uint32_t t = 0; t < num_threads; ++t) {
            ws.emplace_back([&]() {
                for (;;) {
                    uint64_t lo = next.fetch_add(CHUNK, std::memory_order_relaxed);
                    if (lo >= n) break;
                    uint64_t hi = std::min(lo + CHUNK, n);
                    for (uint64_t i = lo; i < hi; ++i) seed_one(i);
                }
            });
        }
        for (auto& w : ws) w.join();
    }
    store.advance();
}

// Doubling rounds on id-tigs. Each round pairs co-bucketed tigs presenting the
// same canonical boundary k-mer with compatible orientation, joins them, and
// re-routes survivors; closed tigs are emitted as chains. Terminates when a
// round makes zero joins, then flushes the remaining open tigs as terminal
// chains. Buckets within a round are independent (own waiting-map; id-seeded
// rng), so a persistent worker pool fans the per-round bucket loop across
// num_threads -- same structure as ext_run_rounds. `chain_sink` MUST be
// thread-safe when num_threads > 1 (it is called from worker threads).
template <typename Store, typename ChainSink>
inline void id_run_rounds(Store& store, uint32_t k, uint32_t num_buckets, ChainSink&& chain_sink,
                          uint32_t num_threads = 1) {
    if (num_threads == 0) num_threads = 1;
    constexpr uint32_t MAX_ROUNDS = 4096;

    auto process_bucket = [&](uint32_t b) -> uint64_t {
        std::vector<id_tig> tigs = store.take_input_bucket(b);
        if (tigs.empty()) return 0;
        const uint64_t NT = tigs.size();
        uint64_t joined = 0;
        std::vector<ext_end> ends(NT);
        std::vector<uint8_t> has_end(NT, 0), consumed(NT, 0);
        ankerl::unordered_dense::map<kmer_int_t, uint64_t, kmer_hasher> waiting;
        waiting.reserve(NT);

        for (uint64_t i = 0; i < NT; ++i) {
            ext_end e;
            if (id_end(tigs[i], k, e)) {
                e.tig_idx = i;
                ends[i] = e;
                has_end[i] = 1;
            }
        }
        for (uint64_t i = 0; i < NT; ++i) {
            if (consumed[i] || !has_end[i]) continue;
            kmer_int_t j = ends[i].junction;
            auto it = waiting.find(j);
            if (it == waiting.end()) {
                waiting.emplace(j, i);
                continue;
            }
            uint64_t a_idx = it->second;
            if (consumed[a_idx]) {
                it->second = i;
                continue;
            }
            if (!ext_pair_compatible(ends[a_idx], ends[i])) continue;  // real terminus
            id_tig merged = id_join(tigs[a_idx], ends[a_idx].side, tigs[i], ends[i].side, k);
            consumed[a_idx] = 1;
            consumed[i] = 1;
            waiting.erase(it);
            ++joined;
            if (merged.open_flags == 0)
                chain_sink(id_chain{std::move(merged.entries), 0});
            else
                id_route(std::move(merged), k, store);
        }
        for (uint64_t i = 0; i < NT; ++i) {
            if (consumed[i]) continue;
            if (tigs[i].open_flags == 0) {
                chain_sink(id_chain{std::move(tigs[i].entries), 0});
                continue;
            }
            id_tig t = std::move(tigs[i]);
            bool both = (t.open_flags & UNITIG_OPEN_LEFT) && (t.open_flags & UNITIG_OPEN_RIGHT);
            if (both) ext_rng_next(t.rng);  // re-roll presented end
            id_route(std::move(t), k, store);
        }
        return joined;
    };

    // Persistent worker pool, spawned once and reused across rounds (round count
    // is O(log L); per-round re-spawn adds measurable barrier overhead). Mirrors
    // ext_run_rounds: the main thread publishes a fresh cursor + generation and
    // signals go; workers drain buckets via the shared next_b cursor; the last
    // to finish signals done. A generation counter avoids lost/duplicate wakeups.
    std::mutex pool_mu;
    std::condition_variable cv_go, cv_done;
    std::atomic<uint32_t> next_b{0};
    std::atomic<uint64_t> joined_this_round{0};
    uint32_t generation = 0;
    uint32_t active = 0;
    bool pool_stop = false;
    const bool parallel = (num_threads > 1);

    std::vector<std::thread> pool;
    if (parallel) {
        pool.reserve(num_threads);
        for (uint32_t t = 0; t < num_threads; ++t) {
            pool.emplace_back([&]() {
                uint32_t seen = 0;
                for (;;) {
                    {
                        std::unique_lock<std::mutex> lk(pool_mu);
                        cv_go.wait(lk, [&] { return pool_stop || generation != seen; });
                        if (pool_stop) return;
                        seen = generation;
                    }
                    uint64_t local = 0;
                    for (;;) {
                        uint32_t b = next_b.fetch_add(1, std::memory_order_relaxed);
                        if (b >= num_buckets) break;
                        local += process_bucket(b);
                    }
                    joined_this_round.fetch_add(local, std::memory_order_relaxed);
                    {
                        std::unique_lock<std::mutex> lk(pool_mu);
                        if (--active == 0) cv_done.notify_one();
                    }
                }
            });
        }
    }
    struct pool_guard {
        std::mutex& mu;
        std::condition_variable& cv;
        bool& stop;
        std::vector<std::thread>& pool;
        ~pool_guard() {
            {
                std::lock_guard<std::mutex> lk(mu);
                stop = true;
            }
            cv.notify_all();
            for (auto& w : pool) w.join();
        }
    } guard{pool_mu, cv_go, pool_stop, pool};

    auto t_rounds0 = std::chrono::steady_clock::now();
    uint32_t rounds_run = 0;

    for (uint32_t round_no = 0;; ++round_no) {
        joined_this_round.store(0, std::memory_order_relaxed);
        if (!parallel) {
            for (uint32_t b = 0; b < num_buckets; ++b)
                joined_this_round.fetch_add(process_bucket(b), std::memory_order_relaxed);
        } else {
            std::unique_lock<std::mutex> lk(pool_mu);
            next_b.store(0, std::memory_order_relaxed);
            active = num_threads;
            ++generation;
            cv_go.notify_all();
            cv_done.wait(lk, [&] { return active == 0; });
        }
        store.advance();
        ++rounds_run;
        const uint64_t joined = joined_this_round.load(std::memory_order_relaxed);
        if (joined == 0 || round_no >= MAX_ROUNDS) {
            if (round_no >= MAX_ROUNDS)
                std::cerr << "[id-stitch] WARNING: round cap " << MAX_ROUNDS << " hit\n";
            // No further join is possible: flush every remaining open tig as a
            // terminal chain (its open end has no partner in the graph).
            for (uint32_t b = 0; b < num_buckets; ++b) {
                std::vector<id_tig> rem = store.take_input_bucket(b);
                for (auto& t : rem) chain_sink(id_chain{std::move(t.entries), t.open_flags});
            }
            break;
        }
    }
    double secs =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t_rounds0).count();
    std::cout << "  [id-stitch] " << rounds_run << " rounds, " << secs << "s\n";
}

// Split an assembled tig into monochromatic unitigs. Kept separate so the
// file/parallel path can reuse it. Builds a scratch ext_tig only to feed
// ext_split_monochromatic (which takes one).
template <typename Sink>
inline void id_tig_assembled_split(std::string& seq, std::vector<color_run>& runs,
                                   uint8_t open_flags, uint32_t k, Sink&& sink) {
    ext_tig t;
    t.seq = std::move(seq);
    t.runs = std::move(runs);
    t.open_flags = open_flags;
    ext_split_monochromatic(t, k, open_flags, [&](stitchable_unitig&& u) { sink(std::move(u)); });
}

// Assemble one chain's bases+colors and split into monochromatic unitigs. This
// is a fold of the SAME ext_concat_runs + seq-append (drop the k-base overlap)
// the rounds' ext_join performs, so the result is identical to what the
// base-carrying stitch produces for the same chain -- then the same
// ext_split_monochromatic cut. `frag` provides seq_view(id)/runs(id).
template <typename Source, typename Sink>
inline void id_assemble_chain(id_chain const& c, Source& frag, uint32_t k, Sink&& sink) {
    std::string seq;
    std::vector<color_run> runs;
    for (size_t idx = 0; idx < c.entries.size(); ++idx) {
        const uint64_t fid = id_entry_frag(c.entries[idx]);
        const bool rc = id_entry_rc(c.entries[idx]);
        std::string_view sv = frag.seq_view(fid);
        std::string s(sv.data(), sv.size());
        std::vector<color_run> r = frag.runs(fid);
        if (rc) {
            s = revcomp_string(s);
            ext_reverse_runs(r);
        }
        if (idx == 0) {
            seq = std::move(s);
            runs = std::move(r);
        } else {
            ext_concat_runs(runs, std::move(r));    // reconcile + drop shared X unit
            seq.append(s.begin() + k, s.end());      // drop the k-base overlap
        }
    }
    id_tig_assembled_split(seq, runs, c.open_flags, k, sink);
}

// --- Step B.2: scalable disk-based assembly ---------------------------------
// A one-shot bucketed byte spill: append opaque records into num_buckets files
// (LZ4-framed, per-bucket write batch), then read each bucket's decoded bytes
// back once. Same frame format as the round store; unlike it there are no rounds
// (write phase, then read phase). Records are parsed by the caller (it knows the
// fixed/var layout), so this only moves bytes.
class frame_spill_writer {
public:
    frame_spill_writer(std::string dir, std::string prefix, uint32_t num_buckets)
        : m_dir(std::move(dir)), m_prefix(std::move(prefix)), m_num_buckets(num_buckets),
          m_batch(num_buckets), m_files(num_buckets, nullptr), m_locks(num_buckets) {}

    ~frame_spill_writer() {
        for (auto* f : m_files)
            if (f) std::fclose(f);
    }

    frame_spill_writer(frame_spill_writer const&) = delete;
    frame_spill_writer& operator=(frame_spill_writer const&) = delete;

    uint32_t num_buckets() const { return m_num_buckets; }

    std::string bucket_path(uint32_t b) const {
        return m_dir + "/" + m_prefix + std::to_string(b) + ".bin";
    }

    void put(uint32_t b, void const* data, size_t n) {
        std::lock_guard<std::mutex> lk(m_locks[b]);
        uint8_t const* p = (uint8_t const*)data;
        m_batch[b].insert(m_batch[b].end(), p, p + n);
        if (m_batch[b].size() >= FRAME_BUDGET) flush_frame(b);
    }

    // Flush + EOF-mark + close all buckets (call once, after all puts).
    void close() {
        for (uint32_t b = 0; b < m_num_buckets; ++b) {
            if (!m_batch[b].empty()) flush_frame(b);
            if (m_files[b]) {
                uint32_t eof = 0;
                std::fwrite(&eof, sizeof(eof), 1, m_files[b]);
                std::fclose(m_files[b]);
                m_files[b] = nullptr;
            }
        }
    }

private:
    static constexpr size_t FRAME_BUDGET = 4u * 1024 * 1024;

    void flush_frame(uint32_t b) {  // caller holds m_locks[b]
        auto& batch = m_batch[b];
        if (batch.empty()) return;
        if (!m_files[b]) {
            std::string path = bucket_path(b);
            m_files[b] = std::fopen(path.c_str(), "wb");
            if (!m_files[b])
                throw std::runtime_error("cannot open spill file: " + path + ": " +
                                         std::strerror(errno));
        }
        int src = (int)batch.size();
        std::vector<uint8_t> scratch((size_t)LZ4_compressBound(src));
        int comp = LZ4_compress_default((char const*)batch.data(), (char*)scratch.data(), src,
                                        (int)scratch.size());
        if (comp <= 0) throw std::runtime_error("spill LZ4 compress failed");
        uint32_t u = (uint32_t)src, c = (uint32_t)comp;
        std::FILE* f = m_files[b];
        if (std::fwrite(&u, sizeof(u), 1, f) != 1 || std::fwrite(&c, sizeof(c), 1, f) != 1 ||
            std::fwrite(scratch.data(), 1, (size_t)comp, f) != (size_t)comp)
            throw std::runtime_error("short write to spill file");
        batch.clear();
    }

    std::string m_dir, m_prefix;
    uint32_t m_num_buckets;
    std::vector<std::vector<uint8_t>> m_batch;
    std::vector<std::FILE*> m_files;
    std::vector<std::mutex> m_locks;
};

// Decode every frame of one spill bucket file into a flat byte vector, then
// delete the file. Returns empty if the bucket file does not exist (no record
// was ever routed there). Peak extra RAM is this one bucket's decoded bytes --
// the bounded working set of the assembly.
inline std::vector<uint8_t> read_spill_bucket_decoded(std::string const& path) {
    std::vector<uint8_t> raw;
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return raw;
    std::vector<uint8_t> comp;
    for (;;) {
        uint32_t u = 0;
        if (std::fread(&u, sizeof(u), 1, f) != 1) break;
        if (u == 0) break;
        uint32_t c = 0;
        if (std::fread(&c, sizeof(c), 1, f) != 1) break;
        if (comp.size() < c) comp.resize(c);
        if (std::fread(comp.data(), 1, c, f) != c) break;
        size_t base = raw.size();
        raw.resize(base + u);
        int decoded = LZ4_decompress_safe((char const*)comp.data(), (char*)raw.data() + base,
                                          (int)c, (int)u);
        if (decoded < 0 || (uint32_t)decoded != u)
            throw std::runtime_error("spill bucket LZ4 decode failed: " + path);
    }
    std::fclose(f);
    std::remove(path.c_str());
    return raw;
}

}  // namespace detail

// In-RAM id-only compaction stitch (STEP 1): same Source/Sink contract as
// stitch_unitigs_extmem, for validation against it in test_stitch. Seeds + runs
// the id-only doubling collecting closed chains, then assembles + splits each.
template <typename Source, typename Sink>
inline void compact_stitch_mem(Source& frag, uint32_t k, Sink&& sink, uint32_t num_buckets = 0) {
    if (num_buckets == 0) num_buckets = 256;
    detail::id_round_store_mem store(num_buckets);
    std::vector<detail::id_chain> chains;
    auto chain_sink = [&](detail::id_chain&& c) { chains.push_back(std::move(c)); };
    detail::id_seed(frag, k, store, chain_sink);
    detail::id_run_rounds(store, k, num_buckets, chain_sink);
    for (auto const& c : chains) detail::id_assemble_chain(c, frag, k, sink);
}

// File-backed id-only compaction stitch (STEP B.1): same as compact_stitch_mem
// but the doubling round store lives on disk (one bucket resident), validating
// the id_tig codec + file store. Chains and assembly are still in RAM here; the
// scalable disk-based chain spill + re-bucket assembly come in B.2.
template <typename Source, typename Sink>
inline void compact_stitch_file(Source& frag, uint32_t k, std::string const& tmp_dir, Sink&& sink,
                                uint32_t num_buckets = 0) {
    if (num_buckets == 0) num_buckets = 1024;
    detail::id_round_store_file store(tmp_dir, num_buckets);
    std::vector<detail::id_chain> chains;
    auto chain_sink = [&](detail::id_chain&& c) { chains.push_back(std::move(c)); };
    detail::id_seed(frag, k, store, chain_sink);
    detail::id_run_rounds(store, k, num_buckets, chain_sink);
    for (auto const& c : chains) detail::id_assemble_chain(c, frag, k, sink);
}

// Compact in-RAM fragment store: a flat base-byte arena + offset array (NOT
// 92.5M std::strings, which would cost ~5 GiB of object overhead alone), a flat
// color-run arena + offsets, and a flags byte per fragment. Satisfies the
// Source interface (size/open_flags/seq_view/runs) with O(1) random access, so
// the id-only stitch can seed AND assemble straight from RAM -- no second pass
// over the frag spill and no disk re-bucket. Resident bytes ~= frag spill size
// plus ~16 B/fragment of offsets; bounded by the caller's budget check.
class arena_frag_source {
public:
    void reserve(uint64_t n_frags, uint64_t seq_bytes, uint64_t run_count) {
        m_flags.reserve(n_frags);
        m_seq.reserve(seq_bytes);
        m_seq_off.reserve(n_frags + 1);
        m_runs.reserve(run_count);
        m_runs_off.reserve(n_frags + 1);
    }
    // Append one fragment in frag_id order (frag_id == call index).
    void append(uint8_t flags, std::string_view seq, std::vector<color_run> const& runs) {
        m_flags.push_back(flags);
        m_seq.insert(m_seq.end(), seq.begin(), seq.end());
        m_seq_off.push_back(m_seq.size());
        m_runs.insert(m_runs.end(), runs.begin(), runs.end());
        m_runs_off.push_back(m_runs.size());
    }
    size_t size() const { return m_flags.size(); }
    uint8_t open_flags(size_t i) const { return m_flags[i]; }
    std::string_view seq_view(size_t i) const {
        return std::string_view(m_seq.data() + m_seq_off[i],
                                (size_t)(m_seq_off[i + 1] - m_seq_off[i]));
    }
    std::vector<color_run> runs(size_t i) const {
        return std::vector<color_run>(m_runs.begin() + (ptrdiff_t)m_runs_off[i],
                                      m_runs.begin() + (ptrdiff_t)m_runs_off[i + 1]);
    }
    uint64_t seq_bytes() const { return m_seq.size(); }

private:
    std::vector<char> m_seq;             // all fragment bases concatenated
    std::vector<uint64_t> m_seq_off{0};  // seq[i] = [off[i], off[i+1])
    std::vector<color_run> m_runs;       // all color runs concatenated
    std::vector<uint64_t> m_runs_off{0};
    std::vector<uint8_t> m_flags;
};

// In-RAM id-only compaction stitch. The whole frag store is resident (arena),
// so seed and assembly read bases straight from RAM: there is NO second spill
// read and NO disk re-bucket (the only reason compact_stitch_scalable spills is
// to get frag_id->chain locality on disk -- with random access that vanishes).
// The id-only doubling round store stays on disk (small, keeps RAM for the
// arena). Rounds and the final assembly both run across num_threads. Use this
// when the arena fits the RAM budget; else use compact_stitch_scalable.
template <typename Source, typename Sink>
inline void compact_stitch_inram(Source& frag, uint32_t k, std::string const& tmp_dir, Sink&& sink,
                                 uint32_t num_buckets = 0, uint32_t num_threads = 1) {
    using namespace detail;
    if (num_buckets == 0) num_buckets = 1024;
    id_round_store_file store(tmp_dir, num_buckets);
    std::vector<id_chain> chains;
    std::mutex chains_mu;
    // chain_sink is called from the parallel round workers, so it is guarded.
    auto chain_sink = [&](id_chain&& c) {
        std::lock_guard<std::mutex> lk(chains_mu);
        chains.push_back(std::move(c));
    };

    auto t_seed = std::chrono::steady_clock::now();
    id_seed(frag, k, store, chain_sink, num_threads);
    std::cout << "  [id-stitch] seed "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t_seed).count()
              << "s\n";
    id_run_rounds(store, k, num_buckets, chain_sink, num_threads);  // prints its own timing

    // Assemble: chains are independent, so fan across num_threads, reading each
    // fragment's bases directly from the arena. The sink is not thread-safe, so
    // each worker batches and drains under sink_mu.
    auto t_a = std::chrono::steady_clock::now();
    std::mutex sink_mu;
    constexpr size_t SINK_BATCH_BYTES = 2u * 1024 * 1024;
    auto flush_batch = [&](std::vector<stitchable_unitig>& batch) {
        if (batch.empty()) return;
        std::lock_guard<std::mutex> lk(sink_mu);
        for (auto& u : batch) sink(std::move(u));
        batch.clear();
    };
    std::atomic<size_t> next{0};
    const size_t nc = chains.size();
    auto work = [&]() {
        std::vector<stitchable_unitig> batch;
        size_t bytes = 0;
        for (;;) {
            size_t i = next.fetch_add(1, std::memory_order_relaxed);
            if (i >= nc) break;
            id_assemble_chain(chains[i], frag, k, [&](stitchable_unitig&& u) {
                bytes += u.seq.size();
                batch.push_back(std::move(u));
                if (bytes >= SINK_BATCH_BYTES) {
                    flush_batch(batch);
                    bytes = 0;
                }
            });
        }
        flush_batch(batch);
    };
    if (num_threads <= 1) {
        work();
    } else {
        std::vector<std::thread> ws;
        ws.reserve(num_threads);
        for (uint32_t t = 0; t < num_threads; ++t) ws.emplace_back(work);
        for (auto& w : ws) w.join();
    }
    std::cout << "  [id-stitch] assemble "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t_a).count()
              << "s\n";
}
// RAM: the only resident structures are one round-store bucket, one frag-id
// range's member array, and one chain bucket's bases -- so peak RAM is
// independent of fragment/chain COUNT (the extmem invariant).
//
// `for_each_frag(fn)` must call fn(uint8_t open_flags, std::vector<color_run>&
// runs, std::string& seq) for every fragment in frag_id order, and be
// re-invocable (it is called twice: seed, then phase-2 base attach). `n_frags`
// is the total fragment count (known by the producer). `dir` is scratch.
//
// Pipeline:
//   seed+rounds : id-only doubling -> per-chain-member records (frag_id,
//                 chain_id, pos, rc, open_flags) bucketed by frag_id RANGE.
//   phase 2     : stream the fragments once; for each, look up its member
//                 (one range's array resident at a time), orient by rc, and
//                 re-bucket its bases+runs by chain_id % C.
//   phase 3     : per chain bucket, group by chain_id, sort by pos, fold-
//                 assemble (same ext_concat_runs + seq-append + split) -> sink.
template <typename ForEachFrag, typename Sink>
inline void compact_stitch_scalable(ForEachFrag&& for_each_frag, uint64_t n_frags, uint32_t k,
                                    std::string const& dir, Sink&& sink, uint32_t num_buckets = 0,
                                    uint32_t frag_ranges = 0, uint32_t chain_buckets = 0,
                                    uint32_t num_threads = 1) {
    using namespace detail;
    if (num_buckets == 0) num_buckets = 1024;
    if (frag_ranges == 0) frag_ranges = 64;
    // More chain buckets than threads so phase 3 load-balances and only a few
    // small buckets are co-resident per worker.
    if (chain_buckets == 0) chain_buckets = std::max<uint32_t>(64, num_threads * 8);
    if (n_frags == 0) return;
    const uint64_t range_size = (n_frags + frag_ranges - 1) / frag_ranges;

    // --- seed + doubling: emit member records bucketed by frag_id range ------
    id_round_store_file store(dir, num_buckets);
    frame_spill_writer members(dir, "idmem_", frag_ranges);
    // Atomic: chain_sink is called from the round driver's worker threads. Only
    // uniqueness matters (chains bucket by chain_id % C), not order.
    std::atomic<uint64_t> chain_counter{0};
    auto chain_sink = [&](id_chain&& c) {
        const uint64_t cid = chain_counter.fetch_add(1, std::memory_order_relaxed);
        const uint8_t of = c.open_flags;
        for (uint32_t pos = 0; pos < c.entries.size(); ++pos) {
            const uint64_t fid = id_entry_frag(c.entries[pos]);
            const uint8_t pack = (uint8_t)((id_entry_rc(c.entries[pos]) ? 1u : 0u) | (of << 1));
            uint8_t rec[8 + 8 + 4 + 1];
            std::memcpy(rec, &fid, 8);
            std::memcpy(rec + 8, &cid, 8);
            std::memcpy(rec + 16, &pos, 4);
            rec[20] = pack;
            members.put((uint32_t)(fid / range_size), rec, sizeof(rec));
        }
    };

    auto t_seed = std::chrono::steady_clock::now();
    {
        // Seed reads each fragment's open_flags + boundary k-mers (it ignores
        // runs, but the reader fills them). frag_id = iteration order.
        uint64_t fid = 0;
        std::vector<color_run> runs;
        std::string seq;
        for_each_frag([&](uint8_t of, std::vector<color_run>& r, std::string& s) {
            (void)r;
            if (k < 2 || of == 0) {
                chain_sink(id_chain{std::vector<uint64_t>{id_entry(fid, false)}, 0});
            } else {
                id_tig t;
                t.open_flags = of;
                t.rng = id_rng_seed(fid);
                t.entries.push_back(id_entry(fid, false));
                std::string_view sv(s.data(), s.size());
                if (of & UNITIG_OPEN_LEFT) t.kl = id_fwd_kmer(sv, k, SIDE_LEFT);
                if (of & UNITIG_OPEN_RIGHT) t.kr = id_fwd_kmer(sv, k, SIDE_RIGHT);
                id_route(std::move(t), k, store);
            }
            ++fid;
        });
        store.advance();
    }
    std::cout << "  [id-stitch] seed " << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - t_seed).count() << "s\n";
    id_run_rounds(store, k, num_buckets, chain_sink, num_threads);
    members.close();

    // --- phase 2: attach bases, re-bucket by chain_id % C --------------------
    // Single-threaded and pipelined: the frag spill is one stream so the read is
    // serial (read-bound on spinning disk), and pipelining the read with the
    // orient+serialize+route keeps the reader's seq/runs buffers reused (no
    // per-fragment allocation). A range-buffered parallel variant was tried and
    // REVERTED -- buffering stole the reader's buffers (forcing ~num_frags
    // reallocations) and, since phase 2 is read-bound here, the parallelism
    // could not pay that back. Peak extra RAM is one frag-id range's member array.
    auto t_p2 = std::chrono::steady_clock::now();
    frame_spill_writer bases(dir, "idbase_", chain_buckets);
    {
        struct Slot {
            uint64_t cid;
            uint32_t pos;
            uint8_t pack;
        };
        std::vector<Slot> arr;     // current range's members, indexed by (fid - lo)
        uint64_t cur_range = UINT64_MAX, lo = 0;
        uint64_t fid = 0;
        std::vector<uint8_t> rbuf;  // reused base-record buffer
        for_each_frag([&](uint8_t of, std::vector<color_run>& runs, std::string& seq) {
            (void)of;
            const uint64_t range = fid / range_size;
            if (range != cur_range) {
                std::vector<uint8_t> raw = read_spill_bucket_decoded(members.bucket_path(range));
                lo = range * range_size;
                const uint64_t hi = std::min(lo + range_size, n_frags);
                arr.assign((size_t)(hi - lo), Slot{0, 0, 0xFF});
                for (size_t p = 0; p + 21 <= raw.size(); p += 21) {
                    uint64_t rfid, rcid;
                    uint32_t rpos;
                    std::memcpy(&rfid, raw.data() + p, 8);
                    std::memcpy(&rcid, raw.data() + p + 8, 8);
                    std::memcpy(&rpos, raw.data() + p + 16, 4);
                    arr[(size_t)(rfid - lo)] = Slot{rcid, rpos, raw[p + 20]};
                }
                cur_range = range;
            }
            const Slot& m = arr[(size_t)(fid - lo)];
            const bool rc = (m.pack & 1u) != 0;
            const uint8_t cof = (uint8_t)(m.pack >> 1);
            std::string s = seq;
            std::vector<color_run> r = runs;
            if (rc) {
                s = revcomp_string(s);
                ext_reverse_runs(r);
            }
            // base record: [u64 cid][u32 pos][u8 open_flags][u32 seq_len][seq]
            //              [u32 nruns][nruns x (u64 cid + u32 num_kmers)]
            rbuf.clear();
            auto put = [&](void const* d, size_t n) {
                uint8_t const* b = (uint8_t const*)d;
                rbuf.insert(rbuf.end(), b, b + n);
            };
            uint32_t slen = (uint32_t)s.size(), nruns = (uint32_t)r.size();
            put(&m.cid, 8);
            put(&m.pos, 4);
            put(&cof, 1);
            put(&slen, 4);
            put(s.data(), slen);
            put(&nruns, 4);
            for (auto const& cr : r) {
                put(&cr.cid, sizeof(cr.cid));
                put(&cr.num_kmers, sizeof(cr.num_kmers));
            }
            bases.put((uint32_t)(m.cid % chain_buckets), rbuf.data(), rbuf.size());
            ++fid;
        });
    }
    bases.close();
    std::cout << "  [id-stitch] phase2 (attach+rebucket) "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t_p2).count()
              << "s\n";

    // --- phase 3: per chain bucket, group + sort + fold-assemble + split -----
    // Chain buckets are fully independent (a chain lives in exactly one), so the
    // work fans across num_threads. The sink is not thread-safe, so each worker
    // batches its split unitigs and drains them under sink_mu (same scheme as
    // ext_run_rounds). Peak extra RAM is num_threads co-resident chain buckets.
    auto t_p3 = std::chrono::steady_clock::now();
    struct BaseRec {
        uint32_t pos;
        uint8_t open_flags;
        std::string seq;
        std::vector<color_run> runs;
    };
    std::mutex sink_mu;
    constexpr size_t SINK_BATCH_BYTES = 2u * 1024 * 1024;
    auto flush_batch = [&](std::vector<stitchable_unitig>& batch) {
        if (batch.empty()) return;
        std::lock_guard<std::mutex> lk(sink_mu);
        for (auto& u : batch) sink(std::move(u));
        batch.clear();
    };
    auto process_cb = [&](uint32_t cb) {
        std::vector<uint8_t> raw = read_spill_bucket_decoded(bases.bucket_path(cb));
        if (raw.empty()) return;
        ankerl::unordered_dense::map<uint64_t, std::vector<BaseRec>> groups;
        size_t p = 0;
        while (p + 17 <= raw.size()) {
            uint64_t cid;
            uint32_t pos, slen;
            std::memcpy(&cid, raw.data() + p, 8);
            std::memcpy(&pos, raw.data() + p + 8, 4);
            uint8_t of = raw[p + 12];
            std::memcpy(&slen, raw.data() + p + 13, 4);
            p += 17;
            if (p + slen + 4 > raw.size()) break;
            BaseRec br;
            br.pos = pos;
            br.open_flags = of;
            br.seq.assign((char const*)raw.data() + p, slen);
            p += slen;
            uint32_t nruns;
            std::memcpy(&nruns, raw.data() + p, 4);
            p += 4;
            if (p + (size_t)nruns * 12 > raw.size()) break;
            br.runs.resize(nruns);
            for (uint32_t i = 0; i < nruns; ++i) {
                std::memcpy(&br.runs[i].cid, raw.data() + p, 8);
                std::memcpy(&br.runs[i].num_kmers, raw.data() + p + 8, 4);
                p += 12;
            }
            groups[cid].push_back(std::move(br));
        }
        std::vector<stitchable_unitig> batch;
        size_t batch_bytes = 0;
        for (auto& kv : groups) {
            std::vector<BaseRec>& recs = kv.second;
            std::sort(recs.begin(), recs.end(),
                      [](BaseRec const& a, BaseRec const& b) { return a.pos < b.pos; });
            std::string seq = std::move(recs[0].seq);
            std::vector<color_run> runs = std::move(recs[0].runs);
            const uint8_t of = recs[0].open_flags;
            for (size_t i = 1; i < recs.size(); ++i) {
                ext_concat_runs(runs, std::move(recs[i].runs));
                seq.append(recs[i].seq.begin() + k, recs[i].seq.end());
            }
            id_tig_assembled_split(seq, runs, of, k, [&](stitchable_unitig&& u) {
                batch_bytes += u.seq.size();
                batch.push_back(std::move(u));
                if (batch_bytes >= SINK_BATCH_BYTES) {
                    flush_batch(batch);
                    batch_bytes = 0;
                }
            });
        }
        flush_batch(batch);
    };
    if (num_threads <= 1) {
        for (uint32_t cb = 0; cb < chain_buckets; ++cb) process_cb(cb);
    } else {
        std::atomic<uint32_t> next_cb{0};
        std::vector<std::thread> workers;
        workers.reserve(num_threads);
        for (uint32_t t = 0; t < num_threads; ++t) {
            workers.emplace_back([&]() {
                for (;;) {
                    uint32_t cb = next_cb.fetch_add(1, std::memory_order_relaxed);
                    if (cb >= chain_buckets) break;
                    process_cb(cb);
                }
            });
        }
        for (auto& w : workers) w.join();
    }
    std::cout << "  [id-stitch] phase3 (assemble+split) "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t_p3).count()
              << "s\n";
}

}  // namespace cdbg
