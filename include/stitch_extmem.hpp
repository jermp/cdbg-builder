#pragma once

// External-memory cross-bucket stitch via GGCAT-style hash-bucketed
// iterative doubling.
//
// Replaces the in-RAM stitch (stitch.hpp), whose by_junction map + adj
// array + visited array are all O(num_fragments) resident -- hundreds of
// GB on the Blackwell 661k pangenome (5.7e9 fragments). The doubling
// approach holds only one hash-bucket's worth of tigs in RAM at a time,
// so peak is bounded by the bucket count, not the fragment count.
//
// ALGORITHM (faithful to GGCAT v2 extend_unitigs, adapted to our model):
//
//   A "tig" is a (partially) merged fragment: { cid, open_flags, seq }.
//   An OPEN side means the unitig continues into another bucket (degree
//   exactly 1, neighbour elsewhere); by the per-bucket walk's OPEN
//   invariant a (k-1) junction has AT MOST 2 open ends globally.
//
//   Up front, compute the JOINABLE junction set: the (k-1) junctions
//   with exactly 2 open ends, equal cid, compatible orientation -- the
//   in-RAM stitch's exact join predicate. This is invariant under
//   merging, so it is computed once. An open end NOT at a joinable
//   junction is a terminal unitig boundary.
//
//   Round loop (until the store drains):
//     1. Each non-terminal tig presents ONE chosen open side and is
//        routed to bucket H(junction) % NUM_BUCKETS. A both-open tig
//        picks its presented end from a per-tig RNG, re-rolled each
//        round it survives, so a chain of L fragments collapses in
//        O(log L) expected rounds.
//     2. Within each bucket, two tigs presenting the same joinable
//        junction are concatenated (k-1 overlap dropped); the merged
//        tig is sunk if now terminal, else carried to the next round.
//     3. A tig all of whose open ends are terminal (non-joinable) is
//        sunk immediately as a complete unitig. This is what guarantees
//        termination: non-joinable ends never become joinable, so no
//        tig dangles forever.
//
// DIFFERENCES from GGCAT (intentional, for our data model):
//   - Key is the canonical (k-1)-mer junction, not the full boundary
//     k-mer (our fragments overlap by k-1). GGCAT's full-k-mer key makes
//     dBG branches land in different buckets for free; our (k-1) key
//     collides them, so we exclude non-joinable junctions explicitly via
//     the precomputed set.
//   - Join is gated on equal cid (our unitigs are monochromatic;
//     GGCAT enforces this differently). GGCAT joins on k-mer adjacency
//     alone.
//
// Two stores share one round driver: round_store_mem (in-RAM rounds;
// fast, used by tests/reference) and round_store_file (LZ4-framed files,
// one bucket resident at a time; the production external-memory path).

#include <array>
#include <atomic>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include <lz4.h>
#include <unordered_dense/unordered_dense.h>

#include "bucket_walker.hpp"  // stitchable_unitig, UNITIG_OPEN_*
#include "kmer.hpp"
#include "stitch.hpp"  // detail::side_junction_canonical, SIDE_*, revcomp_string

namespace cdgb {

namespace detail {

// Set of (k-1) junctions that are JOINABLE: exactly 2 open ends share
// them, those ends have equal cid, and their side/orientation is
// compatible -- i.e. precisely the in-RAM stitch's join predicate. Any
// open end whose junction is NOT in this set is a terminal unitig
// boundary (branch with 3+ ends, singleton, or 2-but-incompatible).
//
// All three properties are INVARIANT under merging: a junction's two
// ends are always the same two original-fragment ends (merges only
// relocate a tig's OTHER, outer end), so the predicate is computable
// once over the initial fragments and never changes. This is the single
// source of truth that fixes both failure modes seen earlier:
//   - false joins (joining a count==2-but-incompatible junction), and
//   - the hang (carrying forever a tig whose only open end is at a
//     non-joinable junction -- it's actually terminal and must be sunk).
//
// NOTE: in-RAM for now -- the one remaining unbounded structure in the
// otherwise external-memory stitch; bucketing it to disk is a
// documented follow-up. Correctness first.
using joinable_set_t = ankerl::unordered_dense::set<kmer_int_t, kmer_hasher>;

// A tig in the doubling loop. `rng` is a per-tig random state used to
// pick which open end to present when both are open. It MUST vary
// between tigs (we seed it from sequence content) so two both-open tigs
// don't flip in lockstep -- if they did, a pair needing to meet on a
// specific shared end would never align and the loop would not
// terminate. Randomized selection is what gives O(log L) expected
// rounds (GGCAT does the same).
struct ext_tig {
    uint64_t cid;
    uint8_t open_flags;
    uint64_t rng;
    std::string seq;
};

// xorshift step; returns the next state.
inline uint64_t ext_rng_next(uint64_t& s) {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return s;
}

// Seed a tig's rng deterministically from its sequence so runs are
// reproducible but distinct tigs get distinct streams.
inline uint64_t ext_seed_from_seq(std::string const& seq) {
    uint64_t h = 1469598103934665603ull;
    for (unsigned char c : seq) {
        h ^= c;
        h *= 1099511628211ull;
    }
    return h ? h : 0x9e3779b97f4a7c15ull;  // never zero (xorshift fixed point)
}

// One chosen open end of a tig, in the current round's bucket.
struct ext_end {
    uint64_t tig_idx;       // index into the round-bucket's tig vector
    kmer_int_t junction;    // canonical (k-1)-mer
    uint8_t side;           // SIDE_LEFT / SIDE_RIGHT of the keyed end
    bool is_canonical_fwd;  // junction-in-own-frame == canonical
};

// Pick the open side a tig presents this round. For a both-open tig the
// side is read from the current rng state (NOT advanced here, so repeated
// queries within one round are stable); the caller advances rng when it
// carries the tig to the next round. Returns false if fully closed.
inline bool ext_choose_side(ext_tig const& t, uint8_t& out_side) {
    bool l = (t.open_flags & UNITIG_OPEN_LEFT) != 0;
    bool r = (t.open_flags & UNITIG_OPEN_RIGHT) != 0;
    if (l && r) {
        out_side = (t.rng & 1) ? SIDE_RIGHT : SIDE_LEFT;
        return true;
    }
    if (l) {
        out_side = SIDE_LEFT;
        return true;
    }
    if (r) {
        out_side = SIDE_RIGHT;
        return true;
    }
    return false;
}

// Compatibility of two ends sharing a junction (same rule as the in-RAM
// stitch's pair_compatible).
inline bool ext_pair_compatible(ext_end const& a, ext_end const& b) {
    if (a.side != b.side) return a.is_canonical_fwd == b.is_canonical_fwd;
    return a.is_canonical_fwd != b.is_canonical_fwd;
}

// Join two tigs at a shared, compatible junction. Orient `a` so its
// keyed end (a_side) is at the RIGHT and `b` so its keyed end (b_side)
// is at the LEFT, then merged = a_oriented + b_oriented[k-1..]. The two
// (k-1) overlaps are guaranteed equal by ext_pair_compatible (verified
// for all four side/orientation cases). Returns the merged tig.
inline ext_tig ext_join(ext_tig const& a, uint8_t a_side, ext_tig const& b, uint8_t b_side,
                        uint32_t k) {
    // a oriented with keyed end at RIGHT.
    std::string as = (a_side == SIDE_RIGHT) ? a.seq : revcomp_string(a.seq);
    // b oriented with keyed end at LEFT.
    std::string bs = (b_side == SIDE_LEFT) ? b.seq : revcomp_string(b.seq);

    ext_tig m;
    m.cid = a.cid;
    m.seq.reserve(as.size() + bs.size() - (k - 1));
    m.seq = as;
    m.seq.append(bs.begin() + (k - 1), bs.end());
    // Mix both parents' rng so the merged tig gets a fresh stream.
    m.rng = a.rng ^ (b.rng * 0x9e3779b97f4a7c15ull);
    if (m.rng == 0) m.rng = 0x9e3779b97f4a7c15ull;

    // Merged LEFT = a's other end (the non-keyed end, after orienting a
    // with keyed end at RIGHT, that's a's LEFT in the `as` frame).
    //   a_side==RIGHT (as=a):            merged-left open == a OPEN_LEFT
    //   a_side==LEFT  (as=revcomp(a)):   merged-left open == a OPEN_RIGHT
    bool left_open = (a_side == SIDE_RIGHT) ? (a.open_flags & UNITIG_OPEN_LEFT)
                                            : (a.open_flags & UNITIG_OPEN_RIGHT);
    // Merged RIGHT = b's other end.
    //   b_side==LEFT  (bs=b):            merged-right open == b OPEN_RIGHT
    //   b_side==RIGHT (bs=revcomp(b)):   merged-right open == b OPEN_LEFT
    bool right_open = (b_side == SIDE_LEFT) ? (b.open_flags & UNITIG_OPEN_RIGHT)
                                            : (b.open_flags & UNITIG_OPEN_LEFT);
    m.open_flags = (uint8_t)((left_open ? UNITIG_OPEN_LEFT : 0) |
                             (right_open ? UNITIG_OPEN_RIGHT : 0));
    return m;
}

// In-memory round store: NUM_BUCKETS vectors of tigs. The file-backed
// store (added next) has the same interface: emit(bucket, tig),
// for_each_bucket(fn), swap_to_input(), clear_output(). Holding the
// whole round in RAM is only for algorithm validation; the file store
// holds one bucket at a time.
struct round_store_mem {
    explicit round_store_mem(uint32_t num_buckets)
        : m_in(num_buckets), m_out(num_buckets), m_num_buckets(num_buckets) {}

    uint32_t num_buckets() const { return m_num_buckets; }

    void emit(uint32_t b, ext_tig&& t) { m_out[b].push_back(std::move(t)); }

    // Move and return input bucket b (caller consumes it; storage is
    // freed). Matches the file store's take-by-value semantics so the
    // driver is store-agnostic.
    std::vector<ext_tig> take_input_bucket(uint32_t b) {
        std::vector<ext_tig> v = std::move(m_in[b]);
        m_in[b].clear();
        return v;
    }

    // Move this round's output to be the next round's input.
    void advance() {
        m_in.swap(m_out);
        for (auto& v : m_out) v.clear();
    }

private:
    std::vector<std::vector<ext_tig>> m_in;
    std::vector<std::vector<ext_tig>> m_out;
    uint32_t m_num_buckets;
};

// ext_tig (de)serialization for the file-backed round store. Layout:
//   [u64 cid][u8 open_flags][u64 rng][u32 seq_len][seq_len bytes ACGT]
inline void ext_tig_serialize(ext_tig const& t, std::vector<uint8_t>& out) {
    auto put = [&](void const* p, size_t n) {
        uint8_t const* b = (uint8_t const*)p;
        out.insert(out.end(), b, b + n);
    };
    put(&t.cid, sizeof(t.cid));
    put(&t.open_flags, sizeof(t.open_flags));
    put(&t.rng, sizeof(t.rng));
    uint32_t len = (uint32_t)t.seq.size();
    put(&len, sizeof(len));
    put(t.seq.data(), len);
}

// Returns bytes consumed, 0 on malformed/EOF.
inline size_t ext_tig_deserialize(uint8_t const* buf, size_t buf_len, ext_tig& t) {
    constexpr size_t HDR = sizeof(uint64_t) + 1 + sizeof(uint64_t) + sizeof(uint32_t);
    if (buf_len < HDR) return 0;
    size_t p = 0;
    std::memcpy(&t.cid, buf + p, sizeof(t.cid));
    p += sizeof(t.cid);
    t.open_flags = buf[p];
    p += 1;
    std::memcpy(&t.rng, buf + p, sizeof(t.rng));
    p += sizeof(t.rng);
    uint32_t len;
    std::memcpy(&len, buf + p, sizeof(len));
    p += sizeof(len);
    if (buf_len < p + len) return 0;
    t.seq.assign((char const*)buf + p, len);
    p += len;
    return p;
}

// File-backed round store: one LZ4-framed file per bucket per round.
// emit() appends to per-bucket write batches (flushed as frames at a
// size threshold). take_input_bucket() decodes one bucket's file fully
// into a vector, then deletes it -- so peak RAM is ONE bucket's tigs
// plus the per-bucket write batches, not the whole round. Frame format
// matches bucket_io.hpp: [u32 uncompressed][u32 compressed][bytes],
// terminated by [u32 0].
class round_store_file {
public:
    round_store_file(std::string dir, uint32_t num_buckets)
        : m_dir(std::move(dir)), m_num_buckets(num_buckets), m_batch(num_buckets),
          m_files(num_buckets, nullptr) {
        open_round_files(/*round=*/0);
    }

    ~round_store_file() {
        for (auto* f : m_files)
            if (f) std::fclose(f);
    }

    round_store_file(round_store_file const&) = delete;
    round_store_file& operator=(round_store_file const&) = delete;

    uint32_t num_buckets() const { return m_num_buckets; }

    void emit(uint32_t b, ext_tig&& t) {
        ext_tig_serialize(t, m_batch[b]);
        if (m_batch[b].size() >= FRAME_BUDGET) flush_frame(b, /*output=*/true);
    }

    // Decode input bucket b fully into a vector; delete its file.
    std::vector<ext_tig> take_input_bucket(uint32_t b) {
        std::vector<ext_tig> out;
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
                throw std::runtime_error("ext-stitch round bucket LZ4 decode failed: " + path);
        }
        std::fclose(f);
        std::remove(path.c_str());
        size_t pos = 0;
        while (pos < raw.size()) {
            ext_tig t;
            size_t got = ext_tig_deserialize(raw.data() + pos, raw.size() - pos, t);
            if (got == 0) break;
            pos += got;
            out.push_back(std::move(t));
        }
        return out;
    }

    // Finish the output round (flush + EOF-mark all bucket files), then
    // make it the next round's input and open fresh output files.
    void advance() {
        for (uint32_t b = 0; b < m_num_buckets; ++b) {
            if (!m_batch[b].empty()) flush_frame(b, /*output=*/true);
            if (m_files[b]) {
                uint32_t eof = 0;
                std::fwrite(&eof, sizeof(eof), 1, m_files[b]);
                std::fclose(m_files[b]);
                m_files[b] = nullptr;
            }
        }
        m_in_round = m_out_round;
        m_out_round = m_in_round + 1;
        open_round_files(m_out_round);
    }

private:
    static constexpr size_t FRAME_BUDGET = 4u * 1024 * 1024;

    std::string bucket_path(uint32_t round, uint32_t b) const {
        return m_dir + "/stitch_r" + std::to_string(round) + "_b" + std::to_string(b) + ".bin";
    }

    void open_round_files(uint32_t round) {
        m_out_round = round;
        for (uint32_t b = 0; b < m_num_buckets; ++b) m_files[b] = nullptr;
        // Lazily open on first emit to avoid creating num_buckets empty
        // files every round (most buckets are empty in late rounds).
    }

    void flush_frame(uint32_t b, bool /*output*/) {
        auto& batch = m_batch[b];
        if (batch.empty()) return;
        if (!m_files[b]) {
            std::string path = bucket_path(m_out_round, b);
            m_files[b] = std::fopen(path.c_str(), "wb");
            if (!m_files[b])
                throw std::runtime_error("cannot open stitch round file: " + path + ": " +
                                         std::strerror(errno));
        }
        int src = (int)batch.size();
        int bound = LZ4_compressBound(src);
        if (m_out.size() < (size_t)bound) m_out.resize((size_t)bound);
        int comp = LZ4_compress_default((char const*)batch.data(), (char*)m_out.data(), src,
                                        (int)m_out.size());
        if (comp <= 0) throw std::runtime_error("ext-stitch LZ4 compress failed");
        uint32_t u = (uint32_t)src, c = (uint32_t)comp;
        std::FILE* f = m_files[b];
        if (std::fwrite(&u, sizeof(u), 1, f) != 1 || std::fwrite(&c, sizeof(c), 1, f) != 1 ||
            std::fwrite(m_out.data(), 1, (size_t)comp, f) != (size_t)comp)
            throw std::runtime_error("short write to stitch round file");
        batch.clear();
    }

    std::string m_dir;
    uint32_t m_num_buckets;
    std::vector<std::vector<uint8_t>> m_batch;  // per-bucket pending output bytes
    std::vector<std::FILE*> m_files;            // per-bucket OUTPUT file handles
    std::vector<uint8_t> m_out;                 // LZ4 scratch
    uint32_t m_in_round = 0;
    uint32_t m_out_round = 0;
};

// Route a tig to a round-output bucket by its chosen open end's
// junction, recording the keyed side in toggle-independent form. Returns
// false (and does nothing) if the tig is fully closed -- caller sinks it.
template <typename Store>
inline bool ext_route(ext_tig&& t, uint32_t k, Store& store, uint64_t& joined_counter) {
    (void)joined_counter;
    uint8_t side;
    if (!ext_choose_side(t, side)) return false;  // fully closed
    bool fwd;
    kmer_int_t j = side_junction_canonical(std::string_view(t.seq.data(), t.seq.size()), k, side,
                                           fwd);
    uint32_t b = (uint32_t)((kmer_hasher{}(j) >> 1) % store.num_buckets());
    store.emit(b, std::move(t));
    return true;
}

// Store-templated round loop. `store` is seeded (round-0 input ready);
// runs doubling rounds until none join, sinking finished/maximal/
// circular unitigs. Works for both round_store_mem and
// round_store_file (same emit / take_input_bucket / advance interface).
template <typename Store, typename Sink>
inline void ext_run_rounds(Store& store, uint32_t k, uint32_t num_buckets,
                           joinable_set_t const& joinable, Sink&& sink,
                           std::atomic<uint64_t>* done) {
    // Whether `side` of tig t sits at a joinable junction. A tig is
    // TERMINAL (a complete unitig) iff none of its open ends is
    // joinable; such a tig is sunk immediately. Because joinability is a
    // fixed global property, a non-joinable open end never becomes
    // joinable, so nothing dangles -- this is what makes the drain
    // terminate (fixes the earlier hang).
    auto end_joinable = [&](ext_tig const& t, uint8_t side) -> bool {
        bool fwd;
        kmer_int_t j = detail::side_junction_canonical(
            std::string_view(t.seq.data(), t.seq.size()), k, side, fwd);
        return joinable.find(j) != joinable.end();
    };
    auto is_terminal = [&](ext_tig const& t) -> bool {
        if ((t.open_flags & UNITIG_OPEN_LEFT) && end_joinable(t, SIDE_LEFT)) return false;
        if ((t.open_flags & UNITIG_OPEN_RIGHT) && end_joinable(t, SIDE_RIGHT)) return false;
        return true;
    };
    auto sink_unitig = [&](ext_tig& t) {
        if (done) done->fetch_add(1, std::memory_order_relaxed);
        stitchable_unitig u;
        u.cid = t.cid;
        u.open_flags = 0;  // complete unitig
        u.seq = std::move(t.seq);
        sink(std::move(u));
    };

    // Loop until the store drains. Each round joins co-bucketed pairs at
    // joinable junctions, sinks terminal tigs, and carries the rest
    // (re-rolling which end a both-open tig presents). Termination: no
    // tig is carried (store empty after advance) -- guaranteed because
    // every join reduces the live tig count and terminal tigs leave the
    // loop, while non-joinable ends are classified terminal up front.
    constexpr uint32_t MAX_ROUNDS = 4096;  // safety net; O(log L) expected
    for (uint32_t round_no = 0;; ++round_no) {
        uint64_t carried_this_round = 0;

        for (uint32_t b = 0; b < num_buckets; ++b) {
            std::vector<ext_tig> tigs = store.take_input_bucket(b);
            if (tigs.empty()) continue;

            const uint64_t NT = tigs.size();
            std::vector<ext_end> ends(NT);
            std::vector<uint8_t> has_end(NT, 0);
            std::vector<uint8_t> consumed(NT, 0);
            // Map each presented joinable-end junction to the first tig
            // that presented it this round; the second arrival joins.
            ankerl::unordered_dense::map<kmer_int_t, uint64_t, kmer_hasher> waiting;
            waiting.reserve(NT);

            for (uint64_t i = 0; i < NT; ++i) {
                uint8_t side;
                if (!ext_choose_side(tigs[i], side)) continue;  // closed (shouldn't reach)
                bool fwd;
                kmer_int_t j = detail::side_junction_canonical(
                    std::string_view(tigs[i].seq.data(), tigs[i].seq.size()), k, side, fwd);
                ends[i] = ext_end{i, j, side, fwd};
                has_end[i] = 1;
            }

            // Pairing pass: a tig presenting a joinable junction either
            // becomes the waiter or joins the existing waiter. The
            // global `joinable` set already guarantees exactly-2-ends +
            // equal-cid + compatible-orientation, so a present pair at a
            // joinable junction is always a valid join.
            for (uint64_t i = 0; i < NT; ++i) {
                if (consumed[i] || !has_end[i]) continue;
                kmer_int_t j = ends[i].junction;
                if (joinable.find(j) == joinable.end()) continue;  // terminal end
                auto it = waiting.find(j);
                if (it == waiting.end()) {
                    waiting.emplace(j, i);
                    continue;
                }
                uint64_t a_idx = it->second;
                if (consumed[a_idx]) {  // stale; replace waiter
                    it->second = i;
                    continue;
                }
                ext_tig merged = ext_join(tigs[a_idx], ends[a_idx].side, tigs[i], ends[i].side, k);
                consumed[a_idx] = 1;
                consumed[i] = 1;
                waiting.erase(it);
                if (is_terminal(merged)) {
                    sink_unitig(merged);
                } else {
                    uint64_t dummy = 0;
                    detail::ext_route(std::move(merged), k, store, dummy);
                    ++carried_this_round;
                }
            }

            // Survivors: sink terminal tigs; carry the rest (re-roll the
            // presented end for both-open tigs so each end eventually
            // meets its partner).
            for (uint64_t i = 0; i < NT; ++i) {
                if (consumed[i]) continue;
                if (is_terminal(tigs[i])) {
                    sink_unitig(tigs[i]);
                    continue;
                }
                ext_tig t = std::move(tigs[i]);
                bool both = (t.open_flags & UNITIG_OPEN_LEFT) &&
                            (t.open_flags & UNITIG_OPEN_RIGHT);
                if (both) detail::ext_rng_next(t.rng);  // re-roll presented end
                uint64_t dummy = 0;
                detail::ext_route(std::move(t), k, store, dummy);
                ++carried_this_round;
            }
        }

        store.advance();
        if (carried_this_round == 0) break;  // store drained
        if (round_no >= MAX_ROUNDS) {
            std::cerr << "[ext-stitch] WARNING: round cap " << MAX_ROUNDS << " hit; flushing "
                      << carried_this_round << " survivors\n";
            for (uint32_t b = 0; b < num_buckets; ++b) {
                std::vector<ext_tig> rem = store.take_input_bucket(b);
                for (auto& t : rem) sink_unitig(t);
            }
            break;
        }
    }
}

// Seed round 0 from the fragment source into `store`, sinking already-
// closed fragments. Shared by the in-RAM and file-backed entry points.
template <typename Source, typename Store, typename Sink>
inline void ext_seed_round0(Source& frag, uint32_t k, Store& store, Sink&& sink,
                            joinable_set_t& joinable, std::atomic<uint64_t>* done) {
    const uint64_t n = (uint64_t)frag.size();

    // Compute the JOINABLE junction set (the in-RAM stitch's exact join
    // predicate, precomputed once). For each (k-1) junction, track its
    // open ends: count (capped at 3), and -- for the first 2 -- their
    // cid and side/orientation. A junction is joinable iff it has
    // exactly 2 ends, equal cid, and ext_pair_compatible holds. Every
    // other open end is a terminal unitig boundary.
    //
    // Invariant under merging: a junction's two ends are always the same
    // two original-fragment ends (merges relocate only a tig's OTHER,
    // outer end), so this set is fixed for the whole run.
    {
        struct jinfo {
            uint8_t count = 0;
            uint64_t cid[2] = {0, 0};
            uint8_t side[2] = {0, 0};
            bool fwd[2] = {false, false};
        };
        ankerl::unordered_dense::map<kmer_int_t, jinfo, kmer_hasher> jmap;
        for (uint64_t i = 0; i < n; ++i) {
            uint8_t fl = frag.open_flags(i);
            if (fl == 0) continue;
            std::string_view sv = frag.seq_view(i);
            if ((uint32_t)sv.size() < k) continue;
            uint64_t cid = frag.cid(i);
            for (uint8_t side : {SIDE_LEFT, SIDE_RIGHT}) {
                uint8_t bit = (side == SIDE_LEFT) ? UNITIG_OPEN_LEFT : UNITIG_OPEN_RIGHT;
                if (!(fl & bit)) continue;
                bool fwd;
                kmer_int_t j = side_junction_canonical(sv, k, side, fwd);
                jinfo& info = jmap[j];
                if (info.count < 2) {
                    info.cid[info.count] = cid;
                    info.side[info.count] = side;
                    info.fwd[info.count] = fwd;
                }
                if (info.count < 3) ++info.count;
            }
        }
        for (auto const& kv : jmap) {
            jinfo const& info = kv.second;
            if (info.count != 2) continue;
            if (info.cid[0] != info.cid[1]) continue;
            ext_end ea{0, kv.first, info.side[0], info.fwd[0]};
            ext_end eb{0, kv.first, info.side[1], info.fwd[1]};
            if (!ext_pair_compatible(ea, eb)) continue;
            joinable.insert(kv.first);
        }
    }

    for (uint64_t i = 0; i < n; ++i) {
        ext_tig t;
        t.cid = frag.cid(i);
        t.open_flags = frag.open_flags(i);
        std::string_view sv = frag.seq_view(i);
        t.seq.assign(sv.data(), sv.size());
        t.rng = ext_seed_from_seq(t.seq);
        if (k < 2 || t.open_flags == 0) {
            if (done) done->fetch_add(1, std::memory_order_relaxed);
            stitchable_unitig u;
            u.cid = t.cid;
            u.open_flags = t.open_flags;
            u.seq = std::move(t.seq);
            sink(std::move(u));
            continue;
        }
        uint64_t dummy = 0;
        ext_route(std::move(t), k, store, dummy);
    }
    store.advance();  // round-0 emissions become round-1 input
}

}  // namespace detail

// External-memory stitch with IN-RAM round storage. Fast; used by tests
// and as the reference. Same Source/Sink interface as
// stitch_unitigs_streaming. `num_buckets` controls hash fan-out; 0 picks
// a default.
template <typename Source, typename Sink>
inline void stitch_unitigs_extmem(Source& frag, uint32_t k, Sink&& sink,
                                  uint32_t num_buckets = 0,
                                  std::atomic<uint64_t>* done = nullptr) {
    if (num_buckets == 0) num_buckets = 256;
    detail::round_store_mem store(num_buckets);
    detail::joinable_set_t joinable;
    detail::ext_seed_round0(frag, k, store, sink, joinable, done);
    detail::ext_run_rounds(store, k, num_buckets, joinable, sink, done);
}

// External-memory stitch with FILE-backed round storage. Peak RAM is one
// bucket's tigs at a time (plus per-bucket write batches), independent
// of fragment count. Round files live under `tmp_dir`. This is the
// production path for inputs too large for the in-RAM stitch.
template <typename Source, typename Sink>
inline void stitch_unitigs_extmem_file(Source& frag, uint32_t k, std::string const& tmp_dir,
                                       Sink&& sink, uint32_t num_buckets,
                                       std::atomic<uint64_t>* done = nullptr) {
    if (num_buckets == 0) num_buckets = 1024;
    detail::round_store_file store(tmp_dir, num_buckets);
    detail::joinable_set_t joinable;
    detail::ext_seed_round0(frag, k, store, sink, joinable, done);
    detail::ext_run_rounds(store, k, num_buckets, joinable, sink, done);
}

}  // namespace cdgb
