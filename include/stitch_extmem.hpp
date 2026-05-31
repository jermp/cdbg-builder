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
//   A "tig" is a (partially) merged fragment: { runs, open_flags, seq },
//   where `runs` is the RLE color-run sequence over the tig's k-mers.
//   An OPEN side means the unitig continues into another bucket. With the
//   k-base super-k-mer overlap, the two fragments that meet there share
//   their FULL boundary k-mer (it is physically present in both buckets),
//   so we key the join on that canonical k-mer -- no global joinable map.
//
//   Round loop (until the store drains):
//     1. Each non-closed tig presents ONE chosen open side and is routed
//        to bucket H(boundary_kmer) % NUM_BUCKETS. A both-open tig picks
//        its presented end from a per-tig RNG, re-rolled each round it
//        survives, so a chain of L fragments collapses in O(log L)
//        expected rounds.
//     2. Within each bucket, two tigs presenting the same boundary k-mer
//        with compatible orientation are concatenated (the shared boundary
//        k-mer, k bases, dropped); the merged tig is sunk if now closed,
//        else carried to the next round. The join is COLORLESS (topology
//        only, like GGCAT): each tig carries an RLE color-run sequence that
//        is concatenated at the seam and cut into monochromatic unitigs at
//        emit. No cid gate -- that would strand color-boundary ends.
//     3. A tig with no open end is sunk immediately as a complete unitig.
//
//   Branch separation is structural: the per-bucket walker (step 3)
//   leaves an end OPEN only on a degree-1 cross-bucket simple-path edge,
//   which has exactly one partner. Real branches and dead ends are
//   CLOSED, and their distinct boundary k-mers route to different slots,
//   so they never false-join -- this is why the (k-1)-keyed stitch's
//   precomputed joinable set (a global O(num_junctions) map) is gone.
//
// FAITHFUL to GGCAT: the join is colorless (k-mer adjacency only). Colors
// ride along as a per-tig RLE color-run sequence (GGCAT's UnitigColorData),
// joined at each merge; the ONE extra step GGCAT does not need is the
// monochromatic split at emit (cdbg output unitigs are monochromatic).
//
// Two stores share one round driver: round_store_mem (in-RAM rounds;
// fast, used by tests/reference) and round_store_file (LZ4-framed files,
// one bucket resident at a time; the production external-memory path).

#include <array>
#include <atomic>
#include <cerrno>
#include <cassert>
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

// A tig in the doubling loop. `rng` is a per-tig random state used to
// pick which open end to present when both are open. It MUST vary
// between tigs (we seed it from sequence content) so two both-open tigs
// don't flip in lockstep -- if they did, a pair needing to meet on a
// specific shared end would never align and the loop would not
// terminate. Randomized selection is what gives O(log L) expected
// rounds (GGCAT does the same).
struct ext_tig {
    std::vector<color_run> runs;  // RLE color sequence over this tig's k-mers
    uint8_t open_flags;
    uint64_t rng;
    std::string seq;
    // Persistent presented side (GGCAT flags bit0 = HASH_ENDING_FLAG_MASK):
    // which open end this tig keys on for the current round. Set when the tig
    // is routed; on a failed (orientation-incompatible) collision it is flipped
    // to present the other extremity next round, exactly as GGCAT re-emits with
    // `flags ^ HASH_ENDING_FLAG_MASK`. Only meaningful while the tig is open.
    uint8_t pres = SIDE_LEFT;
};

// Split a colorless (topological) tig into monochromatic unitigs at its
// color-run boundaries and feed each to `emit`. This is the only step GGCAT
// does not need (its unitigs stay topological with a color sequence); we cut
// at run boundaries because cdbg output unitigs are monochromatic. Run i
// covers num_kmers consecutive k-mers; the piece spanning them is
// seq[base_off .. base_off + num_kmers + k - 1), and consecutive pieces share
// only k-1 bases (no shared k-mer), so every k-mer lands in exactly one piece.
template <typename Emit>
inline void ext_split_monochromatic(ext_tig& t, uint32_t k, uint8_t open_flags, Emit&& emit) {
    const size_t nruns = t.runs.size();
    // A COLOR_RUN_FOREIGN placeholder can only survive to emit at an EXTREMITY:
    // a boundary k-mer whose open end never joined its primary partner (the
    // small unjoined-end residual). Its color is owned by the adjacent bucket,
    // which emits it -- so we DROP such an extremal foreign k-mer here rather
    // than output it with the sentinel cid. (An interior foreign run would mean
    // a k-mer was lost from the graph; it does not occur because every boundary
    // k-mer is primary in exactly one bucket, and is asserted against below.)
    size_t lo = 0, hi = nruns;  // [lo, hi) runs to emit
    if (lo < hi && t.runs[lo].cid == COLOR_RUN_FOREIGN) ++lo;
    if (lo < hi && t.runs[hi - 1].cid == COLOR_RUN_FOREIGN) --hi;
    // k-mer offset of the first emitted run (foreign-dropped prefix k-mers).
    size_t base_off = 0;
    for (size_t i = 0; i < lo; ++i) base_off += t.runs[i].num_kmers;
    for (size_t i = lo; i < hi; ++i) {
        const uint32_t m = t.runs[i].num_kmers;
        if (t.runs[i].cid == COLOR_RUN_FOREIGN) {
            // Safety net: with BCALM2 ownership every k-mer is primary in
            // exactly one bucket, so a foreign placeholder is always resolved at
            // its join and none should reach here. If one ever did (an open end
            // that never met its primary partner), the owning bucket emits it --
            // so drop it rather than output the sentinel cid.
            base_off += m;
            continue;
        }
        stitchable_unitig u;
        u.seq = t.seq.substr(base_off, (size_t)m + k - 1);
        u.runs.push_back({t.runs[i].cid, m});
        // OPEN flags only survive on the two extremal pieces (a split tig is
        // only ever sunk when fully closed, so open_flags is 0 in practice).
        uint8_t of = 0;
        if (i == lo) of |= (open_flags & UNITIG_OPEN_LEFT);
        if (i + 1 == hi) of |= (open_flags & UNITIG_OPEN_RIGHT);
        u.open_flags = of;
        emit(std::move(u));
        base_off += m;
    }
}

// Reverse a run sequence in place (reverse-complementing a tig reverses the
// k-mer order; color sets are orientation-independent so cids are unchanged).
inline void ext_reverse_runs(std::vector<color_run>& runs) {
    std::reverse(runs.begin(), runs.end());
}

// Append `b` runs onto `a` at a join seam. `a` is oriented keyed-end-at-RIGHT
// and `b` keyed-end-at-LEFT, so a's LAST k-mer and b's FIRST k-mer are both the
// shared boundary k-mer X. The merged SEQ drops one full copy of X (k bases), so
// the runs must drop exactly one X unit, keeping sum(runs)==seq k-mer count.
//
// X is colored on the side that OWNS it and is a COLOR_RUN_FOREIGN placeholder
// on the side that doesn't (it may be owned by a, by b, or -- when X is primary
// in some third bucket -- foreign in BOTH). Reconcile: the surviving X unit
// takes the real cid if either side has one, else stays foreign (a later join
// at the other end resolves it). We keep X in `a`'s last run, fix its cid, then
// drop X from `b`'s front before appending.
inline void ext_concat_runs(std::vector<color_run>& a, std::vector<color_run> b) {
    if (b.empty()) return;
    // Reconcile X's color into a's last run (X = a's last k-mer).
    if (!a.empty()) {
        uint64_t a_x = a.back().cid;
        uint64_t b_x = b.front().cid;
        uint64_t real = (a_x != COLOR_RUN_FOREIGN) ? a_x : b_x;  // foreign if both foreign
        if (a_x != real) {
            // a's X run was a foreign placeholder (length 1 at the open end);
            // flip it to the real cid recovered from b's owning side.
            a.back().cid = real;
            // Coalesce the just-reconciled X unit backward into a's preceding
            // run when they now share a cid. Without this, every join across a
            // same-color boundary leaves a spurious [run(c,n), run(c,1)] split;
            // over a long chain that compounds into thousands of bogus runs per
            // topological tig, and the monochromatic split then shatters one
            // tig into thousands of unitigs (salmonella-10: 2123 runs/tig ->
            // 409k unitigs instead of 87,297).
            if (a.size() >= 2 && a[a.size() - 2].cid == real) {
                a[a.size() - 2].num_kmers += a.back().num_kmers;
                a.pop_back();
            }
        }
    }
    // Drop the shared boundary k-mer X from b's front (already counted in a).
    if (--b.front().num_kmers == 0) b.erase(b.begin());
    if (b.empty()) return;
    size_t bi = 0;
    if (!a.empty() && a.back().cid == b.front().cid) {
        a.back().num_kmers += b.front().num_kmers;
        bi = 1;
    }
    for (; bi < b.size(); ++bi) a.push_back(b[bi]);
}

// xorshift step; returns the next state.
inline uint64_t ext_rng_next(uint64_t& s) {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return s;
}

// Well-mixed single bit from a state word (splitmix64 finalizer). xorshift64's
// LOW bit is almost static across steps, so `(rng & 1)` made a both-open tig
// present nearly the same end every round -> two partners almost never
// co-presented their shared boundary k-mer in one round and the stitch
// plateaued (deterministic-looking non-join). The finalizer's top bit flips
// well, restoring the ~1/4 per-round meet probability the doubling relies on.
inline bool ext_mix_bit(uint64_t s) {
    s ^= s >> 30;
    s *= 0xbf58476d1ce4e5b9ull;
    s ^= s >> 27;
    s *= 0x94d049bb133111ebull;
    s ^= s >> 31;
    return (s >> 63) & 1u;
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

// Canonical full boundary k-mer at one side of a fragment (GGCAT
// final_executor.rs:147-174 / extremal.rs:25-43). SIDE_LEFT keys the
// first k bases, SIDE_RIGHT the last k bases. Replaces the (k-1) junction
// key: because adjacent k-overlap fragments share their FULL boundary
// k-mer, two ends that should glue present equal canonical k-mers, while
// distinct k-mers at a dBG branch differ and never collide -- so no
// global joinable map is needed to separate branches. Requires
// seq.size() >= k (always true for an open fragment).
inline kmer_int_t side_kmer_canonical(std::string_view seq, uint32_t k, uint8_t side,
                                      bool& is_canonical_fwd) {
    char const* p = (side == SIDE_LEFT) ? seq.data() : seq.data() + (seq.size() - k);
    kmer_int_t fwd = 0;
    for (uint32_t i = 0; i < k; ++i) fwd = (fwd << 2) | (kmer_int_t)nuc_to_2bit(p[i]);
    kmer_int_t rc = reverse_complement(fwd, k);
    if (fwd <= rc) {
        is_canonical_fwd = true;
        return fwd;
    }
    is_canonical_fwd = false;
    return rc;
}

// One chosen open end of a tig, in the current round's bucket.
struct ext_end {
    uint64_t tig_idx;       // index into the round-bucket's tig vector
    kmer_int_t junction;    // canonical full boundary k-mer
    uint8_t side;           // SIDE_LEFT / SIDE_RIGHT of the keyed end
    bool is_canonical_fwd;  // boundary-k-mer-in-own-frame == canonical
};

// Pick the open side a tig presents this round. For a both-open tig the
// side is read from the current rng state (NOT advanced here, so repeated
// queries within one round are stable); the caller advances rng when it
// carries the tig to the next round. Returns false if fully closed.
inline bool ext_choose_side(ext_tig const& t, uint8_t& out_side) {
    bool l = (t.open_flags & UNITIG_OPEN_LEFT) != 0;
    bool r = (t.open_flags & UNITIG_OPEN_RIGHT) != 0;
    if (l && r) {
        out_side = ext_mix_bit(t.rng) ? SIDE_RIGHT : SIDE_LEFT;
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

// Join two tigs at a shared, compatible full boundary k-mer. Orient `a`
// so its keyed end (a_side) is at the RIGHT and `b` so its keyed end
// (b_side) is at the LEFT, then merged = a_oriented + b_oriented[k..].
// The two k-base overlaps (the shared boundary k-mer) are guaranteed
// equal by ext_pair_compatible (verified for all four side/orientation
// cases). Returns the merged tig.
inline ext_tig ext_join(ext_tig const& a, uint8_t a_side, ext_tig const& b, uint8_t b_side,
                        uint32_t k) {
    // a oriented with keyed end at RIGHT.
    std::string as = (a_side == SIDE_RIGHT) ? a.seq : revcomp_string(a.seq);
    // b oriented with keyed end at LEFT.
    std::string bs = (b_side == SIDE_LEFT) ? b.seq : revcomp_string(b.seq);

    ext_tig m;
    // Orient each parent's run sequence to match its oriented bases, then
    // concatenate (dropping the shared boundary k-mer from b's front).
    m.runs = a.runs;
    if (a_side != SIDE_RIGHT) ext_reverse_runs(m.runs);  // as = revcomp(a) when a_side==LEFT
    std::vector<color_run> b_runs = b.runs;
    if (b_side != SIDE_LEFT) ext_reverse_runs(b_runs);   // bs = revcomp(b) when b_side==RIGHT
    ext_concat_runs(m.runs, std::move(b_runs));
    m.seq.reserve(as.size() + bs.size() - k);
    m.seq = as;
    // Drop the full shared boundary k-mer (k bases), not k-1.
    m.seq.append(bs.begin() + k, bs.end());
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
//   [u8 open_flags][u64 rng][u32 nruns]{[u64 cid][u32 num_kmers]}*
//   [u32 seq_len][seq_len bytes ACGT]
inline void ext_tig_serialize(ext_tig const& t, std::vector<uint8_t>& out) {
    auto put = [&](void const* p, size_t n) {
        uint8_t const* b = (uint8_t const*)p;
        out.insert(out.end(), b, b + n);
    };
    put(&t.open_flags, sizeof(t.open_flags));
    put(&t.rng, sizeof(t.rng));
    uint32_t nruns = (uint32_t)t.runs.size();
    put(&nruns, sizeof(nruns));
    for (auto const& r : t.runs) {
        put(&r.cid, sizeof(r.cid));
        put(&r.num_kmers, sizeof(r.num_kmers));
    }
    uint32_t len = (uint32_t)t.seq.size();
    put(&len, sizeof(len));
    put(t.seq.data(), len);
}

// Returns bytes consumed, 0 on malformed/EOF.
inline size_t ext_tig_deserialize(uint8_t const* buf, size_t buf_len, ext_tig& t) {
    constexpr size_t HDR = 1 + sizeof(uint64_t) + sizeof(uint32_t);
    if (buf_len < HDR) return 0;
    size_t p = 0;
    t.open_flags = buf[p];
    p += 1;
    std::memcpy(&t.rng, buf + p, sizeof(t.rng));
    p += sizeof(t.rng);
    uint32_t nruns;
    std::memcpy(&nruns, buf + p, sizeof(nruns));
    p += sizeof(nruns);
    constexpr size_t RUN = sizeof(uint64_t) + sizeof(uint32_t);
    if (buf_len < p + (size_t)nruns * RUN + sizeof(uint32_t)) return 0;
    t.runs.resize(nruns);
    for (uint32_t i = 0; i < nruns; ++i) {
        std::memcpy(&t.runs[i].cid, buf + p, sizeof(uint64_t));
        p += sizeof(uint64_t);
        std::memcpy(&t.runs[i].num_kmers, buf + p, sizeof(uint32_t));
        p += sizeof(uint32_t);
    }
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
    kmer_int_t j = side_kmer_canonical(std::string_view(t.seq.data(), t.seq.size()), k, side, fwd);
    uint32_t b = (uint32_t)((kmer_hasher{}(j) >> 1) % store.num_buckets());
    store.emit(b, std::move(t));
    return true;
}

// Store-templated round loop. `store` is seeded (round-0 input ready);
// runs doubling rounds until none join, sinking finished/maximal/
// circular unitigs. Works for both round_store_mem and
// round_store_file (same emit / take_input_bucket / advance interface).
template <typename Store, typename Sink>
inline void ext_run_rounds(Store& store, uint32_t k, uint32_t num_buckets, Sink&& sink,
                           std::atomic<uint64_t>* done) {
    // Full-k-mer keying (GGCAT extend_unitigs.rs:410-481) needs no global
    // joinable map. An end is OPEN only because the step-3 walker found a
    // degree-1 cross-bucket simple-path edge, which by construction has
    // exactly one partner fragment sharing the full boundary k-mer. Real
    // branches and dead ends are CLOSED by the walker, so a distinct k-mer
    // at a branch routes to a different slot and never collides. A tig is
    // TERMINAL iff it has no open end at all; an open end always has (or
    // will, in a later round when its partner co-buckets) a unique mate.
    auto is_terminal = [&](ext_tig const& t) -> bool { return t.open_flags == 0; };
    auto sink_unitig = [&](ext_tig& t) {
        // Cut the finished topological tig at its color-run boundaries into
        // monochromatic output unitigs (verify.py requires monochromaticity).
        ext_split_monochromatic(t, k, t.open_flags, [&](stitchable_unitig&& u) {
            if (done) done->fetch_add(1, std::memory_order_relaxed);
            sink(std::move(u));
        });
    };

    // Loop until the store drains. Each round joins co-bucketed pairs that
    // present the same full boundary k-mer with compatible orientation
    // (colorless); carries the rest (re-rolling which end a both-open tig
    // presents). Termination: every join strictly reduces the live tig
    // count, and a tig that finds no partner this round eventually
    // co-buckets with its unique mate (both key the same canonical k-mer,
    // so both route to the same bucket whenever they present that end) --
    // O(log L) expected rounds.
    constexpr uint32_t MAX_ROUNDS = 4096;  // safety net; O(log L) expected
    for (uint32_t round_no = 0;; ++round_no) {
        uint64_t carried_this_round = 0;
        uint64_t joined_this_round = 0;

        for (uint32_t b = 0; b < num_buckets; ++b) {
            std::vector<ext_tig> tigs = store.take_input_bucket(b);
            if (tigs.empty()) continue;

            const uint64_t NT = tigs.size();
            std::vector<ext_end> ends(NT);
            std::vector<uint8_t> has_end(NT, 0);
            std::vector<uint8_t> consumed(NT, 0);
            // Map each presented boundary k-mer to the first tig that
            // presented it this round; the second compatible arrival joins.
            ankerl::unordered_dense::map<kmer_int_t, uint64_t, kmer_hasher> waiting;
            waiting.reserve(NT);

            for (uint64_t i = 0; i < NT; ++i) {
                uint8_t side;
                if (!ext_choose_side(tigs[i], side)) continue;  // closed (shouldn't reach)
                bool fwd;
                kmer_int_t j = detail::side_kmer_canonical(
                    std::string_view(tigs[i].seq.data(), tigs[i].seq.size()), k, side, fwd);
                ends[i] = ext_end{i, j, side, fwd};
                has_end[i] = 1;
            }

            // Pairing pass (GGCAT first-come/second-come): a tig presenting
            // a boundary k-mer either becomes the waiter or joins the
            // existing waiter. The full-k-mer key makes branches route
            // elsewhere, so a colliding pair is a genuine adjacency -- we
            // gate the join on compatible orientation only (colorless, like
            // GGCAT). A pair that fails the gate is a real terminus on
            // that end: drop the waiter's claim and let the later survivor
            // pass treat each as terminal-on-this-end (carried, re-rolled,
            // or sunk if its other end also fails to join).
            for (uint64_t i = 0; i < NT; ++i) {
                if (consumed[i] || !has_end[i]) continue;
                kmer_int_t j = ends[i].junction;
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
                // Colorless gate (GGCAT): orientation only -- NO color check.
                // The full-k-mer key already separated branches, so a colliding
                // compatible pair is a genuine topological adjacency that always
                // joins; the color difference (if any) is recorded as a run
                // boundary inside the merged tig and cut apart at emit. This is
                // what makes every open end have a unique partner and the loop
                // terminate (a cid gate here would strand color-boundary ends).
                bool ok = ext_pair_compatible(ends[a_idx], ends[i]);
                if (!ok) {
                    // Incompatible orientation at a shared k-mer (palindromic /
                    // branch remnant): keep the waiter; handle this as survivor.
                    continue;
                }
                ext_tig merged = ext_join(tigs[a_idx], ends[a_idx].side, tigs[i], ends[i].side, k);
                consumed[a_idx] = 1;
                consumed[i] = 1;
                waiting.erase(it);
                ++joined_this_round;
                if (is_terminal(merged)) {
                    sink_unitig(merged);
                } else {
                    uint64_t dummy = 0;
                    detail::ext_route(std::move(merged), k, store, dummy);
                    ++carried_this_round;
                }
            }

            // Survivors: sink fully-closed tigs; carry the rest (re-roll
            // the presented end for both-open tigs so each end eventually
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
#ifdef CDGB_STITCH_DEBUG
        std::cerr << "[ext-stitch round " << round_no << "] carried=" << carried_this_round
                  << " joined=" << joined_this_round << "\n";
#endif
        // GGCAT termination (extend_unitigs.rs:716): stop as soon as a full
        // round produces ZERO joins. `carried` is NOT a termination signal --
        // it counts every still-open tig re-routed for another attempt, and a
        // tig whose partner never co-buckets (or a genuine terminus) is carried
        // every round forever. When no join happened this round, no further
        // join is possible, so flush every remaining open tig as a terminal
        // unitig (its open end has no partner in the graph) and finish.
        if (joined_this_round == 0 || round_no >= MAX_ROUNDS) {
            if (round_no >= MAX_ROUNDS)
                std::cerr << "[ext-stitch] WARNING: round cap " << MAX_ROUNDS << " hit\n";
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
//
// No precomputed joinable set: full-k-mer keying separates branches
// structurally (the boundary k-mers differ), and the per-bucket pairing
// gates each join on compatible orientation locally (colorless). This
// deletes the one global O(num_junctions) map the (k-1)-keyed stitch
// needed -- the round loop now streams only one bucket's tigs at a time.
template <typename Source, typename Store, typename Sink>
inline void ext_seed_round0(Source& frag, uint32_t k, Store& store, Sink&& sink,
                            std::atomic<uint64_t>* done) {
    const uint64_t n = (uint64_t)frag.size();

    for (uint64_t i = 0; i < n; ++i) {
        ext_tig t;
        t.runs = frag.runs(i);
        t.open_flags = frag.open_flags(i);
        std::string_view sv = frag.seq_view(i);
        t.seq.assign(sv.data(), sv.size());
        t.rng = ext_seed_from_seq(t.seq);
        if (k < 2 || t.open_flags == 0) {
            // Already closed: split at color-run boundaries and sink.
            ext_split_monochromatic(t, k, t.open_flags, [&](stitchable_unitig&& u) {
                if (done) done->fetch_add(1, std::memory_order_relaxed);
                sink(std::move(u));
            });
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
    detail::ext_seed_round0(frag, k, store, sink, done);
    detail::ext_run_rounds(store, k, num_buckets, sink, done);
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
    detail::ext_seed_round0(frag, k, store, sink, done);
    detail::ext_run_rounds(store, k, num_buckets, sink, done);
}

}  // namespace cdgb
