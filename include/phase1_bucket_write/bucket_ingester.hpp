#pragma once

// Minimizer-based bucketing ingest. Each input file is one color.
//
// For every ACGT-only run in every sequence, we slide a window of (k-m+1)
// m-mers and find the canonical-ntHash minimum. Whenever the window's minimum
// changes, the current super-k-mer ends and a new one begins. Each super-k-mer
// is appended to the bucket file selected by `(min_hash >> 1) % B`,
// matching GGCAT's bucket-id derivation in cn_nthash.rs.
//
// The output is a directory of B bucket files, each holding a stream of
// super-k-mer records (see super_kmer.hpp). Thread safety: each ingest worker
// holds a private per_thread_bucket_buffers; concurrency on the shared bucket
// files is one mutex per bucket.

#include <atomic>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

#include "phase1_bucket_write/bucket_io.hpp"
#include "kmer.hpp"
#include "phase1_bucket_write/minimizer.hpp"
#include "phase1_bucket_write/seq_reader.hpp"
#include "phase1_bucket_write/super_kmer.hpp"

namespace cdbg {

namespace detail {

// Emit super-k-mers for a single ACGT-only run of length L (>= k).
// `bases` holds 2-bit-encoded bases. Records are appended via `sink` to the
// appropriate bucket files.
inline void emit_super_kmers(uint8_t const* bases, uint32_t L, uint32_t k, uint32_t m,
                             uint32_t color, uint32_t num_buckets,
                             per_thread_bucket_buffers& sink)  //
{
    assert(L >= k);                // caller only emits runs with run_len >= k
    const uint32_t K = L - k + 1;  // number of k-mers; (k-1)-mers indexed by 0..K

    // Minimize over (k-1)-mers, NOT k-mers. A (k-1)-mer has (k-1)-m+1 = k-m m-mers, so the
    // window is W = k - m. This is what co-locates branches: every dBG edge
    // X->Y shares its junction (k-1)-mer S = suffix(X) = prefix(Y), and both
    // endpoints route to bucket(min(S)). A branch X->{Y1,Y2} shares one S, so
    // X, Y1 and Y2 all land in bucket(min(S)) -- the branch is locally visible
    // in one bucket. (k-mer minimizers, W = k-m+1, scattered the two arms into
    // different buckets so no bucket saw the branch -> strict-topology
    // internal-branch failures.) It also makes each k-mer's primary bucket
    // intrinsic: X is primary in bucket(min(prefix(X))), foreign elsewhere.
    const uint32_t W = k - m;  // window size, in m-mer indices (a (k-1)-mer)

    // Skip the bottom bit to match GGCAT's "uniqueness-flag" reservation,
    // then map the minimizer hash into [0, num_buckets) with a modulo.
    // num_buckets need not be a power of two (it's sized from the RAM model),
    // and the divide is negligible here: once per super-k-mer, dwarfed by the
    // per-base ntHash/minimizer scan.
    //
    // A multiply-shift range reduction (h * B >> 64) was tried and REVERTED:
    // it maps proportionally to the hash VALUE, so it needs h uniform over
    // [0, 2^64). But canonical_mhash = min(fwd, rc) is biased toward small
    // values (the min of two uniforms), and the >>1 caps it at 2^63 -- so the
    // reduction piled almost everything into a few low buckets (2353/11915
    // nonempty, max bucket 243 MiB), collapsing dedup. The modulo keys on the
    // LOW bits, which stay uniform even when the magnitude is skewed.
    auto bucket_of = [&](uint64_t h) -> uint32_t { return (uint32_t)((h >> 1) % num_buckets); };

    // O(1)-space sliding-window minimizer: it owns the ntHash roll and re-scans
    // on expiry (no ring buffer, no modulo). init() covers (k-1)-mer 0's window.
    windowed_min mq;
    if (!mq.init(bases, m, (int32_t)W)) return;  // run is ACGT-only, shouldn't happen

    uint32_t run_a = 0;  // first (k-1)-mer index of the running minimizer block
    bool is_first_super = true;
    bool cur_owns_first = true;        // first super's idx0 is the ACGT-begin (not a boundary)
    uint64_t cur_min = mq.min_hash();  // minimizer of (k-1)-mer 0
    uint32_t cur_bucket = bucket_of(cur_min);
    // Base offset (in `bases`) of the m-mer achieving the running minimizer. The
    // whole run shares one minimizer value, so this position is constant across
    // the run; captured whenever the running minimizer changes. Recorded per
    // super-k-mer as the §3.8 normalization anchor (relative to the super start).
    int32_t cur_min_mmer_pos = mq.cur_min_pos;

    auto emit_super = [&](uint32_t first_kmer, uint32_t last_kmer, uint32_t bucket,
                          bool is_run_begin, bool is_run_end, bool owns_first, bool owns_last,
                          int32_t min_mmer_pos) {
        uint32_t base_start = first_kmer;
        uint32_t base_len = (last_kmer - first_kmer) + k;
        uint8_t flags = 0;
        if (is_run_begin) flags |= SK_FLAG_IS_ACGT_BEGIN;
        if (is_run_end) flags |= SK_FLAG_IS_ACGT_END;
        if (owns_first) flags |= SK_FLAG_OWNS_FIRST;
        if (owns_last) flags |= SK_FLAG_OWNS_LAST;
        // Minimizer position relative to this super-k-mer's first base. The
        // minimizer m-mer lies within the run's (k-1)-mer windows, so it is
        // always inside [base_start, base_start + base_len - m].
        int32_t min_pos = min_mmer_pos - (int32_t)base_start;
        if (min_pos < 0) min_pos = 0;  // defensive; should not happen
        sink.append(bucket, flags, color, bases + base_start, base_len, (uint32_t)min_pos);
    };

    // Walk the (k-1)-mers j = 1..K. A maximal run [run_a..b] of equal-minimizer
    // (k-1)-mers becomes a super-k-mer spanning k-mers [run_a-1 .. b] (clamped
    // at the read start): the k-mers whose prefix OR suffix junction is in the
    // run. Consecutive supers therefore share exactly their boundary k-mer
    // X_b (k-overlap): it is the last k-mer of one super (its prefix junction
    // P_b is in this run) and the begin-ignored first k-mer of the next (its
    // suffix junction P_{b+1} starts the next run). bucket(min P_b) != bucket(min
    // P_{b+1}) is exactly when that k-mer lands in two buckets -- the BCALM2
    // "k-mer in two buckets iff its left and right minimizers differ" rule.
    //
    // Ownership (which bucket COLORS the shared k-mer): X = A.last = B.first is
    // colored in bucket(min(mA, mB)). At the split mA = cur_min, mB = new_min:
    // the ending super A owns its last X iff mA < mB; the starting super B owns
    // its first X iff mB < mA. (mA != mB at a split, so exactly one owns it.)
    for (uint32_t j = 1; j <= K; ++j) {
        mq.advance();  // slide to (k-1)-mer j's window (rightmost m-mer j+k-m-1)
        uint64_t new_min = mq.min_hash();  // minimizer of (k-1)-mer j
        if (new_min != cur_min) {
            uint32_t first_kmer = (run_a == 0) ? 0 : run_a - 1;
            bool owns_last = (cur_min < new_min);
            emit_super(first_kmer, j - 1, cur_bucket, is_first_super, /*is_run_end=*/false,
                       cur_owns_first, owns_last, cur_min_mmer_pos);
            is_first_super = false;
            cur_owns_first = (new_min < cur_min);  // B owns its first X iff mB < mA
            cur_min = new_min;
            cur_bucket = bucket_of(new_min);
            cur_min_mmer_pos = mq.cur_min_pos;  // minimizer m-mer position of the new run
            run_a = j;
        }
    }
    uint32_t first_kmer = (run_a == 0) ? 0 : run_a - 1;
    // Run-end super: its last k-mer is the ACGT-end terminus (not a cross-bucket
    // boundary), so it always owns it; its first k-mer ownership was decided at
    // the preceding split (cur_owns_first).
    emit_super(first_kmer, K - 1, cur_bucket, is_first_super, /*is_run_end=*/true, cur_owns_first,
               /*owns_last=*/true, cur_min_mmer_pos);
}

inline void ingest_file_bucketed(std::string const& path, uint32_t k, uint32_t m,
                                 uint32_t num_buckets, uint32_t color,
                                 per_thread_bucket_buffers& sink) {
    auto& prof = bucket_prof();
    seq_reader r(path);
    char const* s = nullptr;
    size_t l = 0;
    std::vector<uint8_t> bases_buf;
    // Accumulate the scan/decode timer in a local and flush it to the shared
    // atomic once per file (below), instead of fetch_add-ing the contended
    // global counter once per ACGT run in the hot loop.
    uint64_t scan_ns = 0;
    for (;;) {
        auto t_read = bucket_write_prof::clock::now();
        bool ok = r.next(s, l);
        prof.ns_seq_read.fetch_add(bucket_write_prof::since(t_read), std::memory_order_relaxed);
        if (!ok) break;
        auto t_body = bucket_write_prof::clock::now();
        // One scratch buffer per sequence, grown (never per-run resized) to the
        // sequence length; runs are written into its front [0, run_len).
        if (bases_buf.size() < l) bases_buf.resize(l);
        uint8_t* buf = bases_buf.data();
        size_t pos = 0;
        while (pos < l) {
            // Single pass: scan an ACGT-only run starting at pos, decoding each
            // base to 2-bit straight into `buf` and stopping at the first
            // non-ACGT char. (Previously a scan-for-end loop and a separate
            // decode loop each called nuc_to_2bit over the run -- two passes;
            // this is one.) This scan + decode is what a SIMD FASTX parser
            // replaces; time it separately from emit_super_kmers.
            auto t_scan = bucket_write_prof::clock::now();
            size_t run_len = 0;
            while (pos < l) {
                uint8_t b = nuc_to_2bit(s[pos]);
                if (b == 0xff) break;
                buf[run_len++] = b;
                ++pos;
            }
            scan_ns += bucket_write_prof::since(t_scan);
            if (run_len >= k) {
                emit_super_kmers(buf, (uint32_t)run_len, k, m, color, num_buckets, sink);
            }
            while (pos < l and nuc_to_2bit(s[pos]) == 0xff) ++pos;
        }
        prof.ns_loop_body.fetch_add(bucket_write_prof::since(t_body), std::memory_order_relaxed);
    }
    prof.ns_scan2bit.fetch_add(scan_ns, std::memory_order_relaxed);
    prof.n_files.fetch_add(1, std::memory_order_relaxed);
}

}  // namespace detail

// Parallel driver. Spawns `num_threads` workers, each pulling files from a
// shared queue. `done` (if non-null) is incremented after each file finishes.
inline void ingest_bucketed(std::vector<std::string> const& files, uint32_t k, uint32_t m,
                            uint32_t num_buckets, bucket_writer& writer, uint32_t num_threads,
                            std::atomic<uint64_t>* done = nullptr,
                            std::vector<uint64_t> const* file_sizes = nullptr,
                            std::atomic<uint64_t>* done_bytes = nullptr) {
    if (num_threads == 0) num_threads = 1;
    std::atomic<size_t> next{0};
    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    auto run = [&]() {
        per_thread_bucket_buffers bufs(writer);
        for (;;) {
            size_t i = next.fetch_add(1);
            if (i >= files.size()) break;
            try {
                detail::ingest_file_bucketed(files[i], k, m, num_buckets, (uint32_t)i, bufs);
            } catch (std::exception& e) {
                std::cerr << "error ingesting " << files[i] << ": " << e.what() << '\n';
            }
            if (done) done->fetch_add(1, std::memory_order_relaxed);
            // Advance the byte-based progress by this file's on-disk size, so
            // the bar tracks work (bytes) rather than file count -- files vary
            // ~3x in size, so file-% badly misreports true progress.
            if (done_bytes and file_sizes and i < file_sizes->size())
                done_bytes->fetch_add((*file_sizes)[i], std::memory_order_relaxed);
        }
        bufs.flush_all();
    };

    for (uint32_t t = 0; t < num_threads; ++t) workers.emplace_back(run);
    for (auto& w : workers) w.join();
}

}  // namespace cdbg
