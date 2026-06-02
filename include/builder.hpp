#pragma once

// Public API for building a colored compacted dBG.
//
// Wraps the full pipeline (minimizer-bucketed ingest -> per-bucket dBG
// build + global color interning -> stitch -> emit FASTA + .color_sets) so
// that downstream tools can construct a `build_config`, instantiate a
// `builder`, and call `build()`. Mirrors the call shape used by
// Fulgor's `index<ColorSets>::builder` for its ccdBG dependency
// (https://github.com/jermp/fulgor/blob/main/include/builders/builder.hpp).
//
// Typical use:
//
//   cdbg::build_config cfg;
//   cfg.filenames_list = "...";
//   cfg.out_basename   = "...";
//   cfg.k = 31;
//   cfg.num_threads = 8;
//
//   cdbg::builder b(cfg);
//   b.build();
//   // After build():
//   //   b.num_colors()         -- one per input file
//   //   b.num_unitigs()        -- count of stitched unitigs in <basename>.fa
//   //   b.num_color_classes()  -- count of distinct color sets in <basename>.color_sets
//
// build() throws std::runtime_error on configuration errors or I/O
// failures. On success, three artifacts are written:
//   <basename>.fa          colored unitigs in FASTA, headers = cid
//   <basename>.u2c         unitig-to-color-set bit_vector (run-end
//                          marker, popcount = num_color_classes)
//   <basename>.color_sets  hybrid-encoded color sets + EF offsets
// The scratch directory (cfg.tmp_dir or an mkdtemp'd one) is removed.

#include <atomic>
#include <cerrno>
#include <chrono>
#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include <unistd.h>

#if defined(__GLIBC__)
#include <malloc.h>
#endif

#include <bit_vector.hpp>
#include <essentials.hpp>

#include "bucket_io.hpp"
#include "bucket_ingester.hpp"
#include "bucket_walker.hpp"
#include "minimizer.hpp"
#include "stitch.hpp"
#include "stitch_extmem.hpp"
#include "streaming_color_set_dict.hpp"
#include "unitig_spill.hpp"
#include "util.hpp"

namespace cdbg {

struct builder {
    builder() = default;
    explicit builder(build_config const& cfg) : m_cfg(cfg) {}

    // Run the full pipeline. Throws std::runtime_error on bad config or
    // I/O error. On success, writes m_cfg.out_basename + {".fa", ".color_sets"}
    // and removes the scratch directory.
    void build() {
        validate_and_resolve_config();

        auto files = read_filenames(m_cfg.filenames_list);
        if (files.empty()) throw std::runtime_error("no input files");
        if (files.size() > (uint64_t)UINT32_MAX)
            throw std::runtime_error("too many colors (max 2^32 - 1)");
        m_num_colors = (uint32_t)files.size();

        uint32_t const num_buckets = 1u << m_cfg.bucket_log2;
        std::cout << "k = " << m_cfg.k << ", m = " << m_cfg.m << ", num_colors = " << m_num_colors
                  << ", num_threads = " << m_cfg.num_threads << ", num_buckets = " << num_buckets
                  << "\n";
        std::cout << "  bucket-write tuning: flush_bases=" << format_bytes(m_flush_bases)
                  << ", spill_bytes=" << format_bytes(m_spill_bytes);
        if (m_cfg.max_ram_gb > 0 && PLATFORM_RAM_OVERHEAD > 1.0) {
            std::cout << " (platform RAM overhead " << PLATFORM_RAM_OVERHEAD << "x)";
        }
        std::cout << "\n";

        auto const t_start = std::chrono::steady_clock::now();
        auto print_total = [&] {
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                          std::chrono::steady_clock::now() - t_start)
                          .count();
            std::cout << "[total construction time] " << (ms / 1000.0) << " s\n";
        };

        std::string const tmp_dir = resolve_tmp_dir();
        std::cout << "  tmp_dir = " << tmp_dir << "\n";

        auto writer = std::make_unique<bucket_writer>(tmp_dir, num_buckets, m_flush_bases,
                                                      m_spill_bytes);
        // When -g is set, arm a background RSS watcher with
        // hysteresis. Bucket-write must leave room for what comes
        // after: bucket-process adds ~1 GiB on top on multi-thousand-
        // genome inputs, stitch adds ~300 MiB. We reserve those by
        // budgeting bucket-write at 60% of -g (HIGH) with a
        // release point at 45% (LOW). Once live RSS crosses HIGH, the
        // watcher trips pressure and sweeps all compactors once;
        // ingest threads then spill every batch until live RSS falls
        // below LOW, at which point pressure clears and normal
        // m_bytes-based spilling resumes. The hysteresis avoids the
        // death-spiral of a sticky flag (every batch spills forever).
        //
        // If /proc/self/status isn't readable we fall back to the
        // monotonic getrusage peak; in that case there's no release,
        // pressure is sticky, and we do over-spill -- that's the
        // safer-but-slower path on platforms without VmRSS.
        if (m_cfg.max_ram_gb > 0) {
            uint64_t budget_bytes =
                (uint64_t)(m_cfg.max_ram_gb * 1024.0 * 1024.0 * 1024.0);
            uint64_t high_threshold_bytes = (uint64_t)(0.60 * (double)budget_bytes);
            uint64_t low_threshold_bytes = (uint64_t)(0.45 * (double)budget_bytes);
            writer->start_rss_watcher(high_threshold_bytes, low_threshold_bytes);
        }
        {
            phase_rss_marker rss("bucket-write");
            {
                timer _("bucket-write");
                std::atomic<uint64_t> done{0};
                progress prog("bucket-write", done, files.size());
                ingest_bucketed(files, m_cfg.k, m_cfg.m, m_cfg.bucket_log2, *writer,
                                m_cfg.num_threads, &done);
                prog.stop();
            }
            // Snapshot before writer.close() frees the compactor
            // hashmaps -- their memory is part of the bucket-write peak.
            rss.stop();
        }
        // Stop the RSS watcher before close(): we don't want a stray
        // try_spill firing during close()'s final spill+fclose loop.
        writer->stop_rss_watcher();
        if (writer->pressure_was_engaged()) {
            std::cout << "  bucket-write: RSS pressure engaged; observed live-RSS high "
                      << format_bytes(writer->observed_rss_high()) << "\n";
        }
        writer->close();
        std::cout << "  bucket bytes written: " << writer->total_bytes() << " (compressed; "
                  << writer->total_uncompressed_bytes() << " uncompressed)\n";
        // Per-bucket size distribution: tells us whether the bucket-process
        // working set is driven by a few outliers (resplit fixes it) or is
        // broad (need more buckets / concurrency cap). num_threads = how many
        // buckets are resident at once.
        writer->report_bucket_size_distribution(m_cfg.num_threads);
        bucket_prof().print(m_cfg.num_threads);

        // Bucket processing emits stitchable fragments AND merges
        // per-bucket color sets into the shared global dict on the
        // fly. By the time we exit this block, every fragment already
        // carries a global cid.
        //
        // Fragments are NOT accumulated in an in-RAM vector -- on
        // dense pangenome inputs that vector grew to multi-GB during
        // bucket-process and was the dominant peak contributor.
        // Instead each fragment is streamed through a disk-backed
        // frag_unitig_writer to a single tmp file. Stitch then reads
        // that file back with a streaming frag_unitig_stream_reader that
        // holds one fragment at a time, so no full
        // vector<stitchable_unitig> -- and no O(num_fragments) index --
        // is ever materialised.
        //
        // Streaming dict: encodes each new color set into its bvb at
        // intern() time and immediately flushes complete 64-bit words
        // to the final <basename>.color_sets file. Per-class memory
        // is just 16 bytes of metadata; the compressed bit_vector
        // never sits in RAM. EF offsets are appended to the file at
        // finalize().
        frag_unitig_writer frag_sink(tmp_dir + "/frag_unitigs.bin");
        streaming_color_set_dict global_dict(m_num_colors,
                                             m_cfg.out_basename + ".color_sets");
        std::mutex global_mu;
        {
            phase_rss_marker rss("bucket-process");
            {
                timer _("bucket-process");
                std::atomic<uint64_t> done{0};
                progress prog("bucket-process", done, num_buckets);
                // bucket-process holds one bucket's kmer_info resident per
                // in-flight bucket. Rather than cap the user's -t, we cap the
                // total kmer_info RAM via a memory-admission gate: workers
                // (all num_threads of them) wait to load a bucket until it fits
                // the budget. Under a tight -g this loads FEWER large buckets
                // at once (and more small ones) -- "load less in RAM at a time"
                // -- without ever reducing -t. Budget = share of -g remaining
                // after the live carry-in.
                uint64_t bp_budget = 0;
                if (m_cfg.max_ram_gb > 0) {
                    const uint64_t total =
                        (uint64_t)(m_cfg.max_ram_gb * 1024.0 * 1024.0 * 1024.0);
                    const uint64_t carry = current_rss_bytes();
                    const uint64_t avail = total > carry ? (total - carry) : total;
                    bp_budget = (uint64_t)((double)avail * 0.80);  // bucket-process share
                }
                process_buckets(*writer, m_cfg.k, m_num_colors, m_cfg.num_threads,
                                std::ref(frag_sink), global_dict, global_mu, &done,
                                bp_budget);
                prog.stop();
                std::cout << "  bucket fragments: " << frag_sink.count() << "\n";
                std::cout << "  distinct color classes: " << global_dict.size() << "\n";
                std::cout << "  global color dict resident: ~"
                          << format_bytes(global_dict.resident_bytes())
                          << " (stays in RAM through stitch + emit)\n";
            }
            rss.stop();
        }
        process_prof().print(m_cfg.num_threads);
        m_num_color_classes = global_dict.size();
        // Interning is done. Free the dict's dedup index + per-class hash
        // vector NOW (finalize only needs the class count + the on-disk
        // bits/offsets). At high class counts this was the dominant cross-
        // phase carry-in -- it stayed resident through stitch + emit and
        // pushed stitch over budget (100K -g16: ~15 GiB carried into stitch).
        global_dict.release_index();
        frag_sink.close_for_writing();
        // bucket_writer's per-bucket compactor state (m_dict_classes
        // and other reuse-friendly buffers) is alive at high-water
        // until the writer is destroyed -- ~300 MB on 25K scale that
        // would otherwise carry into stitch. Drop it now; the bucket
        // *files* on disk are still there (cleanup_tmp_dir removes
        // them at end of build), and process_buckets has already read
        // them.
        writer.reset();
        // glibc holds free'd allocations in per-thread arenas across
        // phase boundaries; on a 16-thread bucket-process this can
        // be 500 MB - 1 GB of "free but not returned to OS" memory
        // that still counts toward RSS during stitch. Force release
        // back to the kernel before stitch starts.
        release_free_heap_to_os_();

        // Stitch streams each finished unitig into a cid-range
        // unitig_bucket_writer. The K bucket count auto-scales so
        // per-bucket peak at emit stays under ~10% of -g. We size K
        // up-front from the frag writer's own seq-byte counter (no walk
        // over the spill).
        std::unique_ptr<unitig_bucket_writer> uwriter_ptr;
        {
            // The stitch reads the frag spill with a STREAMING reader that
            // holds one fragment at a time -- no O(num_fragments) index. The
            // total seq bytes (for bucket sizing) and fragment count (for the
            // progress bar) come from the writer's own counters, so we never
            // walk an index. emit only needs uwriter (cid-bucketed unitig
            // spill) and global_dict (streaming color sets).
            phase_rss_marker rss("stitch");
            {
                timer _("stitch");
                const uint64_t total_frag_seq_bytes = frag_sink.total_seq_bytes();
                const uint64_t n_frags = frag_sink.count();

                const uint32_t unitig_bucket_count = pick_unitig_bucket_count_(
                    m_num_color_classes, total_frag_seq_bytes, m_cfg.max_ram_gb);
                uwriter_ptr = std::make_unique<unitig_bucket_writer>(
                    tmp_dir, m_num_color_classes, unitig_bucket_count);

                frag_unitig_stream_reader frag_reader(frag_sink.path());
                std::atomic<uint64_t> done{0};
                progress prog("stitch", done, n_frags);
                // External-memory iterative-doubling stitch: per-round
                // tigs are bucketed to LZ4-framed files under tmp_dir, and the
                // parallel round loop holds num_threads buckets resident at
                // once. We DIMENSION THE BUCKET COUNT (never the user's thread
                // count) so num_threads resident buckets fit stitch's share of
                // -g: more, smaller buckets. The user's -t is honored as-is.
                const uint32_t stitch_buckets = pick_stitch_buckets_(
                    total_frag_seq_bytes, m_cfg.max_ram_gb, m_cfg.num_threads);
                std::cout << "  stitch buckets: " << stitch_buckets
                          << ", threads: " << m_cfg.num_threads << "\n";
                // Parallel over the per-round bucket loop: buckets are
                // independent within a round (own waiting-map; content-seeded
                // RNG), so this fans out cleanly. Output is unchanged --
                // emit_fasta sorts each cid-bucket, so .fa is identical
                // regardless of which thread emitted which unitig.
                stitch_unitigs_extmem_file_stream(frag_reader, m_cfg.k, tmp_dir,
                                                  std::ref(*uwriter_ptr), stitch_buckets, &done,
                                                  m_cfg.num_threads);
                prog.stop();
                std::cout << "  unitigs after stitching: " << uwriter_ptr->total_unitigs() << "\n";
                // frag_reader destroyed here -- file handle closed.
            }
            rss.stop();
        }
        // Spill file no longer needed; safe to unlink.
        frag_sink.unlink();
        m_num_unitigs = uwriter_ptr->total_unitigs();

        {
            phase_rss_marker rss("emit-fasta");
            emit_fasta(*uwriter_ptr);
            rss.stop();
        }
        {
            phase_rss_marker rss("emit-colors");
            emit_colors(global_dict);
            rss.stop();
        }

        cleanup_tmp_dir(tmp_dir);

        // Peak resident set size across the whole build, as tracked by
        // the kernel (VmHWM in /proc/self/status). One read at the end;
        // zero hot-path cost.
        m_peak_rss_bytes = process_peak_rss_bytes();
        if (m_peak_rss_bytes) {
            std::cout << "[peak resident memory] " << format_bytes(m_peak_rss_bytes);
            if (m_cfg.max_ram_gb > 0) {
                uint64_t budget = (uint64_t)(m_cfg.max_ram_gb * 1024.0 * 1024.0 * 1024.0);
                if (m_peak_rss_bytes <= budget) {
                    std::cout << "  (within budget of " << format_bytes(budget) << ")";
                } else {
                    std::cout << "  (OVER budget of " << format_bytes(budget) << " by "
                              << format_bytes(m_peak_rss_bytes - budget) << ")";
                }
            }
            std::cout << "\n";
        }

        std::cout << "done. wrote " << m_cfg.out_basename << ".fa, " << m_cfg.out_basename
                  << ".u2c, and " << m_cfg.out_basename << ".color_sets\n";
        print_total();
    }

    // Stats populated by build(); zero before build() runs.
    uint32_t num_colors() const { return m_num_colors; }
    uint64_t num_unitigs() const { return m_num_unitigs; }
    uint64_t num_color_classes() const { return m_num_color_classes; }
    uint64_t peak_rss_bytes() const { return m_peak_rss_bytes; }
    build_config const& config() const { return m_cfg; }

    // --- emit-phase microbenchmark entry point (src/emit_bench.cpp) ------
    // Runs ONLY the emit phase against a pre-populated unitig_bucket_writer
    // and streaming_color_set_dict, using the SAME emit_fasta / emit_colors
    // the real pipeline uses, so the RSS measured is faithful. Lets us
    // reproduce the emit peak at 661k proportions without the ~14 h
    // ingest/process/stitch upstream. `num_unitigs_total` must equal the
    // count pushed into `uwriter` (sizes the u2c bit_vector, as build() does).
    // `out_basename` and `max_ram_gb` are taken from the config passed at
    // construction.
    void run_emit_only(unitig_bucket_writer& uwriter, streaming_color_set_dict& dict,
                       uint64_t num_unitigs_total) {
        m_num_unitigs = num_unitigs_total;
        {
            phase_rss_marker rss("emit-fasta");
            emit_fasta(uwriter);
            rss.stop();
        }
        {
            phase_rss_marker rss("emit-colors");
            emit_colors(dict);
            rss.stop();
        }
        m_peak_rss_bytes = process_peak_rss_bytes();
        std::cout << "[emit-only peak resident memory] " << format_bytes(m_peak_rss_bytes);
        if (m_cfg.max_ram_gb > 0) {
            uint64_t budget = (uint64_t)(m_cfg.max_ram_gb * 1024.0 * 1024.0 * 1024.0);
            std::cout << (m_peak_rss_bytes <= budget ? "  (within budget of "
                                                     : "  (OVER budget of ")
                      << format_bytes(budget) << ")";
        }
        std::cout << "\n";
    }

    // Public forwarder so the emit-bench harness sizes its unitig bucket
    // count EXACTLY as build() does (else the emit-fasta read_bucket peak,
    // which scales with records-per-bucket = num_unitigs / K, is wrong).
    static uint32_t pick_unitig_bucket_count(uint64_t num_color_classes,
                                             uint64_t total_seq_bytes_estimate,
                                             double max_ram_gb) {
        return pick_unitig_bucket_count_(num_color_classes, total_seq_bytes_estimate, max_ram_gb);
    }

private:
    // Soft-cap policy. Reads m_cfg.max_ram_gb (0 = unset) and decides:
    //   - bucket_log2 (more buckets -> smaller per-bucket data structures)
    //   - color-bvb spill threshold (in bytes; 0 = never spill)
    // Tighter budgets push bucket_log2 toward 13 (8192 buckets, GGCAT's
    // upper end) and shrink the bvb spill threshold proportionally.
    // The user-facing CLI override (-b) takes precedence
    // over this auto-tune.
    static constexpr uint32_t MIN_BUCKETS_LOG2 = 10;
    static constexpr uint32_t MAX_BUCKETS_LOG2 = 13;

    void validate_and_resolve_config() {
        if (m_cfg.filenames_list.empty())
            throw std::runtime_error("build_config::filenames_list is empty");
        if (m_cfg.out_basename.empty())
            throw std::runtime_error("build_config::out_basename is empty");
        if (m_cfg.k == 0 || m_cfg.k > MAX_K)
            throw std::runtime_error("k must satisfy 1 <= k <= " + std::to_string(MAX_K));
        if (m_cfg.num_threads == 0) m_cfg.num_threads = 1;
        if (m_cfg.m == 0) m_cfg.m = compute_best_m(m_cfg.k);
        if (m_cfg.m < 2 || m_cfg.m > m_cfg.k)
            throw std::runtime_error("invalid m=" + std::to_string(m_cfg.m) +
                                     " (need 2 <= m <= k)");
        if (m_cfg.bucket_log2 == 0) {
            // User didn't pin -b. Auto-pick.
            m_cfg.bucket_log2 = auto_bucket_log2();
        } else if (m_cfg.bucket_log2 < MIN_BUCKETS_LOG2 || m_cfg.bucket_log2 > MAX_BUCKETS_LOG2) {
            throw std::runtime_error("-b must be in [" +
                                     std::to_string(MIN_BUCKETS_LOG2) + ", " +
                                     std::to_string(MAX_BUCKETS_LOG2) + "]");
        }
        // bucket_writer opens one FILE per bucket (LZ4-framed, raw
        // stdio). Raise the soft FD limit if needed; clamp bucket_log2
        // if even the hard limit isn't enough.
        m_cfg.bucket_log2 = ensure_fd_capacity_for_buckets_(m_cfg.bucket_log2);
        auto_tune_bucket_write_params();
    }

    // Try to raise RLIMIT_NOFILE so we can open `1 << log2` bucket
    // files plus a small headroom for stdin/stdout/sidecar/inputs. If
    // even the hard limit is too small, clamp log2 down and warn.
    static uint32_t ensure_fd_capacity_for_buckets_(uint32_t log2) {
        constexpr uint32_t HEADROOM = 64;
        struct rlimit r;
        if (::getrlimit(RLIMIT_NOFILE, &r) != 0) return log2;  // best effort
        uint32_t needed = (1u << log2) + HEADROOM;
        if (r.rlim_cur >= needed) return log2;
        rlim_t target = std::min<rlim_t>(needed, r.rlim_max);
        struct rlimit nr = r;
        nr.rlim_cur = target;
        ::setrlimit(RLIMIT_NOFILE, &nr);
        ::getrlimit(RLIMIT_NOFILE, &nr);
        if (nr.rlim_cur >= needed) return log2;
        // Hard limit too small. Clamp log2 to what we can actually open.
        uint32_t avail = (uint32_t)nr.rlim_cur > HEADROOM ? (uint32_t)nr.rlim_cur - HEADROOM : 0;
        uint32_t fitted = MIN_BUCKETS_LOG2;
        while (fitted < log2 && (1u << (fitted + 1)) <= avail) ++fitted;
        if (fitted < log2) {
            std::cerr << "warning: RLIMIT_NOFILE hard limit " << nr.rlim_max
                      << " can't accommodate 2^" << log2 << " buckets; clamping bucket_log2 to "
                      << fitted << ". Raise the hard limit (e.g. ulimit -Hn) for tighter budgets.\n";
        }
        return fitted;
    }

    // Per-platform RAM-overhead multiplier. The auto-tune below models
    // bucket-write peak using overhead constants calibrated on macOS,
    // where libsystem_malloc has heavier per-allocation bookkeeping and
    // RSS accounting includes pages glibc would reclaim. On Linux the
    // same workload typically sits ~50% lower in RSS for the same
    // logical state, so we'd over-tighten the auto-tune if we used the
    // macOS calibration there.
    //
    // We compensate by dividing the effective budget that the auto-tune
    // sees by this constant. A value of 2.0 on macOS means the auto-
    // tune behaves as if the user passed half their -g (smaller
    // flush_bases / spill_bytes / fewer buckets) so the resulting
    // structural footprint actually fits under the original cap once
    // platform overhead is added back. On Linux the constant is 1.0
    // (no derating) because the existing constants already match.
    //
    // The watcher's HIGH/LOW thresholds are NOT divided by this
    // multiplier -- the watcher measures real RSS, which is what we
    // want to cap.
    static constexpr double PLATFORM_RAM_OVERHEAD =
#if defined(__APPLE__)
        2.0;
#else
        1.0;
#endif

    // Fraction of -g the bucket-write auto-tune is allowed to
    // plan for (per-thread buffers + compactor + compressor). The
    // remainder absorbs the bucket-process working set, which on Mac
    // expands by ~1.75 GiB on the 4546-genome workload (k-mer rsids
    // map + per-bucket compact_color_set_dict + libsystem_malloc bookkeeping).
    // We hand bucket-write a smaller slice on Mac so the cap holds
    // once bucket-process expands on top.
    static constexpr double BUCKET_WRITE_SHARE =
#if defined(__APPLE__)
        0.40;
#else
        0.50;
#endif

    double effective_max_ram_gb() const {
        return m_cfg.max_ram_gb / PLATFORM_RAM_OVERHEAD;
    }

    // Pick num_buckets so that the auto-tune can give us spill_bytes
    // at least TARGET_SPILL_BYTES. Otherwise the bucket files explode
    // with redundant records:
    //
    //   spill_bytes ≈ (SHARE × effective_ram / 2) / (B × COMPACTOR_OVERHEAD)
    //
    // is the value the compactor-half of the bucket-write share lets
    // us afford. Each spill clears the per-bucket dedup state, so a
    // popular super-k-mer that's hit again after a spill is emitted
    // as a *separate* on-disk record. The smaller spill_bytes is, the
    // more often this happens, and on the 4546-genome workload
    // pushing spill_bytes to its 16 KiB floor produced 1.96 M spills
    // and 10.6 GB of uncompressed bucket bytes -- vs ~313 K spills
    // and ~3 GB at "reasonable" spill_bytes. The 3.5x bloat then
    // overwhelmed any savings from the new compression.
    //
    // Strategy: solve for the largest B in [MIN, MAX] such that
    // B × TARGET_SPILL_BYTES × COMPACTOR_OVERHEAD ≤ (SHARE × eff)/2.
    // Bigger budgets get more buckets (which helps bucket-process
    // parallelism); tighter budgets get fewer buckets (so spill_bytes
    // can stay at the target and each spill carries enough data that
    // dedup pays off).
    //
    // No -g set -> historical 1024 buckets.
    uint32_t auto_bucket_log2() const {
        if (m_cfg.max_ram_gb <= 0) return MIN_BUCKETS_LOG2;
        constexpr double SHARE = BUCKET_WRITE_SHARE;
        constexpr double COMPACTOR_OVERHEAD = 7.0;
        constexpr size_t TARGET_SPILL_BYTES = 64 * 1024;  // good for LZ4 dedup

        const double compactor_share_bytes =
            effective_max_ram_gb() * 1024.0 * 1024.0 * 1024.0 * (SHARE / 2.0);
        const double max_b =
            compactor_share_bytes / ((double)TARGET_SPILL_BYTES * COMPACTOR_OVERHEAD);

        uint32_t log2 = MAX_BUCKETS_LOG2;
        while (log2 > MIN_BUCKETS_LOG2 && (double)((uint64_t)1 << log2) > max_b) {
            --log2;
        }
        return log2;
    }

    // Joint auto-tune of (flush_bases, spill_bytes) against a fixed
    // share of -g, given (num_threads, num_buckets).
    //
    // Bucket-write peak comes from three terms that we model as:
    //
    //   per-thread buffers   ≈ num_threads × num_buckets × flush_bases × BUFFER_OVERHEAD
    //   compactor footprint  ≈ num_buckets × spill_bytes × COMPACTOR_OVERHEAD
    //   compressor state     ≈ num_buckets × COMPRESSOR_BYTES_PER_BUCKET
    //
    // The overhead factors are empirical: nominal storage substantially
    // undercounts what the process actually resident-pages because of
    //   (a) std::vector capacity slack (per-thread buffers retain their
    //       high-water capacity for the whole phase),
    //   (b) std::string heap allocations for compactor keys (each ~80 B
    //       header + bytes + malloc header beyond the 15 B SSO),
    //   (c) glibc/libsystem_malloc fragmentation across millions of
    //       small allocs churned by repeated spill cycles,
    //   (d) ankerl::unordered_dense map structure overhead per entry.
    //
    // The compressor term used to dominate (zlib level-1 at ~256-448
    // KiB per gzFile). With LZ4's block API (LZ4_compress_default is
    // stateless from the caller's view) there is no persistent
    // compressor state per bucket -- only the m_batch_buf and
    // m_out_buf scratch buffers, which are accounted for inside
    // COMPACTOR_OVERHEAD (raised to 7x). COMPRESSOR_BYTES_PER_BUCKET
    // here is just a small constant for unmodelled per-bucket
    // bookkeeping.
    //
    // Strategy: reserve a target share (default 50%) of -g for
    // bucket-write. Subtract the (small) compressor floor. Split what
    // remains evenly between per-thread buffers and compactor data,
    // then solve for flush_bases and spill_bytes.
    //
    // No -g set -> keep historical defaults (this keeps the
    // small-input dev path identical and avoids surprising regressions
    // for users who don't care about a budget).
    void auto_tune_bucket_write_params() {
        if (m_cfg.max_ram_gb <= 0) {
            m_flush_bases = 64 * 1024;
            m_spill_bytes = DEFAULT_COMPACTOR_SPILL_BYTES;  // 256 KiB
            return;
        }
        constexpr double SHARE = BUCKET_WRITE_SHARE;  // share of budget for bucket-write
        constexpr double COMPACTOR_OVERHEAD = 7.0;    // compactor: structure + key allocs + frag
        constexpr double BUFFER_OVERHEAD = 2.0;     // per-thread: capacity slack + recs vec
        // LZ4 block API has no persistent compressor state -- the
        // batch and out buffers we hold between spills are accounted
        // for by COMPACTOR_OVERHEAD (raised from 5x to 7x to cover
        // them). Keep a small constant here for unmodelled per-bucket
        // bookkeeping (FILE handle, mutex padding, etc.).
        constexpr size_t COMPRESSOR_BYTES_PER_BUCKET = 4 * 1024;
        constexpr size_t MIN_FLUSH_BASES = 4 * 1024;
        constexpr size_t MIN_SPILL_BYTES = 16 * 1024;
        constexpr size_t MAX_FLUSH_BASES = 64 * 1024;   // historical default
        constexpr size_t MAX_SPILL_BYTES = 256 * 1024;  // historical default

        double const effective_ram_gb = effective_max_ram_gb();
        const uint64_t budget_bytes =
            (uint64_t)(effective_ram_gb * 1024.0 * 1024.0 * 1024.0 * SHARE);
        const uint64_t B = 1ull << m_cfg.bucket_log2;
        const uint64_t T = std::max<uint64_t>(1, m_cfg.num_threads);

        const uint64_t compressor_total = B * (uint64_t)COMPRESSOR_BYTES_PER_BUCKET;
        uint64_t remainder = budget_bytes > compressor_total ? budget_bytes - compressor_total : 0;
        uint64_t half = remainder / 2;

        double flush_d =
            T * B == 0 ? (double)MAX_FLUSH_BASES : (double)half / (T * B * BUFFER_OVERHEAD);
        size_t flush = (size_t)flush_d;
        if (flush < MIN_FLUSH_BASES) flush = MIN_FLUSH_BASES;
        if (flush > MAX_FLUSH_BASES) flush = MAX_FLUSH_BASES;

        double spill_d =
            B == 0 ? (double)MAX_SPILL_BYTES : (double)half / (B * COMPACTOR_OVERHEAD);
        size_t spill = (size_t)spill_d;
        if (spill < MIN_SPILL_BYTES) spill = MIN_SPILL_BYTES;
        if (spill > MAX_SPILL_BYTES) spill = MAX_SPILL_BYTES;

        m_flush_bases = flush;
        m_spill_bytes = spill;
    }

    static std::vector<std::string> read_filenames(std::string const& path) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("cannot open filenames list: " + path);
        std::vector<std::string> v;
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty()) v.push_back(line);
        }
        return v;
    }

    // If m_cfg.tmp_dir is empty, mkdtemp under $TMPDIR. Otherwise use it
    // (creating it if missing); if it already exists it must be an empty
    // directory, since cleanup_tmp_dir wipes the whole thing.
    std::string resolve_tmp_dir() {
        if (m_cfg.tmp_dir.empty()) {
            std::string tmpl =
                (std::filesystem::temp_directory_path() / "cdbg_buckets_XXXXXX").string();
            std::vector<char> buf(tmpl.begin(), tmpl.end());
            buf.push_back('\0');
            if (mkdtemp(buf.data()) == nullptr)
                throw std::runtime_error(std::string("mkdtemp failed: ") + std::strerror(errno));
            return std::string(buf.data());
        }
        std::error_code ec;
        if (std::filesystem::exists(m_cfg.tmp_dir, ec)) {
            if (!std::filesystem::is_directory(m_cfg.tmp_dir, ec))
                throw std::runtime_error("-d " + m_cfg.tmp_dir +
                                         " exists but is not a directory");
            if (!std::filesystem::is_empty(m_cfg.tmp_dir, ec))
                throw std::runtime_error("-d " + m_cfg.tmp_dir +
                                         " is not empty (the tool will remove the directory on"
                                         " exit, so it must start empty)");
        } else {
            std::filesystem::create_directories(m_cfg.tmp_dir, ec);
            if (ec)
                throw std::runtime_error("cannot create -d " + m_cfg.tmp_dir + ": " +
                                         ec.message());
        }
        return m_cfg.tmp_dir;
    }

    static void cleanup_tmp_dir(std::string const& tmp_dir) {
        timer _("removing tmp files");
        std::error_code ec;
        std::filesystem::remove_all(tmp_dir, ec);
    }

    // Pick the cid-range bucket count for the unitig spill, so the
    // peak in-RAM seq footprint at emit time -- one bucket loaded
    // and sorted -- stays under a target fraction of -g.
    //
    // Each bucket's records: roughly total_seq_bytes / K plus
    // ~2x overhead (std::string capacity slack + per-record book-
    // keeping in the sort vector). We size K so that target
    // per-bucket footprint is <= 10% of -g, leaving the
    // remaining 90% for whatever else is resident at emit time
    // (the streaming color-set dict's metadata, the FILE buffer,
    // the u2c bit_vector, etc.). Floor at MIN_K so we don't end
    // up with one giant bucket when -g is unset or huge,
    // ceiling at MAX_K to stay within RLIMIT_NOFILE headroom
    // (bucket-write's FDs are already closed when stitch starts,
    // so we have ~1024 FDs available).
    // Best-effort release of free'd glibc-arena memory back to the
    // kernel. No-op on non-glibc allocators (musl, jemalloc). Worth
    // calling at phase boundaries because on 16-thread workloads the
    // per-thread arena freelists can hold 1 GB+ that still counts in
    // RSS even after all containers were destroyed.
    static void release_free_heap_to_os_() {
#if defined(__GLIBC__)
        ::malloc_trim(0);
#endif
    }

    static uint32_t pick_unitig_bucket_count_(uint64_t num_color_classes,
                                              uint64_t total_seq_bytes_estimate,
                                              double max_ram_gb) {
        if (num_color_classes == 0) return 1;
        constexpr uint32_t MIN_K = 16;
        constexpr uint32_t MAX_K = 1024;
        constexpr uint32_t DEFAULT_K = 64;
        constexpr double OVERHEAD = 2.0;
        constexpr double SHARE = 0.10;

        if (max_ram_gb <= 0 || total_seq_bytes_estimate == 0) {
            return (uint32_t)std::min<uint64_t>(num_color_classes, DEFAULT_K);
        }
        const uint64_t budget_bytes =
            (uint64_t)(max_ram_gb * 1024.0 * 1024.0 * 1024.0 * SHARE);
        if (budget_bytes == 0) {
            return (uint32_t)std::min<uint64_t>(num_color_classes, DEFAULT_K);
        }
        const uint64_t needed = (uint64_t)((double)total_seq_bytes_estimate * OVERHEAD);
        uint64_t k = (needed + budget_bytes - 1) / budget_bytes;
        if (k < MIN_K) k = MIN_K;
        if (k > MAX_K) k = MAX_K;
        if (k > num_color_classes) k = num_color_classes;
        return (uint32_t)k;
    }

    // Pick the hash fan-out for the external-memory stitch. Peak RAM is
    // ~one round-bucket's tigs at a time, so we want each bucket to hold
    // at most a share of -g worth of (sequence + per-tig overhead)
    // bytes. More buckets -> smaller per-bucket resident set but more
    // round files; fewer buckets -> larger resident set. With no -g a
    // fixed default keeps small inputs in a handful of rounds.
    // Pick the stitch round-store bucket count AND the effective stitch
    // concurrency together, so peak RSS stays within budget.
    //
    // The parallel stitch holds ONE bucket resident per in-flight thread, so
    // peak ~= effective_threads * (one bucket's tigs). To keep that within
    // stitch's SHARE of -g we size buckets for `num_threads` resident copies:
    // count >= needed * num_threads / budget. If that hits the MAX_BUCKETS
    // (file-count) ceiling, we instead CAP the thread count so
    // effective_threads * one_bucket still fits. Returns {buckets, threads}.
    // Dimension the stitch round-store bucket count so `num_threads` buckets
    // (the parallel round loop holds one decoded bucket per in-flight thread)
    // fit stitch's share of -g. We size the BUCKET COUNT, never the user's
    // thread count -- more, smaller buckets keep full -t concurrency.
    static uint32_t pick_stitch_buckets_(uint64_t total_frag_seq_bytes, double max_ram_gb,
                                         uint32_t num_threads) {
        constexpr uint32_t MIN_BUCKETS = 64;
        constexpr uint32_t MAX_BUCKETS = 1u << 20;  // 1M files cap (fd headroom)
        constexpr uint32_t DEFAULT_BUCKETS = 1024;
        // Per-resident-bucket RAM is OVERHEAD x its on-disk seq bytes. Measured
        // on 100K -g16: a 3.0x model bounded stitch to 8 GiB but it actually
        // spiked to ~17 GiB across 48 resident round-0 buckets -> real ~6.4x.
        // The extra is ext_tig's per-fragment FIXED overhead
        // (std::vector<color_run> + std::string headers, ~120 B/tig) that
        // dominates for short fragments. 7.0x with margin (validated: 100K
        // -g16 stitch peak +1.16 GiB, whole run 15.77 < 16).
        constexpr double OVERHEAD = 7.0;
        constexpr double SHARE = 0.50;  // stitch's share of -g, across ALL resident buckets
        if (num_threads == 0) num_threads = 1;
        if (max_ram_gb <= 0 || total_frag_seq_bytes == 0) return DEFAULT_BUCKETS;

        // Share of what REMAINS after the live carry-in (the color dict is
        // freed before stitch via release_index, so this is small now).
        const uint64_t carry_in = current_rss_bytes();
        const uint64_t budget_total = (uint64_t)(max_ram_gb * 1024.0 * 1024.0 * 1024.0);
        const uint64_t avail =
            budget_total > carry_in ? (budget_total - carry_in) : budget_total;
        const uint64_t budget = (uint64_t)((double)avail * SHARE);
        if (budget == 0) return DEFAULT_BUCKETS;

        const uint64_t needed = (uint64_t)((double)total_frag_seq_bytes * OVERHEAD);
        // count s.t. needed/count * num_threads <= budget (num_threads resident).
        uint64_t count = ((needed * num_threads) + budget - 1) / budget;
        if (count < MIN_BUCKETS) count = MIN_BUCKETS;
        if (count > MAX_BUCKETS) count = MAX_BUCKETS;  // hard fd ceiling; see note below
        return (uint32_t)count;
    }

    // FASTA emit. Hand-rolled 1 MiB buffer + std::to_chars for the
    // integer header + memcpy for the sequence body. Significantly
    // faster than std::ofstream's default 8 KiB buffer + stream
    // operators on millions of small records.
    //
    // Reads unitigs from the disk-backed bucket spill in bucket order
    // (bucket 0 = lowest cids, ..., bucket K-1 = highest). Within a
    // bucket records arrive in stitch order, so we sort by cid to
    // make the .fa output strictly cid-ascending. Per-bucket peak
    // in RAM is one bucket's seqs (~few MB on the 4546-genome
    // workload) plus the per-record vector.
    //
    // While we already have the unitigs in their final cid-ascending
    // emission order, also build the unitig-to-color-set "u2c"
    // bit_vector and serialize it to <basename>.u2c. Bit i is set
    // iff unitig i (in .fa emission order) is the last unitig of a
    // color-set run. Length = num_unitigs, popcount =
    // num_color_classes. Downstream consumers (Fulgor) recover the
    // per-unitig color-set id via rank1(unitig_id) using a rank9
    // index built at load time. Matches the bit_vector layout in
    // https://github.com/jermp/fulgor/blob/main/include/index.hpp .
    void emit_fasta(unitig_bucket_writer& uwriter) const {
        timer _("emit fasta");

        FILE* fa = std::fopen((m_cfg.out_basename + ".fa").c_str(), "wb");
        if (!fa) throw std::runtime_error("cannot open " + m_cfg.out_basename + ".fa");
        constexpr size_t BUF_BYTES = 1 << 20;
        std::vector<char> buf(BUF_BYTES);
        size_t pos = 0;
        auto flush_buf = [&] {
            if (pos == 0) return;
            if (std::fwrite(buf.data(), 1, pos, fa) != pos) {
                std::fclose(fa);
                throw std::runtime_error("short write to " + m_cfg.out_basename + ".fa");
            }
            pos = 0;
        };
        auto reserve = [&](size_t n) {
            if (pos + n > BUF_BYTES) flush_buf();
        };

        bits::bit_vector::builder u2c_bvb((uint64_t)m_num_unitigs, /*init=*/false);
        size_t emitted = 0;
        uint64_t prev_cid = 0;

        // Per-bucket RAM cap for the cid-sort. read_bucket_sorted sorts in RAM
        // when a bucket fits this, else external merge-sorts -- so the
        // emit-fasta peak is bounded by this cap REGARDLESS of cid skew (which
        // otherwise made one low-cid bucket dominate; 661k: +36 GiB). Use a
        // modest share of -g; fall back to a fixed cap with no -g.
        const uint64_t emit_mem_cap =
            m_cfg.max_ram_gb > 0
                ? (uint64_t)(m_cfg.max_ram_gb * 1024.0 * 1024.0 * 1024.0 * 0.10)
                : (uint64_t)(1ull << 30);  // 1 GiB default

        auto emit_record = [&](uint64_t cid, std::string_view seq) {
            // Mark the final unitig of the previous run. cid is globally
            // non-decreasing (buckets ascending by cid range, sorted within),
            // so cid != prev_cid is exactly a color-set group boundary.
            if (emitted > 0 && cid != prev_cid) u2c_bvb.set(emitted - 1, 1);
            prev_cid = cid;
            ++emitted;

            reserve(22);  // '>' + up to 20 digits (uint64_t) + '\n'
            buf[pos++] = '>';
            auto rr = std::to_chars(buf.data() + pos, buf.data() + pos + 20, cid);
            pos = (size_t)(rr.ptr - buf.data());
            buf[pos++] = '\n';
            size_t s_pos = 0;
            while (s_pos < seq.size()) {
                if (pos == BUF_BYTES) flush_buf();
                size_t take = std::min(BUF_BYTES - pos, seq.size() - s_pos);
                std::memcpy(buf.data() + pos, seq.data() + s_pos, take);
                pos += take;
                s_pos += take;
            }
            reserve(1);
            buf[pos++] = '\n';
        };

        for (uint32_t b = 0; b < uwriter.num_buckets(); ++b) {
            uwriter.read_bucket_sorted(b, emit_mem_cap, emit_record);
        }
        flush_buf();
        std::fclose(fa);

        std::cout << "  [emit-fasta] " << uwriter.num_buckets()
                  << " buckets, cid-sort RAM cap " << format_bytes(emit_mem_cap)
                  << " (external merge-sort if a bucket exceeds it)\n";

        // Close out the very last run.
        if (emitted > 0) u2c_bvb.set(emitted - 1, 1);
        bits::bit_vector u2c;
        u2c_bvb.build(u2c);
        essentials::save(u2c, (m_cfg.out_basename + ".u2c").c_str());

        uwriter.close_and_unlink();
    }

    void emit_colors(streaming_color_set_dict& global_dict) const {
        timer _("emit color_sets");
        // Encoding already happened during bucket-process via
        // global_dict.intern(); each intern flushed complete 64-bit
        // words to the final file. finalize() flushes the trailing
        // partial word, builds + appends the elias_fano over per-
        // class bit-offsets, then fseeks back to write the now-known
        // header totals.
        global_dict.finalize();
    }

    build_config m_cfg;
    uint32_t m_num_colors = 0;
    uint64_t m_num_unitigs = 0;
    uint64_t m_num_color_classes = 0;
    uint64_t m_peak_rss_bytes = 0;
    // bucket-write knobs picked by auto_tune_bucket_write_params().
    size_t m_flush_bases = 0;
    size_t m_spill_bytes = 0;
};

}  // namespace cdbg
