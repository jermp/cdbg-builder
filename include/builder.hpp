#pragma once

// Public API for building a colored compacted dBG.
//
// Wraps the full pipeline so
// that downstream tools can construct a `build_config`, instantiate a
// `builder`, and call `build()`.
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

#include "phase1_bucket_write/bucket_io.hpp"
#include "phase1_bucket_write/bucket_ingester.hpp"
#include "phase2_bucket_process/bucket_walker.hpp"
#include "phase1_bucket_write/minimizer.hpp"
#include "phase3_stitch/stitch.hpp"
#include "phase3_stitch/stitch_extmem.hpp"
#include "phase2_bucket_process/streaming_color_set_dict.hpp"
#include "phase2_bucket_process/unitig_spill.hpp"
#include "phase4_emit/emit.hpp"
#include "ram_governor.hpp"
#include "util.hpp"

namespace cdbg {

struct builder {
    // `explicit` to avoid an implicit build_config -> builder conversion.
    explicit builder(build_config const& cfg) : m_cfg(cfg) {}

    // Run the full pipeline. Throws std::runtime_error on bad config or
    // I/O error. On success, writes m_cfg.out_basename + {".fa", ".color_sets"}
    // and removes the scratch directory.
    void build() {
        validate_and_resolve_config();

        // RAM governor: the measured-RSS backstop that reclaims glibc-retained
        // pages under pressure and spills registered participants (the unitig
        // writer during stitch -- the buffer the per-phase watcher couldn't
        // reach). It is constructed here but STARTED only at stitch (below):
        // bucket-write and bucket-process have their own RAM controls and the
        // governor never fires there (RSS stays well under its high-watermark),
        // so polling during them is pure overhead for no benefit -- it must not
        // slow the fast phases. budget 0 => disabled. Trip a touch below -g
        // (0.80/0.62) so the spill+trim lag on a spinning disk has room before
        // the hard cap. See algorithm.md §5.5.
        ram_governor governor(m_cfg.max_ram_gb > 0
                                  ? (uint64_t)(m_cfg.max_ram_gb * 1024.0 * 1024.0 * 1024.0)
                                  : 0,
                              /*high_frac=*/0.80, /*low_frac=*/0.62);

        // Read the filenames list: one input path per line, blanks skipped.
        // File at line i is color i.
        std::vector<std::string> files;
        {
            std::ifstream in(m_cfg.filenames_list);
            if (!in) {
                throw std::runtime_error("cannot open filenames list: " + m_cfg.filenames_list);
            }
            std::string line;
            while (std::getline(in, line)) {
                if (!line.empty()) files.push_back(line);
            }
        }
        if (files.empty()) throw std::runtime_error("no input files");
        if (files.size() > (uint64_t)UINT32_MAX) {
            throw std::runtime_error("too many colors (max 2^32 - 1)");
        }
        m_num_colors = files.size();

        uint32_t const num_buckets = m_cfg.num_buckets;
        std::cout << "k = " << m_cfg.k << ", m = " << m_cfg.m << ", num_colors = " << m_num_colors
                  << ", num_threads = " << m_cfg.num_threads << ", num_buckets = " << num_buckets
                  << "\n";
        std::cout << "  bucket-write model: flush_bases=" << format_bytes(m_flush_bases)
                  << ", spill_bytes=" << format_bytes(m_spill_bytes) << ", alpha=" << m_cfg.alpha
                  << ", beta=" << m_cfg.beta;
        if (m_cfg.max_ram_gb > 0 and PLATFORM_RAM_OVERHEAD > 1.0) {
            std::cout << " (platform RAM overhead " << PLATFORM_RAM_OVERHEAD << "x)";
        }
        std::cout << "\n";

        seconds_timer build_timer;
        build_timer.start();

        std::string const tmp_dir = resolve_tmp_dir();
        std::cout << "  tmp_dir = " << tmp_dir << "\n";

        auto writer =
            std::make_unique<bucket_writer>(tmp_dir, num_buckets, m_flush_bases, m_spill_bytes);

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
            uint64_t budget_bytes = m_cfg.max_ram_gb * 1024.0 * 1024.0 * 1024.0;
            uint64_t high_threshold_bytes = 0.60 * budget_bytes;
            uint64_t low_threshold_bytes = 0.45 * budget_bytes;
            writer->start_rss_watcher(high_threshold_bytes, low_threshold_bytes);
        }

        {
            phase_rss_marker rss("bucket-write");
            {
                timer _("bucket-write");
                std::atomic<uint64_t> done{0};

                std::vector<uint64_t> file_sizes(files.size(), 0);
                uint64_t total_bytes = 0;
                for (uint64_t i = 0; i < files.size(); ++i) {
                    std::error_code ec;
                    auto sz = std::filesystem::file_size(files[i], ec);
                    file_sizes[i] = ec ? 0 : (uint64_t)sz;
                    total_bytes += file_sizes[i];
                }
                std::atomic<uint64_t> done_bytes{0};
                std::unique_ptr<progress> prog =
                    total_bytes > 0
                        ? std::make_unique<progress>("bucket-write", done_bytes, total_bytes,
                                                     /*render_bytes=*/true, &done,
                                                     (uint64_t)files.size(), "files")
                        : std::make_unique<progress>("bucket-write", done, (uint64_t)files.size());
                ingest_bucketed(files, m_cfg.k, m_cfg.m, m_cfg.num_buckets, *writer,
                                m_cfg.num_threads, &done, &file_sizes, &done_bytes);
                prog->stop();
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
        // Per-bucket size distribution: with the bucket-process admission gate,
        // this previews how many buckets fit in -g at once (the effective
        // concurrency) and whether a few large outliers dominate the working
        // set. num_threads = the max buckets that can be resident at once.
        writer->report_bucket_size_distribution(m_cfg.num_threads);
        bucket_prof().print(m_cfg.num_threads);

        // Free the per-bucket compactors (dedup hashmaps + LZ4 scratch buffers)
        // now that bucket files are written -- they are not needed by
        // bucket-process (which reads the files) and otherwise stay alive until
        // writer.reset() AFTER the phase, pinning glibc's churned arenas.
        writer->release_compactors();
        // Return bucket-write's freed memory to the OS BEFORE bucket-process.
        // The per-thread buffer pool and the just-freed compactors are freed but
        // glibc holds them in its arenas, so VmRSS still reflects the bucket-write
        // peak here. bucket-process's admission gate reads current_rss() as its
        // baseline carry; if that carry still includes the freed-but-unreturned
        // bucket-write memory, avail collapses and the gate admits ~one bucket at
        // a time (serial!). Trimming now drops carry to the true live baseline
        // so the gate admits the full set of buckets and runs -t-parallel.
        release_free_heap_to_os_();

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
        // holds one fragment at a time.
        //
        // Streaming dict: encodes each new color set into its bvb at
        // intern() time and immediately flushes complete 64-bit words
        // to the final <basename>.color_sets file. Per-class memory
        // is just 16 bytes of metadata; the compressed bit_vector
        // never sits in RAM. EF offsets are appended to the file at
        // finalize().
        frag_unitig_writer frag_sink(tmp_dir + "/frag_unitigs.bin");
        streaming_color_set_dict global_dict(m_num_colors, m_cfg.out_basename + ".color_sets");
        std::mutex global_mu;
        {
            phase_rss_marker rss("bucket-process");
            {
                timer _("bucket-process");
                std::atomic<uint64_t> done{0};
                progress prog("bucket-process", done, num_buckets);
                // bucket-process holds one bucket's kmer_info resident per
                // in-flight bucket, PLUS the global color dict, which grows
                // monotonically through the phase. The admission gate keeps the
                // whole phase under a total-RSS target: it admits a bucket only
                // while  carry + reserved_kmer_info + live_dict <= target,
                // re-reading the live dict each time. So as the dict grows,
                // fewer buckets load at once -- the dict is reserved for
                // dynamically, not left to overflow on top of a static
                // kmer_info budget.
                // All num_threads threads stay alive; -t is never reduced.
                uint64_t bp_budget = 0;
                uint64_t dedup_budget = 0;
                if (m_cfg.max_ram_gb > 0) {
                    const uint64_t total = m_cfg.max_ram_gb * 1024.0 * 1024.0 * 1024.0;
                    bp_budget = BUCKET_PROCESS_BUDGET_FRAC * total;
                    // Share of -g the color-sets-dedup-map may use before its
                    // overflow path engages. MEASURE-ONLY today (reported below,
                    // not yet enforced) -- see colorset-dedup-externalization.md.
                    dedup_budget = COLORSET_DEDUP_BUDGET_FRAC * total;
                }
                global_dict.set_dedup_budget(dedup_budget);
                process_buckets(*writer, m_cfg.k, m_num_colors, m_cfg.num_threads, frag_sink,
                                global_dict, global_mu, &done, bp_budget,
                                /*delete_consumed_buckets=*/!m_cfg.keep_tmp);
                prog.stop();
                std::cout << "  bucket fragments: " << frag_sink.count() << "\n";
                std::cout << "  distinct color classes: " << global_dict.size() << "\n";
                // color-sets-dedup-map occupancy vs its -g budget. Peak == final,
                // since the index only grows until release_index(). MEASURE-ONLY:
                // shows whether/by how much the externalization's overflow path
                // would have engaged, before any spill machinery exists.
                {
                    const uint64_t resident = global_dict.resident_bytes();
                    std::cout << "  color-sets-dedup-map: ~" << format_bytes(resident)
                              << " resident";
                    if (dedup_budget > 0) {
                        std::cout << " / " << format_bytes(dedup_budget) << " budget";
                        if (resident <= dedup_budget) {
                            std::cout << " (within -> no spill)";
                        } else {
                            std::cout << " (OVER by " << format_bytes(resident - dedup_budget)
                                      << " -> overflow path would spill)";
                        }
                    }
                    std::cout << "\n";
                }
            }
            rss.stop();
        }
        process_prof().print(m_cfg.num_threads);
        m_num_color_classes = global_dict.size();

        // Interning is done. Free the dict's dedup index + per-class hash
        // vector NOW (finalize only needs the class count + the on-disk
        // bits/offsets).
        global_dict.release_index();
        frag_sink.close_for_writing();
        writer.reset();
        release_free_heap_to_os_();

        // Start the governor now: stitch (and emit) is where it earns its keep
        // (the unitig writer spill + glibc reclaim near -g). Bucket-write/process
        // ran without it. It stops when `governor` goes out of scope at build end.
        governor.start();

        // Stitch streams each finished unitig into a cid-range
        // unitig_bucket_writer. The K bucket count auto-scales so
        // per-bucket peak at emit stays under ~10% of -g. We size K
        // up-front from the frag writer's own seq-byte counter (no walk
        // over the spill).
        std::unique_ptr<unitig_bucket_writer> uwriter_ptr;
        {
            // The stitch reads the frag spill with a STREAMING reader that
            // holds one fragment at a time. The
            // total seq bytes (for bucket sizing) and fragment count (for the
            // progress bar) come from the writer's own counters, so we never
            // walk an index. emit only needs uwriter (cid-bucketed unitig
            // spill) and global_dict (streaming color sets).
            phase_rss_marker rss("stitch");
            {
                timer _("stitch");
                const uint64_t total_frag_seq_bytes = frag_sink.total_seq_bytes();
                const uint64_t n_frags = frag_sink.count();

                const uint64_t unitig_bucket_count = pick_unitig_bucket_count_(
                    m_num_color_classes, total_frag_seq_bytes, m_cfg.max_ram_gb);

                // The write-through unitig_bucket_writer only needs to know
                // "finite (=> write-through, under -g)" vs "SIZE_MAX (=> keep in
                // RAM, no -g)". No -g => SIZE_MAX (keep every unitig in RAM, skip
                // the temp round-trip); under -g => a finite budget so the writer
                // streams every unitig straight to its bucket file.
                size_t unitig_ram_budget = SIZE_MAX;
                if (m_cfg.max_ram_gb > 0)
                    unitig_ram_budget = (size_t)(m_cfg.max_ram_gb * 1024.0 * 1024.0 * 1024.0 *
                                                 STITCH_BUDGET_FRAC);
                uwriter_ptr = std::make_unique<unitig_bucket_writer>(
                    tmp_dir, m_num_color_classes, unitig_bucket_count, unitig_ram_budget);
                // Govern the unitig writer for the stitch's duration: under -g
                // pressure the governor spills its largest in-RAM cid-buckets to
                // disk. Unregistered when this block exits (before emit, which
                // reads it single-threaded and so runs lock-free).
                auto uwriter_reg = governor.add(uwriter_ptr.get());

                std::atomic<uint64_t> done{0};
                progress prog("stitch", done, n_frags);
                // External-memory iterative-doubling stitch: per-round
                // tigs are bucketed to LZ4-framed files under tmp_dir, and the
                // parallel round loop holds num_threads buckets resident at
                // once. We DIMENSION THE BUCKET COUNT (never the user's thread
                // count) so num_threads resident buckets fit stitch's share of
                // -g: more, smaller buckets. The user's -t is honored as-is.
                const uint64_t stitch_buckets =
                    pick_stitch_buckets_(total_frag_seq_bytes, m_cfg.max_ram_gb, m_cfg.num_threads);
                std::cout << "  stitch buckets: " << stitch_buckets
                          << ", threads: " << m_cfg.num_threads << "\n";
                // Parallel over the per-round bucket loop: buckets are
                // independent within a round (own waiting-map; content-seeded
                // RNG), so this fans out cleanly. Output is unchanged --
                // emit_fasta sorts each cid-bucket, so .fa is identical
                // regardless of which thread emitted which unitig.
                // Yield each fragment's RAW record bytes via a block-buffered
                // reader (one fread per block, no per-field fread); the seed
                // parses them on worker threads. `done` (the progress bar) is
                // ticked once per record here on the reader thread.
                auto for_each_raw = [&](auto&& fn) {
                    frag_unitig_block_reader rd(frag_sink.path());
                    uint8_t const* rec;
                    uint32_t len;
                    while (rd.next_raw(rec, len)) {
                        fn(rec, len);
                        done.fetch_add(1, std::memory_order_relaxed);
                    }
                };
                // Base-carrying: assemble unitig sequences IN the doubling rounds
                // from a single frag-spill read, through an always-on-disk round
                // store (one bucket per thread resident). The proven README path:
                // 661k in ~3.5h at ~41 GiB, within -g. -g is held by keeping little
                // in RAM + the bucket-count sizing; the governor still backstops the
                // unitig writer. `done` is driven by for_each_raw above, so the
                // stitch takes nullptr.
                stitch_unitigs_extmem_file_streaming(for_each_raw, frag_record_parse, m_cfg.k,
                                                     tmp_dir, std::ref(*uwriter_ptr), stitch_buckets,
                                                     /*done=*/nullptr, m_cfg.num_threads);
                prog.stop();
                std::cout << "  unitigs after stitching: " << uwriter_ptr->total_unitigs() << "\n";
            }
            rss.stop();
        }
        // Spill file no longer needed; safe to unlink (kept under --keep-tmp for
        // debugging).
        if (!m_cfg.keep_tmp) frag_sink.unlink();
        m_num_unitigs = uwriter_ptr->total_unitigs();

        // Stitch allocates large transient buffers that are freed when it
        // returns, but glibc keeps the freed pages in its per-thread arenas
        // (RSS stays high), starving emit's bucket I/O of page cache. Hand them
        // back to the OS before emit.
        release_free_heap_to_os_();

        // emit streams the u2c bit_vector straight to <out>.u2c (its run-end
        // bits are produced in ascending unitig order), so it never holds the
        // num_unitigs-bit bitmap in RAM -- no -g floor here.
        {
            phase_rss_marker rss("emit-fasta");
            emit_fasta(*uwriter_ptr, m_cfg.out_basename, m_num_unitigs, m_cfg.max_ram_gb);
            rss.stop();
        }
        {
            phase_rss_marker rss("emit-colors");
            emit_colors(global_dict);
            rss.stop();
        }

        // Remove the scratch dir and everything under it (unless --keep-tmp,
        // which preserves the scratch dir for debugging).
        if (!m_cfg.keep_tmp) {
            timer _("removing tmp files");
            std::error_code ec;
            std::filesystem::remove_all(tmp_dir, ec);
        } else {
            std::cout << "  --keep-tmp: leaving scratch dir " << tmp_dir << " in place\n";
        }

        // Peak resident set size across the whole build, as tracked by
        // the kernel (VmHWM in /proc/self/status). One read at the end;
        // zero hot-path cost.
        m_peak_rss_bytes = process_peak_rss_bytes();
        if (m_peak_rss_bytes) {
            std::cout << "[peak resident memory] " << format_bytes(m_peak_rss_bytes);
            if (m_cfg.max_ram_gb > 0) {
                uint64_t budget = m_cfg.max_ram_gb * 1024.0 * 1024.0 * 1024.0;
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
        build_timer.stop();
        std::cout << "[total construction time] " << build_timer.elapsed() << " s\n";
    }

    // Stats populated by build(); zero before build() runs.
    uint64_t num_colors() const { return m_num_colors; }
    uint64_t num_unitigs() const { return m_num_unitigs; }
    uint64_t num_color_classes() const { return m_num_color_classes; }
    uint64_t peak_rss_bytes() const { return m_peak_rss_bytes; }
    build_config const& config() const { return m_cfg; }

private:
    // ---- RAM budget (-g) enforcement: honored AT ALL COSTS ------------------
    // -g is a HARD limit. Most structures spill to disk to stay under it, and a
    // few that don't spill are instead SCALED to fit: the stitch's per-range
    // member array shrinks by raising frag_ranges, and emit's u2c bit_vector is
    // streamed straight to disk. The genuinely NON-SPILLABLE floors that remain
    // (bucket-write per-thread buffers, a single bucket's kmer_info) ABORT with a
    // precise message when their budget share is exceeded rather than silently
    // blow -g. And -g below MIN_RAM_GB is rejected outright, since below that the
    // fixed per-thread/working floors cannot be honored.
    static constexpr double MIN_RAM_GB = 4.0;

    // Uniform "cannot honor -g" abort. `need`/`avail` are bytes.
    [[noreturn]] static void fail_ram_budget_(std::string const& what, uint64_t need,
                                              uint64_t avail, std::string const& remedy) {
        throw std::runtime_error("-g cannot be honored: " + what + " needs " + format_bytes(need) +
                                 " of NON-SPILLABLE RAM but only " + format_bytes(avail) +
                                 " is available for it. " + remedy +
                                 " (-g is a hard limit; aborting rather than exceeding it.)");
    }


    //   - bucket_log2 (more buckets -> smaller per-bucket data structures)
    //   - color-bvb spill threshold (in bytes; 0 = never spill)
    // Tighter budgets push bucket_log2 toward 16 and shrink the bvb spill threshold proportionally.
    // The user-facing CLI override (-b) takes precedence
    // over this auto-tune.
    // Bounds for the -b power-of-two override only. The auto path computes an
    // arbitrary (non-power-of-two) count from the RAM model and is bounded by
    // MIN_AUTO_BUCKETS below and the OS fd limit above.
    static constexpr uint32_t MIN_BUCKETS_LOG2 = 10;
    static constexpr uint32_t MAX_BUCKETS_LOG2 = 16;

    // Floor for the auto-derived bucket count (avoid degenerate single-bucket
    // runs at very tight -g; the model can still land here when g is tiny).
    static constexpr uint32_t MIN_AUTO_BUCKETS = 64;

    // Fraction of -g the bucket-write phase is sized against. 0.50 matches the
    // long-standing BUCKET_WRITE_SHARE that kept bucket-write lean. Phase 1
    // doesn't need most of -g -- only enough that B = M/(alpha*T*flush+beta*
    // spill) is large enough for a T-way bucket-process.
    static constexpr double BUCKET_WRITE_BUDGET_FRAC = 0.50;

    // Good default per-bucket batching payloads. B is derived as B = frac*g /
    // (alpha*T*flush + beta*spill), so bucket-write RAM stays at frac*g
    // regardless of these -- they trade B (parallelism) for dedup. 64K is the
    // sweet spot: a bigger spill (128K) dedups a bit better (faster load) but
    // makes B smaller, so each compactor's hashmap is larger and bucket-write's
    // hashmap/sort-unique cost dominates -- net SLOWER overall (measured: 100K
    // g16 went 3883s -> 4100s at 128K) and a higher bucket-process peak.
    // Overridable via --flush / --spill.
    static constexpr uint64_t DEFAULT_FLUSH_BASES = 4 * 1024;
    static constexpr uint64_t DEFAULT_SPILL_BYTES = 64 * 1024;

    // Total-RSS target for bucket-process, as a fraction of -g. The admission
    // gate keeps  carry + reserved working set + live color dict  under this,
    // with a hard live-RSS ceiling as backstop; the remaining ~18% is headroom
    // for in-flight load lag, the frag sink, and glibc fragmentation.
    static constexpr double BUCKET_PROCESS_BUDGET_FRAC = 0.82;

    // Share of -g the color-sets-dedup-map (the dominant non-spillable structure
    // in bucket-process) may use before its overflow path engages. It coexists
    // with the per-bucket walk sets under BUCKET_PROCESS_BUDGET_FRAC, so it gets
    // a sub-share. MEASURE-ONLY today (reported, not enforced) -- a starting
    // point to validate against real 100k/661k numbers before the overflow
    // machinery lands. See colorset-dedup-externalization.md.
    static constexpr double COLORSET_DEDUP_BUDGET_FRAC = 0.50;

    // Stitch's whole working-RAM share of -g (unitig writer + RAM-first stores +
    // phase-3 assemble transient all sized to sum within it). The remaining
    // (1-frac)*g is headroom for glibc arena retention + scratch, which scale with
    // thread churn rather than -g, so the fraction is deliberately conservative
    // (measured: at -t32 on 100k the stitch glibc retention alone was ~0.29*g).
    // The measured-RSS watcher (trip 0.70*g) is the backstop if the model still
    // undercounts; this fraction is what keeps the watcher from having to fire.
    static constexpr double STITCH_BUDGET_FRAC = 0.60;

    void validate_and_resolve_config() {
        if (m_cfg.filenames_list.empty()) {
            throw std::runtime_error("build_config::filenames_list is empty");
        }
        if (m_cfg.out_basename.empty()) {
            throw std::runtime_error("build_config::out_basename is empty");
        }
        if (m_cfg.k == 0 or m_cfg.k > MAX_K) {
            throw std::runtime_error("k must satisfy 1 <= k <= " + std::to_string(MAX_K));
        }
        if (m_cfg.num_threads == 0) m_cfg.num_threads = 1;
        // -g floor: below MIN_RAM_GB the non-spillable working set can't be held.
        if (m_cfg.max_ram_gb > 0 and m_cfg.max_ram_gb < MIN_RAM_GB) {
            throw std::runtime_error(
                "-g must be at least " + std::to_string((int)MIN_RAM_GB) + " GiB (got " +
                format_bytes((uint64_t)(m_cfg.max_ram_gb * 1024.0 * 1024.0 * 1024.0)) +
                "): below that the non-spillable per-thread and per-phase working structures "
                "cannot be honored. Raise -g, or omit it entirely for unbounded RAM.");
        }
        if (m_cfg.m == 0) m_cfg.m = compute_best_m(m_cfg.k);
        if (m_cfg.m < 2 or m_cfg.m > m_cfg.k) {
            throw std::runtime_error("invalid m=" + std::to_string(m_cfg.m) +
                                     " (need 2 <= m <= k)");
        }
        // Per-bucket batching knobs (CLI-overridable). The RAM model sizes the
        // bucket COUNT against these, so smaller values -> larger B -> smaller
        // buckets (faster bucket-process) while bucket-write stays in budget.
        m_flush_bases = m_cfg.flush_bases ? m_cfg.flush_bases : DEFAULT_FLUSH_BASES;
        m_spill_bytes = m_cfg.spill_bytes ? m_cfg.spill_bytes : DEFAULT_SPILL_BYTES;

        // bucket-write non-spillable floor (the "high thread count" case). The
        // RSS watcher can spill the compactor hashmaps but NOT the per-thread
        // keys/recs buffers, whose floor is num_buckets*(alpha*T*flush +
        // beta*spill). The model sizes num_buckets to fit BUCKET_WRITE_BUDGET_FRAC
        // *g but cannot go below MIN_AUTO_BUCKETS; if even that many buckets'
        // buffers exceed the share -- which happens as -t (T) inflates the
        // per-bucket term -- bucket-write would exceed -g unrecoverably. Abort,
        // and report the largest -t that fits. (Skipped under -b, which the user
        // fixes explicitly.)
        if (m_cfg.max_ram_gb > 0 and m_cfg.bucket_log2 == 0) {
            const double g = m_cfg.max_ram_gb * 1024.0 * 1024.0 * 1024.0;
            const double M = BUCKET_WRITE_BUDGET_FRAC * g / PLATFORM_RAM_OVERHEAD;
            const double per_bucket =
                m_cfg.alpha * m_cfg.num_threads * m_flush_bases + m_cfg.beta * m_spill_bytes;
            const double floor_min = (double)MIN_AUTO_BUCKETS * per_bucket;
            if (floor_min > M) {
                const double t_room = M / (double)MIN_AUTO_BUCKETS - m_cfg.beta * m_spill_bytes;
                const long t_max =
                    t_room > 0 ? (long)(t_room / (m_cfg.alpha * m_flush_bases)) : 0;
                std::string remedy =
                    t_max >= 1
                        ? "Reduce -t to <= " + std::to_string(t_max) + ", or raise -g."
                        : "Raise -g (even -t 1 does not fit; spill_bytes/flush_bases too large).";
                fail_ram_budget_("bucket-write per-thread buffers at -t " +
                                     std::to_string(m_cfg.num_threads),
                                 (uint64_t)floor_min, (uint64_t)M, remedy);
            }
        }

        // Resolve the bucket COUNT (need not be a power of two).
        if (m_cfg.bucket_log2 != 0) {
            // -b override: forces a power-of-two count, range-checked.
            if (m_cfg.bucket_log2 < MIN_BUCKETS_LOG2 or m_cfg.bucket_log2 > MAX_BUCKETS_LOG2)
                throw std::runtime_error("-b must be in [" + std::to_string(MIN_BUCKETS_LOG2) +
                                         ", " + std::to_string(MAX_BUCKETS_LOG2) + "]");
            m_cfg.num_buckets = 1u << m_cfg.bucket_log2;
        } else {
            m_cfg.num_buckets = auto_bucket_count_();
        }
        // bucket_writer opens one FILE per bucket (LZ4-framed, raw stdio).
        // Raise the soft FD limit if needed; clamp the count if even the
        // hard limit isn't enough.
        m_cfg.num_buckets = ensure_fd_capacity_for_buckets_(m_cfg.num_buckets);
    }

    // Size the bucket COUNT to the bucket-write budget so both of its pools fit:
    //
    //   W = alpha*T*B*flush + beta*B*spill = B*(alpha*T*flush + beta*spill) <= frac*g
    //   => B = frac*g / (alpha*T*flush + beta*spill)
    //
    // This is U-INDEPENDENT (bucket-write RAM doesn't depend on input size), so
    // there is no up-front estimate to get wrong. bucket-process then gets the
    // most buckets bucket-write can afford -- the best B for it regardless of U
    // -- and its runtime admission gate handles the (measured) per-bucket
    // kmer_info. No -g set -> historical default count.
    uint64_t auto_bucket_count_() const {
        if (m_cfg.max_ram_gb <= 0) return 1u << MIN_BUCKETS_LOG2;
        const double g = m_cfg.max_ram_gb * 1024.0 * 1024.0 * 1024.0;
        const double M = BUCKET_WRITE_BUDGET_FRAC * g / PLATFORM_RAM_OVERHEAD;
        const double T = std::max<uint32_t>(1, m_cfg.num_threads);
        const double per_bucket = m_cfg.alpha * T * m_flush_bases + m_cfg.beta * m_spill_bytes;
        const double b = M / per_bucket;
        // Floor at a small sane minimum (avoid degenerate single-bucket runs at
        // tiny -g); the fd clamp handles the ceiling. If the model lands below
        // num_threads, bucket-process simply won't engage every thread -- safe,
        // just less parallel, and surfaced by the bucket size distribution.
        uint64_t count = b < 1.0 ? 1 : b;
        if (count < MIN_AUTO_BUCKETS) count = MIN_AUTO_BUCKETS;
        if (count > UINT32_MAX) count = UINT32_MAX;
        return count;
    }

    // Try to raise RLIMIT_NOFILE so we can open `count` bucket files plus a
    // small headroom for stdin/stdout/sidecar/inputs. If even the hard limit
    // is too small, clamp the count down and warn.
    static uint32_t ensure_fd_capacity_for_buckets_(uint32_t count) {
        constexpr uint32_t HEADROOM = 64;
        struct rlimit r;
        if (::getrlimit(RLIMIT_NOFILE, &r) != 0) return count;  // best effort
        uint64_t needed = (uint64_t)count + HEADROOM;
        if (r.rlim_cur >= needed) return count;
        rlim_t target = std::min<rlim_t>((rlim_t)needed, r.rlim_max);
        struct rlimit nr = r;
        nr.rlim_cur = target;
        ::setrlimit(RLIMIT_NOFILE, &nr);
        ::getrlimit(RLIMIT_NOFILE, &nr);
        if (nr.rlim_cur >= needed) return count;
        // Hard limit too small. Clamp the count to what we can actually open.
        uint32_t avail = (uint32_t)nr.rlim_cur > HEADROOM ? (uint32_t)nr.rlim_cur - HEADROOM : 1;
        if (avail < count) {
            std::cerr << "warning: RLIMIT_NOFILE hard limit " << nr.rlim_max
                      << " can't accommodate " << count << " buckets; clamping to " << avail
                      << ". Raise the hard limit (e.g. ulimit -Hn) for tighter budgets.\n";
            return avail;
        }
        return count;
    }

    // Per-platform RAM-overhead multiplier. The bucket-count model uses
    // overhead constants (alpha/beta) whose calibration differs by allocator:
    // macOS libsystem_malloc has heavier per-allocation bookkeeping and counts
    // pages glibc would reclaim, so the same logical state sits higher in RSS.
    // We compensate by dividing the budget the model sizes against by this
    // constant (2.0 on macOS, 1.0 on Linux), so the bucket count comes out
    // smaller on macOS and the footprint still fits once platform overhead is
    // added back. The watcher's thresholds are NOT divided by it -- the watcher
    // measures real RSS, which is what we want to cap.
    static constexpr double PLATFORM_RAM_OVERHEAD =
#if defined(__APPLE__)
        2.0;
#else
        1.0;
#endif

    // If m_cfg.tmp_dir is empty, mkdtemp under $TMPDIR. Otherwise use it
    // (creating it if missing); if it already exists it must be an empty
    // directory, since cleanup_tmp_dir wipes the whole thing.
    std::string resolve_tmp_dir() {
        if (m_cfg.tmp_dir.empty()) {
            std::string tmpl =
                (std::filesystem::temp_directory_path() / "cdbg_buckets_XXXXXX").string();
            std::vector<char> buf(tmpl.begin(), tmpl.end());
            buf.push_back('\0');
            if (mkdtemp(buf.data()) == nullptr) {
                throw std::runtime_error(std::string("mkdtemp failed: ") + std::strerror(errno));
            }
            return std::string(buf.data());
        }
        std::error_code ec;
        if (std::filesystem::exists(m_cfg.tmp_dir, ec)) {
            if (!std::filesystem::is_directory(m_cfg.tmp_dir, ec)) {
                throw std::runtime_error("-d " + m_cfg.tmp_dir + " exists but is not a directory");
            }
            if (!std::filesystem::is_empty(m_cfg.tmp_dir, ec)) {
                throw std::runtime_error("-d " + m_cfg.tmp_dir +
                                         " is not empty (the tool will remove the directory on"
                                         " exit, so it must start empty)");
            }
        } else {
            std::filesystem::create_directories(m_cfg.tmp_dir, ec);
            if (ec) {
                throw std::runtime_error("cannot create -d " + m_cfg.tmp_dir + ": " + ec.message());
            }
        }
        return m_cfg.tmp_dir;
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

    static uint64_t pick_unitig_bucket_count_(uint64_t num_color_classes,
                                              uint64_t total_seq_bytes_estimate,
                                              double max_ram_gb) {
        if (num_color_classes == 0) return 1;
        constexpr uint64_t MIN_K = 16;
        constexpr uint64_t MAX_K = 1024;
        constexpr uint64_t DEFAULT_K = 64;
        constexpr double OVERHEAD = 2.0;
        constexpr double SHARE = 0.10;

        if (max_ram_gb <= 0 or total_seq_bytes_estimate == 0) {
            return std::min<uint64_t>(num_color_classes, DEFAULT_K);
        }
        const uint64_t budget_bytes = max_ram_gb * 1024.0 * 1024.0 * 1024.0 * SHARE;
        if (budget_bytes == 0) { return std::min<uint64_t>(num_color_classes, DEFAULT_K); }
        const uint64_t needed = total_seq_bytes_estimate * OVERHEAD;
        uint64_t k = (needed + budget_bytes - 1) / budget_bytes;
        if (k < MIN_K) k = MIN_K;
        if (k > MAX_K) k = MAX_K;
        if (k > num_color_classes) k = num_color_classes;
        return k;
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
    static uint64_t pick_stitch_buckets_(uint64_t total_frag_seq_bytes, double max_ram_gb,
                                         uint32_t num_threads) {
        constexpr uint64_t MIN_BUCKETS = 64;
        constexpr uint64_t MAX_BUCKETS = 1u << 20;  // 1M files cap (fd headroom)
        constexpr uint64_t DEFAULT_BUCKETS = 1024;
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
        if (max_ram_gb <= 0 or total_frag_seq_bytes == 0) return DEFAULT_BUCKETS;

        // Share of what REMAINS after the live carry-in (the color dict is
        // freed before stitch via release_index, so this is small now).
        const uint64_t carry_in = current_rss_bytes();
        const uint64_t budget_total = max_ram_gb * 1024.0 * 1024.0 * 1024.0;
        const uint64_t avail = budget_total > carry_in ? (budget_total - carry_in) : budget_total;
        const uint64_t budget = avail * SHARE;
        if (budget == 0) return DEFAULT_BUCKETS;

        const uint64_t needed = total_frag_seq_bytes * OVERHEAD;
        // count s.t. needed/count * num_threads <= budget (num_threads resident).
        uint64_t count = ((needed * num_threads) + budget - 1) / budget;
        if (count < MIN_BUCKETS) count = MIN_BUCKETS;
        if (count > MAX_BUCKETS) count = MAX_BUCKETS;  // hard fd ceiling; see note below
        return count;
    }

    build_config m_cfg;
    // num_colors is enforced < 2^32 (one per input file) but kept 64-bit so it
    // never needs a narrowing cast as it flows through the pipeline. It is
    // serialized as a u32 only at the .color_sets header boundary (the dicts).
    uint64_t m_num_colors = 0;
    uint64_t m_num_unitigs = 0;
    uint64_t m_num_color_classes = 0;
    uint64_t m_peak_rss_bytes = 0;
    // bucket-write batching payload, fixed at good defaults in
    // validate_and_resolve_config(); the bucket COUNT is sized against these.
    uint64_t m_flush_bases = 0;
    uint64_t m_spill_bytes = 0;
};

}  // namespace cdbg
