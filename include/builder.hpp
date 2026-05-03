#pragma once

// Public API for building a colored compacted dBG.
//
// Wraps the full pipeline (minimizer-bucketed ingest -> per-bucket dBG
// build + global color interning -> stitch -> emit FASTA + .colors) so
// that downstream tools can construct a `build_config`, instantiate a
// `builder`, and call `build()`. Mirrors the call shape used by
// Fulgor's `index<ColorSets>::builder` for its ccdBG dependency
// (https://github.com/jermp/fulgor/blob/main/include/builders/builder.hpp).
//
// Typical use:
//
//   cdgb::build_config cfg;
//   cfg.filenames_list = "...";
//   cfg.out_basename   = "...";
//   cfg.k = 31;
//   cfg.num_threads = 8;
//
//   cdgb::builder b(cfg);
//   b.build();
//   // After build():
//   //   b.num_colors()         -- one per input file
//   //   b.num_unitigs()        -- count of stitched unitigs in <basename>.fa
//   //   b.num_color_classes()  -- count of distinct color sets in <basename>.colors
//
// build() throws std::runtime_error on configuration errors or I/O
// failures. On success, <basename>.fa and <basename>.colors are written
// and the scratch directory (cfg.tmp_dir or an mkdtemp'd one) is removed.

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

#include <essentials.hpp>

#include "bucket_io.hpp"
#include "bucket_ingester.hpp"
#include "bucket_walker.hpp"
#include "color_set_dict.hpp"
#include "hybrid_color_sets.hpp"
#include "minimizer.hpp"
#include "stitch.hpp"
#include "streaming_color_set_dict.hpp"
#include "util.hpp"

namespace cdgb {

struct builder {
    builder() = default;
    explicit builder(build_config const& cfg) : m_cfg(cfg) {}

    // Run the full pipeline. Throws std::runtime_error on bad config or
    // I/O error. On success, writes m_cfg.out_basename + {".fa", ".colors"}
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

        bucket_writer writer(tmp_dir, num_buckets, m_flush_bases, m_spill_bytes);
        // When --max-ram is set, arm a background RSS watcher with
        // hysteresis. Bucket-write must leave room for what comes
        // after: bucket-process adds ~1 GiB on top on multi-thousand-
        // genome inputs, stitch adds ~300 MiB. We reserve those by
        // budgeting bucket-write at 60% of --max-ram (HIGH) with a
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
            writer.start_rss_watcher(high_threshold_bytes, low_threshold_bytes);
        }
        {
            phase_rss_marker rss("bucket-write");
            {
                timer _("bucket-write");
                std::atomic<uint64_t> done{0};
                progress prog("bucket-write", done, files.size());
                ingest_bucketed(files, m_cfg.k, m_cfg.m, m_cfg.bucket_log2, writer,
                                m_cfg.num_threads, &done);
                prog.stop();
            }
            // Snapshot before writer.close() frees the compactor
            // hashmaps -- their memory is part of the bucket-write peak.
            rss.stop();
        }
        // Stop the RSS watcher before close(): we don't want a stray
        // try_spill firing during close()'s final spill+gzclose loop.
        writer.stop_rss_watcher();
        if (writer.pressure_was_engaged()) {
            std::cout << "  bucket-write: RSS pressure engaged; observed live-RSS high "
                      << format_bytes(writer.observed_rss_high()) << "\n";
        }
        writer.close();
        std::cout << "  bucket bytes written: " << writer.total_bytes() << " (compressed; "
                  << writer.total_uncompressed_bytes() << " uncompressed)\n";
        bucket_prof().print(m_cfg.num_threads);

        // Bucket processing emits stitchable fragments AND merges per-bucket
        // color sets into the shared global dict on the fly. By the time we
        // exit this block, every fragment already carries a global cid.
        std::vector<stitchable_unitig> frag_unitigs;
        std::mutex out_mu;
        // Streaming dict: encodes each new color set into its bvb at
        // intern() time, holding only metadata (32 B/class) plus the
        // compressed bits. Compared to the previous in-RAM-vectors
        // dict, this caps peak RAM during the dominant bucket-process
        // phase by the *compressed* color-set size, not the sum of
        // class sizes.
        streaming_color_set_dict global_dict(m_num_colors);
        if (m_color_bvb_spill_bytes > 0) {
            global_dict.enable_spill(tmp_dir + "/colors.bits", m_color_bvb_spill_bytes);
        }
        std::mutex global_mu;
        {
            phase_rss_marker rss("bucket-process");
            {
                timer _("bucket-process");
                std::atomic<uint64_t> done{0};
                progress prog("bucket-process", done, num_buckets);
                process_buckets(writer, m_cfg.k, m_cfg.num_threads, frag_unitigs, out_mu,
                                global_dict, global_mu, &done);
                prog.stop();
                std::cout << "  bucket fragments: " << frag_unitigs.size() << "\n";
                std::cout << "  distinct color classes: " << global_dict.size() << "\n";
            }
            rss.stop();
        }
        m_num_color_classes = global_dict.size();

        std::vector<stitchable_unitig> all_unitigs;
        {
            phase_rss_marker rss("stitch");
            {
                timer _("stitch");
                std::atomic<uint64_t> done{0};
                progress prog("stitch", done, frag_unitigs.size());
                stitch_unitigs(frag_unitigs, m_cfg.k, all_unitigs, &done);
                prog.stop();
                frag_unitigs = {};
                std::cout << "  unitigs after stitching: " << all_unitigs.size() << "\n";
            }
            rss.stop();
        }
        m_num_unitigs = all_unitigs.size();

        {
            phase_rss_marker rss("emit");
            emit_fasta(all_unitigs, global_dict.size());
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

        std::cout << "done. wrote " << m_cfg.out_basename << ".fa and " << m_cfg.out_basename
                  << ".colors\n";
        print_total();
    }

    // Stats populated by build(); zero before build() runs.
    uint32_t num_colors() const { return m_num_colors; }
    uint64_t num_unitigs() const { return m_num_unitigs; }
    uint64_t num_color_classes() const { return m_num_color_classes; }
    uint64_t peak_rss_bytes() const { return m_peak_rss_bytes; }
    build_config const& config() const { return m_cfg; }

private:
    // Soft-cap policy. Reads m_cfg.max_ram_gb (0 = unset) and decides:
    //   - bucket_log2 (more buckets -> smaller per-bucket data structures)
    //   - color-bvb spill threshold (in bytes; 0 = never spill)
    // Tighter budgets push bucket_log2 toward 13 (8192 buckets, GGCAT's
    // upper end) and shrink the bvb spill threshold proportionally.
    // The user-facing CLI override (--buckets-log2) takes precedence
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
            // User didn't pin --buckets-log2. Auto-pick.
            m_cfg.bucket_log2 = auto_bucket_log2();
        } else if (m_cfg.bucket_log2 < MIN_BUCKETS_LOG2 || m_cfg.bucket_log2 > MAX_BUCKETS_LOG2) {
            throw std::runtime_error("--buckets-log2 must be in [" +
                                     std::to_string(MIN_BUCKETS_LOG2) + ", " +
                                     std::to_string(MAX_BUCKETS_LOG2) + "]");
        }
        // bucket_writer opens one gzFile per bucket. Raise the soft FD
        // limit if needed; clamp bucket_log2 if even the hard limit
        // isn't enough.
        m_cfg.bucket_log2 = ensure_fd_capacity_for_buckets_(m_cfg.bucket_log2);
        m_color_bvb_spill_bytes = auto_color_bvb_spill_bytes();
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

    // Heuristic, balancing two competing pressures:
    //   - bucket-process per-bucket data shrinks with more buckets.
    //   - bucket-write zlib state grows ~448 KiB per bucket (un-tunable
    //     internal state + gzbuffer), which becomes the dominant peak
    //     contributor on tight budgets.
    //
    // Strategy: pick the largest bucket_log2 such that the zlib floor
    // num_buckets * ~448 KiB stays under ~30% of --max-ram. That keeps
    // bucket-write's unspillable structural floor bounded; the
    // remaining 70% absorbs per-thread buffers, compactor data, and
    // the bucket-process working set.
    //
    // No --max-ram set -> historical 1024 buckets.
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
    // tune behaves as if the user passed half their --max-ram (smaller
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

    double effective_max_ram_gb() const {
        return m_cfg.max_ram_gb / PLATFORM_RAM_OVERHEAD;
    }

    // Pick the bucket count by enforcing that the *full* bucket-write
    // structural floor at MIN values fits in the bucket-write share of
    // the (platform-derated) budget. The structural floor at the
    // floor-tunable values is:
    //
    //   B * (compressor_per_bucket
    //        + T * MIN_FLUSH_BASES * BUFFER_OVERHEAD
    //        + MIN_SPILL_BYTES * COMPACTOR_OVERHEAD)
    //
    // i.e. the per-bucket cost when flush_bases / spill_bytes have
    // already collapsed to their minimums. Picking log2 such that this
    // fits in SHARE * effective_ram leaves room for the auto-tune to
    // actually choose values above the floor (and so to spill less
    // often, write fewer disk records, and stay within real RSS).
    //
    // The earlier "compressor-only" check picked too high a log2 once
    // LZ4 made the compressor cheap: at --max-ram 4 / 8 threads it
    // accepted 8192 buckets, which collapsed flush_bases and
    // spill_bytes to floors and pushed the per-thread + compactor
    // contributions far past budget.
    //
    // No --max-ram set -> historical 1024 buckets.
    uint32_t auto_bucket_log2() const {
        if (m_cfg.max_ram_gb <= 0) return MIN_BUCKETS_LOG2;
        // Match the constants used by auto_tune_bucket_write_params --
        // they have to be consistent for this estimate to be right.
        constexpr double SHARE = 0.50;
        constexpr double COMPACTOR_OVERHEAD = 5.0;
        constexpr double BUFFER_OVERHEAD = 2.0;
        constexpr uint64_t COMPRESSOR_BYTES_PER_BUCKET = 64ull * 1024ull;
        constexpr uint64_t MIN_FLUSH_BASES = 4ull * 1024ull;
        constexpr uint64_t MIN_SPILL_BYTES = 16ull * 1024ull;

        const uint64_t T = std::max<uint64_t>(1, m_cfg.num_threads);
        const uint64_t per_bucket_floor =
            COMPRESSOR_BYTES_PER_BUCKET +
            (uint64_t)(T * MIN_FLUSH_BASES * BUFFER_OVERHEAD) +
            (uint64_t)(MIN_SPILL_BYTES * COMPACTOR_OVERHEAD);
        const double share_bytes =
            effective_max_ram_gb() * 1024.0 * 1024.0 * 1024.0 * SHARE;

        uint32_t log2 = MAX_BUCKETS_LOG2;
        while (log2 > MIN_BUCKETS_LOG2 &&
               (double)((uint64_t)1 << log2) * (double)per_bucket_floor > share_bytes) {
            --log2;
        }
        return log2;
    }

    // Heuristic. Cap the streaming-color bvb at ~1/4 of the budget, so
    // it spills before crowding out the rest.
    uint64_t auto_color_bvb_spill_bytes() const {
        if (m_cfg.max_ram_gb <= 0) return 0;  // never spill
        double bytes = m_cfg.max_ram_gb * (1024.0 * 1024.0 * 1024.0) * 0.25;
        if (bytes < (double)(64 * 1024 * 1024)) bytes = 64 * 1024 * 1024;  // floor at 64 MiB
        return (uint64_t)bytes;
    }

    // Joint auto-tune of (flush_bases, spill_bytes) against a fixed
    // share of --max-ram, given (num_threads, num_buckets).
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
    // KiB per gzFile) but with LZ4 frame (max64KB block,
    // blockIndependent, autoFlush=1) the per-bucket cctx + small
    // output scratch is ~40-64 KiB. We bake that conservatively into
    // COMPRESSOR_BYTES_PER_BUCKET = 64 KiB.
    //
    // Strategy: reserve a target share (default 50%) of --max-ram for
    // bucket-write. Subtract the (small) compressor floor. Split what
    // remains evenly between per-thread buffers and compactor data,
    // then solve for flush_bases and spill_bytes.
    //
    // No --max-ram set -> keep historical defaults (this keeps the
    // small-input dev path identical and avoids surprising regressions
    // for users who don't care about a budget).
    void auto_tune_bucket_write_params() {
        if (m_cfg.max_ram_gb <= 0) {
            m_flush_bases = 64 * 1024;
            m_spill_bytes = DEFAULT_COMPACTOR_SPILL_BYTES;  // 256 KiB
            return;
        }
        constexpr double SHARE = 0.50;              // half the budget for bucket-write
        constexpr double COMPACTOR_OVERHEAD = 5.0;  // compactor: structure + key allocs + frag
        constexpr double BUFFER_OVERHEAD = 2.0;     // per-thread: capacity slack + recs vec
        constexpr size_t COMPRESSOR_BYTES_PER_BUCKET = 64 * 1024;  // LZ4 cctx + out_buf
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
                (std::filesystem::temp_directory_path() / "cdgb_buckets_XXXXXX").string();
            std::vector<char> buf(tmpl.begin(), tmpl.end());
            buf.push_back('\0');
            if (mkdtemp(buf.data()) == nullptr)
                throw std::runtime_error(std::string("mkdtemp failed: ") + std::strerror(errno));
            return std::string(buf.data());
        }
        std::error_code ec;
        if (std::filesystem::exists(m_cfg.tmp_dir, ec)) {
            if (!std::filesystem::is_directory(m_cfg.tmp_dir, ec))
                throw std::runtime_error("--tmp-dir " + m_cfg.tmp_dir +
                                         " exists but is not a directory");
            if (!std::filesystem::is_empty(m_cfg.tmp_dir, ec))
                throw std::runtime_error("--tmp-dir " + m_cfg.tmp_dir +
                                         " is not empty (the tool will remove the directory on"
                                         " exit, so it must start empty)");
        } else {
            std::filesystem::create_directories(m_cfg.tmp_dir, ec);
            if (ec)
                throw std::runtime_error("cannot create --tmp-dir " + m_cfg.tmp_dir + ": " +
                                         ec.message());
        }
        return m_cfg.tmp_dir;
    }

    static void cleanup_tmp_dir(std::string const& tmp_dir) {
        timer _("removing tmp files");
        std::error_code ec;
        std::filesystem::remove_all(tmp_dir, ec);
    }

    // FASTA emit. Hand-rolled 1 MiB buffer + std::to_chars for the integer
    // header + memcpy for the sequence body. Significantly faster than
    // std::ofstream's default 8 KiB buffer + stream operators on millions
    // of small records.
    void emit_fasta(std::vector<stitchable_unitig> const& all_unitigs,
                    uint32_t num_color_classes) const {
        timer _("emit fasta");
        std::vector<std::vector<size_t>> by_class(num_color_classes);
        for (size_t i = 0; i < all_unitigs.size(); ++i) {
            by_class[all_unitigs[i].cid].push_back(i);
        }

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
        for (uint32_t cid = 0; cid < num_color_classes; ++cid) {
            for (size_t idx : by_class[cid]) {
                std::string const& seq = all_unitigs[idx].seq;
                // Header: '>' + decimal cid + '\n' fits in 12 chars for any uint32_t.
                reserve(12);
                buf[pos++] = '>';
                auto r = std::to_chars(buf.data() + pos, buf.data() + pos + 11, cid);
                pos = (size_t)(r.ptr - buf.data());
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
            }
        }
        flush_buf();
        std::fclose(fa);
    }

    void emit_colors(streaming_color_set_dict& global_dict) const {
        timer _("emit colors");
        // Encoding already happened during bucket-process via
        // global_dict.intern(); finalize() just builds the hybrid
        // wrapper (moves the bit_vector out of the bvb, builds an
        // elias_fano of the per-class bit_offsets) so we can
        // essentials::save it.
        hybrid h;
        global_dict.finalize(h);
        essentials::save(h, (m_cfg.out_basename + ".colors").c_str());
    }

    build_config m_cfg;
    uint32_t m_num_colors = 0;
    uint64_t m_num_unitigs = 0;
    uint64_t m_num_color_classes = 0;
    uint64_t m_peak_rss_bytes = 0;
    uint64_t m_color_bvb_spill_bytes = 0;  // 0 = never spill
    // bucket-write knobs picked by auto_tune_bucket_write_params().
    size_t m_flush_bases = 0;
    size_t m_spill_bytes = 0;
};

}  // namespace cdgb
