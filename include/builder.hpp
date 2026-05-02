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

        auto const t_start = std::chrono::steady_clock::now();
        auto print_total = [&] {
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                          std::chrono::steady_clock::now() - t_start)
                          .count();
            std::cout << "[total construction time] " << (ms / 1000.0) << " s\n";
        };

        std::string const tmp_dir = resolve_tmp_dir();
        std::cout << "  tmp_dir = " << tmp_dir << "\n";

        bucket_writer writer(tmp_dir, num_buckets);
        {
            timer _("bucket-write");
            std::atomic<uint64_t> done{0};
            progress prog("bucket-write", done, files.size());
            ingest_bucketed(files, m_cfg.k, m_cfg.m, m_cfg.bucket_log2, writer, m_cfg.num_threads,
                            &done);
            prog.stop();
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
        color_set_dict global_dict;
        std::mutex global_mu;
        {
            timer _("bucket-process");
            std::atomic<uint64_t> done{0};
            progress prog("bucket-process", done, num_buckets);
            process_buckets(writer, m_cfg.k, m_cfg.num_threads, frag_unitigs, out_mu, global_dict,
                            global_mu, &done);
            prog.stop();
            std::cout << "  bucket fragments: " << frag_unitigs.size() << "\n";
            std::cout << "  distinct color classes: " << global_dict.size() << "\n";
        }
        m_num_color_classes = global_dict.size();

        std::vector<stitchable_unitig> all_unitigs;
        {
            timer _("stitch");
            std::atomic<uint64_t> done{0};
            progress prog("stitch", done, frag_unitigs.size());
            stitch_unitigs(frag_unitigs, m_cfg.k, all_unitigs, &done);
            prog.stop();
            frag_unitigs = {};
            std::cout << "  unitigs after stitching: " << all_unitigs.size() << "\n";
        }
        m_num_unitigs = all_unitigs.size();

        emit_fasta(all_unitigs, global_dict.size());
        emit_colors(global_dict);

        cleanup_tmp_dir(tmp_dir);

        // Peak resident set size across the whole build, as tracked by
        // the kernel (VmHWM in /proc/self/status). One read at the end;
        // zero hot-path cost.
        m_peak_rss_bytes = process_peak_rss_bytes();
        if (m_peak_rss_bytes) {
            std::cout << "[peak resident memory] " << format_bytes(m_peak_rss_bytes) << "\n";
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

    void emit_colors(color_set_dict const& global_dict) const {
        timer _("emit colors");
        hybrid_builder hb(m_num_colors);
        hb.encode_parallel(global_dict, m_cfg.num_threads);
        hybrid h;
        hb.build(h);
        essentials::save(h, (m_cfg.out_basename + ".colors").c_str());
    }

    build_config m_cfg;
    uint32_t m_num_colors = 0;
    uint64_t m_num_unitigs = 0;
    uint64_t m_num_color_classes = 0;
    uint64_t m_peak_rss_bytes = 0;
};

}  // namespace cdgb
