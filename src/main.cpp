#include <algorithm>
#include <atomic>
#include <cerrno>
#include <charconv>
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

#include <essentials.hpp>
#include <parser.hpp>

#include "bucket_io.hpp"
#include "bucket_ingester.hpp"
#include "bucket_walker.hpp"
#include "color_set_dict.hpp"
#include "hybrid_color_sets.hpp"
#include "minimizer.hpp"
#include "stitch.hpp"
#include "util.hpp"

namespace {
class timer {
public:
    timer(const char* label) : m_label(label), m_t0(std::chrono::steady_clock::now()) {}
    ~timer() {
        auto t1 = std::chrono::steady_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - m_t0).count();
        std::cout << "[" << m_label << "] " << (ms / 1000.0) << " s\n";
    }

private:
    const char* m_label;
    std::chrono::steady_clock::time_point m_t0;
};

bool parse_args(int argc, char** argv, cdgb::build_config& cfg) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("filenames_list",
               "Text file with one input path per line. The file at line i has color i.", "-i",
               true);
    parser.add("out_basename", "Output basename. Produces <basename>.fa and <basename>.colors.",
               "-o", true);
    parser.add("k", "K-mer length (must be <= " + std::to_string(cdgb::MAX_K) + ").", "-k", true);
    parser.add("num_threads", "Number of worker threads (default 1).", "-t", false);
    parser.add("m", "Minimizer length (default: auto, derived from k).", "-m", false);
    parser.add("buckets_log2", "log2 of the bucket count (default 10 -> 1024).", "--buckets-log2",
               false);
    parser.add("tmp_dir", "Scratch directory for bucket files (default: mkdtemp under $TMPDIR).",
               "--tmp-dir", false);
    parser.add("verbose", "Verbose output.", "--verbose", false, true);

    if (!parser.parse()) return false;

    cfg.filenames_list = parser.get<std::string>("filenames_list");
    cfg.out_basename = parser.get<std::string>("out_basename");
    cfg.k = parser.get<uint32_t>("k");
    if (parser.parsed("num_threads")) cfg.num_threads = parser.get<uint32_t>("num_threads");
    if (parser.parsed("m")) cfg.m = parser.get<uint32_t>("m");
    if (parser.parsed("buckets_log2")) cfg.bucket_log2 = parser.get<uint32_t>("buckets_log2");
    if (parser.parsed("tmp_dir")) cfg.tmp_dir = parser.get<std::string>("tmp_dir");
    if (parser.parsed("verbose")) cfg.verbose = parser.get<bool>("verbose");

    if (cfg.k == 0 || cfg.k > cdgb::MAX_K) {
        std::cerr << "error: k must satisfy 1 <= k <= " << cdgb::MAX_K << "\n";
        return false;
    }
    return true;
}

std::vector<std::string> read_filenames(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open filenames list: " + path);
    std::vector<std::string> v;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty()) v.push_back(line);
    }
    return v;
}

}  // namespace

int main(int argc, char** argv) {
    cdgb::build_config cfg;
    if (!parse_args(argc, argv, cfg)) return 1;

    auto files = read_filenames(cfg.filenames_list);
    if (files.empty()) {
        std::cerr << "no input files\n";
        return 1;
    }
    if (files.size() > (uint64_t)UINT32_MAX) {
        std::cerr << "too many colors (max 2^32 - 1)\n";
        return 1;
    }

    if (cfg.m == 0) cfg.m = cdgb::compute_best_m(cfg.k);
    if (cfg.m < 2 || cfg.m > cfg.k) {
        std::cerr << "invalid m=" << cfg.m << " (need 2 <= m <= k)\n";
        return 1;
    }
    const uint32_t num_buckets = 1u << cfg.bucket_log2;

    std::cout << "k = " << cfg.k << ", m = " << cfg.m << ", num_colors = " << files.size()
              << ", num_threads = " << cfg.num_threads << ", num_buckets = " << num_buckets << "\n";

    auto t_start = std::chrono::steady_clock::now();
    auto print_total = [&] {
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                      std::chrono::steady_clock::now() - t_start)
                      .count();
        std::cout << "[total construction time] " << (ms / 1000.0) << " s\n";
    };

    // Resolve a scratch directory. If the user supplied --tmp-dir, we
    // create it if missing. Once construction starts, the tool wipes
    // the entire directory at the end -- so if the path already exists,
    // it must be a directory and it must be empty (we won't clobber
    // unrelated user files).
    std::string tmp_dir = cfg.tmp_dir;
    if (tmp_dir.empty()) {
        std::string tmpl =
            (std::filesystem::temp_directory_path() / "cdgb_buckets_XXXXXX").string();
        std::vector<char> buf(tmpl.begin(), tmpl.end());
        buf.push_back('\0');
        if (mkdtemp(buf.data()) == nullptr) {
            std::cerr << "mkdtemp failed: " << std::strerror(errno) << "\n";
            return 1;
        }
        tmp_dir = buf.data();
    } else {
        std::error_code ec;
        if (std::filesystem::exists(tmp_dir, ec)) {
            if (!std::filesystem::is_directory(tmp_dir, ec)) {
                std::cerr << "error: --tmp-dir " << tmp_dir << " exists but is not a directory\n";
                return 1;
            }
            if (!std::filesystem::is_empty(tmp_dir, ec)) {
                std::cerr << "error: --tmp-dir " << tmp_dir
                          << " is not empty (the tool will remove the directory on exit, so it"
                             " must start empty)\n";
                return 1;
            }
        } else {
            std::filesystem::create_directories(tmp_dir, ec);
            if (ec) {
                std::cerr << "error: cannot create --tmp-dir " << tmp_dir << ": " << ec.message()
                          << "\n";
                return 1;
            }
        }
    }
    std::cout << "  tmp_dir = " << tmp_dir << "\n";

    cdgb::bucket_writer writer(tmp_dir, num_buckets);
    {
        timer _("bucket-write");
        std::atomic<uint64_t> done{0};
        cdgb::progress prog("bucket-write", done, files.size());
        cdgb::ingest_bucketed(files, cfg.k, cfg.m, cfg.bucket_log2, writer, cfg.num_threads, &done);
        prog.stop();
    }
    writer.close();
    std::cout << "  bucket bytes written: " << writer.total_bytes() << " (compressed; "
              << writer.total_uncompressed_bytes() << " uncompressed)\n";
    cdgb::bucket_prof().print(cfg.num_threads);

    // Bucket processing emits stitchable fragments AND merges per-bucket
    // color sets into the shared global dict on the fly. By the time we
    // exit this block, every fragment already carries a global cid; no
    // separate single-threaded intern pass is needed.
    std::vector<cdgb::stitchable_unitig> frag_unitigs;
    std::mutex out_mu;
    cdgb::color_set_dict global_dict;
    std::mutex global_mu;
    {
        timer _("bucket-process");
        std::atomic<uint64_t> done{0};
        cdgb::progress prog("bucket-process", done, num_buckets);
        cdgb::process_buckets(writer, cfg.k, cfg.num_threads, frag_unitigs, out_mu, global_dict,
                              global_mu, &done);
        prog.stop();
        std::cout << "  bucket fragments: " << frag_unitigs.size() << "\n";
        std::cout << "  distinct color classes: " << global_dict.size() << "\n";
    }

    std::vector<cdgb::stitchable_unitig> all_unitigs;
    {
        timer _("stitch");
        std::atomic<uint64_t> done{0};
        cdgb::progress prog("stitch", done, frag_unitigs.size());
        cdgb::stitch_unitigs(frag_unitigs, cfg.k, all_unitigs, &done);
        prog.stop();
        frag_unitigs = {};
        std::cout << "  unitigs after stitching: " << all_unitigs.size() << "\n";
    }

    // Emit FASTA + .colors. Both files use the cid that's already on
    // each unitig (assigned during process_buckets when each bucket's
    // local color_set_dict was merged into the global one).
    {
        std::vector<std::vector<size_t>> by_class(global_dict.size());
        for (size_t i = 0; i < all_unitigs.size(); ++i) {
            by_class[all_unitigs[i].cid].push_back(i);
        }

        // FASTA write. std::ofstream's default ~8 KiB buffer + per-token
        // formatting via `<<` is slow on millions of tiny records. Use a
        // hand-rolled 1 MiB buffer fed by std::to_chars (for the integer
        // header) + std::memcpy (for the sequence) and a single fwrite
        // when the buffer fills. About 3-5x faster than the stream layer
        // on the 1.9M-unitig output.
        {
            timer _("emit fasta");
            FILE* fa = std::fopen((cfg.out_basename + ".fa").c_str(), "wb");
            if (!fa) throw std::runtime_error("cannot open " + cfg.out_basename + ".fa");
            constexpr size_t BUF_BYTES = 1 << 20;
            std::vector<char> buf(BUF_BYTES);
            size_t pos = 0;
            auto flush_buf = [&] {
                if (pos == 0) return;
                if (std::fwrite(buf.data(), 1, pos, fa) != pos) {
                    std::fclose(fa);
                    throw std::runtime_error("short write to " + cfg.out_basename + ".fa");
                }
                pos = 0;
            };
            auto reserve = [&](size_t n) {
                if (pos + n > BUF_BYTES) flush_buf();
            };
            for (uint32_t cid = 0; cid < global_dict.size(); ++cid) {
                for (size_t idx : by_class[cid]) {
                    const std::string& seq = all_unitigs[idx].seq;
                    // Header: '>' + decimal cid + '\n' is at most 12 chars
                    // for any uint32_t. The sequence write below handles
                    // arbitrary lengths via chunking.
                    reserve(12);
                    buf[pos++] = '>';
                    auto r = std::to_chars(buf.data() + pos, buf.data() + pos + 11, cid);
                    pos = (size_t)(r.ptr - buf.data());
                    buf[pos++] = '\n';
                    // Sequence: memcpy into the buffer, flushing when full.
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

        {
            timer _("emit colors");
            cdgb::hybrid_builder hb((uint32_t)files.size());
            hb.encode_parallel(global_dict, cfg.num_threads);
            cdgb::hybrid h;
            hb.build(h);
            essentials::save(h, (cfg.out_basename + ".colors").c_str());
        }
    }

    // Clean up the scratch directory. We've enforced at startup that we
    // own its contents (either we mkdtemp'd it, or the user passed an
    // empty directory), so a single recursive remove is safe.
    {
        timer _("removing tmp files");
        std::error_code ec;
        std::filesystem::remove_all(tmp_dir, ec);
    }

    std::cout << "done. wrote " << cfg.out_basename << ".fa and " << cfg.out_basename
              << ".colors\n";
    print_total();
    return 0;
}
