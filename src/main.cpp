#include <atomic>
#include <cerrno>
#include <chrono>
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
#include "build_config.hpp"
#include "color_set_dict.hpp"
#include "hybrid_color_sets.hpp"
#include "minimizer.hpp"
#include "prof.hpp"
#include "progress.hpp"
#include "stitch.hpp"

namespace {
class Timer {
public:
    Timer(const char* label) : m_label(label), m_t0(std::chrono::steady_clock::now()) {}
    ~Timer() {
        auto t1 = std::chrono::steady_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - m_t0).count();
        std::cout << "[" << m_label << "] " << (ms / 1000.0) << " s\n";
    }

private:
    const char* m_label;
    std::chrono::steady_clock::time_point m_t0;
};

bool parse_args(int argc, char** argv, cdgb::BuildConfig& cfg) {
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
    cdgb::BuildConfig cfg;
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

    // Resolve a scratch directory.
    std::string tmp_dir = cfg.tmp_dir;
    bool tmp_owned = false;
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
        tmp_owned = true;
    }
    std::cout << "  tmp_dir = " << tmp_dir << "\n";

    cdgb::BucketWriter writer(tmp_dir, num_buckets);
    {
        Timer _("bucket-write");
        std::atomic<uint64_t> done{0};
        cdgb::Progress prog("bucket-write", done, files.size());
        cdgb::ingest_bucketed(files, cfg.k, cfg.m, cfg.bucket_log2, writer, cfg.num_threads, &done);
        prog.stop();
    }
    writer.close();
    std::cout << "  bucket bytes written: " << writer.total_bytes() << " (compressed; "
              << writer.total_uncompressed_bytes() << " uncompressed)\n";
    cdgb::bucket_prof().print(cfg.num_threads);

    std::vector<cdgb::StitchableUnitig> frag_unitigs;
    std::mutex out_mu;
    {
        Timer _("bucket-process");
        std::atomic<uint64_t> done{0};
        cdgb::Progress prog("bucket-process", done, num_buckets);
        cdgb::process_buckets(writer, cfg.k, cfg.num_threads, frag_unitigs, out_mu, &done);
        prog.stop();
        std::cout << "  bucket fragments: " << frag_unitigs.size() << "\n";
    }

    std::vector<cdgb::StitchableUnitig> all_unitigs;
    {
        Timer _("stitch");
        cdgb::stitch_unitigs(frag_unitigs, cfg.k, all_unitigs);
        frag_unitigs = {};
        std::cout << "  unitigs after stitching: " << all_unitigs.size() << "\n";
    }

    // Globally intern color sets and write FASTA + .colors.
    cdgb::ColorSetDict global_dict;
    std::vector<uint32_t> unitig_cid(all_unitigs.size());
    {
        Timer _("emit fasta + colors");
        for (size_t i = 0; i < all_unitigs.size(); ++i) {
            unitig_cid[i] = global_dict.intern(std::move(all_unitigs[i].colors));
        }
        std::cout << "  distinct color classes: " << global_dict.size() << "\n";

        std::vector<std::vector<size_t>> by_class(global_dict.size());
        for (size_t i = 0; i < all_unitigs.size(); ++i) { by_class[unitig_cid[i]].push_back(i); }

        std::ofstream fa(cfg.out_basename + ".fa");
        if (!fa) throw std::runtime_error("cannot open " + cfg.out_basename + ".fa");
        for (uint32_t cid = 0; cid < global_dict.size(); ++cid) {
            for (size_t idx : by_class[cid]) {
                fa << '>' << cid << '\n' << all_unitigs[idx].seq << '\n';
            }
        }
        fa.close();

        cdgb::HybridBuilder hb((uint32_t)files.size());
        for (uint32_t cid = 0; cid < global_dict.size(); ++cid) {
            const auto& cs = global_dict.at(cid);
            hb.encode_color_set(cs.data(), cs.size());
        }
        cdgb::Hybrid h;
        hb.build(h);
        essentials::save(h, (cfg.out_basename + ".colors").c_str());
    }

    // Clean up bucket files. If we created the directory ourselves, remove it
    // entirely; if the user supplied --tmp-dir, only remove the bucket files
    // we wrote (leave the directory and any other contents alone).
    {
        std::error_code ec;
        if (tmp_owned) {
            std::filesystem::remove_all(tmp_dir, ec);
        } else {
            for (uint32_t b = 0; b < num_buckets; ++b) {
                std::filesystem::remove(writer.bucket_path(b), ec);
            }
        }
    }

    std::cout << "done. wrote " << cfg.out_basename << ".fa and " << cfg.out_basename
              << ".colors\n";
    print_total();
    return 0;
}
