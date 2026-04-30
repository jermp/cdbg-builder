#include <atomic>
#include <chrono>
#include <fstream>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#include <essentials.hpp>

#include "build_config.hpp"
#include "color_set_dict.hpp"
#include "concurrent_kmer_map.hpp"
#include "finalize.hpp"
#include "hybrid_color_sets.hpp"
#include "ingester.hpp"
#include "progress.hpp"
#include "unitig_walker.hpp"

namespace {

void print_usage(const char* argv0) {
    std::cerr <<
        "cdgb-build: build a colored compacted de Bruijn graph and emit\n"
        "  - <out>.fa      FASTA of monochromatic colored unitigs (header = color_set_id)\n"
        "  - <out>.colors  binary in Fulgor's `hybrid` color-set format\n"
        "\n"
        "usage:\n"
        "  " << argv0 << " -i <filenames_list> -k <k> -o <out_basename> [-t <num_threads>]\n"
        "\n"
        "<filenames_list> is a text file containing one input path per line.\n"
        "Each input file is one color, in line order (file at line i has color i).\n"
        "Inputs may be FASTA, FASTQ, or gzipped variants of either.\n";
}

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
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const char* opt) {
            if (i + 1 >= argc) throw std::runtime_error(std::string("missing value for ") + opt);
            return std::string(argv[++i]);
        };
        if (a == "-i") cfg.filenames_list = need("-i");
        else if (a == "-k") cfg.k = (uint32_t)std::stoul(need("-k"));
        else if (a == "-o") cfg.out_basename = need("-o");
        else if (a == "-t") cfg.num_threads = (uint32_t)std::stoul(need("-t"));
        else if (a == "--shards") cfg.shard_log2 = (uint32_t)std::stoul(need("--shards"));
        else if (a == "-v" || a == "--verbose") cfg.verbose = true;
        else if (a == "-h" || a == "--help") return false;
        else { std::cerr << "unknown arg: " << a << "\n"; return false; }
    }
    if (cfg.filenames_list.empty() || cfg.out_basename.empty() || cfg.k == 0) return false;
    if (cfg.k > cdgb::MAX_K) throw std::runtime_error("k must be <= 63");
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
    try {
        if (!parse_args(argc, argv, cfg)) { print_usage(argv[0]); return 1; }
    } catch (std::exception& e) {
        std::cerr << "argument error: " << e.what() << "\n";
        print_usage(argv[0]);
        return 1;
    }

    auto files = read_filenames(cfg.filenames_list);
    if (files.empty()) { std::cerr << "no input files\n"; return 1; }
    if (files.size() > (uint64_t)UINT32_MAX) {
        std::cerr << "too many colors (max 2^32 - 1)\n"; return 1;
    }
    std::cout << "k = " << cfg.k << ", num_colors = " << files.size()
              << ", num_threads = " << cfg.num_threads << "\n";

    const uint64_t num_shards = uint64_t(1) << cfg.shard_log2;

    cdgb::ConcurrentKmerMap raw_map(cfg.shard_log2);
    {
        Timer _("ingest");
        std::atomic<uint64_t> done{0};
        cdgb::Progress prog("ingest", done, files.size());
        cdgb::ingest_parallel(files, cfg.k, raw_map, cfg.num_threads, &done);
        prog.stop();
        std::cout << "  k-mers ingested: " << raw_map.num_kmers() << "\n";
    }

    cdgb::FinalKmerMap fkm(cfg.shard_log2);
    cdgb::ColorSetDict global_dict;
    {
        Timer _("finalize");
        std::atomic<uint64_t> done{0};
        cdgb::Progress prog("finalize", done, num_shards);
        cdgb::finalize(raw_map, fkm, global_dict, cfg.num_threads, &done);
        prog.stop();
        std::cout << "  distinct color classes: " << global_dict.size() << "\n";
        std::cout << "  k-mers retained: " << fkm.size() << "\n";
    }

    std::vector<cdgb::Unitig> unitigs;
    {
        Timer _("unitig walk");
        std::mutex mu;
        cdgb::UnitigWalker w(fkm, cfg.k);
        std::atomic<uint64_t> done{0};
        cdgb::Progress prog("unitig walk", done, num_shards);
        w.walk_all(cfg.num_threads, unitigs, mu, &done);
        prog.stop();
        std::cout << "  unitigs: " << unitigs.size() << "\n";
    }

    // Group unitigs by color_set_id so the FASTA is written in color-set order.
    {
        Timer _("emit fasta + colors");
        std::vector<std::vector<size_t>> by_class(global_dict.size());
        for (size_t i = 0; i < unitigs.size(); ++i) {
            by_class[unitigs[i].color_set_id].push_back(i);
        }

        // FASTA: header = color_set_id, body = unitig sequence.
        std::ofstream fa(cfg.out_basename + ".fa");
        if (!fa) throw std::runtime_error("cannot open " + cfg.out_basename + ".fa");
        for (uint32_t cid = 0; cid < global_dict.size(); ++cid) {
            for (size_t idx : by_class[cid]) {
                fa << '>' << cid << '\n' << unitigs[idx].seq << '\n';
            }
        }
        fa.close();

        // Hybrid color sets in canonical order (id 0, 1, ...).
        cdgb::HybridBuilder hb((uint32_t)files.size());
        for (uint32_t cid = 0; cid < global_dict.size(); ++cid) {
            const auto& cs = global_dict.at(cid);
            hb.encode_color_set(cs.data(), cs.size());
        }
        cdgb::Hybrid h;
        hb.build(h);
        essentials::save(h, (cfg.out_basename + ".colors").c_str());
    }

    std::cout << "done. wrote " << cfg.out_basename << ".fa and "
              << cfg.out_basename << ".colors\n";
    return 0;
}
