#include <exception>
#include <iostream>
#include <string>

#include <parser.hpp>

#include "cdbg_builder.hpp"
#include "kmer.hpp"
#include "util.hpp"

namespace {

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
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    cdgb::build_config cfg;
    if (!parse_args(argc, argv, cfg)) return 1;

    try {
        cdgb::cdbg_builder builder(cfg);
        builder.build();
    } catch (std::exception const& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
