#include <exception>
#include <iostream>
#include <string>

#include <parser.hpp>

#include "builder.hpp"
#include "kmer.hpp"
#include "util.hpp"

namespace {

bool parse_args(int argc, char** argv, cdbg::build_config& cfg) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("filenames_list",
               "Text file with one input path per line. The file at line i has color i.", "-i",
               true);
    parser.add("out_basename",
               "Output basename. Produces <basename>.fa, <basename>.u2c, and "
               "<basename>.color_sets.",
               "-o", true);
    parser.add("k", "K-mer length (must be <= " + std::to_string(cdbg::MAX_K) + ").", "-k", true);
    parser.add("num_threads", "Number of worker threads (default 1).", "-t", false);
    parser.add("m", "Minimizer length (default: auto, derived from k).", "-m", false);
    parser.add("buckets_log2",
               "log2 of the bucket count (override; default: auto -- the bucket COUNT is"
               " derived from the RAM model when -g is set, else 1024).",
               "-b", false);
    parser.add("tmp_dir", "Scratch directory for bucket files (default: mkdtemp under $TMPDIR).",
               "-d", false);
    parser.add("max_ram_gb",
               "Soft RAM budget in GiB. The builder sizes the bucket count from the model"
               " B = 0.50*g / (alpha*T*flush + beta*spill) so bucket-write fits, and spills"
               " the color bit_vector to disk; the actual peak RSS is reported at the end"
               " (no hard kill if exceeded).",
               "-g", false);
    parser.add("alpha",
               "RAM-model per-thread buffer overhead multiplier (default 2.0). Tunable knob"
               " for re-calibrating the bucket-count model per machine/allocator.",
               "--alpha", false);
    parser.add("beta",
               "RAM-model per-bucket compactor overhead multiplier (default 7.0). Tunable"
               " knob for re-calibrating the bucket-count model per machine/allocator.",
               "--beta", false);
    parser.add("flush_bases",
               "Per-thread->compactor handoff size in bytes (default 4096). Smaller ->"
               " larger derived bucket count B -> smaller, faster buckets.",
               "--flush", false);
    parser.add("spill_bytes",
               "Compactor dedup window in bytes before a disk frame (default 65536)."
               " Smaller -> larger B but weaker dedup / bigger bucket files.",
               "--spill", false);
    parser.add("verbose", "Verbose output.", "--verbose", false, true);
    parser.add("resume_stitch",
               "BENCHMARK/RECOVERY: skip bucket-write+process and run stitch (+emit_fasta) from an "
               "existing -d <tmp_dir>/frag_unitigs.bin left by a killed run. Times the stitch; does "
               "NOT finalize <out>.color_sets. Requires -d; -i is ignored.",
               "--resume-stitch", false, true);
    parser.add("keep_tmp", "Do not delete -d <tmp_dir> on success (keeps frag spill for --resume-stitch).",
               "--keep-tmp", false, true);
    parser.add("base_stitch",
               "Use the base-carrying stitch (assembles unitig sequences in the doubling "
               "rounds): faster on small/medium inputs that fit RAM, but it moves ALL bases "
               "through every round -- multi-TB disk and a hard-to-bound peak at large scale. "
               "Default is the id-only stitch, which HARD-enforces -g and is disk-frugal.",
               "--base-stitch", false, true);

    if (!parser.parse()) return false;

    cfg.filenames_list = parser.get<std::string>("filenames_list");
    cfg.out_basename = parser.get<std::string>("out_basename");
    cfg.k = parser.get<uint32_t>("k");
    if (parser.parsed("num_threads")) cfg.num_threads = parser.get<uint32_t>("num_threads");
    if (parser.parsed("m")) cfg.m = parser.get<uint32_t>("m");
    if (parser.parsed("buckets_log2")) cfg.bucket_log2 = parser.get<uint32_t>("buckets_log2");
    if (parser.parsed("tmp_dir")) cfg.tmp_dir = parser.get<std::string>("tmp_dir");
    if (parser.parsed("max_ram_gb")) cfg.max_ram_gb = parser.get<double>("max_ram_gb");
    if (parser.parsed("alpha")) cfg.alpha = parser.get<double>("alpha");
    if (parser.parsed("beta")) cfg.beta = parser.get<double>("beta");
    if (parser.parsed("flush_bases")) cfg.flush_bases = parser.get<uint64_t>("flush_bases");
    if (parser.parsed("spill_bytes")) cfg.spill_bytes = parser.get<uint64_t>("spill_bytes");
    if (parser.parsed("verbose")) cfg.verbose = parser.get<bool>("verbose");
    if (parser.parsed("resume_stitch")) cfg.resume_stitch = parser.get<bool>("resume_stitch");
    if (parser.parsed("keep_tmp")) cfg.keep_tmp = parser.get<bool>("keep_tmp");
    if (parser.parsed("base_stitch")) cfg.base_stitch = parser.get<bool>("base_stitch");
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    cdbg::build_config cfg;
    if (!parse_args(argc, argv, cfg)) return 1;

    try {
        cdbg::builder builder(cfg);
        builder.build();
    } catch (std::exception const& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
