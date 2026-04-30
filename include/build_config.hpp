#pragma once

#include <cstdint>
#include <string>

namespace cdgb {

struct BuildConfig {
    std::string filenames_list;     // text file: one input path per line
    std::string out_basename;       // produces <basename>.fa and <basename>.colors
    uint32_t k = 31;
    uint32_t num_threads = 1;
    uint32_t shard_log2 = 8;        // 2^shard_log2 hash shards (in-memory path)
    bool verbose = false;

    // Minimizer-bucketed disk-based ingest path (--bucketed).
    bool bucketed = false;
    uint32_t m = 0;                 // 0 = auto (compute_best_m(k))
    uint32_t bucket_log2 = 10;      // 2^10 = 1024 buckets
    std::string tmp_dir;            // empty -> mkdtemp under $TMPDIR
};

}  // namespace cdgb
