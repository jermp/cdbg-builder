#pragma once

#include <cstdint>
#include <string>

namespace cdgb {

struct BuildConfig {
    std::string filenames_list;     // text file: one input path per line
    std::string out_basename;       // produces <basename>.fa and <basename>.colors
    uint32_t k = 31;
    uint32_t num_threads = 1;
    uint32_t shard_log2 = 8;        // 2^shard_log2 hash shards
    bool verbose = false;
};

}  // namespace cdgb
