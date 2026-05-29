#pragma once

// Synthetic pangenome FASTA generator for the per-phase tests.
//
// Produces a set of related "genomes" (one file = one color) that share
// a common backbone with per-genome point mutations and optional
// accessory segments. High backbone sharing is what drives the dedup
// and color-set machinery, so this faithfully exercises the real
// bucket-write / bucket-process code paths -- a set of unrelated random
// genomes would dedup to almost nothing and miss the interesting
// behavior.
//
// Files are written as plain (uncompressed) FASTA; seq_reader handles
// plain and .gz transparently, and plain keeps the generator free of a
// gzip-compression dependency (our libdeflate build is decompress-only).
//
// All counts/indices are uint64_t. Determinism: everything is driven by
// a single seed.

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace cdgb_test {

struct pangenome_params {
    uint64_t num_genomes = 10;     // = number of colors
    uint64_t backbone_len = 50000; // shared core length in bases
    double mutation_rate = 0.01;   // per-base substitution probability per genome
    uint64_t accessory_len = 0;    // extra random bases appended per genome (0 = none)
    uint64_t seed = 1;
    uint64_t line_width = 80;      // FASTA wrap width
};

// Convenience size tiers. Tune freely; "large" is meant to be run with
// a real RAM cap to observe scaling, not in CI.
inline pangenome_params tier_small(uint64_t seed = 1) {
    return {10, 50000, 0.01, 2000, seed, 80};
}
inline pangenome_params tier_medium(uint64_t seed = 1) {
    return {1000, 200000, 0.02, 20000, seed, 80};
}
inline pangenome_params tier_large(uint64_t seed = 1) {
    return {50000, 2000000, 0.03, 100000, seed, 80};
}

namespace detail {

inline char rand_base(std::mt19937_64& rng) {
    static const char b[4] = {'A', 'C', 'G', 'T'};
    return b[rng() & 3u];
}

inline void write_fasta(std::string const& path, std::string const& header,
                        std::string const& seq, uint64_t line_width) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) throw std::runtime_error("cannot open synth FASTA for write: " + path);
    std::fputc('>', f);
    std::fwrite(header.data(), 1, header.size(), f);
    std::fputc('\n', f);
    for (uint64_t off = 0; off < seq.size(); off += line_width) {
        uint64_t n = std::min<uint64_t>(line_width, seq.size() - off);
        std::fwrite(seq.data() + off, 1, n, f);
        std::fputc('\n', f);
    }
    std::fclose(f);
}

}  // namespace detail

// Generate `params.num_genomes` FASTA files under `dir` (created if
// needed). Returns the list of file paths in color order (file i =
// color i). Also writes `<dir>/filenames.txt` (one path per line) so a
// test can hand it straight to the builder.
//
// genome[i] = backbone with each base independently substituted with
// probability mutation_rate, then `accessory_len` random bases appended
// (the accessory tail is unique per genome, producing color-1
// super-k-mers; the mutated backbone produces the shared, multi-color
// super-k-mers that dominate dedup).
inline std::vector<std::string> generate_pangenome(pangenome_params const& params,
                                                   std::string const& dir) {
    std::filesystem::create_directories(dir);
    std::mt19937_64 rng(params.seed);

    // Shared backbone.
    std::string backbone(params.backbone_len, 'A');
    for (uint64_t i = 0; i < params.backbone_len; ++i) backbone[i] = detail::rand_base(rng);

    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<std::string> paths;
    paths.reserve(params.num_genomes);

    for (uint64_t g = 0; g < params.num_genomes; ++g) {
        std::string seq = backbone;
        for (uint64_t i = 0; i < seq.size(); ++i) {
            if (unit(rng) < params.mutation_rate) {
                // Substitute with a different base (so it's a real SNP).
                char c = seq[i];
                char n;
                do {
                    n = detail::rand_base(rng);
                } while (n == c);
                seq[i] = n;
            }
        }
        for (uint64_t i = 0; i < params.accessory_len; ++i) seq.push_back(detail::rand_base(rng));

        std::string path = dir + "/genome_" + std::to_string(g) + ".fasta";
        detail::write_fasta(path, "genome_" + std::to_string(g), seq, params.line_width);
        paths.push_back(path);
    }

    std::string list = dir + "/filenames.txt";
    std::FILE* lf = std::fopen(list.c_str(), "wb");
    if (!lf) throw std::runtime_error("cannot open filenames list for write: " + list);
    for (auto const& p : paths) {
        std::fwrite(p.data(), 1, p.size(), lf);
        std::fputc('\n', lf);
    }
    std::fclose(lf);

    return paths;
}

}  // namespace cdgb_test
