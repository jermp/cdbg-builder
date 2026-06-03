// Emit-phase microbenchmark.
//
// Reproduces the builder's EMIT phase (emit_fasta + emit_colors) at arbitrary
// scale WITHOUT the ~14 h ingest/process/stitch upstream, so we can measure the
// emit peak RSS at 661k proportions (389M color classes, 1.31B unitigs,
// 1.48e12 total color integers) in minutes.
//
// It synthesizes the exact two inputs emit consumes, using the REAL classes:
//   1. streaming_color_set_dict: intern `num_classes` synthetic color sets
//      whose sizes sum to ~`total_integers` (so m_classes, m_index, the
//      sidecar offsets, and the finalize() EF build all match production).
//   2. unitig_bucket_writer: push `num_unitigs` synthetic unitigs whose seq
//      lengths sum to ~`total_bases`, each tagged with a cid in [0,num_classes)
//      (so read_bucket + the per-bucket in-RAM sort match production).
// Then builder::run_emit_only() runs the real emit_fasta + emit_colors with the
// same phase_rss_marker instrumentation.
//
// Usage:
//   emit-bench --out PFX -g RAM_GB \
//     --num-colors C --num-classes NC --total-integers TI \
//     --num-unitigs NU --total-bases TB \
//     [--unitig-buckets K] [--threads T] [--seed S]
//
// All synthetic content is random but size-faithful; emit RAM depends on the
// COUNTS and SIZES, not the actual color/base values, so this is representative.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "builder.hpp"
#include "streaming_color_set_dict.hpp"
#include "unitig_spill.hpp"

using namespace cdbg;

namespace {

uint64_t arg_u64(int argc, char** argv, char const* key, uint64_t dflt) {
    for (int i = 1; i + 1 < argc; ++i)
        if (std::strcmp(argv[i], key) == 0) return std::strtoull(argv[i + 1], nullptr, 10);
    return dflt;
}
double arg_f64(int argc, char** argv, char const* key, double dflt) {
    for (int i = 1; i + 1 < argc; ++i)
        if (std::strcmp(argv[i], key) == 0) return std::strtod(argv[i + 1], nullptr);
    return dflt;
}
std::string arg_str(int argc, char** argv, char const* key, std::string dflt) {
    for (int i = 1; i + 1 < argc; ++i)
        if (std::strcmp(argv[i], key) == 0) return argv[i + 1];
    return dflt;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string out = arg_str(argc, argv, "--out", "/tmp/emit_bench");
    const double ram_gb = arg_f64(argc, argv, "-g", 0.0);
    const uint32_t num_colors = (uint32_t)arg_u64(argc, argv, "--num-colors", 661405);
    const uint64_t num_classes = arg_u64(argc, argv, "--num-classes", 1000000);
    const uint64_t total_integers = arg_u64(argc, argv, "--total-integers", 0);
    const uint64_t num_unitigs = arg_u64(argc, argv, "--num-unitigs", 1000000);
    const uint64_t total_bases = arg_u64(argc, argv, "--total-bases", 0);
    const uint32_t unitig_buckets = (uint32_t)arg_u64(argc, argv, "--unitig-buckets", 0);
    const uint32_t threads = (uint32_t)arg_u64(argc, argv, "--threads", 1);
    const uint64_t seed = arg_u64(argc, argv, "--seed", 42);
    const uint32_t k = (uint32_t)arg_u64(argc, argv, "-k", 31);
    // cid skew: 0 = uniform; >0 biases unitigs toward LOW cids (small color
    // classes / common k-mers), as in real pangenomes. cid = floor(num_classes
    // * u^(1+skew)), u uniform in [0,1). skew=1 ~ quadratic bias; higher =
    // fatter low-cid head. Reproduces the emit-fasta read_bucket skew that
    // uniform cids miss (661k: ~9x fattest-vs-average).
    const double cid_skew = arg_f64(argc, argv, "--cid-skew", 0.0);

    // Mean per-class color-set size (>=1). If total_integers unset, use 1.
    const double mean_class_size =
        total_integers > 0 ? (double)total_integers / (double)num_classes : 1.0;
    // Mean per-unitig seq length (>= k). If total_bases unset, use 2k.
    const double mean_unitig_len =
        total_bases > 0 ? (double)total_bases / (double)num_unitigs : (double)(2 * k);

    std::printf(
        "emit-bench: num_colors=%u num_classes=%llu total_integers=%llu "
        "(mean class %.1f)\n            num_unitigs=%llu total_bases=%llu (mean len %.1f) "
        "-g %.2f threads=%u\n",
        num_colors, (unsigned long long)num_classes, (unsigned long long)total_integers,
        mean_class_size, (unsigned long long)num_unitigs, (unsigned long long)total_bases,
        mean_unitig_len, ram_gb, threads);

    build_config cfg;
    cfg.out_basename = out;
    cfg.k = k;
    cfg.num_threads = threads;
    cfg.max_ram_gb = ram_gb;
    builder b(cfg);

    std::string tmp_dir = out + ".tmp";
    std::filesystem::create_directories(tmp_dir);

    // ---- 1) Populate the streaming color-set dict (drives emit_colors) ----
    {
        streaming_color_set_dict dict(num_colors, out + ".color_sets");
        std::mt19937_64 rng(seed);
        // Class sizes ~ geometric around the mean, clamped to [1, num_colors].
        // Distinct content per class: a sorted random subset of colors. We make
        // each class's content unique (salt by id) so dedup never collapses
        // them -- we want exactly num_classes distinct interned classes.
        std::vector<uint32_t> colors;
        colors.reserve(num_colors);
        std::geometric_distribution<uint64_t> gsize(
            mean_class_size > 1.0 ? 1.0 / mean_class_size : 0.9999);
        // Each class's content is a per-class random subset of [0,num_colors),
        // seeded by the class id, so distinct classes get distinct content
        // (collision over the run is ~num_classes^2 / 2^(bits of content),
        // negligible). This guarantees dict.size() == num_classes -- a faithful
        // benchmark needs exactly the requested class count, not a deduped
        // subset. The first color is forced to (c % num_colors) and, for sz==1,
        // that single value rotates so even singletons stay distinct up to
        // num_colors classes; larger classes are distinct by their random body.
        for (uint64_t c = 0; c < num_classes; ++c) {
            uint64_t sz = gsize(rng) + 1;
            if (sz > num_colors) sz = num_colors;
            colors.clear();
            std::mt19937_64 crng(seed * 1000003ull + c);  // per-class stream
            if (sz >= num_colors) {
                for (uint32_t x = 0; x < num_colors; ++x) colors.push_back(x);
            } else {
                // Sample `sz` distinct colors via a small hash set.
                ankerl::unordered_dense::set<uint32_t> pick;
                pick.reserve(sz);
                pick.insert((uint32_t)(c % num_colors));  // anchor for singleton distinctness
                std::uniform_int_distribution<uint32_t> ud(0, num_colors - 1);
                while ((uint64_t)pick.size() < sz) pick.insert(ud(crng));
                colors.assign(pick.begin(), pick.end());
                std::sort(colors.begin(), colors.end());
            }
            std::vector<uint32_t> cp = colors;
            dict.intern(std::move(cp));
            if (((c + 1) % 20000000ull) == 0)
                std::printf("  interned %llu / %llu classes\n", (unsigned long long)(c + 1),
                            (unsigned long long)num_classes);
        }
        std::printf("  dict: %llu classes (requested %llu), %llu total integers, %llu bits\n",
                    (unsigned long long)dict.size(), (unsigned long long)num_classes,
                    (unsigned long long)dict.total_integers(),
                    (unsigned long long)dict.total_bits());
        if (dict.size() != num_classes)
            std::printf("  WARNING: interned class count != requested (dedup collisions); "
                        "emit-colors RAM will be sized to the smaller count\n");

        // ---- 2) Populate the unitig bucket writer (drives emit_fasta) ----
        // emit_fasta's read_bucket loads ONE cid-bucket's records into RAM and
        // sorts them, so its peak scales with records-per-bucket =
        // num_unitigs / K. K MUST match what build() picks, or the measured
        // emit-fasta RAM is off by the K ratio. Default (--unitig-buckets 0):
        // derive K exactly as the builder does from total_bases and -g.
        uint32_t ub = unitig_buckets;
        if (ub == 0)
            ub = builder::pick_unitig_bucket_count(num_classes, total_bases, ram_gb);
        std::printf("  unitig buckets K = %u (records/bucket ~= %.1f M)\n", ub,
                    (double)num_unitigs / ub / 1e6);
        unitig_bucket_writer uwriter(tmp_dir, num_classes, ub);
        std::mt19937_64 urng(seed ^ 0x9e3779b97f4a7c15ull);
        std::geometric_distribution<uint64_t> glen(
            mean_unitig_len > (double)k ? 1.0 / (mean_unitig_len - (double)k + 1.0) : 0.5);
        std::uniform_int_distribution<uint64_t> cid_pick(0, num_classes - 1);
        std::uniform_real_distribution<double> u01(0.0, 1.0);
        const double skew_exp = 1.0 + cid_skew;  // cid = num_classes * u^skew_exp
        const char* B = "ACGT";
        std::string seq;
        for (uint64_t u = 0; u < num_unitigs; ++u) {
            uint64_t len = (uint64_t)k + glen(urng);
            seq.resize(len);
            uint64_t bits = urng();
            for (uint64_t i = 0; i < len; ++i) {
                if ((i & 31) == 0) bits = urng();
                seq[i] = B[(bits >> ((i & 31) * 2)) & 3];
            }
            uint64_t cid;
            if (cid_skew <= 0.0) {
                cid = cid_pick(urng);
            } else {
                double uu = u01(urng);
                cid = (uint64_t)((double)num_classes * std::pow(uu, skew_exp));
                if (cid >= num_classes) cid = num_classes - 1;
            }
            stitchable_unitig su;
            su.seq.swap(seq);
            su.set_mono(cid, k);
            uwriter(std::move(su));
            seq.clear();
            if (((u + 1) % 100000000ull) == 0)
                std::printf("  pushed %llu / %llu unitigs\n", (unsigned long long)(u + 1),
                            (unsigned long long)num_unitigs);
        }
        std::printf("  uwriter: %llu unitigs in %u buckets\n",
                    (unsigned long long)uwriter.total_unitigs(), uwriter.num_buckets());

        // ---- 3) Run the REAL emit phase, measuring RSS per sub-phase ----
        b.run_emit_only(uwriter, dict, num_unitigs);
    }

    std::filesystem::remove_all(tmp_dir);
    return 0;
}
