// Correctness + scaling test for the stitch phase.
//
// Oracle (generate-by-construction): build K random "true" unitigs,
// each with a distinct color-class id. Split each into overlapping
// fragments (overlap k-1; internal boundaries OPEN, endpoints CLOSED),
// randomly reverse-complement individual fragments, shuffle them all,
// and feed them to stitch. A correct stitch must reconstruct exactly
// the K true unitigs (modulo orientation), each with its original cid
// and both ends closed.
//
// Distinct cids guarantee no cross-unitig gluing. Junction (k-1)-mers
// are random 30-mers, so collisions between distinct unitigs are
// astronomically unlikely (~n^2 / 4^30) and don't occur in practice
// for the sizes here.
//
// Exit code 0 = all pass; non-zero = at least one failure.

#include <atomic>
#include <mutex>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <map>
#include <random>
#include <string>
#include <vector>

#include <filesystem>
#include <unistd.h>

#include "compact_extmem.hpp"  // GGCAT-style id-only compaction (parallel-stitch)
#include "gen.hpp"
#include "stitch.hpp"
#include "stitch_extmem.hpp"
#include "unitig_spill.hpp"  // frag_unitig_writer (+ companion links spill) + reader

namespace {
uint64_t g_tmp_counter = 0;
}

using cdbg::stitchable_unitig;
using cdbg_test::canonical;
using cdbg_test::gen_shared_cid_branch;
using cdbg_test::random_dna;
using cdbg_test::split_unitig;

namespace {

// (canonical-seq, cid) multiset, as a map<key,count>.
using multiset_t = std::map<std::pair<std::string, uint32_t>, uint64_t>;

multiset_t to_multiset(std::vector<stitchable_unitig> const& v) {
    multiset_t m;
    for (auto const& u : v) ++m[{canonical(u.seq), (uint32_t)u.mono_cid()}];
    return m;
}

enum class which_stitch {
    ext_mem,
    ext_file,
    ext_file_mt,
    compact_mem,
    compact_file,
    compact_scalable,
    compact_scalable_links,
    compact_inram
};

char const* stitch_name(which_stitch w) {
    switch (w) {
        case which_stitch::ext_mem:
            return "ext_mem";
        case which_stitch::ext_file:
            return "ext_file";
        case which_stitch::ext_file_mt:
            return "ext_file_mt";
        case which_stitch::compact_mem:
            return "compact_mem";
        case which_stitch::compact_file:
            return "compact_file";
        case which_stitch::compact_scalable:
            return "compact_scalable";
        case which_stitch::compact_scalable_links:
            return "compact_scalable_links";
        default:
            return "compact_inram";
    }
}

// Run the selected stitch implementation on `frags`, collecting output.
void run_stitch(which_stitch w, std::vector<stitchable_unitig>& frags, uint32_t k,
                std::vector<stitchable_unitig>& out) {
    out.clear();
    // ext_file_mt runs the parallel round driver, whose sink is called from
    // multiple threads -- guard the test accumulator with a mutex.
    std::mutex out_mu;
    auto sink = [&](stitchable_unitig&& u) {
        std::lock_guard<std::mutex> lk(out_mu);
        out.push_back(std::move(u));
    };
    if (w == which_stitch::ext_mem) {
        cdbg::vector_frag_source src(frags);
        // Small bucket count exercises multi-round doubling on tiny
        // inputs; 0 would pick the production default.
        cdbg::stitch_unitigs_extmem(src, k, sink, /*num_buckets=*/16);
    } else if (w == which_stitch::compact_mem) {
        // GGCAT-style id-only compaction + fold assembly (parallel-stitch Step 1).
        // Must produce the same unitig multiset as ext_mem.
        cdbg::vector_frag_source src(frags);
        cdbg::compact_stitch_mem(src, k, sink, /*num_buckets=*/16);
    } else if (w == which_stitch::compact_file) {
        // id-only compaction with the file-backed round store (parallel-stitch
        // Step B.1) -- validates the id_tig codec + on-disk doubling.
        cdbg::vector_frag_source src(frags);
        std::string dir = std::filesystem::temp_directory_path().string() +
                          "/cdbg_stitch_test_" + std::to_string(::getpid()) + "_" +
                          std::to_string(g_tmp_counter++);
        std::filesystem::create_directories(dir);
        cdbg::compact_stitch_file(src, k, dir, sink, /*num_buckets=*/16);
        std::filesystem::remove_all(dir);
    } else if (w == which_stitch::compact_scalable) {
        // Full RAM-bounded path (parallel-stitch Step B.2): id-only doubling +
        // disk-based re-bucket assembly. Small frag_ranges/chain_buckets force
        // the multi-bucket join + group/sort logic even on tiny inputs.
        cdbg::vector_frag_source src(frags);
        std::string dir = std::filesystem::temp_directory_path().string() +
                          "/cdbg_stitch_test_" + std::to_string(::getpid()) + "_" +
                          std::to_string(g_tmp_counter++);
        std::filesystem::create_directories(dir);
        auto for_each = [&](auto&& fn) {
            std::vector<cdbg::color_run> runs;
            std::string seq;
            for (uint64_t i = 0; i < src.size(); ++i) {
                std::string_view sv = src.seq_view(i);
                seq.assign(sv.data(), sv.size());
                runs = src.runs(i);
                fn(src.open_flags(i), runs, seq);
            }
        };
        // 4 threads exercises the parallel round driver + thread-safe chain
        // sink; output is a thread-count-independent multiset (compared canonical).
        cdbg::compact_stitch_scalable(for_each, src.size(), k, dir, sink, /*num_buckets=*/16,
                                      /*frag_ranges=*/4, /*chain_buckets=*/4, /*num_threads=*/4);
        std::filesystem::remove_all(dir);
    } else if (w == which_stitch::compact_scalable_links) {
        // Same as compact_scalable, but the seed reads the companion LINKS spill
        // written by the real frag_unitig_writer (the production producer) rather
        // than recomputing boundary k-mers from the bases. Routes the frags
        // through the writer (frag spill + links), then runs the scalable stitch
        // with links_path set. Must match the oracle exactly.
        std::string dir = std::filesystem::temp_directory_path().string() +
                          "/cdbg_stitch_test_" + std::to_string(::getpid()) + "_" +
                          std::to_string(g_tmp_counter++);
        std::filesystem::create_directories(dir);
        std::string spill = dir + "/frag_unitigs.bin";
        cdbg::frag_unitig_writer fw(spill, k);
        for (auto const& u : frags) {
            stitchable_unitig copy = u;  // writer moves the seq out
            fw(std::move(copy));
        }
        fw.close_for_writing();
        const std::string links = fw.links_path();
        const uint64_t n = fw.count();
        auto for_each = [&](auto&& fn) {
            cdbg::frag_unitig_stream_reader rd(spill);
            uint8_t of;
            std::vector<cdbg::color_run> runs;
            std::string seq;
            while (rd.next(of, runs, seq)) fn(of, runs, seq);
        };
        cdbg::compact_stitch_scalable(for_each, n, k, dir, sink, /*num_buckets=*/16,
                                      /*frag_ranges=*/4, /*chain_buckets=*/4, /*num_threads=*/4,
                                      links);
        std::filesystem::remove_all(dir);
    } else if (w == which_stitch::compact_inram) {
        // In-RAM arena path: load fragments into the compact arena, then seed +
        // assemble straight from RAM (disk round store only). Must match oracle.
        cdbg::vector_frag_source src(frags);
        std::string dir = std::filesystem::temp_directory_path().string() +
                          "/cdbg_stitch_test_" + std::to_string(::getpid()) + "_" +
                          std::to_string(g_tmp_counter++);
        std::filesystem::create_directories(dir);
        cdbg::arena_frag_source arena;
        arena.reserve(src.size(), 0, 0);
        for (uint64_t i = 0; i < src.size(); ++i)
            arena.append(src.open_flags(i), src.seq_view(i), src.runs(i));
        cdbg::compact_stitch_inram(arena, k, dir, sink, /*num_buckets=*/16, /*num_threads=*/4);
        std::filesystem::remove_all(dir);
    } else {
        cdbg::vector_frag_source src(frags);
        std::string dir = std::filesystem::temp_directory_path().string() +
                          "/cdbg_stitch_test_" + std::to_string(::getpid()) + "_" +
                          std::to_string(g_tmp_counter++);
        std::filesystem::create_directories(dir);
        // ext_file single-threaded; ext_file_mt with 4 threads exercises the
        // parallel per-round bucket loop + per-output-bucket locking.
        uint32_t threads = (w == which_stitch::ext_file_mt) ? 4u : 1u;
        cdbg::stitch_unitigs_extmem_file(src, k, dir, sink, /*num_buckets=*/16,
                                         /*done=*/nullptr, threads);
        std::filesystem::remove_all(dir);
    }
}

// Run one randomized correctness case against the chosen stitcher.
bool run_case_impl(which_stitch w, uint64_t seed, uint32_t k, uint64_t num_unitigs,
                   uint64_t max_len, uint64_t max_frags) {
    std::mt19937_64 rng(seed);

    std::vector<stitchable_unitig> truth;  // the K true unitigs (closed both ends)
    std::vector<stitchable_unitig> frags;  // shuffled fragments fed to stitch
    truth.reserve(num_unitigs);

    for (uint64_t u = 0; u < num_unitigs; ++u) {
        const uint64_t len = k + (rng() % (max_len - k + 1));
        std::string S = random_dna(len, rng);
        const uint32_t cid = (uint32_t)u;  // distinct -> no cross-unitig gluing
        const uint64_t nf = 1 + (rng() % max_frags);
        auto parts = split_unitig(S, k, cid, nf, rng);
        for (auto& p : parts) frags.push_back(std::move(p));

        stitchable_unitig t;
        t.seq = std::move(S);
        t.set_mono(cid, k);
        t.open_flags = 0;
        truth.push_back(std::move(t));
    }

    std::shuffle(frags.begin(), frags.end(), rng);

    std::vector<stitchable_unitig> out;
    run_stitch(w, frags, k, out);

    // 1) Multiset of (canonical seq, cid) must match the truth exactly.
    multiset_t want = to_multiset(truth);
    multiset_t got = to_multiset(out);
    bool ok = (want == got);

    // 2) Every reconstructed unitig should be fully closed (both true
    //    endpoints were closed).
    uint64_t open_outputs = 0;
    for (auto const& o : out)
        if (o.open_flags != 0) ++open_outputs;

    if (!ok || open_outputs != 0) {
        std::fprintf(stderr,
                     "[FAIL %s] seed=%llu k=%u unitigs=%llu: want %zu distinct, got %zu "
                     "distinct, out_count=%zu open_outputs=%llu\n",
                     stitch_name(w),
                     (unsigned long long)seed, k, (unsigned long long)num_unitigs, want.size(),
                     got.size(), out.size(), (unsigned long long)open_outputs);
        return false;
    }
    return true;
}

// Run a case against BOTH stitchers.
bool run_case(uint64_t seed, uint32_t k, uint64_t num_unitigs, uint64_t max_len,
              uint64_t max_frags) {
    bool b = run_case_impl(which_stitch::ext_mem, seed, k, num_unitigs, max_len, max_frags);
    bool c = run_case_impl(which_stitch::ext_file, seed, k, num_unitigs, max_len, max_frags);
    bool d = run_case_impl(which_stitch::ext_file_mt, seed, k, num_unitigs, max_len, max_frags);
    bool e = run_case_impl(which_stitch::compact_mem, seed, k, num_unitigs, max_len, max_frags);
    bool f = run_case_impl(which_stitch::compact_file, seed, k, num_unitigs, max_len, max_frags);
    bool g = run_case_impl(which_stitch::compact_scalable, seed, k, num_unitigs, max_len, max_frags);
    bool gl =
        run_case_impl(which_stitch::compact_scalable_links, seed, k, num_unitigs, max_len, max_frags);
    bool h = run_case_impl(which_stitch::compact_inram, seed, k, num_unitigs, max_len, max_frags);
    return b && c && d && e && f && g && gl && h;
}

// Shared-cid branchy correctness case. Builds `num_clusters` branch
// clusters; each cluster is 3 true unitigs that share ONE cid and meet at
// a dBG branch (see gen_shared_cid_branch). A correct stitch must
// reproduce every true unitig exactly -- it must NOT glue same-cid
// unitigs across the branch, and must keep different-cid clusters apart.
bool run_branch_case_impl(which_stitch w, uint64_t seed, uint32_t k, uint64_t num_clusters) {
    std::mt19937_64 rng(seed);
    std::vector<stitchable_unitig> truth;  // 3 unitigs per cluster
    std::vector<stitchable_unitig> frags;  // shuffled fragments fed to stitch
    for (uint64_t c = 0; c < num_clusters; ++c) {
        gen_shared_cid_branch(k, (uint32_t)c, rng, frags, truth);
    }
    std::shuffle(frags.begin(), frags.end(), rng);

    std::vector<stitchable_unitig> out;
    run_stitch(w, frags, k, out);

    multiset_t want = to_multiset(truth);
    multiset_t got = to_multiset(out);
    bool ok = (want == got);
    uint64_t open_outputs = 0;
    for (auto const& o : out)
        if (o.open_flags != 0) ++open_outputs;

    if (!ok || open_outputs != 0) {
        std::fprintf(stderr,
                     "[FAIL branch %s] seed=%llu k=%u clusters=%llu: want %zu distinct, got %zu "
                     "distinct, out_count=%zu open_outputs=%llu\n",
                     stitch_name(w), (unsigned long long)seed, k,
                     (unsigned long long)num_clusters, want.size(), got.size(), out.size(),
                     (unsigned long long)open_outputs);
        return false;
    }
    return true;
}

bool run_branch_case(uint64_t seed, uint32_t k, uint64_t num_clusters) {
    bool b = run_branch_case_impl(which_stitch::ext_mem, seed, k, num_clusters);
    bool c = run_branch_case_impl(which_stitch::ext_file, seed, k, num_clusters);
    bool d = run_branch_case_impl(which_stitch::ext_file_mt, seed, k, num_clusters);
    bool e = run_branch_case_impl(which_stitch::compact_mem, seed, k, num_clusters);
    bool f = run_branch_case_impl(which_stitch::compact_file, seed, k, num_clusters);
    bool g = run_branch_case_impl(which_stitch::compact_scalable, seed, k, num_clusters);
    bool gl = run_branch_case_impl(which_stitch::compact_scalable_links, seed, k, num_clusters);
    bool h = run_branch_case_impl(which_stitch::compact_inram, seed, k, num_clusters);
    return b && c && d && e && f && g && gl && h;
}

// Larger case: report timing + fragment/unitig counts, plus the same
// correctness check across the ext_mem and ext_file stitch paths.
bool run_scale(uint64_t seed, uint32_t k, uint64_t num_unitigs) {
    std::mt19937_64 rng(seed);
    std::vector<stitchable_unitig> truth;
    std::vector<stitchable_unitig> frags;
    for (uint64_t u = 0; u < num_unitigs; ++u) {
        const uint64_t len = k + (rng() % 2000);
        std::string S = random_dna(len, rng);
        auto parts = split_unitig(S, k, (uint32_t)u, 1 + (rng() % 16), rng);
        for (auto& p : parts) frags.push_back(std::move(p));
        stitchable_unitig t;
        t.seq = std::move(S);
        t.set_mono((uint32_t)u, k);
        truth.push_back(std::move(t));
    }
    std::shuffle(frags.begin(), frags.end(), rng);
    const uint64_t n_frags = frags.size();
    const multiset_t want = to_multiset(truth);

    std::vector<stitchable_unitig> out;
    bool ok = true;
    for (which_stitch w :
         {which_stitch::ext_mem, which_stitch::ext_file, which_stitch::ext_file_mt,
          which_stitch::compact_mem, which_stitch::compact_file, which_stitch::compact_scalable,
          which_stitch::compact_scalable_links, which_stitch::compact_inram}) {
        auto t0 = std::chrono::steady_clock::now();
        run_stitch(w, frags, k, out);
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        bool this_ok = (want == to_multiset(out));
        ok = ok && this_ok;
        std::fprintf(stderr, "[scale %s] unitigs=%llu frags=%llu out=%zu time=%.0fms %s\n",
                     stitch_name(w),
                     (unsigned long long)num_unitigs, (unsigned long long)n_frags, out.size(),
                     ms, this_ok ? "OK" : "MISMATCH");
    }
    return ok;
}

}  // namespace

int main() {
    const uint32_t k = 31;
    int failures = 0;

    // 50 randomized correctness cases, varied sizes.
    for (uint64_t seed = 1; seed <= 50; ++seed) {
        std::mt19937_64 meta(seed * 2654435761ull);
        const uint64_t num_unitigs = 1 + (meta() % 300);
        const uint64_t max_len = 40 + (meta() % 800);
        const uint64_t max_frags = 1 + (meta() % 12);
        if (!run_case(seed, k, num_unitigs, max_len, max_frags)) ++failures;
    }

    // A couple of fixed edge cases.
    if (!run_case(1001, k, 1, k, 1)) ++failures;       // single k-length unitig, 1 frag
    if (!run_case(1002, k, 1, k + 5, 6)) ++failures;   // 1 unitig, many frags
    if (!run_case(1003, k, 500, 60, 10)) ++failures;   // many short unitigs

    // Shared-cid branchy cases: catch same-cid false joins at dBG
    // branches (the distinct-cid oracle above cannot).
    for (uint64_t seed = 2001; seed <= 2010; ++seed) {
        const uint64_t clusters = 1 + (seed % 8);
        if (!run_branch_case(seed, k, clusters)) ++failures;
    }

    // Medium scale (correctness + timing).
    if (!run_scale(9001, k, 50000)) ++failures;

    if (failures == 0) {
        std::fprintf(stderr, "[test_stitch] ALL PASS\n");
        return 0;
    }
    std::fprintf(stderr, "[test_stitch] %d FAILURE(S)\n", failures);
    return 1;
}
