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
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "gen.hpp"
#include "stitch.hpp"
#include "stitch_extmem.hpp"

using cdgb::stitchable_unitig;
using cdgb_test::canonical;
using cdgb_test::random_dna;
using cdgb_test::split_unitig;

namespace {

// (canonical-seq, cid) multiset, as a map<key,count>.
using multiset_t = std::map<std::pair<std::string, uint32_t>, uint64_t>;

multiset_t to_multiset(std::vector<stitchable_unitig> const& v) {
    multiset_t m;
    for (auto const& u : v) ++m[{canonical(u.seq), u.cid}];
    return m;
}

enum class which_stitch { in_ram, ext_mem };

// Run the selected stitch implementation on `frags`, collecting output.
void run_stitch(which_stitch w, std::vector<stitchable_unitig>& frags, uint32_t k,
                std::vector<stitchable_unitig>& out) {
    out.clear();
    if (w == which_stitch::in_ram) {
        cdgb::stitch_unitigs(frags, k, out);
    } else {
        auto sink = [&](stitchable_unitig&& u) { out.push_back(std::move(u)); };
        cdgb::vector_frag_source src(frags);
        // Small bucket count exercises multi-round doubling on tiny
        // inputs; 0 would pick the production default.
        cdgb::stitch_unitigs_extmem(src, k, sink, /*num_buckets=*/16);
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
        t.cid = cid;
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
                     w == which_stitch::in_ram ? "in_ram" : "ext_mem",
                     (unsigned long long)seed, k, (unsigned long long)num_unitigs, want.size(),
                     got.size(), out.size(), (unsigned long long)open_outputs);
        return false;
    }
    return true;
}

// Run a case against BOTH stitchers.
bool run_case(uint64_t seed, uint32_t k, uint64_t num_unitigs, uint64_t max_len,
              uint64_t max_frags) {
    bool a = run_case_impl(which_stitch::in_ram, seed, k, num_unitigs, max_len, max_frags);
    bool b = run_case_impl(which_stitch::ext_mem, seed, k, num_unitigs, max_len, max_frags);
    return a && b;
}

// Larger case: report timing + fragment/unitig counts, plus the same
// correctness check. Single-thread in-RAM stitch; this is the baseline
// the external-memory redesign must match.
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
        t.cid = (uint32_t)u;
        truth.push_back(std::move(t));
    }
    std::shuffle(frags.begin(), frags.end(), rng);
    const uint64_t n_frags = frags.size();
    const multiset_t want = to_multiset(truth);

    std::vector<stitchable_unitig> out;
    bool ok = true;
    for (which_stitch w : {which_stitch::in_ram, which_stitch::ext_mem}) {
        auto t0 = std::chrono::steady_clock::now();
        run_stitch(w, frags, k, out);
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        bool this_ok = (want == to_multiset(out));
        ok = ok && this_ok;
        std::fprintf(stderr, "[scale %s] unitigs=%llu frags=%llu out=%zu time=%.0fms %s\n",
                     w == which_stitch::in_ram ? "in_ram" : "ext_mem",
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

    // Medium scale (correctness + timing).
    if (!run_scale(9001, k, 50000)) ++failures;

    if (failures == 0) {
        std::fprintf(stderr, "[test_stitch] ALL PASS\n");
        return 0;
    }
    std::fprintf(stderr, "[test_stitch] %d FAILURE(S)\n", failures);
    return 1;
}
