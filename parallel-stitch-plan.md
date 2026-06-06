# parallel-stitch: GGCAT-style id-only compaction + single base/color assembly

## Why (the finding)

Our `stitch_extmem` is faithful to GGCAT's `links_compaction` in every structural
respect (hash-bucketed iterative doubling, full-boundary-k-mer keying, colorless
join, structural branch separation, terminate-on-zero-joins, per-bucket parallel).
**The one decisive difference:** GGCAT doubles *id-only* records (segment-id list +
orientation flags, no bases, no colors) and assembles bases/colors exactly once in a
final linear pass (`reorganize_reads` -> `build_unitigs`). We carry the full `seq` +
color `runs` inside every `ext_tig` and LZ4-recompress the growing sequences into the
round store *every round*. That is the whole ~5 s (GGCAT) vs ~72 s (us) gap, plus our
seed reads all 92.5 M fragments' bases (33 s) just to compute boundary k-mers.

Item counts are comparable (our 92.5 M fragments ≈ their boundary-edge links, same
minimizer bucketing). The win is removing bases/colors from the loop, not parallelism.

## Target pipeline

```
bucket-process                 (mostly unchanged; emits 2 streams now)
  ├─ reads spill   : per fragment {frag_id, bases, runs}      (= today's frag spill, now id-addressable)
  └─ links spill   : per OPEN end {boundary_kmer, frag_id, side, canon}  -> hash-bucketed

compaction  (NEW, id-only)     replaces today's seed+rounds
  - seed   : read links, one id-tig per fragment, route by H(boundary_kmer)
  - rounds : doubling on id-tigs (clone of ext_choose_side/ext_pair_compatible/ext_join
             but payload = entries:[(frag_id,rc)] + left/right junction, NO bases/runs)
  - output : per CLOSED chain, the ordered (frag_id, rc) list

assembly + emit  (NEW)         replaces today's emit-fasta read-back
  - gather each chain's fragment bases+runs (see "assembly strategy" below)
  - concat (drop k-1 overlap, RC per flag), join runs, split monochromatic -> cid bucket
  - emit-fasta: cid-sort + write .fa/.u2c (today's back end, parallelizable)
```

## What we reuse verbatim (base-agnostic already)
- `ext_choose_side`, `ext_pair_compatible`, `ext_mix_bit` — operate on open_flags / side /
  canonical / rng only.
- The round driver (`ext_run_rounds`): persistent pool, per-round barrier, route/pair/
  terminate. Templated on the tig payload, so it adapts.
- `round_store_file` (LZ4-framed buckets) — now stores tiny id-tigs.

## What changes
- `ext_tig` -> `id_tig { uint64_t jl, jr; uint8_t canon, open_flags; uint64_t rng;
  std::vector<uint64_t> entries; uint8_t pres; }` (entry = frag_id<<1 | rc).
- `ext_join` -> id-only: concat/reverse `entries` (flip rc on reverse), recompute jl/jr from
  the parents' outer ends. No seq/runs.
- rng seed: from `frag_id` (was from seq) — just needs to vary per tig.
- boundary k-mer (`side_kmer_canonical`) computed in **bucket-process** (it has the bases),
  written into the link record, NOT recomputed in the loop.
- `ext_split_monochromatic` moves to the **assembly** stage.

## Assembly strategy (the one real design fork — RAM scalability)
Fetching each chain's fragment bases by id without an O(num_frags) RAM index (the 137 GB
index the extmem design exists to avoid). GGCAT does this with `reorganize_reads`
(re-bucket read sequences by final-unitig, then assemble per-bucket).

**Plan: stage it.**
- **Step A (correctness):** simplest possible assembler — for TEST sizes, an in-RAM
  `frag_id -> (offset|bases)` map — to validate the id-only compaction produces exactly the
  same unitig set as today's stitch (oracle = existing `stitch_unitigs_extmem` in
  `test_stitch`). Proves the hard/novel part (id doubling + orientation + color join order)
  before optimizing.
- **Step B (scalability):** replace with GGCAT-style sort/bucket reorganize — bounded RAM,
  no per-fragment resident index:
  1. compaction emits per chain member `(frag_id, chain_id, pos, rc)`, bucketed by frag_id range
  2. sequential pass over reads spill (frag_id order) merge-joins to attach bases+runs ->
     emit `(chain_id, pos, bases, runs)` bucketed by chain_id range
  3. per chain-bucket: sort by (chain_id,pos), assemble, split monochromatic -> cid bucket
  Each base byte is moved ~2-3x total (vs ~several× through today's rounds).
- **Step C:** parallelize assembly + emit-fasta (independent across buckets), integrate into
  `builder.hpp`, benchmark vs current.

## Rollout / test anchors (keep today's stitch as the reference the whole time)
1. New header `compact_extmem.hpp`; id-only structs + adapted doubling. Unit-test the
   doubling in isolation (chains of known fragments -> expected grouping/orientation).
2. Step-A assembler; wire a `which_stitch::compact` case into `test_stitch` comparing the
   assembled unitig multiset against `ext_mem`. Must match on all randomized + branch cases.
3. Step-B scalable reorganize; same test must still pass; add a RAM-bound assertion.
4. Step-C parallel + builder integration; bw20k benchmark.

## Expected payoff
seed 33 s -> ~3 s (links only); rounds 38 s -> ~8 s (id-only); + assembly ~20 s (was folded
into the 38 s + 16 s emit). Rough target: stitch+emit ~88 s -> ~40 s, i.e. total ~250 s ->
~205 s, near/below GGCAT, with assembly+emit parallelism as further headroom.

## Risk
Largest change to the most correctness-sensitive code. Mitigation: keep the existing stitch
intact as the oracle; land id-only compaction behind a test harness first (Step A) and only
swap the production path after the scalable assembler (Step B) matches output byte-for-byte.
