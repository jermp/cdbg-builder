# Branch B — deterministic memory-bounded stitch (§3.9 + §3.9.1): implementation map

Working design for porting GGCAT v2's deterministic min-index merge into cdbg's
phase 3, replacing the randomized iterative-doubling stitch. Grounded in the
current cdbg stitch (`include/phase3_stitch/stitch_extmem.hpp`) and the GGCAT v2
source at commit `65749e68cc30e1fde2377d3b8d61033b8a885eb9`
(`crates/assembler_pipeline/src/extend_unitigs.rs`,
`crates/sequence_output/src/sequences_joiner.rs`,
`crates/io/src/partial_unitigs_extra_data.rs`).

Not part of the repo build; a reference for the implementer.

---

## 0. Baseline to beat (main, this branch's fork point)

verify.py green both datasets:
- salmonella_10: 86,630 unitigs / 171 color sets; stitch 28 rounds, 2.5 s, peak 326.6 MiB.
- se-50: 274,241 unitigs / 9,540 color sets; stitch 29 rounds, 5.5 s, peak 387.2 MiB.

Real-scale poles (661k, from README): stitch is `+1.87 GiB` (sets the 36.52 GiB
overall peak) and ~3,400 s of doubling rounds within a 7,524 s phase.

---

## 1. What is being replaced, and what is reused

`stitch_unitigs_extmem_file_streaming` (stitch_extmem.hpp:1120) =
`round_store_file` + `ext_seed_round0_raw_par` (random-side seed) +
`ext_run_rounds` (RNG doubling: each round processes ALL buckets in parallel,
each tig presents a RANDOM open end routed to `H(kmer) % B`, joins co-presented
compatible pairs, re-routes survivors; O(log L) expected rounds).

REUSE unchanged (these are correct and orthogonal to the merge policy):
- `ext_tig` {seq, runs, open_flags, rng}, `ext_tig_serialize/deserialize`.
- `side_kmer_canonical(seq,k,side,&fwd)` → canonical boundary k-mer + is_fwd.
- `ext_pair_compatible(a,b)` — join-legality (side + orientation).
- `ext_join(a,a_side,b,b_side,k)` — orient, drop shared boundary k-mer, concat
  seq + `ext_concat_runs` (KEEP the all-foreign-bridge reconciliation).
- `ext_split_monochromatic` — cut a finished tig at color-run boundaries at sink.
- `ext_end`, `SIDE_LEFT/SIDE_RIGHT`.

REPLACE:
- routing policy (random side + `H(kmer)%B`) → **min-index** (present the open
  end whose bucket is smaller; deterministic tie-break when equal).
- driver (`ext_run_rounds` doubling) → **smallest-index-first sweep** with
  in-bucket sub-rounds and forward-deferral.

Gate behind `build_config::deterministic_stitch` (`--deterministic-stitch`),
selected in builder.hpp where `stitch_unitigs_extmem_file_streaming` is called.

---

## 2. The deterministic merge (§3.9), cdbg form

Each partial unitig has ≤2 open extremity k-mers; each maps to a bucket via
`(H(kmer) >> 1) % B`. Two fragments that must glue share their FULL boundary
k-mer, so they map to the SAME bucket for that end (as today — full-k-mer key).

Routing (`find_best_orientation_and_flags`, extend_unitigs.rs:111-267):
- One open end → present it; bucket = that end's bucket.
- Both open → present the end with the SMALLER bucket; route there.
  - If both ends' buckets are EQUAL → **presence-map tie-break** (§2.1).
- Record the presented side on the tig (repurpose `ext_tig::rng`, unused here:
  store SIDE_LEFT/SIDE_RIGHT). RC/`last_align` handling is unnecessary in cdbg
  (join re-derives orientation from `side_kmer_canonical`).

Sweep (extend_unitigs.rs outer loop @410 + ProcessSubpartition @596-925):
- Process buckets **b = 0..B-1 in order**.
- For each b, sub-round loop until b is exhausted:
  - drain b's current tigs;
  - local pairing hashmap keyed by presented boundary k-mer:
    first arrival waits; second compatible arrival → `ext_join`;
  - joined+completed (both ends sealed) → sink (split monochromatic);
  - joined+open → route via min-index: min-bucket == b → **re-queue into b**
    (next sub-round); min-bucket > b → **defer** into that later bucket;
  - drained unmatched waiters → circular (both extremity k-mers equal) → emit
    circular (trim one base); else re-queue into b.
- **No RNG, no doubling.**

Termination / min-index invariant: joining consumes the presented (min) end, so
both remaining ends have bucket ≥ b; the merged tig routes to ≥ b (re-queue or
forward). Forward progress on b + shrinking per-bucket tig count ⇒ terminates.
Every open end has exactly one partner (degree-1 cross-bucket edge from phase 2),
and that partner's min-bucket is ≤ this end's bucket, so by the time bucket b is
processed all of b's partners are present (initial + deferrals from b' < b).

### 2.1 Presence-map tie-break (REQUIRED for stage 1 determinism)

When a both-open tig has bucket_L == bucket_R (== b), which end it presents is
ambiguous — and two partners sharing end X must present X *together* or they
never pair (random re-roll would fix this but reintroduces RNG; the whole point
is determinism). GGCAT resolves it with per-(sub)partition presence maps
(extend_unitigs.rs:181-214): a `HashMap<kmer, bool>` per bucket recording, for
an inserted extremity, whether the inserter kept it. The partner consults it and
matches. Port (subpartitions_count = 1 ⇒ one presence map per bucket):
- lock/consult `presence[b]` for both extremity hashes;
- if X already present with flag f → change_extremity = !f (match the inserter);
  else if Y present with flag g → change_extremity = g;
- insert X→!change_extremity and Y→change_extremity.
- Also carry `is_circular = (hash_L == hash_R)`.
This is the deterministic analogue of the doubling RNG and CANNOT be dropped in
stage 1.

---

## 3. Staging (each step: verify.py + monochromatic decode green on BOTH datasets)

**Stage 1 — deterministic core, correctness prototype.**
- New `include/phase3_stitch/stitch_deterministic.hpp`: min-index routing +
  presence-map tie-break + smallest-first sweep, ONE hashmap per bucket
  (subpartitions_count = 1), sequences inline, **single-threaded**, in-memory
  bucket vectors (`std::vector<std::vector<ext_tig>>`) for minimal bug surface.
- Flag-gate in builder; `--deterministic-stitch`.
- Expect: green on both; SLOW at scale (no parallelism yet) — correctness only,
  NOT a measurement point. Determinism: byte-identical `.fa` across two runs.

**Stage 2 — parallelism + memory bound.**
- Split each outer bucket into subpartitions (config `MAX_SUBPARTITIONS_COUNT =
  2^11`, `MAX_SUBSUBPARTITION_SIZE = 128 KiB`); thread-parallel across a bucket's
  subpartitions (GGCAT parallelizes WITHIN the outer bucket, not across — outer
  buckets stay sequential). Presence maps become per-subpartition; saturation
  cap `MAX_EXTREMITIES_HASHMAP_SIZE = 524288` → fall back to a coin flip (the ONE
  place a tie can go non-deterministic under memory pressure; acceptable, matches
  GGCAT). Re-tune the outer/subpartition counts against cdbg `-g` (not v2's
  file-size heuristic).
- Spill bucket/subpartition files to disk with LZ4 framing (reuse
  ext_tig_serialize + round_store_file's frame I/O) so only one subpartition's
  hashmap + one bucket's records are resident — the structural memory bound that
  lets the RAM-governor backstop be retired for phase 3.
- This is where a fair wall-time vs. the doubling stitch is measurable.

**Stage 3 — indirect long unitigs (§3.9.1), sequences_joiner.rs.**
- Threshold `T = max(4k, 4096)` (`MAX_INLINE_UNITIG_SIZE = 4096`). A partial
  unitig keeps inline prefix+suffix; interior beyond T spills to a shared
  append-only `oversize-unitigs.dat`, referenced by
  `IndirectReadInfo{file_offset, extra_length, sequence_length<<1|rc_bit}`
  (partial_unitigs_extra_data.rs:105-140).
- Join cases (sequences_joiner.rs `append_sequence`/`splice_join_operation`):
  inline+inline promote-on-overflow; incoming-indirect append the ref list.
  **On RC: reverse the ref order AND flip each ref's rc bit AND transform
  `indirection_start`** — all three, or reconstruction corrupts (the fiddliest
  part; unit-test all join cases with RC on/off BEFORE wiring in).
- Full seq + colors read back from the oversize file only at final emit.
- Confirm `-g` holds on se-50 under multithread.

**Stage 4 — measure phase-3 peak RSS + wall vs the doubling stitch; retire RNG
path behind the flag.**

---

## 4. Risks / open questions

1. **Parallelism.** §3.9 is smallest-first SEQUENTIAL at the outer level; the
   parallelism cdbg's doubling gets for free (all buckets per round × 32-48
   threads) is recovered only WITHIN an outer bucket via subpartitions (stage 2).
   Net wall-time vs. the doubling is the open question — could be a wash. The
   guaranteed wins are determinism + the structural memory bound (retire the RAM
   governor for phase 3). Measure before claiming a speedup.
2. **Tie-break correctness** (§2.1): if two partners fail to co-present, the join
   stalls and unitigs fragment. Test on se-50 and high-branching inputs.
3. **Indirect RC handling** (stage 3): threefold transform; unit-test first.
4. **Boundary color reconciliation**: KEEP `ext_concat_runs`'s all-foreign-bridge
   fix across the new join path; re-run the exact se-50 case.
5. **Circular unitigs**: emitted at the drain step (both extremity k-mers equal),
   not in a doubling loop — ensure cdbg's circular handling maps onto that.

---

## 5. Correctness gate

`test_data/run_both.sh <dir> '<glob>' 31 <t> --deterministic-stitch` → verify.py
== naive ground truth, on salmonella_10 AND se-50, after every stage. Plus a
determinism check: two runs produce byte-identical `.fa` (sort-free).
