# bucket-process RAM cap via outlier resplit — design & plan

Branch: `claude/bucket-process-ramcap` (off `claude/stitch-fully-external`).

## STATUS UPDATE after the full 661k run (-g 64, 48 threads)

The completed 661k run changed the picture. Recorded peaks:

| phase | peak RSS | increment | vs 64 GiB |
|---|---|---|---|
| bucket-write | 30.48 GiB | +30.40 | ok |
| bucket-process | 81.29 GiB | +50.80 | **OVER** |
| stitch | 81.29 GiB | **+0** | ok ✅ (streaming stitch validated at 6.08e9 frags) |
| **emit** | **117.32 GiB** | **+36.03** | **OVER — the actual peak** |

Two corrections to the earlier scorecard:
1. **The stitch is fully validated**: +0 over prior peak at 6.08 billion
   fragments. The streaming-frag-reader work did its job at full scale.
2. **emit is the REAL peak (117 GiB), not bucket-process (81 GiB).** emit
   had been +0 at 4546 genomes so it was assumed capped, but it had never
   been stressed at 389 MILLION color classes. So there are TWO uncapped
   stages, and **emit is the larger violator** — it must be addressed
   first (it sets the 117 GiB peak; fixing bucket-process alone would only
   drop the peak from 117 to ... still 117, because emit dominates).

### emit: what allocates +36 GiB (NOT yet fully explained — instrument first)

emit = emit_fasta + emit_colors. Static sizing from the run's counters
(389e6 classes, 1.48e12 total ints, 1.31e9 unitigs) accounts for only
~17–21 GiB:
  - emit_colors dedup hashtable: 16 B/class metadata + ankerl table
    overhead ≈ **~13 GiB** (this persists from bucket-process, not new).
  - emit_fasta read_bucket: one unitig bucket's records in RAM + sort
    copy ≈ **~4–8 GiB** (MAX_K=1024 cap on emit buckets means a bucket can
    still be large at this scale).
  - EF build at finalize ≈ <1 GiB.
There is a ~15 GiB gap I cannot explain from static reasoning, so the
split `emit-fasta` / `emit-colors` RSS markers were added (this commit)
to localize it on the next run. DO NOT design the fix until that run says
which half owns the +36 GiB.

Likely fixes once localized (to be confirmed):
  - emit_fasta: lower the MAX_K cap so a unitig bucket fits the budget,
    OR stream-sort instead of in-RAM sort.
  - emit_colors: the dedup hashtable is O(num_classes) and unavoidable in
    its current form; may need to drop/spill it before finalize (it isn't
    needed once interning is done), or shrink the per-class metadata.

### bucket-process: broad, not outlier-dominated (pending distribution)

+50.80 GiB / 48 threads ≈ 1.06 GiB per in-flight bucket; 187.6e9 kmers /
8192 buckets ≈ 22.9M kmers/bucket average → ~1 GiB average kmer_info.
This points to the overshoot being BROAD (concurrency × average), so
outlier-only resplit may be necessary-but-not-sufficient; a concurrency
cap and/or more primary buckets may be needed too. The
`report_bucket_size_distribution` instrumentation (committed) will
confirm outlier-vs-broad on the next run.

### Revised priority order
1. **emit** (sets the 117 GiB peak) — instrument (done), then fix.
2. **bucket-process** (81 GiB) — distribution, then resplit ± concurrency.

---

(Original resplit design for bucket-process below; still valid for that
phase once emit is handled.)

---


Goal: make **bucket-process** honor the `-g` budget. It is the last
uncapped phase. On the Blackwell 661k pangenome (`-g 64`) it peaked
~75.9 GiB (~12 GiB over) because a few dense minimizer buckets hold an
O(bucket) `kmer_info` map that doesn't fit a per-thread share of `-g`.
Everything else (bucket-write, stitch, emit) is already capped.

This does NOT change output. Resplit only re-partitions whole
super-k-mer records across more, smaller bucket files; the cross-bucket
stitch reconnects the new boundaries exactly as it does the primary
minimizer boundaries. Validated by `verify.py` (exact ground-truth set
equality) staying green on salmonella-10.

---

## 1. Mechanism (GGCAT kmers-transform, adapted)

A bucket file is an **outlier** if its on-disk size exceeds a
budget-derived cap. An outlier is re-bucketed into `M` sub-buckets by a
**fresh, longer secondary minimizer** `m2 = max(2, k-10)` (GGCAT's
`k-10`). A longer window discriminates among the super-k-mers that all
shared one short primary minimizer, spreading them across the M files.

```
M = clamp(next_pow2(ceil(size / cap)), 4, 1024)      # GGCAT clamp
sub(h) = (h >> 1) & (M-1)                            # primary low-bit convention
```

Sub-bucket files reuse the LZ4-framed bucket-file format, so
`bucket_reader` / `load_bucket` / `process_bucket` consume them
unchanged — process_bucket doesn't even know it's reading a sub-bucket.

## 2. THE compatibility constraint (why the old `bucket_resplit.hpp` can't be reused as-is)

The prior attempt (`claude/bucket-process-resplit`, commit c330f25) was
written against the **old (k−1)-overlap, no-ownership** model. Its
`rewindow_super_kmer` used:
  - `W = k - m2 + 1`  → **k-mer** minimizers (WRONG now)
  - no `SK_FLAG_OWNS_FIRST/LAST`

The current pipeline (this branch's base) uses BCALM2 **(k−1)-mer**
minimizers + **boundary ownership flags**. Resplit MUST mirror the
current `emit_super_kmers` exactly, only swapping `m`→`m2` and
`bucket`→`sub`:
  - `W = k - m2`  (minimize over (k−1)-mers)
  - run-based super-k-mer split (run of equal-minimizer (k−1)-mers →
    super spanning k-mers `[run_a-1 .. b]`)
  - **ownership**: at an m2-minimizer split (mA=cur_min, mB=new_min),
    ending super owns its last k-mer iff `mA < mB`; starting super owns
    its first iff `mB < mA`. Run-end super owns its last (ACGT terminus).
  - BEGIN/END: only the record's true BEGIN/END k-mers keep those flags;
    internal m2 boundaries are open ends (no BEGIN/END) the stitch joins.

**Correctness invariant (the over-fragmentation trap):** routing must be
a deterministic function of each *k-mer's* (k−1)-junction minimizer, NOT
the record's overall minimizer. Identical k-mers in different records
must get identical m2-minimizers → same sub-bucket, or a shared boundary
k-mer is emitted into two sub-buckets and duplicated. Per-k-mer
re-windowing guarantees this. (The doc's note that the old branch
"over-fragments" is most likely the ownership flags being dropped, so the
walker mis-colored boundary k-mers — the exact class of bug we fixed in
the stitch-fully-external work.)

## 3. Where it hooks in (builder bucket-process loop)

`bucket_walker.hpp::process_buckets`, per bucket `b` pulled from the
atomic counter:

```
path = writer.bucket_path(b)
if file_size(path) > cap:
    M = resplit_factor(size, cap)
    paths = resplit_bucket(path, prefix_b, k, M)        # writes M sub files
    for sp in paths: process_bucket(sp, ...)            # same as a normal bucket
    delete sub files
else:
    process_bucket(path, ...)                            # unchanged
```

Notes:
  - Resplit runs inside the worker thread that owns bucket `b`, so it is
    naturally parallel and adds no new threads.
  - `cap` is derived from the per-thread share of `-g`:
    `cap ≈ SHARE * eff_ram / num_threads / OVERHEAD`, where OVERHEAD maps
    on-disk LZ4 bytes → in-RAM `kmer_info` footprint (needs calibration;
    start conservative, measure on a known-dense bucket).
  - Recursive resplit: if a sub-bucket is *still* an outlier (one
    minimizer dominates even under m2), recurse with a different
    hash salt or fall through with a warning. GGCAT recurses; v1 here can
    cap recursion depth at 1–2 and warn, since the 661k overshoot is only
    ~18% (one level of M≈4–8 should suffice).

## 4. Implementation steps (one change → build → test → commit each)

1. Port `bucket_resplit.hpp` but REWRITE `rewindow_super_kmer` to mirror
   the current `emit_super_kmers` ((k−1) minimizers + ownership). Keep
   `resplit_writer`, `resplit_factor` (sound as-is).
2. Add a unit/e2e assertion: resplitting a bucket then processing the sub
   files yields the SAME fragments (multiset) as processing the bucket
   whole. (Force resplit on salmonella-10 with a tiny cap.)
3. Wire into `process_buckets` behind a cap check; add `pick_resplit_cap_()`
   in builder.hpp from `-g`/threads.
4. e2e on salmonella-10 with a cap small enough to FORCE resplit →
   `verify.py` must still report `86630 == ground truth`. This is the
   gate: identical output with resplit active.
5. Measure peak RSS vs cap on salmonella-10 under a tight `-g` to confirm
   bucket-process now tracks budget.
6. (If a large input is available) confirm the 661k-class overshoot is
   gone.

## 5. Guardrails (carried from the stitch work)

- ONE change, build, run BOTH `test_stitch` and the e2e, read numbers,
  commit. Do not batch.
- The e2e (`verify.py` exact ground-truth equality) is the real
  correctness gate — it catches the over-fragmentation / boundary
  mis-coloring this feature is most likely to introduce.
- Resplit is output-invariant by design; any change in the 86630 count
  means a boundary/ownership bug.
- All loop vars / counts `uint64_t`.
