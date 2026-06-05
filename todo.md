# Optimization backlog (TODO)

Forward-looking optimization and hardening items for `cdbg-builder`. The
algorithm itself is described in `algorithm.md`; this file is only the
"what to make faster / more robust next" list.

Priorities are measured against the validated **661K-genome run**
(`-t 48 -g 64`, 389M color classes), whose per-phase split is:

| phase | wall | % | peak RSS |
|---|---|---|---|
| bucket-write | 24,689 s | 46% | 30.81 GiB |
| bucket-process | 14,795 s | 27% | 31.14 GiB |
| stitch | 12,588 s | 23% | **41.03 GiB** (pipeline peak) |
| emit | 1,949 s | 4% | (+0) |
| **total** | **54,067 s (~15 h)** | | **41.03 GiB < 64** |

Correctness gate for *every* change: `salmonella_10` (exact ground-truth
unitig-set equality via `test_data/verify.py`) + `test_stitch`.

---

## 0. Infra first (makes everything below cheap to measure)

- **`--stop-after <phase>`** + a **~20K-genome subset** of the collection
  (warm page cache, fixed `-t`/`-g`) → a ~10-min optimize/measure loop
  instead of the 15 h full run.
- **Split the bucket-write `spill` timer** into `{sort+unique, compress,
  write}`. It currently bundles CPU and HDD I/O, so "spill is all I/O" is
  an unconfirmed assumption.

---

## 1. Container pre-sizing audit (missed `reserve` / realloc)

A pass over every growing map/vector found several that are **not**
pre-sized — they rehash (maps) or realloc-copy (vectors) repeatedly as
they fill, *and* run at a chronically high load factor (longer probe
chains → more cache misses per op).

> **⚠️ Measured caveat — a naive full `reserve` is net-negative on a
> *walked* map.** See the `claude/reserve-presize` study (shelved, not
> merged). Reserving `kmer_info` to ~distinct backfired: `ankerl`'s
> `reserve(n)` rounds the bucket array up to the next power of two, so the
> table ends up ~2× sparser. The find-heavy walk then strides a sparser
> array — more cache misses, and at high thread counts that becomes
> memory-bandwidth contention across co-resident buckets. Same-binary A/B
> on the 20K subset at `-b 12 -t 32` (`--reserve-frac` flip):
>
> | reserve | load | walk | peak RSS | phase |
> |---|---|---|---|---|
> | off (`0`) | 17.2s | **23.5s** | 2.95 GiB | 94.5s |
> | on (`0.6`) | 16.9s | **35.1s** (+49%) | 3.26 GiB | ~100s |
>
> No load win (the per-bucket maps at `-b 12` are too small for growth to
> dominate); a real walk regression that would *worsen* at 661K's 48-way.
> A reserve only paid off at `-b 10` (maps big enough that repeated early
> doublings dominate load, *and* throttled parallelism muted the walk
> contention) — not a default regime. Lesson: the reserve must target the
> *frequent-early-doubling* cost without over-sizing into the
> sparse-table-walk penalty — i.e. a **threshold/capped** reserve, not a
> reserve-to-distinct. See §2.1.

| container | location | grows to | effect |
|---|---|---|---|
| `kmer_info` | `load_bucket` | millions / bucket | bucket-process `load` (but **only big buckets** — see caveat above) |
| global dict `m_classes` (`vector<hash_pair>`) | `streaming_color_set_dict` | ~389M × 16B ≈ 6.2 GB | realloc-copies a multi-GB vector |
| global dict `m_index` (set) | `streaming_color_set_dict` | ~389M | rehashes ~log₂N times |
| local dict `m_index`/`m_classes` | `compact_color_set_dict` | ~tens of K / bucket × 40K | bucket-process `resolve` |
| compactor `m_dedup` | `bucket_io` | ~`spill_bytes`/fill | minor — `clear()` keeps capacity, only the *first* fill rehashes |

(Already reserved, fine: `cid_of`, `visited`, stitch `waiting`,
`m_compactors`.)

**Key connection — this is the same bug as `merge_wait`.** The global
dict's `m_classes.push_back` + `m_index.insert` both run **under
`global_mu`** (the per-bucket→global merge, algorithm.md §4.7). With no
reserve, the last few vector reallocations copy 3–6 GB and the set
rehashes ~389M entries *while holding the lock*, stalling all `T` threads
→ `merge_wait` = 2,726 s (~18% of bucket-process at 661K). Pre-sizing
attacks it at the root, likely cheaper than lock-sharding.

---

## 2. bucket-process (27% of wall)

1. **Threshold-`reserve` `kmer_info` (revised — see §1 caveat).** The naive
   reserve-to-distinct is *net-negative* (slows the walk ~49% via
   pow2-over-sizing + bandwidth contention). The win we actually want is
   bounded: a big 661K bucket grows `0→1→…→N` and pays ~log₂N rehashes,
   each copying the whole table — the early doublings are frequent and the
   late ones copy the most. Goal: skip those **without** over-sizing into
   the sparse-walk penalty. Candidate shapes to try (20K A/B at both `-b 10`
   and `-b 12`, `-t 32`, then confirm direction at 661K `-t 48`):
   - **Floor-only:** `reserve(max(estimate, FLOOR))` with a modest `FLOOR`
     (e.g. 1–4M) — small buckets keep the cheap path; only genuinely large
     buckets get pre-sized, where the load win exists.
   - **Cap the reserve below the pow2 cliff:** size to a target *just under*
     a power-of-two boundary so the rounded array isn't 2× the live set
     (keeps load factor high → walk stays dense). Needs knowing `ankerl`'s
     rounding rule and reserving `≈ 0.9 · 2^j`.
   - **Reserve as a fraction of the *gate's* admitted budget**, not of the
     bucket's own distinct count, so RAM-per-bucket stays bounded.

   Key: measure `walk` *and* `load` *and* peak RSS on every variant — the
   reserve trades all three, and the 20K `-b 12` result proves `load`-only
   measurement is misleading.
2. **Pre-size / de-reallocate the global dict** — segmented `m_classes`
   (so it never realloc-copies) + reserved `m_index`. Cashes out
   `merge_wait` at its root (see §1).
3. **Instrument + cut `resolve`** (33.7s ≈ **36%** of bucket-process at
   20K — the biggest single sub-cost now that the reserve is shelved). It
   braids three costs that we currently can't separate: `record_sets.at`
   color-list **decode**, `local_dict.intern` (**hash + maybe global-dict
   insert**), and the rsid→cid **cache** lookups. *Split-timer first*
   (decode / intern / cache), then attack the dominant one. Note the
   per-bucket `compact_color_set_dict` maps here are **insert-heavy, not
   walked**, so the §1 sparse-walk reserve hazard does *not* apply —
   reserving them should be a clean win (and the same is true of the global
   dict in §2.2).
4. **Sort-for-locality in `load`** — sort each bucket's records by
   canonical k-mer so `out[can]` goes in ascending order (sequential
   memory access, fewer cache misses). Deeper; attacks the residual probe
   floor `reserve` can't remove.
5. **Lock-shard the global merge** — *only if* (2) doesn't already clear
   `merge_wait`.
6. `kmer_info` is the largest per-bucket RSS residual (no easy compression
   — it's actively used during the walk). A streaming-online walker (vs
   load-then-walk) would remove it, but it's a big restructure.

---

## 3. stitch (23% of wall — and the RAM peak)

1. **Live-RSS ceiling / dict-aware budget.** Stitch is the only phase with
   no live-RSS enforcement, and at 661K it was the **pipeline peak**
   (41.03 GiB — the parallel round loop holds `num_threads`=48 decoded
   round-buckets resident at once). It's the phase most likely to bust
   `-g` at a tighter budget or more threads; give it the §4.9 treatment
   (or reserve for the carried ~11 GiB color dict).
2. **Throughput: round 0 dominates** — 6,543 s of 7,700 s is rounds 0–4,
   i.e. the seed pass over 6.08B fragments. Revisit its bucketing /
   parallelism.

---

## 4. bucket-write (46% of wall — biggest single phase)

1. **Is LZ4 worth it?** Only ~1.16× here (432 vs 503 GiB) — the
   2-bit-packed bases + varint-delta colors are near-incompressible. LZ4
   CPU is negligible (~25 s), so it is *not* the slowness, but it adds
   framing complexity for a ~14% I/O saving. Add **`--no-compress`** (raw
   frames) and A/B it on HDD vs. the complexity.
2. **`spill` dominates** (3,914 s of 4,962 s flush). It bundles
   `sort+unique` (CPU — heavy for dense, high-multiplicity color lists) +
   LZ4 + `fwrite` (HDD seeks across 40K files). Split the timer (infra §0)
   before optimizing: if **write**-bound → fewer/bigger files or an
   append-only layout; if **sort**-bound → a smarter color merge.
3. *I/O-bound regime* (cold cache / HDD): prefetch input files in a
   producer pool while compute threads work decompressed buffers; sort the
   file list by inode/disk-position to cut HDD seeks.
4. *CPU-bound regime* (warm cache / SSD): residual hot lines are
   `flush.hashmap` (`m_dedup` find) and `compute` (ntHash + minimizer
   queue). A different flat hashmap (`folly::F14`, `phmap`) is a modest
   lever. (A hand-SIMD ACGT/2-bit branch was tried and dropped — gcc -O3
   already vectorized those loops. A GGCAT-style background compactor was
   tried on `claude/append-only-bucketing` and rejected — it inflated
   per-record encoding ~6× because records must be varint-encoded + 2-bit
   packed even when they'd dedup away.)
5. (minor) reserve compactor `m_dedup` at construction (first-fill only).

---

## 5. emit (4% of wall — effectively done)

Single-threaded by design: `emit_fasta` is a sequential `.fa` write,
`emit_colors` a `finalize()` of the already-streamed file. The old 117 GiB
peak is gone (external merge-sort cap + `release_index`, validated at 389M
classes). Parallelizing the `emit_fasta` cid-bucket loop is possible
(buckets independent; only the final append + u2c ordering need care) but
the payoff is negligible today.

---

## 6. Not on the critical path (known lever, deliberately unbuilt)

**Bucket resplit for fat-tailed minimizer buckets.** The bucket-process
working set is the *sum* of `num_threads` co-resident buckets, not the
single largest, and the measured distribution is **broad** (max ~2.6×
median), not fat-tailed — so the §4.9 admission gate is the right tool;
splitting the largest buckets would not reduce the co-resident sum. Only
useful for a future input with a genuinely fat-tailed bucket.
