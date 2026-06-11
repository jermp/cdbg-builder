# Externalizing the color-sets-dedup-map (design)

Status: **stages 1–2 merged to `main` (encapsulate + measure-only); stages 3–6
NOT built.** The dedup index is still fully in RAM and can overflow `-g` at
scale. See the Staging section for exactly what's done vs planned.

## Problem

Phase-2 (bucket-process) dedups color sets online: each distinct color set gets
a dense **cid**, assigned the first time its 128-bit content hash is seen. The
index that does this — `streaming_color_set_dict::m_shards`, a 256-way sharded
`hash128 → cid` map (the **color-sets-dedup-map**) — is held **entirely in RAM**
and grows ~linearly with the number of distinct color sets:

- ≈ 0.5 GiB @ 20k genomes (≈14.0 M classes)
- ≈ 14 GiB @ 661k genomes

It is the **dominant, fastest-growing, non-spillable** structure in the whole
pipeline (see `algorithm.md` §9.1) and the reason 661k needed `-g 64`. It is
freed (`release_index`) right after bucket-process, so it never reaches stitch;
the problem is purely its **peak during bucket-process**.

## The golden rule (the constraint)

> `-g` must be honored **at all costs**, BUT use the available RAM as much as
> possible — spill to disk **as infrequently as possible**.

So the target is **RAM-first / spill-overflow**, exactly like the stitch stores
(`algorithm.md` §5.5): while the index fits a `-g`-derived budget, behave
**identically to today** (exact online dedup, zero spill, full speed); only when
it would exceed the budget does the overflow go to disk.

## Why this is not a one-liner

The cid is assigned **online** and is **final** the instant a fragment is
written to the frag spill (`bucket_walker.hpp:737` → `r.cid = local_to_global`).
Exact online dedup with *final* cids needs one of:

1. all seen hashes resident in RAM — today's approach, unbounded; or
2. random-access lookups into an on-disk hash — **death on HDD** (millions of
   ~10 ms seeks); violates "spill infrequently"; or
3. **deferred dedup + remap** — assign provisional ids past the budget, dedup
   the overflow at the end, and remap provisional → final.

(1) is what we're replacing, (2) is disqualified by the HDD scratch, so the
design is (3) **layered under** (1): RAM-first exact dedup, with a deferred
overflow path that only engages past the budget.

## Design: RAM-first exact dedup + hash-partitioned overflow + folded remap

### Normal regime (index fits the budget) — unchanged

Carve a dedup budget `D` from `-g` (a share of the bucket-process budget). While
`index.resident_bytes() ≤ D`, `intern_with_hashes` is **exactly today's code**:
probe the in-RAM `hash128 → cid` index, return on hit, else encode + assign the
next cid + insert. Zero spill, zero remap, no format change. Small/medium builds
(the overwhelming common case) never leave this regime.

### Overflow regime (index hit the budget)

The in-RAM index is **frozen** at capacity `M` (cids `[0, M)` are final and stay
resident — every pre-budget class is still deduped against exactly). A candidate
whose hash misses the frozen index is either a genuinely new class or a
duplicate of another *post-budget* class. For these:

- Assign a **provisional id** from a separate space `prov ∈ [0, P)` (a monotone
  counter, distinct from final cids).
- Append `(hash128, prov, encoded_bits)` to one of `R` **hash-partitioned spill
  files**, partition = `hash >> s mod R`. Sequential append → no random I/O.
- Return `prov` as the fragment's cid (so fragments carry provisional ids).

The frozen index still serves all pre-budget classes exactly, so dedup quality
degrades only among the (smaller) tail of post-budget classes — and is made
**exact again** at reconciliation.

### Reconciliation (after bucket-process, before stitch)

For each partition `r` (sized ≤ total_overflow / R, so `R` is chosen to keep one
partition within a RAM share — same "scale the fan-out" trick as `frag_ranges`):

1. Load partition `r`, dedup by `hash128` in a RAM map. First occurrence of each
   distinct hash gets a **final cid** continuing after `M` (so final cids stay
   dense), and its `encoded_bits` are appended to the `.color_sets` stream +
   offset recorded — identical to the normal append path.
2. Emit `(prov → final)` pairs for every provisional id in the partition.

Collect all `(prov → final)` pairs. They form the **overflow remap**, dense over
`prov ∈ [0, P)`. Because pre-budget cids were already final (`prov` only exists
past `M`), the remap is identity for `cid < M` and a lookup for `prov ≥ M`.

### Applying the remap — fold into stitch's frag-spill read (no extra I/O)

Stitch already **streams the frag spill** (seed + attach, `algorithm.md` §5.1/
§5.3). The color runs carry cids; during that read, remap each `cid ≥ M` through
the overflow remap before it flows into assembly. So the remap costs **no extra
pass** — it rides the read stitch already does. (The remap table is dense over
`[M, M+P)`; if it itself exceeds a RAM share it is the one thing that can be
mmap'd / range-loaded, since the access during the frag read is arbitrary — but
`P` ≤ post-budget distinct classes, far smaller than the frag count.)

## Correctness

- **Exactness:** every color set ends with exactly one final cid. Pre-budget:
  deduped by the frozen index. Post-budget: deduped by the per-partition RAM
  map. Identical hashes always collide in the same partition (partition is a
  function of the hash), so no duplicate survives. ⇒ the `.color_sets` file has
  exactly `num_distinct` classes (the 171 / 13.98 M invariants hold).
- **Output identity in the common case:** while under budget, byte-for-byte
  identical to today (same code path, same cid order).
- **cid order at scale:** final cids are `[0, M)` in first-seen order, then the
  partitions in partition order — a *different* permutation than today's pure
  first-seen, but cids are opaque ids, so `.fa`/`.u2c`/`.color_sets` remain
  mutually consistent and verify.py still passes (it checks the *set*, not the
  numbering).

## Honoring `-g`

- The in-RAM index is **capped at `D`** (never grows past it) — hard bound.
- Overflow spill files are sequential, bounded by disk.
- Each reconciliation partition is bounded by choosing `R` so one partition ≤ a
  RAM share (scaled to the measured overflow volume, like `frag_ranges`).
- The remap table is the last resident piece; bounded by `P` and range-loadable.

So the dedup peak is **bounded by `D` regardless of dataset size**, and a build
that fits stays 100 % in RAM at full speed.

## Staging (each stage builds + passes `test_stitch` and the 86630/171 ground truth)

**Status: stages 1–2 merged to `main`; stages 3–6 NOT built.** The dedup index is
still fully in RAM and unbounded — it can still overflow `-g` (≈14 GiB @ 661k).
Stages 1–2 only encapsulate and *measure*; nothing spills yet.

1. ✅ **Encapsulate** the dedup index behind a small `colorset_dedup_index` type
   (sharded `hash128→cid`) — pure refactor, no behavior change. *(the clean seam,
   merged)*
2. ✅ **Budget plumbing (MEASURE-ONLY):** carve `D` from `-g` in the builder
   (`COLORSET_DEDUP_BUDGET_FRAC = 0.50`), pass it down, log resident-vs-budget
   (the `color-sets-dedup-map: ~X / Y budget (within → no spill)` line). **No
   spill** — it only reports whether the not-yet-built overflow path *would*
   engage. (merged)
3. ☐ **Overflow spill + provisional ids** in the dict (partitioned append, frozen
   index past `D`). *(this is where actual externalization begins — not started)*
4. ☐ **Reconciliation** (per-partition dedup → final cids + `.color_sets` append +
   `prov→final`).
5. ☐ **Remap fold** into the stitch frag-spill read.
6. ☐ Bench 100k / 661k at a tight `-g` to confirm the dedup peak is bounded and
   spill is rare under a comfortable `-g`.

## Open questions / knobs

- Dedup budget share `D` of `-g` (coexists with the per-bucket walk set under
  the §4.9 admission gate — the two must jointly fit `0.82·g`).
- Partition count `R` and the remap-table residency at extreme `P`.
- Whether to keep the 128-bit hash or fold in a byte-equality check at
  reconciliation (we discard colors after encoding today, so we'd compare
  encoded bits instead — cheap within a partition).
