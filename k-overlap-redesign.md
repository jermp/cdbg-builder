# k-overlap walker + stitch redesign — design & handoff

Session-resume artifact. Goal: make the cross-bucket stitch phase honor
the `-g` RAM budget at 661k-genome scale by adopting GGCAT's **k-base
super-k-mer overlap** + **full-k-mer stitch keying**, which eliminates
the two unbounded in-RAM structures the current stitch still holds.

A fresh session should read this end-to-end, then implement the
file-by-file plan in §7, validating against the test oracle at each
step. Branch to use: `claude/stitch-fully-external` (off main; currently
holds only the merged ext-mem stitch, no in-progress edits).

> **Correctness note (count is 86630, not 87297).** Earlier revisions of
> this doc cited 87297 unitigs on salmonella_10, taken from the `main`
> branch on the assumption it was ground truth. It is not: `main`
> OVER-splits 626 unitigs into 1293 at provably simple-path points (0
> internal dBG branches, 0 color changes, unique junction k-mers with
> real edges). The correct colored compacted dBG has **86630** unitigs,
> confirmed by an independent naive Python builder (now folded into
> `test_data/verify.py`, which compares by exact unitig-set equality). The
> shipped pipeline produces exactly that. Two walker/stitch bugs that had
> inflated the count to 440k were fixed at the same time (walk every
> k-mer, not just color-owned ones; coalesce reconciled color runs).
> Do NOT "fix" the builder to reproduce 87297 — that count is wrong.

---

## 1. The problem (why this work exists)

On the Blackwell 661k pangenome run (`-g 64`):
- bucket-write peaked 30.6 GiB (under budget — fine).
- bucket-process peaked 75.9 GiB (~12 GiB OVER — separate issue, the
  resplit work; NOT this task).
- The OLD in-RAM stitch **segfaulted**: its `by_junction` map + `adj`
  array + `visited` were all O(num_fragments) resident — hundreds of GB
  at 5.7e9 fragments.

The ext-mem stitch (merged, PR #31, `stitch_extmem.hpp`) fixed the
segfault and is **correct** (exact, deterministic 86630 unitigs on
salmonella_10; all test paths pass). BUT it is **not memory-bounded**:
two structures are still fully in-RAM and blow `-g` at 661k:

- **(A) frag_unitig_reader index**: a `vector<entry>` (24 B/fragment)
  built once over the whole frag spill → **~137 GB** at 5.7e9 fragments.
- **(B) joinable_set / seed-pass jmap**: a global per-junction map built
  in one pass to decide which (k-1) junctions are joinable →
  **tens of GB to 100+** at 661k.

This task removes BOTH by aligning with GGCAT's design.

---

## 2. Root cause: we use k-1 overlap; GGCAT uses k

Our minimizer-bucketing (`bucket_ingester.hpp::emit_super_kmers`) splits
super-k-mers at minimizer changes with the next super starting at k-mer
`i` while the old super ended at k-mer `i-1`. Consecutive supers overlap
by **k-1 bases** — they share only the (k-1)-mer *junction* between the
last k-mer of one and the first k-mer of the next. The two boundary
k-mers (at positions i-1 and i) are DIFFERENT and live in DIFFERENT
buckets.

Consequences that cascade from this single choice:
- The stitch must match adjacent fragments on the **(k-1)-mer junction**.
- A (k-1)-mer junction can be shared by 3+ fragments at a dBG branch
  (predecessor's right + two successors' lefts), so the stitch needs an
  explicit **global "is this junction joinable / a branch" map** (B) to
  avoid false joins. That map is the price of (k-1) keying's lossiness.

**GGCAT splits supers with k-base overlap** (verified — see §3). Two
consecutive supers share their full boundary **k-mer**, which is
physically present in BOTH buckets. The stitch then keys on the full
boundary k-mer, and:
- Branch detection becomes **implicit and local**: distinct k-mers
  around a dBG branch hash to different buckets/slots and never collide,
  so NO global joinable map is needed (B vanishes).
- Each round just streams its bucket files, never a global fragment
  index, so (A) vanishes too.

So the fix is one architectural change — adopt k-overlap — and both (A)
and (B) fall out.

---

## 3. What GGCAT actually does (verified from /tmp/ggcat-v2 source)

Two subagent studies this session established the following, with
file:line cites preserved here for the new session.

### 3.1 k-base overlap is created at the bucketing step
`crates/assembler_minimizer_bucketing/src/lib.rs:240-256`: the next
super-k-mer's `last_index` is set to `index` and the slice is
`(last_index-1)..(index + k - 1)`, so consecutive supers overlap by
exactly **k** bases (one full boundary k-mer). `include_first=false`
after the first split clears `READ_FLAG_INCL_BEGIN` for all but the
first super of a read.

### 3.2 Per-k-mer "ignored boundary" flags
`crates/assembler_kmers_merge/src/unitigs_extender/hashmap.rs:382-398`:
each k-mer gets `begin_ignored` / `end_ignored` flags. The boundary
k-mer is the LAST k-mer of a super in bucket A (→ `end_ignored` in A)
and the FIRST k-mer of a super in bucket B (→ `begin_ignored` in B).
These flags (`READ_FLAG_INCL_BEGIN`/`_END`) tell the walker that k-mer
is the cross-bucket boundary — it is the mechanism that prevents the
duplicated boundary k-mer from being mis-counted as a local branch.

### 3.3 Walker includes the boundary k-mer, then breaks
`hashmap.rs:276-295`: the new base is pushed to `output` BEFORE the
contig-break check; `forward_seq`/`backward_seq` are pre-seeded with the
full k-mer (`hashmap.rs:489-496`). So the emitted fragment's last (or
first) k bases ARE the boundary k-mer. The break is signalled when the
entry's flag == exactly `READ_FLAG_INCL_BEGIN` or `_END`.

### 3.4 Cross-bucket join keys on the full k-mer
`final_executor.rs:147-174`: `get_extremal_hash(out_seq, k,
hash_beginning)` hashes a length-**k** window (`hashes/src/extremal.rs:
25-43`), canonicalises (`should_rc = !hash.is_forward()`), and routes to
a bucket by that hash. `extend_unitigs.rs:369-419`: the join compares
the full k-base slice via `get_extremity` + `equality_compare_start_
zero` — full-k-mer equality, not k-1.

### 3.5 Pairing rule (no global map, no branch count)
`extend_unitigs.rs:410-481`: per-bucket thread-local `HashTable<usize>`
keyed by the boundary k-mer. First arrival stashed (Vacant); second
arrival → `try_join` → on success `entry.remove()` and the joined
fragment is re-routed to next round by its NEW open end's hash. A third
arrival just refills the slot. **No count, no branch flag** — the full
k-mer key makes real branches structurally route elsewhere.

### 3.6 Per-round memory bound
`extend_unitigs.rs:288` buckets_count = `log2(num_inputs)`;
`:317-319` processes one bucket at a time per worker (LPT order,
largest first); `:707-710` clears per-bucket state; input bucket files
unlinked on read (`RemoveFileMode::Remove`). Resident set per round ≈
(#threads) × (one bucket's fragments + its hash table). No global index,
no global junction map.

---

## 4. What we tried this session and what broke

The ext-mem stitch was built with (k-1) keying + a precomputed
`joinable_set` (the in-RAM stitch's exact predicate, invariant under
merging). It is CORRECT and merged (PR #31), but holds (A) and (B) in
RAM.

Attempt to adopt k-overlap with a 1-line change to
`emit_super_kmers` (`super_start = i - 1` instead of `i`): **broke
badly**. e2e on salmonella_10 showed:
- bucket fragments 1986195 vs ~87k (23× over-fragmentation),
- distinct color classes 688 vs 171 (color-class duplication),
- correctness FAIL: one cid mapped to two different color bitmasks.

**Why**: the per-bucket walker (`bucket_walker.hpp`) assumes each
canonical k-mer lives in exactly ONE bucket. With k-overlap, the
boundary k-mer is now in TWO buckets, and:
- the walker's `kmer_info` + phantom-edge logic mis-handle the
  duplicated boundary k-mer (it has a local edge in each bucket that it
  didn't have before, and its successor/predecessor relationships
  change), causing massive fragmentation;
- the same canonical k-mer in two buckets gets interned into two
  INDEPENDENT per-bucket color dicts → two different local cids → global
  dedup only partially reconciles → duplicated/wrong color classes.

So k-overlap is an ARCHITECTURAL change to the walker, not a 1-liner.
The walker needs GGCAT's "ignored boundary" k-mer ownership model
(§3.2): a boundary k-mer is PRIMARY in exactly one bucket (the one whose
minimizer owns it) and IGNORED in the other (present only to provide the
overlap k-mer for stitch matching). The ignored copy must:
- NOT be claimed/walked as a normal node,
- NOT contribute its color to a per-bucket color class (or it must, but
  consistently — see §6 open question),
- only serve as the boundary k-mer the fragment includes at its open
  end.

---

## 5. The plan (architectural)

Adopt GGCAT's model end to end:

1. **bucket-write** (`bucket_ingester.hpp`, `super_kmer.hpp`,
   `bucket_io.hpp`): k-base overlap between consecutive supers. Add a
   per-super-k-mer or per-boundary flag identifying which boundary
   k-mers are "ignored" (primary in the adjacent bucket). GGCAT marks
   the begin/end k-mer of each super; replicate that. The super_kmer
   record format may need an extra flag bit or a length convention.

2. **bucket-walker** (`bucket_walker.hpp`): teach `load_bucket` +
   the unitig walk the ignored-boundary model. A boundary k-mer that is
   IGNORED in this bucket:
   - is included in fragments at OPEN ends so the emitted fragment's
     boundary k bases equal the boundary k-mer (matches the adjacent
     bucket's fragment),
   - is NOT treated as an interior node for branching/degree decisions,
   - its color contribution is handled so the per-bucket color class of
     the OWNING bucket is the authoritative one.
   This is the hard part. Mirror GGCAT's `READ_FLAG_INCL_BEGIN/END`
   handling in `hashmap.rs:276-398`.

3. **stitch** (`stitch_extmem.hpp`): switch keying from the (k-1)
   `side_junction_canonical` to the **full boundary k-mer** (canonical).
   Delete `joinable_set_t` and all branch-set machinery. Pairing becomes
   GGCAT's first-come/second-come on a per-bucket hash table keyed by the
   full-k-mer hash; keep OUR cid gate inside the join (GGCAT has none —
   our unitigs are monochromatic). Concat drops **k** bases (the shared
   boundary k-mer), not k-1.

4. **in-RAM stitch** (`stitch.hpp`): update its concat to drop k and its
   matching to full-k-mer, OR retire it. It currently serves as the
   test oracle's reference; if kept, it must match the new convention.

5. **streaming Source**: replace `frag_unitig_reader`'s eager
   `vector<entry>` index with a forward-only stream (`frag_unitig_
   stream`) yielding (cid, open_flags, seq) one record at a time. After
   this, the only consumer (round 0 of stitch) reads it sequentially, so
   no 137 GB index. Easy once §3-style streaming is the only access.

6. **test oracle** (`test/gen.hpp::split_unitig`,
   `test/test_stitch.cpp`): `split_unitig` must produce k-overlap
   fragments (adjacent fragments share their full boundary k-mer), and
   add a SHARED-CID / branchy synthetic case (the current oracle uses
   distinct cids per unitig, which structurally cannot catch same-cid
   false joins — the e2e caught those instead this session).

---

## 6. Open questions to resolve while implementing

1. **Color of the ignored boundary k-mer.** In GGCAT colors attach to
   k-mers and the ignored boundary is handled by the flag. In OUR design
   colors attach per super-k-mer record (a record carries its color
   list). When the boundary k-mer is duplicated into two buckets, does
   each bucket's record carry the same color list? If the supers come
   from the same input run they should — verify. The per-bucket color
   dict must not produce two different cids for what is one color class.
   This is the source of the "688 vs 171 / cid maps to two bitmasks"
   failure; get it right.

2. **k-mer count accounting.** Boundary k-mers now appear in 2 buckets.
   `n_kmers` profiling double-counts them. The verify.py k-mer-set check
   must still see each canonical k-mer's color set as the UNION across
   its occurrences — confirm the global dict dedup reconciles the two
   bucket-local copies into one global class.

3. **Bucket-size inflation.** k-overlap adds one duplicated k-mer per
   super-k-mer boundary. Quantify the bucket-byte increase (should be
   modest — boundaries are rare vs interior k-mers — but measure on 25K).

4. **resplit interaction.** `bucket_resplit.hpp` re-windows supers with
   m2; it must preserve the k-overlap + ignored-boundary flags. Defer
   until the base k-overlap design works, but keep in mind.

---

## 7. File-by-file change list (implementation order)

Validate `bash test_data/run_test.sh` (expect 86630 unitigs == the
independent ground-truth ccdBG, 171 color sets) AND `build/test_stitch`
after each step that can be tested.

1. `test/gen.hpp` — update `split_unitig` to k-overlap; add a shared-cid
   branchy generator. (Do FIRST so the oracle is ready.)
2. `bucket_ingester.hpp` + `super_kmer.hpp` — k-overlap split + boundary
   flag.
3. `bucket_walker.hpp` — ignored-boundary ownership model in load_bucket
   + walk. (Hardest; expect iteration. Test e2e after.)
4. `stitch_extmem.hpp` — full-k-mer keying, delete joinable_set,
   GGCAT pairing, concat drops k. Test against test_stitch.
5. `stitch.hpp` — full-k-mer concat (oracle) or retire.
6. `frag_unitig_*` in `unitig_spill.hpp` — streaming Source; wire into
   `builder.hpp` stitch call.
7. End-to-end: salmonella_10 matches the ground-truth ccdBG (86630
   unitigs via `verify.py`); then a larger run to confirm peak RAM
   tracks `-g`.

---

## 8. Guardrails (lessons from this session — please honor)

- This area is subtle; EVERY stitch/walker change this session that
  looked like a 1-liner turned out to be architectural. Make ONE change,
  build, run BOTH `test_stitch` and the e2e, read the numbers, commit,
  then next change. Do not batch.
- `test_stitch`'s distinct-cid oracle does NOT catch same-cid false
  joins. The e2e (real shared-cid data) is the real correctness gate.
  Add the shared-cid oracle case (§5.6) early.
- Watch for over-fragmentation (fragment count ≫ expected) and
  color-class duplication (distinct classes ≫ expected) — both are the
  signature of the walker mis-handling duplicated boundary k-mers.
- Never mix file edits and git commands in one parallel tool batch
  (caused repo-state confusion earlier this session).
- Use `uint64_t` for ALL loop variables / indices / counts (project
  rule; cid and k are the only narrower domain types).

---

## 9. Pointers

- This repo: `jermp/cdbg-builder`. Current main has the merged ext-mem
  stitch (`stitch_extmem.hpp`), 64-bit indices, cid=uint64.
- GGCAT v2 reference: clone at `/tmp/ggcat-v2` (may need re-clone:
  `git clone --depth 1 https://github.com/algbio/ggcat.git` +
  `git submodule update --init libs-crates/parallel-processor-rs`).
  Key files: `crates/assembler_minimizer_bucketing/src/lib.rs`,
  `crates/assembler_kmers_merge/src/unitigs_extender/hashmap.rs`,
  `crates/assembler_kmers_merge/src/final_executor.rs`,
  `crates/assembler_pipeline/src/extend_unitigs.rs`,
  `crates/hashes/src/extremal.rs`.
- algorithm.md (this repo) — current pipeline description; §5 (stitch)
  and §4 (bucket-process) will need updating after this redesign.
- Open follow-ups NOT in this task: bucket-process resplit
  (`claude/resplit-on-extmem`, over-fragments — needs its own fix);
  color-sets-partitioner (separate repo, design doc on
  `claude/csp-design`).

---

## 10. Quick-start prompt for the new session

> Read `k-overlap-redesign.md` in full. We're aligning cdbg-builder's
> per-bucket walker and cross-bucket stitch with GGCAT's k-base
> super-k-mer overlap + full-k-mer stitch keying, to drop the two
> in-RAM structures (the 137 GB frag index and the ~100 GB joinable
> map) that currently keep stitch from honoring `-g` at 661k scale.
> Work on branch `claude/stitch-fully-external`. Implement §7 step by
> step, honoring the §8 guardrails: one change at a time, build + run
> BOTH `build/test_stitch` and `bash test_data/run_test.sh`, read the
> numbers, commit, repeat. The GGCAT reference is at /tmp/ggcat-v2
> (§3 has the exact file:line mechanisms). Start with §7 step 1.
