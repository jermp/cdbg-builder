# Algorithm

This document describes how `cdbg-builder` constructs a **colored compacted
de Bruijn graph (ccdBG)** and how it builds and compresses the color sets.


## 1. What the tool produces

Given `N` input files (FASTA / FASTQ, optionally gzipped), each treated as
one **color** (file `i` = color `i`), the tool writes three artifacts under
a user-supplied basename:

- **`<out>.fa`** — colored unitigs in FASTA. Each record is a maximal
  *monochromatic* unitig (a path in the de Bruijn graph whose every
  k-mer has the same color set). Headers are the integer color-set id.
- **`<out>.u2c`** — a `bits::bit_vector` of length `num_unitigs` where
  bit *i* is set iff unitig *i* (in `.fa` emission order) is the last
  unitig of a color-set run. Popcount = `num_color_classes`. Downstream
  consumers (e.g., [Fulgor](https://github.com/jermp/fulgor))
  recover the per-unitig color-set id via
  `rank1(unitig_id)` after building a rank index over this bit_vector.
- **`<out>.color_sets`** — header + hybrid-encoded distinct color sets +
  Elias–Fano of per-class bit-offsets. See §6 for the on-disk layout
  and §5 for the encoding scheme.

A k-mer is **canonical** (the lexicographic min of itself and its
reverse complement). Two adjacent canonical k-mers that share a (k−1)
suffix/prefix in some orientation form a dBG edge. A **unitig** is a
maximal path with internal nodes of in-degree = out-degree = 1.

### What "correct" means (and a count footgun)

Correctness is **exact unitig-set equality** against an independently
built naive ccdBG, checked by `test_data/verify.py` (coverage,
monochromaticity, topology, *and* maximality). On `salmonella_10` the
correct colored compacted dBG has **86,630 unitigs**.



## 2. Pipeline overview

The build runs in four sequential phases: (1) bucket-write, (2) bucket-process, (3) stitch, (4) emit.

The four phases are strictly serial: each consumes the prior phase's
on-disk output and finishes before the next starts. Within a phase,
work is parallelized over input files (bucket-write) or buckets
(bucket-process and stitch). Only emit is single-threaded.

The data flow in one picture:

```
inputs (gzip FASTA, N files)
    │   parse + minimizer + bucket
    ▼
B per-bucket files of LZ4-framed super-k-mer records      ← bucket-write
    │   per-bucket dBG; walk; intern colors
    ▼
single frag_unitigs.bin file (cid + flags + seq) +        ← bucket-process
streaming write to out.color_sets
    │   join open ends across buckets via full boundary k-mer match
    ▼
K cid-range unitig_bucket files                           ← stitch
    │   sort by cid, write FASTA, build u2c
    ▼
out.fa + out.u2c + finalize out.color_sets                ← emit
```

### 2.1 Phase I/O contracts

Each phase communicates with the next **only through on-disk artifacts** —
there is no shared in-RAM state across phases (the one exception is the
streaming color dict, noted below). This is what lets a phase be
optimized, or even swapped, in isolation: as long as it honors the
artifact contract, the rest of the pipeline is unaffected.

**Phase 1 — bucket-write** (`bucket_io.hpp`, `bucket_ingester.hpp`)

- *Reads:* the `N` input gzip-FASTA files (one **color** per file).
- *Produces:* `tmp/bucket_<b>.bin`, `b ∈ [0, B)` — one LZ4-framed file per
  bucket, each a stream of compacted `super_kmer` records
  `{2-bit bases, 4-bit flags, color list}`. Colors are local file indices,
  deduped (union per identical super-k-mer) **within the bucket**.
- *Why it's the natural handoff:* minimizer bucketing (§3.2) routes every
  dBG edge and branch into a single bucket, so each `bucket_<b>.bin` is an
  **independent dBG sub-problem** — phase 2 can process buckets in any
  order, in parallel, with no cross-bucket coordination.

**Phase 2 — bucket-process** (`bucket_walker.hpp`)

- *Reads:* the `B` `tmp/bucket_<b>.bin` files (one bucket per worker).
- *Produces two artifacts:*
  - `tmp/frag_unitigs.bin` — a single stream of **open-ended fragments**
    `{ACGT seq, open_flags (which ends are still extendable), color-run
    sequence}`. The color runs reference **global** cids (the per-bucket
    local→global merge, §4.7, runs before each fragment is written).
  - `<out>.color_sets` — the distinct **global color classes**, interned
    and streamed to disk incrementally during the phase (the only artifact
    written progressively rather than at phase end; finalized by phase 4).
- *Why it's the natural handoff:* fragments are maximal colorless
  (topological) pieces whose open ends each carry their boundary k-mer, so
  phase 3 only has to match boundary k-mers — and global cids decouple the
  (large) color storage from the (small) topology that still needs joining.

**Phase 3 — stitch** (`stitch_extmem.hpp`)

- *Reads:* `tmp/frag_unitigs.bin` (streamed one fragment at a time).
- *Produces:* `tmp/unitig_bucket_<k>.bin`, `k ∈ [0, K)` — **finished,
  monochromatic** unitigs `{ACGT seq, cid}`, partitioned so bucket `k`
  holds every unitig whose `cid` falls in `k`'s range.
- *Why it's the natural handoff:* joining open ends by full boundary-k-mer
  match yields complete unitigs; cid-**range** bucketing means phase 4 can
  emit in strict cid order by walking buckets `0..K-1` (sorting only within
  a bucket), with no global sort over all unitigs.

**Phase 4 — emit** (`builder.hpp`, `unitig_spill.hpp`)

- *Reads:* `tmp/unitig_bucket_<k>.bin` (cid order) **and** `<out>.color_sets`
  (to finalize it).
- *Produces:* `<out>.fa` (cid-ascending FASTA), `<out>.u2c`
  (unitig→color-set bit_vector), and the finalized `<out>.color_sets`
  (bits + Elias–Fano offsets + header). These three are the tool's outputs
  (consumed downstream by Fulgor, §1).

All `tmp/*` artifacts live under the scratch dir and are removed at the
end of `build()`; only the three `<out>.*` files persist.

---

## 3. Phase 1: bucket-write (`include/bucket_io.hpp`, `bucket_ingester.hpp`)

### 3.1 Goal

Read every input file (libdeflate-decompressed if `.gz`), decompose
each ACGT-only run into super-k-mers, compute each super-k-mer's
**canonical minimizer**, and append it to the bucket file selected by
that minimizer. After this phase:

- Each bucket file holds one stream of super-k-mer records.
- Records are *compacted*: a super-k-mer that occurred in many input
  files is stored once with the union of its colors.
- Every k-mer that lives in bucket `b` is reachable from records in
  bucket `b` alone (by construction of the minimizer-bucketing).

### 3.2 Minimizer bucketing ((k−1)-mer minimizers, BCALM2/GGCAT)

Within an ACGT run of length `L ≥ k`, slide a window of size
`W = k − m` m-mers over the run's **(k−1)-mers** (a (k−1)-mer has
`(k−1) − m + 1 = k − m` m-mers). The minimizer of a (k−1)-mer is the
m-mer whose `canonical_mhash` (canonical ntHash; symmetric under RC) is
minimum across its window. Each maximal run of consecutive (k−1)-mers
sharing one minimizer defines a **super-k-mer**; its bucket id is
`(min_hash >> 1) % B` (the low bit is reserved to match GGCAT's
`cn_nthash.rs` convention; `B = num_buckets` need not be a power of two,
so a modulo — not a mask — maps the hash, see §3.7).

Minimizing over **(k−1)-mers**, not k-mers, is what keeps the de Bruijn
graph correct under bucketing (the BCALM2 rule, matching GGCAT's
`BatchMinQueue::new(k - m)`). Every dBG edge `X → Y` shares the junction
(k−1)-mer `S = suffix(X) = prefix(Y)`, so both endpoints route to
`bucket(min(S))`; a branch `X → {Y1, Y2}` shares one `S`, so `X`, `Y1`,
`Y2` all land in `bucket(min(S))` and the branch is locally visible in a
single bucket. (Minimizing over k-mers instead scatters the branch arms
across buckets, so no bucket sees the branch — it would emit false
through-paths the stitch then joins incorrectly.)

A k-mer therefore lands in **two** buckets exactly when its prefix and
suffix (k−1)-mers have different minimizers (`lmin ≠ rmin`); consecutive
super-k-mers share that full boundary k-mer (a **k-base overlap**). The
shared boundary k-mer is **colored in exactly one** bucket —
`bucket(min(lmin, rmin))` — so its color is the full union over its
occurrences there, with no partial-color split (§4).

Each super-k-mer record carries:

- `flags` (4 bits): `IS_ACGT_BEGIN` / `IS_ACGT_END` (the super's first /
  last k-mer has no predecessor / successor in the input run), and
  `OWNS_FIRST` / `OWNS_LAST` — whether this bucket owns (colors) the
  super's first / last boundary k-mer per the `min(lmin, rmin)` rule.
  At a split between an ending super (run minimizer `mA`) and a starting
  super (`mB`), the ending super owns its last shared k-mer iff `mA < mB`
  and the starting super owns its first iff `mB < mA`.
- `colors`: at write time, the single color of this input file.
  Compaction merges color lists across files later in this phase.
- `bases` : 2-bit-packed bases.

### 3.3 Per-thread per-bucket buffers and flush

Each ingest worker thread owns a `per_thread_bucket_buffers` (one
`std::vector<pending_record>` and one `std::vector<uint8_t>` of bases
per bucket). `append()` pushes the record's metadata + appends the
bases bytes + **precomputes a wyhash of the super-k-mer bases** while
they're still hot in L1 from the buffer write — no lock, no hashmap.
The hash is stashed in `pending_record::hash` so the compactor's
later dedup-find doesn't have to re-touch the (typically cold) bases
bytes. When a bucket's bases-buffer crosses `flush_bases` (a fixed
good default of 4 KiB, overridable via `--flush`; §3.7), the thread
acquires that bucket's mutex and calls `bucket_compactor::insert_batch`.

### 3.4 Per-bucket compactor (online dedup)

Each `bucket_compactor` (one per bucket) holds:

- `m_dedup`: an `ankerl::unordered_dense::map<std::string, entry>`
  keyed by the 2-bit-packed super-k-mer bytes, with `entry =
  { colors, flags }`.
- `m_batch_buf`, `m_out_buf`: scratch vectors reused across spills.
- `m_file`: the open bucket file handle.

`insert_batch` walks the incoming batch, finds each super-k-mer in
`m_dedup` via heterogeneous lookup with a `hashed_view{string_view,
uint64_t}` whose hasher returns the precomputed hash directly. The
bases bytes are only re-read on a hash collision (then `string_eq`
runs the byte compare against the stored `std::string`). On hit it
appends the new color to the existing entry's `colors` vector and
AND-shrinks the flags (`e.flags &= r.flags`), so
`IS_ACGT_BEGIN`/`IS_ACGT_END` survive on the deduped record only
when *every* contributor agreed. On miss it constructs the
std::string and emplaces. Sort+unique on `colors` is deferred until
spill.

### 3.5 Spill: serialize, LZ4-compress, write a frame

When `m_dedup`'s tracked-byte counter crosses `spill_bytes`
(a fixed good default of 64 KiB, overridable via `--spill`; §3.7),
`spill_locked` runs:

1. For each entry: sort+unique its `colors`, write a `super_kmer`
   record (varint num_colors + varint color deltas + varint
   length+flags + 2-bit bases) into `m_batch_buf`.
2. `LZ4_compress_default` the whole `m_batch_buf` into `m_out_buf`
   as a single block.
3. Write `[u32 uncompressed_size] [u32 compressed_size] [bytes]` to
   the bucket file. End-of-stream marker is `[u32 0]`.
4. Clear `m_dedup` and `m_batch_buf`.

Per-spill batching is essential: LZ4 has a few-byte per-block framing
overhead; per-record compression of ~100-byte inputs ended up
producing files larger than the input. Per-spill batching gives
proper compression ratio on the joint stream.

### 3.6 RSS pressure watcher

A background thread polls live RSS via `/proc/self/status:VmRSS`
every 100 ms. On HIGH (default 60 % of `-g`), it sweeps every
compactor and force-spills any with pending state. When RSS drops to
LOW (default 45 %), pressure clears and ingest threads return to
their normal `m_bytes`-based spill cadence. Hysteresis prevents the
death-spiral observed with sticky pressure (every batch spilling for
the rest of the phase, exploding spill counts and bucket file sizes).

If `/proc/self/status` is unavailable (macOS, sandbox), the watcher
falls back to monotonic `getrusage` peak; pressure becomes sticky once
tripped (less efficient but still correct).

### 3.7 Sizing the bucket count from the RAM budget

`flush_bases` and `spill_bytes` are **fixed good defaults** (4 KiB /
64 KiB), not auto-tuned; the **bucket count `B` is derived from them**
and the `-g` budget (`builder::auto_bucket_count_`). bucket-write's two
RAM pools are modelled as

```
buffer pool    ≈ alpha * T * B * flush     (per-thread per-bucket buffers)
compactor pool ≈ beta  * B * spill         (per-bucket dedup hashmaps)
```

with empirical overhead factors `alpha = 2` (recs vector + vector
capacity slack) and `beta = 7` (std::string keys + map structure +
glibc fragmentation + LZ4 scratch). Holding the sum to a share `M =
BUCKET_WRITE_BUDGET_FRAC × g` of the budget and solving for B:

```
B = M / (alpha * T * flush + beta * spill)
```

So a **larger** spill (better dedup) yields a **smaller** B, and the
peak stays pinned at `M` *by construction* — the knobs trade B
(bucket-process parallelism, §4.9) for dedup quality (bucket-process
`load` speed), not budget. `BUCKET_WRITE_BUDGET_FRAC = 0.50`; the rest
of `-g` is left for the phases that follow and for the model's ~25 %
real-vs-modelled undercount, which the RSS watcher (§3.6) backs at
runtime. On macOS the budget is additionally divided by
`PLATFORM_RAM_OVERHEAD = 2` (libsystem_malloc retains more pages).

`alpha`, `beta`, `flush`, and `spill` are all overridable on the CLI
(`--alpha/--beta/--flush/--spill`) for per-machine recalibration. `-b`
still forces a power-of-two `B = 2^bucket_log2` override. Because B is
an arbitrary integer (not a power of two), bucketing maps the minimizer
hash with a modulo, `(h >> 1) % B` (§3.2) — the divide is once per
super-k-mer, negligible against the per-base ntHash scan. No `-g` set →
historical `B = 1024`.

### 3.8 Output

Per bucket: a single LZ4-framed file. Frame format:
```
repeat:
  [u32 uncompressed_size]   little-endian; 0 = end-of-stream
  [u32 compressed_size]
  [compressed bytes]        LZ4_compress_default output
```

Inside one decompressed frame: a concatenation of `super_kmer`
records (varint-delta colors + varint length+flags + 2-bit packed
bases). See `super_kmer.hpp`.

---

## 4. Phase 2: bucket-process (`include/bucket_walker.hpp`)

### 4.1 Goal

For each bucket independently:

1. Build the local dBG over canonical k-mers in that bucket, with
   per-k-mer **cross-bucket boundary flags** (GGCAT's begin/end-ignored
   model — no phantom edges, no global joinable map).
2. Walk maximal topological unitig fragments, carrying an RLE color-run
   sequence (colorless extension; the monochromatic split happens at
   emit).
3. Intern color lists into a per-bucket `compact_color_set_dict`, then
   merge per-bucket dicts into a shared `streaming_color_set_dict` whose
   output goes directly to `<out>.color_sets`.

Each bucket's fragments are emitted as records in a single disk-backed
`frag_unitig_writer` for stitch to consume.

### 4.2 Per-bucket walker state

Three structures per in-flight bucket:

- **`kmer_info`** (`bucket_kmer_map`): canonical-k-mer → `{kmer_entry,
  flags, primary}`. `flags` is two bits (`KMER_BOUND_LEFT/RIGHT`)
  marking cross-bucket contig boundaries; `primary` says whether this
  bucket owns (colors) the k-mer. `kmer_entry` is `{first_rsid,
  vector<rest_rsid>}` — almost every k-mer is contributed by a single
  super-k-mer record, so `rest` is empty in the common case.
- **`record_sets`** (`compact_color_set_dict`): interns the sorted-
  deduped color list of every super-k-mer record read from disk.
- **`local_dict`** (`compact_color_set_dict`): per-bucket final color
  classes; the returned local cid rides on each unitig's color runs.

### 4.3 Loading a bucket

`load_bucket`, per super-k-mer record:

1. Intern colors into `record_sets` → get rsid.
2. Roll the canonical k-mer through the record's bases. For each k-mer:
   - **Boundary flags**: the first k-mer gets a boundary bit if the
     record's `IS_ACGT_BEGIN` is clear (it is the begin-ignored overlap
     copy of the previous super); the last gets one if `IS_ACGT_END` is
     clear. The bit is oriented to the k-mer's canonical frame. Flags are
     OR-accumulated over all occurrences. (GGCAT `hashmap.rs:382-398`.)
   - **Ownership / color**: color the k-mer (and mark `primary`) unless
     this bucket does not own it. An interior k-mer is always owned; the
     first / last k-mer is owned iff the record's `OWNS_FIRST` /
     `OWNS_LAST` flag is set (§3.2). So each k-mer is colored in exactly
     one bucket, accumulating the full color union there.

### 4.4 Contig-break flags (no phantom edges)

A k-mer is a **cross-bucket contig boundary** iff exactly one boundary
bit is set (`flags == 1` or `2`): the unitig ends there and continues
into the adjacent bucket, so that end is emitted **OPEN** with the full
boundary k-mer included for the stitch to match on. `flags == 3` (both
sides) is *not* a break — the k-mer is walked through as interior;
`flags == 0` is plain interior (GGCAT `hashmap.rs:291`). To avoid
threading through a real branch whose alternate arm crosses buckets, an
end is only opened when the boundary side's **local degree ≤ 1**;
local degree ≥ 2 is a genuine branch and is CLOSED. Because the (k−1)-mer
bucketing co-locates every branch's arms in one bucket (§3.2), that local
check sees all the arms.

### 4.5 Computing per-k-mer cids and the walk

After load, build `cid_of[k-mer]` for every **primary** (owned) k-mer:
- *Single-rsid* (common): `local_dict.intern(record_sets.at(rsid))`,
  with an `rsid_to_cid` cache for O(1) repeats.
- *Multi-rsid*: sort+unique the rsids; one distinct rsid → cached path;
  else union the referenced color lists and intern that.

After this, `record_sets` and the per-k-mer rsid storage are released;
the walk needs only `kmer_info` (flags + primary) and `cid_of`.

The walk seeds from **every** unvisited k-mer present in the bucket -- owned
*or* foreign -- not just the color-owned ones, and extends both ways along
degree-1 simple paths (GGCAT `hashmap.rs:455-545`), stopping at a contig
break (OPEN, including the boundary k-mer) or a branch/dead-end (CLOSED).
Seeding must NOT be gated on ownership: a short super-k-mer that owns none of
its k-mers (its two boundary k-mers are colored in the two adjacent buckets)
would otherwise never be walked, yet it is exactly the bridge fragment whose
open ends are the stitch partners of those two boundary k-mers -- dropping it
strands both with no partner and the cross-bucket joins fail (this caused
salmonella-10 to emit 440k fragmented unitigs instead of the correct 86,630).
Ownership gates COLORING only. It emits `stitchable_unitig` records `{seq, runs, open_flags}`,
where `runs` is the RLE color-run sequence over the k-mers; a foreign
(unowned) boundary k-mer at an open end carries a `COLOR_RUN_FOREIGN`
placeholder that the stitch reconciles to the owning side's color at the
join. Local cids are remapped to global cids in §4.7.

### 4.6 Streaming the global color-set dict

A single `streaming_color_set_dict` is shared across all bucket
threads, guarded by `global_mu`. The dict exposes two intern
entry points:

- `intern(colors)` — computes the 128-bit content hash inline (wyhash
  primary + fnv1a secondary), then dispatches to ↓.
- `intern_with_hashes(colors, h)` — caller supplies the precomputed
  hash. Used by the batched merge in §4.7 so the (heavy) hash
  compute happens lock-free.

Either way, the hot work under the lock is:

1. Look up in `m_index` (hash-only key, no byte-compare). On hit
   return existing id; collision rate is ~2^-64 per pair, negligible
   over our class counts.
2. On miss: hybrid-encode the new class's bits via
   `hybrid_builder::encode_one` into an in-memory
   `bits::bit_vector::builder`.
3. Append every COMPLETE 64-bit word from the builder to
   `<out>.color_sets` (the file is opened at construction with a
   placeholder header). Only the trailing partial word stays in
   memory.
4. Append the new class's `bit_offset` to a sidecar file
   `<out>.color_sets.tmp_offsets` for the EF later.

Per-class metadata held in RAM: 16 B (primary + secondary hash). The
compressed bit_vector itself is never resident.

### 4.7 Per-bucket → global merge (batched, hash precomputed)

After `process_bucket` returns, the worker thread merges its
`local_dict` into `global_dict`. The naive single-mutex version
(decode + hash + dedup-find + on-miss-encode all under `global_mu`)
serialised 32 threads on multi-thousand-color lists; the
bucket-process diagnostic showed ~80 % of bucket-process wall on
50K was lock contention.

The current path splits per-class work into a lock-free pre-phase
and a much smaller locked phase, in batches of `MERGE_BATCH = 32`:

```cpp
for each batch of up to 32 local classes:
    // Lock-free.
    for i in batch:
        local_dict.at(lc, batch[i].colors);
        batch[i].h = streaming_color_set_dict::compute_hashes(batch[i].colors);
    // Locked.
    {
        lock_guard<mutex> lk(global_mu);
        for i in batch:
            local_to_global[lc] = global_dict.intern_with_hashes(
                std::move(batch[i].colors), batch[i].h);
    }
```

The lock-free pre-phase does the local-dict hybrid decode (~tens of
µs per call on dense pangenome inputs) and the wyhash + fnv1a (~tens
of µs on a multi-thousand-color list). The locked phase is just the
hashmap probe (microseconds) plus, on miss, the hybrid encode +
sidecar offset write. On 50K, this brought `merge_wait` from 936 s
to ~120 s and total bucket-process from 1156 s to 364 s.

After the merge, the worker remaps each unitig fragment's `cid`
local → global and streams the fragment to the disk-backed
`frag_unitig_writer`. Per-fragment: `[u32 cid][u8 open_flags][u32
seq_len][seq bytes]`.

Memory: per-thread merge buffer is `MERGE_BATCH × max-decoded-color-
list`, a few MB at num_colors = 100K.

### 4.8 Bucket-process profiling

`bucket_process_prof` (in `util.hpp`) tracks phase-level wall-equiv
time across worker threads:

- `ns_load` — `bucket_reader` + LZ4 decompress + record intern + the
  per-record k-mer hashmap roll into `kmer_info`.
- `ns_resolve` — the rsid → local cid pass over `kmer_info`.
- `ns_walk` — `classify_left_end` + `extend_and_emit`.
- `ns_pre_decode`, `ns_pre_hash` — the lock-free half of §4.7.
- `ns_merge_lock_wait`, `ns_merge` — wait for vs. work under
  `global_mu`.

Plus counts: `n_buckets`, `n_records`, `n_kmers`, `n_local_classes`,
`n_unitigs`. Overhead is well under 1 % of bucket-process wall and is
left enabled by default.

### 4.9 Memory-admission gate (keeps the phase in `-g`, never caps threads)

bucket-process spawns the user's `T` worker threads; each pops the next
bucket from a shared atomic counter. Rather than reduce `T` when the
working set is large, an **admission gate** bounds how many buckets are
resident *at once* — `T` always stays alive, but under a tight `-g`
fewer (large) buckets load concurrently and more (small) ones do.

Before taking a bucket, a worker waits on a condition variable until
both hold (`process_buckets`):

1. **Reservation** — `carry + Σ reserved_working_set + live_dict +
   this_bucket ≤ BUCKET_PROCESS_BUDGET_FRAC × g` (0.82·g). Each bucket's
   reserved cost is `WORKING_SET_OVERHEAD (24) × its uncompressed bytes`
   — the *full* in-flight footprint (kmer_info + cid_of + visited +
   emitted fragments + local_dict), not just kmer_info. `carry` is the
   RSS measured at phase start and the color dict's resident size is
   re-read **live** each admission, so the monotonically-growing dict is
   reserved for dynamically instead of overflowing on top.
2. **Hard live-RSS ceiling** — `current_rss() + this_bucket ≤ budget`.
   Live RSS already includes carry, the dict, every loaded bucket's real
   footprint, and glibc fragmentation, so this backstops whatever the
   `24×` estimate under-reserves — *workload-independently*, which is why
   bucket-process stays in budget even at scales the constant wasn't
   calibrated on. One bucket is always admitted (forward progress).

For the gate's `carry` baseline to be honest, bucket-write's memory must
be released first: `builder` calls `writer->release_compactors()` (frees
the per-bucket dedup hashmaps + LZ4 buffers; the on-disk files and cached
sizes remain) and `release_free_heap_to_os_()` (`malloc_trim`) **before**
the phase. Skipping either leaves glibc holding bucket-write's pages, so
`current_rss()` reads near the budget and the gate starves to one bucket
at a time (serial). The bucket size distribution printed after
bucket-write (§3.8) previews how many buckets will fit at once.

---

## 5. Phase 3: stitch (`include/stitch_extmem.hpp`)

External-memory, GGCAT-faithful hash-bucketed iterative join
(`extend_unitigs.rs`). It glues the open-ended fragments from §4 into
maximal unitigs while holding only **`num_threads` hash buckets at a
time** in RAM (the per-round bucket loop runs `T`-way parallel; buckets
are independent within a round), so peak memory is bounded by the bucket
count, not the fragment count. The bucket count is dimensioned (not the
thread count) so `T` resident round-buckets fit stitch's share of `-g`.
(The older in-RAM stitch — one global junction map + a per-fragment
adjacency array, O(num_fragments) resident — was retired; `stitch.hpp`
now holds only shared helpers, §1.)

**Keying — the full boundary k-mer.** Two fragments that should glue
share their full boundary **k-mer** (the k-base overlap, §3.2), present
in both. The stitch keys each open end on the canonical full k-mer, so a
dBG branch's distinct k-mers route to different slots and never
false-join — no global "joinable junction" map is needed.

**Round loop.** Each round:

1. Every still-open fragment presents **one** open end (a both-open
   fragment picks which via a per-fragment RNG, re-rolled each round) and
   is routed to bucket `H(boundary_kmer)`.
2. Within a bucket, a hash table keyed by the boundary k-mer pairs
   arrivals first-come/second-come: the first waits, a compatible second
   **joins**. The join concatenates the two oriented sequences dropping
   the shared k-mer (k bases), and concatenates their color-run sequences
   reconciling the shared k-mer's color (a `COLOR_RUN_FOREIGN` placeholder
   from the unowned side takes the owned side's cid). An
   orientation-incompatible collision re-presents the other end next
   round.
3. A joined fragment that still has an open end is re-routed by its **new**
   end's hash for the next round; one with no open end is a finished
   unitig.

The loop ends when a full round produces **zero** joins (GGCAT's
`has_joinable_unitigs`); every remaining open fragment is then flushed as
a terminal unitig. A two-open-ended fragment that re-routes to its own
bucket with both extremal hashes equal is a circular unitig.

**Emit / monochromatic split.** A finished topological tig carries an RLE
color-run sequence; `ext_split_monochromatic` cuts it at run boundaries
into monochromatic output unitigs (cdbg output unitigs are
monochromatic — the one step GGCAT does not need, since it keeps colored
topological unitigs). Any extremal `COLOR_RUN_FOREIGN` placeholder that
never met its primary partner is dropped (the owning bucket emits it).

**Stores.** Two interchangeable round stores share the loop:
`round_store_mem` (in-RAM rounds; tests/dev) and `round_store_file`
(LZ4-framed files; the production external-memory path,
`stitch_unitigs_extmem_file_stream`, which fans the per-round bucket loop
across `num_threads`, holding one decoded bucket per in-flight thread).

---

## 6. Phase 4: emit (`include/builder.hpp`, `include/unitig_spill.hpp`)

### 6.1 Cid-range unitig spill (already done by stitch)

`unitig_bucket_writer` (the stitch sink) partitions finished merged
unitigs into K cid-range buckets: bucket `b` holds every unitig with
`cid ∈ [b * S, (b+1) * S)` where `S = ceil(num_color_classes / K)`.
K auto-scales with `-g` (target: per-bucket peak ~10 % of
budget). Per-bucket file format:
`[u32 cid][u32 seq_len][seq bytes]` repeated.

### 6.2 `emit_fasta` + `u2c`

`emit_fasta` walks K buckets in ascending order. For each bucket:

1. Read all records (one bucket's worth fits in memory by design).
2. Sort by `cid`.
3. For each record: write `>cid\n<seq>\n` to `<out>.fa`.
4. Track the running `cid` and a `bits::bit_vector::builder u2c_bvb`
   constructed pre-sized to `num_unitigs` with `init=false` (one
   `n_unitigs / 8`-byte allocation up front, no growth, no
   zero-fill). On each cid-group boundary, set
   `u2c_bvb.set(emitted - 1, 1)` (the previous unitig was the last
   of its run).

After all buckets: set the very last bit (closes the final run),
build the bit_vector, `essentials::save` to `<out>.u2c`.

`u2c` invariants:

- `num_bits == num_unitigs` (one bit per unitig in `.fa` order).
- `popcount == num_color_classes`.
- The last bit is always set.
- Downstream consumer recovers each unitig's color-set id via
  `rank1(unitig_id)` after building a rank index.

### 6.3 `emit_color_sets`

`global_dict.finalize()` runs:

1. Flush the trailing partial word (zero-padded) to the file.
2. Append the sentinel `total_bit_count` to the offsets sidecar.
3. Rewind the sidecar; feed it through `bits::elias_fano::encode` via
   a single-pass file-backed iterator (`offset_file_iterator`). The
   offsets array never sits in RAM as a vector.
4. Serialize the EF onto the file's tail (essentials::generic_saver
   into a stringstream; one final `fwrite` of the bytes).
5. `fseek(0)` and overwrite the placeholder header with the now-known
   totals.
6. Delete the offsets sidecar.

### 6.4 On-disk layout of `<out>.color_sets`

```
[u32 num_colors]
[u32 sparse_threshold]              = floor(0.25 * num_colors)
[u32 dense_threshold]               = floor(0.75 * num_colors)
[u64 num_color_sets]
[u64 bit_vector_num_bits]
[u64 bit_vector_num_words]          = ceil(bit_vector_num_bits / 64)
[bit_vector_num_words × u64]        color-set bit_vector words
[serialized bits::elias_fano<false,false>]   per-class bit-offsets
                                              (last EF value = bit_vector_num_bits)
```

This is intentionally different from the layout
`essentials::save(hybrid)` produced in older revisions (which
serialized the EF before the bit_vector). The reordering is what
allows incremental writing: the bit_vector is appended word-by-word
during `intern()`, the EF is appended at finalize, and only the
fixed-size header at the front is updated via `fseek`.

---

## 7. Hybrid color-set encoding (`include/hybrid_color_sets.hpp`)

This is the bit-level encoder used by both
`compact_color_set_dict` (in-memory, per bucket) and
`streaming_color_set_dict` (on-disk, global). It encodes a single
sorted-deduped color list (`uint32_t` over alphabet `[0, num_colors)`)
into a `bits::bit_vector::builder`, picking one of three regimes
based on the list's size N relative to the thresholds.

For `sparse_threshold = 0.25 * num_colors` and
`dense_threshold = 0.75 * num_colors`:

### 7.1 Sparse regime (N < sparse_threshold)

```
write_delta(N)
write_delta(set[0])
for i in 1..N:
    write_delta(set[i] - set[i-1] - 1)
```

`write_delta` is Elias delta coding (gamma + binary). Storage scales
with `Σ log₂(gap)`. Best when the color list is short and gaps
between consecutive colors are small or large but compressible.

### 7.2 Dense regime (sparse_threshold ≤ N < dense_threshold)

```
write_delta(N)
plain bitmap of num_colors bits, with set[i] positions set
```

Storage = `num_colors` bits + the `write_delta(N)` overhead. Efficient
in the middle regime where deltas would be small (~`num_colors / N`
≈ 1–4) and Elias-delta on small ints uses ~5–7 bits each.

### 7.3 Very-dense regime (N ≥ dense_threshold)

```
write_delta(N)
delta-gaps over the COMPLEMENT (the (num_colors - N) absent values)
```

When N is close to `num_colors`, the absent set is small and sparse-
codes well. Saves bits vs the dense bitmap.

---

## 8. Output file formats

### 8.1 `<out>.fa`

Standard FASTA. Each record:

```
>cid
<seq bases ACGT>
```

Records are emitted in **strictly cid-ascending** order. Multiple
records can share a `cid` (multiple unitigs in one color class) and
appear consecutively.

### 8.2 `<out>.u2c`

One `bits::bit_vector` serialized with `essentials::save`. See §6.2
for invariants and rank semantics.

### 8.3 `<out>.color_sets`

See §6.4 for the byte layout. Decoding one color set:
1. Read the fixed header (40 bytes).
2. mmap or `fread` the bit_vector words.
3. `essentials::load` the elias_fano from the file tail.
4. To get color set `i`: `bit_offset = ef.access(i)`; build a
   `bits::bit_vector::iterator` at that offset; run the inverse of
   `hybrid_builder::encode_one` against it (regime selection depends
   on the leading `write_delta(N)`).

---

## 9. Memory model

Bucket-write peak (§3.7) is modelled as

```
alpha * T * B * flush_bases    ← per-thread per-bucket buffers (alpha = 2)
+ beta * B * spill_bytes       ← per-bucket dedup hashmaps     (beta  = 7)
```

with `T = num_threads`, `B = num_buckets`. Here `flush_bases` and
`spill_bytes` are **fixed**, and the model is solved instead for **B**
against `BUCKET_WRITE_BUDGET_FRAC × -g` (0.50·g). The RSS pressure
watcher (§3.6) is the runtime safety net for the model's ~25 %
real-vs-modelled undercount.

Independently of the RAM model, bucket-write also holds `B` open
`FILE*` handles concurrently (one per bucket file). This is why
`builder::ensure_fd_capacity_for_buckets_` raises `RLIMIT_NOFILE`
to the largest value the OS allows before bucket-write starts and
clamps the bucket *count* down if the soft+hard limit can't accommodate
`B + slack` descriptors (relevant at scale: ~40 k buckets at 661 k/g64).

Bucket-process peak per in-flight thread:

```
record_sets   (compact, hybrid-encoded)   ~10–30 MB / bucket
local_dict    (compact, hybrid-encoded)   ~10–30 MB / bucket
kmer_info     (per-k-mer rsid+flags+owned) ~hundreds of MB / bucket on dense input
merge batch   (decoded color lists, 32)   few MB at num_colors = 100K
```

`record_sets` and `local_dict` were the dominant contributor with
the live-vector dict (~280 MB each per bucket); the compact
hybrid-encoded variant cut this 10–30×. The current bottleneck per
in-flight thread is `kmer_info` (no easy compression — actively used
during the chain walk). The §4.7 batched merge keeps a small
per-thread buffer (32 decoded classes) so the local-dict decode +
hash work can run lock-free.

The **total** across in-flight buckets — plus the carried color dict
and the per-thread merge batches — is what the §4.9 admission gate holds
to `0.82·g`, enforced by a hard live-RSS ceiling so the phase stays in
budget regardless of how well the `24×` per-bucket estimate matches a
given input. The dict (`streaming_color_set_dict`) carries forward into
stitch + emit, so its index is freed (`release_index`) right after this
phase.

Stitch peak: `num_threads` hash buckets' fragments + their boundary-k-mer
hash tables resident at once (`round_store_file`, the `T`-way parallel
round loop). No global per-fragment index or junction map — round 0 seeds
from the frag spill through a `frag_unitig_stream_reader` that holds one
fragment at a time (not a mmap-backed index), so peak scales with the
bucket count, not the total fragment count. Stitch is sized by model
only (no live-RSS ceiling), so it is the phase most likely to graze `-g`
at scale (it was the pipeline peak — 41 GiB — on the 661K run; see
`todo.md`).

Emit peak: one cid-range bucket of records loaded for sorting (~3 MB
on salmonella-25K), plus the u2c bit_vector builder.

The streaming color_set_dict's RAM cost is bounded by per-class
metadata (16 B + ~12 B index entry per class), independent of
compressed bit count.

---

## 10. Concurrency model

| phase | parallelism | sync |
|---|---|---|
| bucket-write | T worker threads, one input file each at a time; per-thread per-bucket buffers; per-bucket mutex guards `bucket_compactor` | per-bucket mutex + RSS watcher's `try_spill` sweep |
| bucket-process | T worker threads pop next bucket from atomic counter; each owns one bucket end-to-end | per-bucket: thread-local; `global_mu` only during local→global merge; `mem_mu`/cv admission gate (§4.9) bounds resident buckets |
| stitch | T worker threads over the per-round bucket loop (buckets independent within a round) | per-bucket waiting-map is thread-local; round barrier between rounds |
| emit | single-threaded | n/a |

Inter-phase: each phase finishes before the next begins. There is no
overlap between phases (e.g., no GGCAT-style background compactor
running concurrently with bucketing).

---

## 11. File index

| file | role |
|---|---|
| `src/main.cpp`                              | CLI entry point |
| `include/builder.hpp`                       | top-level `build()` orchestrator + auto-tune |
| `include/bucket_io.hpp`                     | per-bucket compactor, LZ4 framing, RSS watcher |
| `include/bucket_ingester.hpp`               | parse + minimizer + bucketing per input file |
| `include/bucket_walker.hpp`                 | per-bucket dBG load + walk; multi-thread driver |
| `include/stitch_extmem.hpp`                 | cross-bucket open-end joining (external-memory) |
| `include/stitch.hpp`                        | shared stitch helpers (side tags, junction, frag source) |
| `include/unitig_spill.hpp`                  | disk-backed frag/unitig sinks + mmap reader |
| `include/super_kmer.hpp`                    | super-k-mer record format (varint + 2-bit) |
| `include/streaming_color_set_dict.hpp`      | global color-set dict; writes `.color_sets` |
| `include/compact_color_set_dict.hpp`        | per-bucket color-set dict (hybrid in-memory) |
| `include/hybrid_color_sets.hpp`             | static `encode_one` (sparse/dense/complementary) |
| `include/kmer.hpp`                          | 2-bit canonical k-mer encoding |
| `include/minimizer.hpp`                     | canonical ntHash + sliding-window minimum |
| `include/seq_reader.hpp`                    | mmap + libdeflate FASTA/FASTQ iterator (kseq over mem_stream) |
| `include/util.hpp`                          | timers, RSS, profiling counters, build_config |
