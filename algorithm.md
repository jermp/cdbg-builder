# Algorithm

This document describes how `cdbg-builder` constructs a **colored compacted
de Bruijn graph (ccdBG)** and how it builds and compresses the color sets.

---

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

---

## 2. Pipeline overview

The build runs in four sequential phases. The wall-time numbers come
from the salmonella-25K benchmark (25 000 input genomes, k = 31, m = 12,
16 threads, --max-ram 4 GB) and are representative of dense pangenome
inputs.

| phase | wall (25K) | what it does | files written |
|---|---|---|---|
| **bucket-write** | ~296 s | partition input k-mers into B disk buckets | `tmp/bucket_*.bin` (LZ4) |
| **bucket-process** | ~419 s | per-bucket dBG + walk; emit fragments | `tmp/frag_unitigs.bin` ; `out.color_sets` (streamed) |
| **stitch** | ~16 s  | join open-ended fragments across buckets | `tmp/unitig_bucket_*.bin` |
| **emit** | ~4 s   | write `.fa` + `.u2c`; finalize `.color_sets` | `out.fa`, `out.u2c`, `out.color_sets` |

The four phases are strictly serial: each consumes the prior phase's
on-disk output and finishes before the next starts. Within a phase,
work is parallelized over either input files (bucket-write) or buckets
(bucket-process). Stitch and emit are single-threaded.

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
    │   join open ends across buckets via (k-1)-mer match
    ▼
K cid-range unitig_bucket files                           ← stitch
    │   sort by cid, write FASTA, build u2c
    ▼
out.fa + out.u2c + finalize out.color_sets                ← emit
```

---

## 3. Phase 1: bucket-write (`include/bucket_io.hpp`, `bucket_ingester.hpp`)

### 3.1 Goal

Read every input file, decompose each ACGT-only run into super-k-mers,
compute each super-k-mer's **canonical minimizer**, and append it to the
bucket file selected by that minimizer. After this phase:

- Each bucket file holds one stream of super-k-mer records.
- Records are *compacted*: a super-k-mer that occurred in many input
  files is stored once with the union of its colors.
- Every k-mer that lives in bucket `b` is reachable from records in
  bucket `b` alone (by construction of the minimizer-bucketing).

### 3.2 Minimizer bucketing

Within an ACGT run of length `L ≥ k`, slide a window of size
`W = k − m + 1` m-mers. The current k-mer's minimizer is the m-mer
whose `canonical_mhash` (canonical ntHash; symmetric under RC) is
minimum across the window. As the window slides, the minimum
sometimes changes; each contiguous span of k-mers sharing one
minimizer is a **super-k-mer**. The super-k-mer's bucket id is
`(min_hash >> 1) & (B − 1)`. The low bit of the hash is reserved
to match GGCAT's `cn_nthash.rs` convention.

Each super-k-mer record carries:

- `flags` : `IS_ACGT_BEGIN` (this super-k-mer's first k-mer has no
  predecessor in input), `IS_ACGT_END` (analogously). The bucket
  walker uses these to insert *phantom* edges across buckets (§4.4).
- `colors`: at write time, the single color of this input file.
  Compaction merges color lists across files later in this phase.
- `bases` : 2-bit-packed bases.

### 3.3 Per-thread per-bucket buffers and flush

Each ingest worker thread owns a `per_thread_bucket_buffers` (one
`std::vector<pending_record>` and one `std::vector<uint8_t>` of bases
per bucket). `append()` pushes the record's metadata + appends the
bases bytes — no lock, no hashmap. When a bucket's bases-buffer
crosses `flush_bases` (auto-tuned, typically 16–64 KiB), the thread
acquires that bucket's mutex and calls
`bucket_compactor::insert_batch`.

### 3.4 Per-bucket compactor (online dedup)

Each `bucket_compactor` (one per bucket) holds:

- `m_dedup`: an `ankerl::unordered_dense::map<std::string, entry>`
  keyed by the 2-bit-packed super-k-mer bytes, with `entry =
  { colors, flags }`.
- `m_batch_buf`, `m_out_buf`: scratch vectors reused across spills.
- `m_file`: the open bucket file handle.

`insert_batch` walks the incoming batch, transparent-finds each
super-k-mer in `m_dedup` (no string allocation on hit), and either
inserts a new entry or appends the new color to an existing one's
`colors` vector. Sort+unique on `colors` is deferred until spill.

### 3.5 Spill: serialize, LZ4-compress, write a frame

When `m_dedup`'s tracked-byte counter crosses `spill_bytes`
(auto-tuned, typically 16–256 KiB), `spill_locked` runs:

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
every 100 ms. On HIGH (default 60 % of `--max-ram`), it sweeps every
compactor and force-spills any with pending state. When RSS drops to
LOW (default 45 %), pressure clears and ingest threads return to
their normal `m_bytes`-based spill cadence. Hysteresis prevents the
death-spiral observed with sticky pressure (every batch spilling for
the rest of the phase, exploding spill counts and bucket file sizes).

If `/proc/self/status` is unavailable (macOS, sandbox), the watcher
falls back to monotonic `getrusage` peak; pressure becomes sticky once
tripped (less efficient but still correct).

### 3.7 Auto-tune

Three knobs (`auto_bucket_log2`, `auto_tune_bucket_write_params`):

- `num_buckets = 2^bucket_log2`. Picked as the largest B such that
  `B × TARGET_SPILL_BYTES × COMPACTOR_OVERHEAD ≤ (BUCKET_WRITE_SHARE
  × max_ram) / 2`. Forces each bucket to afford a spill big enough
  for dedup to pay off. For 4 GB / 16 threads typical = 2048;
  scales with budget.
- `flush_bases`, `spill_bytes`. Joint-tuned against the same share,
  splitting it 50/50 between per-thread buffers and per-bucket
  compactor data. Empirical overhead factors:
  `BUFFER_OVERHEAD = 2x` (vector capacity slack), `COMPACTOR_OVERHEAD
  = 7x` (string allocs + map structure + glibc fragmentation).
- `BUCKET_WRITE_SHARE`: 0.50 on Linux, 0.40 on macOS (libsystem_malloc
  retains pages glibc would reclaim, so the auto-tune budgets less).

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
   awareness of *phantom* edges that cross into neighboring buckets.
2. Walk maximal unitig fragments. Each fragment is monochromatic.
3. Intern its color list into a per-bucket `compact_color_set_dict`,
   then merge per-bucket dicts into a shared
   `streaming_color_set_dict` whose output goes directly to
   `<out>.color_sets`.

Each bucket's fragments are emitted as records in a single
disk-backed `frag_unitig_writer` for stitch to consume.

### 4.2 Per-bucket walker state

Three structures per in-flight bucket:

- **`kmer_info`** (`bucket_kmer_map`): canonical-k-mer → `{kmer_entry,
  phantom flags}`. `kmer_entry` is `{first_rsid, vector<rest_rsid>}`
  — almost every k-mer is contributed by a single super-k-mer record,
  so `rest` is empty in the common case (no heap alloc per k-mer).
- **`record_sets`** (`compact_color_set_dict`): interns the sorted-
  deduped color list of every super-k-mer record read from disk.
  Hybrid-encoded into an in-memory `bit_vector`; the dict gives back
  decoded color lists into a caller-supplied scratch buffer
  (`at(rsid, scratch)`). Memory per class: 24 B metadata +
  hybrid-encoded bits.
- **`local_dict`** (`compact_color_set_dict`): per-bucket final color
  classes. Each unique color list (single-rsid k-mer's colors, or
  the union of multi-rsid k-mer's color lists) is interned here and
  the returned local cid is stored as the unitig fragment's `cid`.

### 4.3 Loading a bucket

`load_bucket`:

1. Open the bucket file, decompress LZ4 frames, parse `super_kmer`
   records.
2. For each record: intern colors into `record_sets` → get rsid.
3. Roll the canonical k-mer through the record's bases. Each k-mer
   gets the rsid added to its `kmer_entry` via `add()`.
4. Set phantom-edge flags on the first/last k-mer if the record's
   `IS_ACGT_BEGIN`/`IS_ACGT_END` bits are NOT set (i.e., the
   super-k-mer continues into another bucket).

### 4.4 Phantom edges

Each k-mer carries two bits (`KMER_PHANTOM_LEFT/RIGHT`) saying
whether a known dBG edge crosses into a different bucket. The walker
treats a phantom edge as having degree 1 on that side: if the only
edge on a side is phantom, the walk stops there and the unitig is
emitted with that end **OPEN**, ready for the stitch phase to glue
it. If the side mixes a local edge and a phantom edge, the global
degree is ≥ 2 and the walk stops as a real branch (CLOSED).

This matches GGCAT's `try_extend_function` logic. Closed unitigs go
straight to the sink. Open unitigs carry their first and last
canonical k-mers and side-flags so stitch can match them.

### 4.5 Computing per-k-mer cids

After load, two passes:

**Pass A**: build `cid_of[k-mer]`. For each k-mer in `kmer_info`:
- *Single-rsid* (common): `local_dict.intern(record_sets.at(rsid))`.
  A small `rsid_to_cid` cache makes this O(1) on subsequent k-mers
  with the same rsid.
- *Multi-rsid*: sort+unique the rsids; if one distinct rsid remains
  use the cached path; else union the referenced color lists into a
  fresh sorted set and intern that.

After this, `record_sets` and the per-k-mer rsid storage are
released; the walker only needs `kmer_info`'s phantom bits + the
`cid_of` map for the chain walk.

**Pass B**: the walk emits `stitchable_unitig` records with
`{seq, cid, open_flags}`. Local cids will be remapped to global cids
in §4.7.

### 4.6 Streaming the global color-set dict

A single `streaming_color_set_dict` is shared across all bucket
threads, guarded by `global_mu`.

`intern(colors)`:

1. Compute primary (`wyhash`) + secondary (`fnv1a`) 128-bit hash.
2. Look up in `m_index` (hash-only key, no byte-compare). On hit
   return existing id; collision rate is ~2^-64 per pair, negligible
   over our class counts.
3. On miss: hybrid-encode the new class's bits via
   `hybrid_builder::encode_one` into an in-memory
   `bits::bit_vector::builder`.
4. Append every COMPLETE 64-bit word from the builder to
   `<out>.color_sets` (the file is opened at construction with a
   placeholder header). Only the trailing partial word stays in
   memory.
5. Append the new class's `bit_offset` to a sidecar file
   `<out>.color_sets.tmp_offsets` for the EF later.

Per-class metadata held in RAM: 16 B (primary + secondary hash). The
compressed bit_vector itself is never resident.

### 4.7 Per-bucket → global merge

After `process_bucket` returns, the worker thread takes `global_mu`
and walks its `local_dict`:

```cpp
for each local class lc:
    local_dict.at(lc, decode_scratch);
    local_to_global[lc] = global_dict.intern(std::move(decode_scratch));
```

Then it remaps each unitig fragment's `cid` from local → global and
streams the fragment to the disk-backed `frag_unitig_writer`. Per-
fragment: `[u32 cid][u8 open_flags][u32 seq_len][seq bytes]`.

### 4.8 Auto-cap on concurrency (none today)

The user explicitly asks for `--threads N`. bucket-process spawns N
worker threads; each pops the next bucket from a shared atomic
counter. There is **no** auto-cap — if the per-thread walker state
won't fit in `--max-ram / N`, the process simply runs over budget. A
clean abort with a "lower --threads or raise --max-ram" message is
listed as future work in §10.

---

## 5. Phase 3: stitch (`include/stitch.hpp`)

Single-threaded. Takes a `frag_unitig_reader` (mmap-backed view over
`tmp/frag_unitigs.bin`) and a `Sink` (the cid-range
`unitig_bucket_writer` for emit).

**Pass 1 — junction indexing**: for every fragment with an open
side, hash the canonical (k−1)-mer at that side. Insert into a
`by_junction` map keyed by junction kmer; the map's value holds up to
two open-end refs per junction (overflow rejected).

**Pass 2 — adjacency build**: for each junction with exactly two
open ends and matching `cid` (gating on the cheap u32 compare we
established at process_bucket time), record the bidirectional link
in `adj` (one packed `uint64_t` per fragment, low/high u32 = LEFT/
RIGHT side encoded as `(other_idx << 1) | other_side`, sentinel
`UINT32_MAX` for "no link"). 8 bytes per fragment.

**Pass 3 — walk chains**: for each unvisited fragment, start a
chain at a free side. Repeatedly hop via `adj`, copying the next
fragment's seq from the mmap (revcomped if the chain orientation
flips). At each hop, append `seq[k-1..]` to the merged unitig (the
shared k-1 prefix is already covered). Emit the merged unitig to the
sink with its propagated `cid`.

mmap'd seq pages are accessed in chain order — random with respect
to the on-disk layout. `MADV_RANDOM` tells the kernel to skip
readahead and reclaim pages aggressively under pressure. On SSD
this is acceptable; on HDD it would be slow.

---

## 6. Phase 4: emit (`include/builder.hpp`, `include/unitig_spill.hpp`)

### 6.1 Cid-range unitig spill (already done by stitch)

`unitig_bucket_writer` (the stitch sink) partitions finished merged
unitigs into K cid-range buckets: bucket `b` holds every unitig with
`cid ∈ [b * S, (b+1) * S)` where `S = ceil(num_color_classes / K)`.
K auto-scales with `--max-ram` (target: per-bucket peak ~10 % of
budget). Per-bucket file format:
`[u32 cid][u32 seq_len][seq bytes]` repeated.

### 6.2 `emit_fasta` + `u2c`

`emit_fasta` walks K buckets in ascending order. For each bucket:

1. Read all records (one bucket's worth fits in memory by design).
2. Sort by `cid`.
3. For each record: write `>cid\n<seq>\n` to `<out>.fa`.
4. Track the running `cid` and a `bits::bit_vector::builder u2c_bvb`:
   on each cid-group boundary, set `u2c_bvb.set(emitted - 1, 1)`
   (the previous unitig was the last of its run).

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

Bucket-write peak (planned by `auto_tune_bucket_write_params`) is:

```
T * B * flush_bases * BUFFER_OVERHEAD       ← per-thread per-bucket buffers
+ B * spill_bytes * COMPACTOR_OVERHEAD      ← per-bucket dedup hashmaps
+ B * COMPRESSOR_BYTES_PER_BUCKET           ← small constant per bucket
```

with `T = num_threads`, `B = num_buckets`, BUFFER_OVERHEAD = 2x,
COMPACTOR_OVERHEAD = 7x. Auto-tune solves this for `flush_bases` and
`spill_bytes` against `BUCKET_WRITE_SHARE × --max-ram`. The RSS
pressure watcher provides a runtime safety net.

Bucket-process peak per in-flight thread:

```
record_sets   (compact, hybrid-encoded)   ~10–30 MB / bucket
local_dict    (compact, hybrid-encoded)   ~10–30 MB / bucket
kmer_info     (per-k-mer rsid + phantom)  ~hundreds of MB / bucket on dense input
```

`record_sets` and `local_dict` were the dominant contributor with
the live-vector dict (~280 MB each per bucket); the compact
hybrid-encoded variant cut this 10–30×. The current bottleneck per
in-flight thread is `kmer_info` (no easy compression — actively used
during the chain walk).

Stitch peak: `frag_unitig_reader` index (~16 B per fragment) +
`by_junction` + `adj` (8 B per fragment, packed). Mmap pages of seq
bytes are accessed with `MADV_RANDOM` and counted in RSS only as
they're touched.

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
| bucket-process | T worker threads pop next bucket from atomic counter; each owns one bucket end-to-end | per-bucket: thread-local; `global_mu` only during local→global merge |
| stitch | single-threaded (the in-memory adjacency and walk are not yet parallelized) | n/a |
| emit | single-threaded | n/a |

Inter-phase: each phase finishes before the next begins. There is no
overlap between phases (e.g., no GGCAT-style background compactor
running concurrently with bucketing).

---

## 11. Known weaknesses and optimization surfaces

### 11.1 bucket-write CPU dominance

bucket-write spends ~205 s of its 296 s wall in `flush` (139 s
hashmap, 64 s spill+LZ4, 2 s lock_wait). At 16 threads, the hashmap
cost is ~95 % of theoretical perfect parallelism, so it's already
well-distributed — the cost is per-record CPU, not coordination.

Realistic levers, ordered by experimental value over effort:
- Try `mimalloc` / `jemalloc` via `LD_PRELOAD`. ~2 B std::string
  allocations through bucket-write. Could buy 5–10 % wall.
- Decouple gz-decode from minimizer-compute (separate thread pool).
  Today gz-decode blocks on the same worker that does compute;
  a producer/consumer split would let the 22 s of gz overlap.
- Try a different hashmap (`folly::F14`, `phmap::flat_hash_map`).
  Modest, ~10–20 % per find.
- Vectorize `pack_2bit` (BMI2 `pext`) and `nthash_roll` (SIMD).
  Modest; pack_2bit currently runs only on dedup'd records so the
  amortized cost is small.

The architectural lever (GGCAT-style background compactor that
removes dedup from the bucketing critical path) was tried on the
`claude/append-only-bucketing` branch and rejected: it inflated the
per-record encoding cost ~6× because raw records have to be
varint-encoded and 2-bit-packed even when they would have deduped
away.

### 11.2 No bucket-process concurrency cap

If `--threads N` × per-bucket walker state exceeds `--max-ram`, the
process runs over budget. The right behavior is a clean abort with a
"lower --threads or raise --max-ram" message; this hasn't been
implemented.

### 11.3 Single-threaded stitch and emit

Stitch is ~16 s on the 25K workload (small relative to the other
phases) but would scale poorly if num_unitigs grew much further.
Emit is ~4 s and IO-bound. Both are unparallelized today; not a
priority.

### 11.4 No SIMD in the hot inner loops

`pack_2bit`, `nthash_roll`, `canonical_mhash` are all scalar.
Realistic ~2–4× speedup with SSE/AVX/BMI2. Worth attempting once
allocator and decoupled gz-decode wins are exhausted.

### 11.5 Per-bucket walker `kmer_info` is the largest residual

`compact_color_set_dict` cut `record_sets` and `local_dict` ~10×.
What remains is `kmer_info`, which is the per-canonical-k-mer
hashtable used during the dBG walk. Compressing it would require
either (a) restructuring the walker to be a streaming-online
algorithm rather than load-then-walk, or (b) accepting some CPU
cost to encode the rsid lists more compactly. Both are substantial.

### 11.6 mmap-backed stitch is SSD-friendly, HDD-painful

`MADV_RANDOM` keeps the resident set bounded but means each
fragment access can be a page fault. On NVMe SSDs ~50 µs per
fragment; on HDD this would be orders of magnitude slower.
Documented but not auto-detected.

### 11.7 Wall-time gap to GGCAT (~17 % on 25K)

We're slower than GGCAT on bucketing (296 s vs 81 s) but **faster**
on bucket-process (419 s vs 532 s, GGCAT's `kmers merge`). Net gap
is ~128 s. Closing it requires the background-compactor architecture
plus the micro-optimizations above; current peak-RSS-and-cap
behavior is already strong.

---

## 12. File index

| file | role |
|---|---|
| `src/main.cpp`                              | CLI entry point |
| `include/builder.hpp`                       | top-level `build()` orchestrator + auto-tune |
| `include/bucket_io.hpp`                     | per-bucket compactor, LZ4 framing, RSS watcher |
| `include/bucket_ingester.hpp`               | parse + minimizer + bucketing per input file |
| `include/bucket_walker.hpp`                 | per-bucket dBG load + walk; multi-thread driver |
| `include/stitch.hpp`                        | cross-bucket open-end joining |
| `include/unitig_spill.hpp`                  | disk-backed frag/unitig sinks + mmap reader |
| `include/super_kmer.hpp`                    | super-k-mer record format (varint + 2-bit) |
| `include/streaming_color_set_dict.hpp`      | global color-set dict; writes `.color_sets` |
| `include/compact_color_set_dict.hpp`        | per-bucket color-set dict (hybrid in-memory) |
| `include/hybrid_color_sets.hpp`             | static `encode_one` (sparse/dense/complementary) |
| `include/kmer.hpp`                          | 2-bit canonical k-mer encoding |
| `include/minimizer.hpp`                     | canonical ntHash + sliding-window minimum |
| `include/seq_reader.hpp`                    | gzip FASTA/FASTQ iterator (kseq.h) |
| `include/util.hpp`                          | timers, RSS, profiling counters, build_config |
