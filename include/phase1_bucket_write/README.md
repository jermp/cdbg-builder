# Phase 1 — bucket-write

Read every input file, decompose each ACGT run into super-k-mers, compute each
super-k-mer's **canonical minimizer**, and append it to the bucket file that
minimizer selects. Minimizer bucketing over **(k−1)-mers** (the BCALM2/GGCAT
rule) guarantees every dBG edge and branch lands co-located in a single bucket,
so each bucket file becomes an independent dBG sub-problem for Phase 2. A
per-bucket online compactor dedups super-k-mers that recur across input files
(unioning their colors) before they hit disk.

See `algorithm.md` §3 for the full description.

## Input
- The `N` input files (FASTA/FASTQ, optionally gzipped). File `i` = **color** `i`.

## Output
- `tmp/bucket_<b>.bin`, `b ∈ [0, B)` — one LZ4-framed file per bucket, each a
  stream of compacted `super_kmer` records `{2-bit bases, 4-bit flags, color
  list}`. Colors are local file indices, deduped within the bucket.

## Files
| file | role |
|---|---|
| `bucket_ingester.hpp` | parse + minimizer + bucketing per input file |
| `bucket_io.hpp` | per-bucket compactor (online dedup), LZ4 framing, RSS-pressure watcher |
| `minimizer.hpp` | canonical ntHash + sliding-window minimum |
| `seq_reader.hpp` | mmap + libdeflate FASTA/FASTQ iterator |
| `super_kmer.hpp` | super-k-mer record format (varint + 2-bit) — this phase's output record |

> The foundational primitives used across all phases — `kmer.hpp` and
> `util.hpp` — live at `include/` root next to the `builder.hpp` orchestrator,
> not in any one phase folder.
