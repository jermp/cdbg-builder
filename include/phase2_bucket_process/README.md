# Phase 2 — bucket-process

For each Phase-1 bucket **independently**: build the local de Bruijn graph over
canonical k-mers (with cross-bucket boundary flags, the GGCAT begin/end-ignored
model), walk maximal topological **fragments** carrying an RLE color-run
sequence, and intern color lists — first into a per-bucket dict, then merged
into the shared global color-set dict streamed straight to `<out>.color_sets`.
Buckets are processed in parallel; a memory-admission gate (hard live-RSS
ceiling) bounds how many load at once so the phase stays within `-g`.

See `algorithm.md` §4 for the full description.

## Input
- `tmp/bucket_<b>.bin`, `b ∈ [0, B)` (one bucket per worker).

## Output
- `tmp/frag_unitigs.bin` — one stream of **open-ended fragments**
  `{ACGT seq, open_flags, color-run sequence}` with **global** cids.
- `tmp/frag_unitigs.bin.links` — companion **links spill**: one fixed-size
  `[u8 open_flags][kbytes kl][kbytes kr]` per fragment, in `frag_id` order, so
  Phase 3 can seed from ~17 B/fragment instead of re-reading the frag spill.
- `<out>.color_sets` — distinct **global color classes**, interned and streamed
  to disk incrementally during this phase (finalized in Phase 4).

## Files
| file | role |
|---|---|
| `bucket_walker.hpp` | per-bucket dBG load + walk; multi-thread driver + admission gate |
| `compact_color_set_dict.hpp` | per-bucket color-set dict (hybrid in-memory) |
| `streaming_color_set_dict.hpp` | global color-set dict; writes `<out>.color_sets` (finalized in Phase 4) |
| `hybrid_color_sets.hpp` | static `encode_one` (sparse/dense/complementary) — shared color encoder |
| `unitig_spill.hpp` | disk-backed frag/unitig sinks (batched write + companion links spill) + mmap reader |

> `unitig_spill.hpp` is the disk-handoff plumbing for the 2→3→4 flow (it also
> holds the `unitig_bucket_writer` Phase 3 fills and Phase 4 reads); it lives
> here as the phase that first produces the frag + links spills.
