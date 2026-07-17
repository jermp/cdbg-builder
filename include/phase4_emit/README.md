# Phase 4 — emit

Turn the stitch's cid-range unitig buckets into the final on-disk outputs.
Walk the K buckets in ascending order; within each, sort by cid (so `<out>.fa`
is strictly cid-ascending) and write `>cid\n<seq>\n`. In the same pass build
the `<out>.u2c` bit_vector (bit *i* set iff unitig *i* is the last of a
color-set run). Finally finalize `<out>.color_sets`, whose bit_vector was
already streamed out word-by-word during Phase 2 — finalize appends the
Elias–Fano of per-class bit-offsets and writes the header totals.

The unitig buckets are **RAM-first**: `unitig_bucket_writer` keeps them in
memory up to a `-g`-derived budget, so emit reads the unitigs straight from
RAM and the sequences never round-trip through disk (no temp-bucket write in
stitch, no read-back here). Only buckets that exceed the budget spill to
`tmp/unitig_bucket_<k>.bin` and are external-merge-sorted at emit.

See `algorithm.md` §6 for the full description.

## Input

- the in-RAM cid-range unitig buckets (the `unitig_bucket_writer` carried over
  from stitch), plus any `tmp/unitig_bucket_<k>.bin` that spilled under `-g`.
- the in-progress `<out>.color_sets` (to finalize) and the running
  `streaming_color_set_dict`.

## Output

- `<out>.fa` — colored unitigs in FASTA, cid-ascending, headers = cid.
- `<out>.u2c` — unitig→color-set `bits::bit_vector` (run-end markers;
  popcount = `num_color_sets`).
- `<out>.color_sets` — finalized: header + hybrid-encoded color sets + EF offsets.

These three `<out>.*` files are the tool's deliverables (consumed downstream by
Fulgor). All `tmp/*` scratch is removed after this phase.

## Files

| file       | role                                                                                                      |
|------------|-----------------------------------------------------------------------------------------------------------|
| `emit.hpp` | `emit_fasta` (FASTA + u2c) and `emit_colors` (finalize `.color_sets`), called by the builder orchestrator |

> The top-level `builder.hpp` orchestrator drives all four phases and lives at
> `include/` root; only the emit step's code lives here.
