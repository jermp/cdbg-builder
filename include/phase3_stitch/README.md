# Phase 3 — stitch

Glue the open-ended fragments from Phase 2 into maximal unitigs across bucket
boundaries, by matching the **full boundary k-mer** the two sides share. Uses
GGCAT-style hash-bucketed iterative doubling in an **id-only** form: the rounds
carry only fragment-**id** chains (no bases, no colors), so they are cheap; the
bases+colors of each closed chain are assembled exactly **once** at the end.
Internal steps (logged as `seed` / `rounds` / `phase2` / `phase3`) are
sub-steps of this single pipeline phase, **not** the four top-level phases. All
working stores are RAM-first with spill-to-disk overflow, so the phase honors
`-g` at all costs.

See `algorithm.md` §5 for the full description.

## Input
- `tmp/frag_unitigs.bin.links` — to seed round 0 (boundary k-mers).
- `tmp/frag_unitigs.bin` — streamed once, to attach bases+colors during assembly.

## Output
- `tmp/unitig_bucket_<k>.bin`, `k ∈ [0, K)` — **finished, monochromatic**
  unitigs `{ACGT seq, cid}`, partitioned so bucket `k` holds every unitig whose
  `cid` falls in `k`'s range (cid-range bucketing lets Phase 4 emit in strict
  cid order without a global sort).

## Files
| file | role |
|---|---|
| `compact_extmem.hpp` | **production stitch**: id-only doubling join + base/color assembly (RAM-first, spill-to-disk) |
| `stitch.hpp` | shared stitch helpers (side tags, junction k-mers, frag source) |
| `stitch_extmem.hpp` | base-carrying external-memory stitch; **retired from the build**, kept as `test_stitch`'s independent reference oracle |
