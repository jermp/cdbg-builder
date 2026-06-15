# Phase 3 — stitch

Glue the open-ended fragments from Phase 2 into maximal unitigs across bucket
boundaries, by matching the **full boundary k-mer** the two sides share. Uses
GGCAT-style hash-bucketed iterative doubling in a **base-carrying** form: each
round carries the growing unitig sequences (and their color runs) and grows them
in place as fragments join. Internal steps (logged as `seed` / `rounds`) are
sub-steps of this single pipeline phase, **not** the four top-level phases. The
round store is always on disk (one bucket's tigs per in-flight thread resident),
so the phase honors `-g` by keeping little in RAM rather than by spilling on
demand.

See `algorithm.md` §5 for the full description.

## Input
- `tmp/frag_unitigs.bin` — streamed once (round-0 seed), carrying bases + colors.

## Output
- `tmp/unitig_bucket_<k>.bin`, `k ∈ [0, K)` — **finished, monochromatic**
  unitigs `{ACGT seq, cid}`, partitioned so bucket `k` holds every unitig whose
  `cid` falls in `k`'s range (cid-range bucketing lets Phase 4 emit in strict
  cid order without a global sort).

## Files
| file | role |
|---|---|
| `stitch_extmem.hpp` | **production stitch**: base-carrying external-memory streaming doubling (`stitch_unitigs_extmem_file_streaming`); the random-access `stitch_unitigs_extmem` / `_file` variants are kept as `test_stitch`'s independent reference oracles |
| `stitch.hpp` | shared stitch helpers (side tags, junction k-mers, frag source) |
