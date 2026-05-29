# Per-phase tests

Each test exercises one pipeline phase in isolation on synthetic,
seeded random data, using a **generate-by-construction** oracle: build
a known-correct answer, derive the phase's input from it, run the
phase, and assert it reconstructs the answer. This verifies both
correctness and (with the larger cases) scaling, independently of the
end-to-end `test_data/run_test.sh` flow.

## Building

Tests build by default with the main project:

```bash
mkdir -p build && cd build && cmake .. && make -j
```

To skip them: `cmake -DCDGB_BUILD_TESTS=OFF ..`.

## Running

```bash
bash test/run_tests.sh
```

Or run an individual test directly: `build/test/test_stitch`. Each
returns 0 on success, non-zero on failure, and prints a `[FAIL] ...`
line per failing case with the seed (reproducible).

## Tests

| executable | phase | oracle |
|---|---|---|
| `test_stitch` | stitch | split K random "true" unitigs (distinct cids) into overlapping fragments with OPEN internal boundaries, shuffle + randomly revcomp, stitch, and assert the K unitigs are reconstructed exactly (modulo orientation), each fully closed |

`gen.hpp` holds the shared generators (`random_dna`, `revcomp`,
`canonical`, `split_unitig`).

## Why this matters for the stitch redesign

`test_stitch` is the regression guard for the upcoming external-memory
stitch (Phase 2). The in-RAM stitch is the baseline; the external-
memory version must produce identical output on every seed. The
`[scale]` case (50k unitigs) also gives a per-phase timing/throughput
number that's isolated from the rest of the pipeline.

## Planned additions

- `test_bucket_write` — synth reads with known minimizers; check
  super-k-mers land in the right buckets and dedup correctly.
- `test_bucket_process` — synth bucket files encoding a known per-bucket
  dBG; check emitted fragments + phantom-edge open flags.
- `test_emit` — synth cid-range unitig buckets; check `.fa` ordering and
  `.u2c` rank semantics.
