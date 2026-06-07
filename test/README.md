# Tests

`test_stitch` exercises the stitch phase in isolation on synthetic,
seeded random data using a **generate-by-construction** oracle: build a
known-correct answer, derive the phase's input from it, run the phase,
and assert it reconstructs the answer. This verifies correctness (and,
with the larger case, scaling) independently of the end-to-end
`test_data/run_test.sh` flow.

## Building

The test builds by default with the main project:

```bash
mkdir -p build && cd build && cmake .. && make -j
```

To skip it: `cmake -DCDGB_BUILD_TESTS=OFF ..`.

## Running

```bash
bash test/run_tests.sh
```

Or run the binary directly: `build/test_stitch`. It returns 0 on
success, non-zero on failure, and prints a `[FAIL] ...` line per failing
case with the seed (reproducible).

## What it checks

| executable | phase | oracle |
|---|---|---|
| `test_stitch` | stitch | split K random "true" unitigs (distinct cids) into overlapping fragments with OPEN internal boundaries, shuffle + randomly revcomp, stitch, and assert the K unitigs are reconstructed exactly (modulo orientation), each fully closed |

It runs every stitch implementation against an independent oracle: the
production id-only compaction stitch (`compact_mem` / `compact_file` /
`compact_scalable` / `compact_scalable_links` / `compact_inram`, in
`compact_extmem.hpp`) is cross-checked against the base-carrying
external-memory reference (`ext_mem` / `ext_file`, in `stitch_extmem.hpp`,
kept solely as this oracle). The `[scale]` case (50k unitigs) also gives
an isolated per-phase timing number.

`gen.hpp` holds the shared generators (`random_dna`, `revcomp`,
`canonical`, `split_unitig`).
