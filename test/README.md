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

It runs every stitch variant against an independent oracle: the production
base-carrying streaming stitch (`ext_stream_ram` / `ext_stream_spill`,
`stitch_unitigs_extmem_file_streaming` in `stitch_extmem.hpp`) is cross-checked
against the random-access external-memory references (`ext_mem` / `ext_file` /
`ext_file_mt`, in `stitch_extmem.hpp`, kept solely as this oracle). The
`[scale]` case (50k unitigs) also gives an isolated per-phase timing number.

`gen.hpp` holds the shared generators (`random_dna`, `revcomp`,
`canonical`, `split_unitig`).
