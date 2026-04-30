# test_data

A tiny end-to-end correctness test for `cdgb-build`.

## Layout

- `inputs/a.fa`, `inputs/b.fa`, `inputs/c.fa` — three reference files used as
  three input colors (color id = position in the filenames list).
- `verify.py` — recomputes the expected canonical-k-mer → {color ids} map
  from the same inputs and validates the produced `<out>.fa` and
  `<out>.colors` files.
- `run_test.sh` — generates a filenames list pointing at `inputs/`, runs
  `cdgb-build`, then runs `verify.py`.

## Running

After building the project (so `build/cdgb-build` exists):

```bash
./test_data/run_test.sh
```

The script exits 0 on success and prints `OK: ...` with summary counts.

## What is checked

For `k = 5` with the bundled inputs the verifier asserts:

- `<out>.fa` and `<out>.colors` exist and `<out>.colors` is non-empty.
- Every canonical input k-mer appears in exactly one output unitig.
- Every k-mer inside a unitig maps to the same input color set.
- All unitigs sharing a `color_set_id` (FASTA header) carry the same
  color set.
- The number of distinct `color_set_id`s equals the number of distinct
  color sets in the inputs, and the sets themselves match.

The expected counts on the bundled data are 4 canonical k-mers across
3 distinct color sets (`{0,1}`, `{1}`, `{2}`).
