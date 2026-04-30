# test_data

End-to-end correctness test for `cdgb-build`.

## Layout

- `salmonella_10/*.fasta.gz` — 10 *Salmonella enterica* genome assemblies.
  These are the same files as in
  [Fulgor's `test_data/salmonella_10/`](https://github.com/jermp/fulgor/tree/main/test_data/salmonella_10),
  used here as 10 input colors (color id = position in the filenames list).
- `verify.py` — recomputes the expected canonical-k-mer → color-bitmask map
  from the inputs and validates the produced `<out>.fa` and `<out>.colors`.
  K-mers are 2-bit packed into Python ints so the verifier scales to the
  ~7 M unique canonical k-mers in this dataset.
- `run_test.sh` — generates a filenames list pointing at the gzipped
  inputs, runs `cdgb-build` with `k = 31`, then runs `verify.py`.

## Running

After building the project (so `build/cdgb-build` exists):

```bash
./test_data/run_test.sh             # uses 4 threads
THREADS=8 ./test_data/run_test.sh   # override
```

The script exits 0 on success and prints `OK: ...` with summary counts.

## What is checked

For every input k-mer the verifier asserts:

- `<out>.fa` and `<out>.colors` exist and `<out>.colors` is non-empty.
- Every canonical input k-mer appears in exactly one output unitig.
- Every k-mer inside a unitig maps to the same input color set.
- All unitigs sharing a `color_set_id` (FASTA header) carry the same
  color set.
- The number of distinct `color_set_id`s equals the number of distinct
  color sets in the inputs, and the sets themselves match.

On the bundled salmonella_10 inputs the build produces around 6.9 M
canonical k-mers across 171 distinct color sets.
