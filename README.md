# cdbg-builder

A small, multi-threaded tool that builds a **colored compacted de Bruijn
graph** from a collection of reference sequences and emits the result in a
format compatible with [Fulgor](https://github.com/jermp/fulgor)'s `hybrid`
color-set encoding.

Given `N` input files (each treated as one *color*), the tool produces:

- `<out>.fa` — a FASTA of *monochromatic colored unitigs*. Each header is the
  integer id of the color set the unitig belongs to.
- `<out>.colors` — the color sets, serialized in Fulgor's `hybrid` format.

Inputs may be FASTA, FASTQ, or gzipped variants of either; the file at line
*i* of the filenames list is assigned color *i*.

## Build

Requirements: a C++17 compiler, CMake ≥ 3.13, and zlib.

```bash
git clone --recursive https://github.com/jermp/cdbg-builder.git
cd cdbg-builder
mkdir build && cd build
cmake ..
make -j
```

If you cloned without `--recursive`, run `git submodule update --init --recursive`
first.

The build produces a single executable, `cdgb-build`, in the `build/` directory.

## Usage

```
cdgb-build -i <filenames_list> -k <k> -o <out_basename> [-t <num_threads>]
```

Options:

| Flag                 | Description                                                       | Default |
|----------------------|-------------------------------------------------------------------|---------|
| `-i PATH`            | Text file with one input path per line (one color each)           | —       |
| `-k INT`             | k-mer length (≤ 63)                                               | —       |
| `-o NAME`            | Output basename; writes `NAME.fa` and `NAME.colors`               | —       |
| `-t INT`             | Number of worker threads                                          | 1       |
| `-m INT`             | Minimizer length used for bucketing                               | auto    |
| `--buckets-log2 INT` | `2^N` minimizer-derived bucket files on disk                      | 10      |
| `--tmp-dir PATH`     | Scratch directory for the bucket files (created if missing; if it already exists it must be empty; removed on success) | mkdtemp |
| `-v`                 | Verbose output                                                    | off     |

The build pipeline is GGCAT-style: stream input → write super-k-mers
into per-minimizer bucket files on disk → walk each bucket independently
→ stitch fragments across buckets via shared (k-1)-mer junctions. RAM
usage stays bounded by the largest single bucket rather than the full
k-mer set.

## Example: *Salmonella enterica* pangenome (4,546 genomes)

This reproduces the dataset used in Fulgor's
[example](https://github.com/jermp/fulgor#indexing-an-example-salmonella-enterica-pangenome).

Download and unpack the genomes:

```bash
wget https://zenodo.org/records/1323684/files/Salmonella_enterica.zip
unzip Salmonella_enterica.zip
```

Build the filenames list (one absolute path per line):

```bash
find $(pwd)/Salmonella_enterica/Genomes/*.fasta > salmonella_4546_filenames.txt
```

Run `cdgb-build` (from the `build/` directory) with `k = 31` and 8 threads:

```bash
./cdgb-build \
    -i ~/salmonella_4546_filenames.txt \
    -o ~/Salmonella_enterica/salmonella_4546 \
    -k 31 \
    -t 8 \
    --verbose
```

This produces:

- `~/Salmonella_enterica/salmonella_4546.fa` — colored unitigs in FASTA form
- `~/Salmonella_enterica/salmonella_4546.colors` — Fulgor-compatible hybrid
  color sets

The two files can then be consumed by downstream tools that accept Fulgor's
hybrid color-set format.
