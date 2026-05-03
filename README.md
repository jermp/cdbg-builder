# cdbg-builder

A small, multi-threaded tool that builds a **colored compacted de Bruijn
graph** from a collection of reference sequences and emits the result in a
format compatible with [Fulgor](https://github.com/jermp/fulgor)'s `hybrid`
color-set encoding.

Given `N` input files (each treated as one *color*), the tool produces:

- `<out>.fa` — a FASTA of *monochromatic colored unitigs*. Each header is the
  integer id of the color set the unitig belongs to.
- `<out>.u2c` — a serialized `bits::bit_vector` of length `num_unitigs` where
  bit *i* is set iff unitig *i* (in the `.fa` emission order) is the last
  unitig of a color-set run. The popcount equals the number of distinct
  color sets and the last bit is always set. Downstream consumers (e.g.
  Fulgor) recover the per-unitig color-set id via `rank1(unitig_id)` after
  building a rank index over this bit_vector.
- `<out>.color_sets` — the color sets, written incrementally during
  bucket-process so the compressed bit_vector never sits whole in RAM.
  Layout: a small fixed header (`num_colors`, sparse / dense thresholds,
  `num_color_sets`, `bit_vector_num_bits`, `bit_vector_num_words`) followed
  by the `num_color_sets` distinct sets encoded with the same hybrid
  sparse / dense / complementary-dense rules as Fulgor's hybrid codec,
  followed by an `elias_fano` of per-class bit-offsets appended at the
  end. (Different on-disk layout from Fulgor's existing `hybrid`
  serialization; the consumer must be updated to match.)

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
| `-o NAME`            | Output basename; writes `NAME.fa`, `NAME.u2c`, and `NAME.color_sets` | —       |
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
- `~/Salmonella_enterica/salmonella_4546.u2c` — unitig-to-color-set
  bit_vector
- `~/Salmonella_enterica/salmonella_4546.color_sets` — Fulgor-compatible hybrid
  color sets

The three files can then be consumed by downstream tools that accept Fulgor's
hybrid color-set format.
