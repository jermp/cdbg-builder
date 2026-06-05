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

Requirements: a C++17 compiler and CMake ≥ 3.13. All other dependencies
(libdeflate, lz4, bits, cmd_line_parser, unordered_dense, kseq) are
vendored as git submodules.

```bash
git clone --recursive https://github.com/jermp/cdbg-builder.git
cd cdbg-builder
mkdir build && cd build
cmake ..
make -j
```

If you cloned without `--recursive`, run `git submodule update --init --recursive`
first.

The build produces a single executable, `cdbg-build`, in the `build/` directory.

## Usage

```
cdbg-build -i <filenames_list> -k <k> -o <out_basename> [-t <num_threads>] [-g <GiB>]
```

Options:

| Flag         | Description                                                       | Default |
|--------------|-------------------------------------------------------------------|---------|
| `-i PATH`    | Text file with one input path per line (one color each)           | —       |
| `-k INT`     | k-mer length (≤ 63)                                               | —       |
| `-o NAME`    | Output basename; writes `NAME.fa`, `NAME.u2c`, and `NAME.color_sets` | —    |
| `-t INT`     | Number of worker threads                                          | 1       |
| `-m INT`     | Minimizer length used for bucketing                               | auto    |
| `-b INT`     | log2 of the bucket count (`2^N` bucket files on disk)             | auto (derived from `-g` if set, else 10) |
| `-d PATH`    | Scratch directory for the bucket files (created if missing; if it already exists it must be empty; removed on success) | mkdtemp |
| `-g FLOAT`   | Soft RAM budget in GiB. Tunes bucket count + spill thresholds and arms a runtime RSS watcher; peak is reported at the end (not a hard cap) | unset |
| `--verbose`  | Verbose output                                                    | off     |

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

Run `cdbg-build` (from the `build/` directory) with `k = 31`, 8 threads,
and a 4 GiB soft RAM budget:

```bash
./cdbg-build \
    -i ~/salmonella_4546_filenames.txt \
    -o ~/Salmonella_enterica/salmonella_4546 \
    -k 31 \
    -t 8 \
    -g 4 \
    --verbose
```

For larger pangenomes, scale `-t` to your core count and `-g` to
roughly half of available memory. The bucket-write phase auto-tunes
the bucket count and per-bucket spill thresholds against `-g` to
keep peak RSS under that budget.

This produces:

- `~/Salmonella_enterica/salmonella_4546.fa` — colored unitigs in FASTA form
- `~/Salmonella_enterica/salmonella_4546.u2c` — unitig-to-color-set
  bit_vector
- `~/Salmonella_enterica/salmonella_4546.color_sets` — Fulgor-compatible hybrid
  color sets

The three files can then be consumed by downstream tools that accept Fulgor's
hybrid color-set format.

## Experiment on Blackwell 661k

	./cdbg-build -i /mnt/hd2/pibiri/DNA/filenames/blackwell-661k_filenames.txt -o /mnt/hd2/pibiri/DNA/bw-661k -k 31 -t 48 --verbose -d /mnt/hd2/pibiri/DNA/tmp_dir/ -g 64
	
	k = 31, m = 12, num_colors = 661405, num_threads = 48, num_buckets = 40329
	  bucket-write model: flush_bases=4.00 KiB, spill_bytes=64.00 KiB, alpha=2, beta=7
	  tmp_dir = /mnt/hd2/pibiri/DNA/tmp_dir/
	[bucket-write] 661405/661405 (100.0%) 24689.1s      
	[bucket-write] 24689.1 s
	  [bucket-write peak RSS] 30.81 GiB (+30.70 GiB) [live in 119.78 MiB -> out 30.66 GiB]
	  bucket bytes written: 464134845529 (compressed; 540328029902 uncompressed)
	[bucket size dist] buckets=40329 nonempty=40329  uncompressed bytes:
	  mean=12.78 MiB  p50=12.26  p90=17.74  p99=24.46  p999=33.59  max=114.36 MiB
	  largest-48-sum=1.94 GiB (worst-case co-resident bucket bytes; kmer_info RAM is a ~constant multiple of this)
	  max-bucket=119910429 bytes  total=503.22 GiB
	[bucket-write profile] (per-thread time, ns/threads -> wall-equiv):
	  seq_read      46.22s   (libdeflate gzip + kseq parsing)
	  compute      1280.42s   (ACGT scan + 2-bit + ntHash + minimizer + per-thread append)
	  flush        4962.68s   (writer.flush total = lock + hashmap + spill)
	    lock_wait  240.85s
	    hashmap    787.38s
	    spill      3914.12s   (sort+unique + write_super_kmer + gzwrite)
	  counts: files=661405 records=275184679404 flushes=2700095542 spills=26990321 inserts=16792621468
	[bucket-process] 40329/40329 (100.0%) 14794.8s      
	  bucket fragments: 6079027539
	  distinct color classes: 389111906
	  global color dict resident: ~11.02 GiB (stays in RAM through stitch + emit)
	[bucket-process] 14794.8 s
	  [bucket-process peak RSS] 31.14 GiB (+332.80 MiB) [live in 748.09 MiB -> out 24.84 GiB]
	[bucket-process profile] (per-thread time, ns/threads -> wall-equiv):
	  load        6153.45s   (bucket_reader + LZ4 + record intern + kmer hashmap roll)
	  resolve     1256.40s   (rsid -> local cid; record_sets.at + local_dict.intern)
	  walk        1675.65s   (classify_left_end + extend_and_emit)
	  pre_decode   139.75s   (local_dict.at, lock-free)
	  pre_hash     117.76s   (wyhash + fnv1a on decoded class, lock-free)
	  merge_wait  2726.36s   (waiting on global_mu)
	  merge        184.44s   (global_dict.intern_with_hashes under global_mu)
	  counts: buckets=40329 records=16792621468 kmers=179239979488 local_classes=2282712132 unitigs=6079027539
	  stitch buckets: 2319, threads: 48
	[stitch] 507235476/6079027539 (8.3%) 4886.9s        [stitch] round-0 seed: 4887.07s
	[stitch] 1313555819/6079027539 (21.6%) 12586.5s        [stitch] 33 rounds, 7699.66s total (6542.66s in rounds 0-4)
	[stitch] 6079027539/6079027539 (100.0%) 12588.0s      
	  unitigs after stitching: 1313598625
	[stitch] 12588 s
	  [stitch peak RSS] 41.03 GiB (+9.89 GiB) [live in 1.51 GiB -> out 9.57 GiB]
	  [emit-fasta] 68 buckets, cid-sort RAM cap 6.40 GiB (external merge-sort if a bucket exceeds it)
	[emit fasta] 1915.19 s
	  [emit-fasta peak RSS] 41.03 GiB (+0) [live in 9.57 GiB -> out 7.07 GiB]
	  num_color_sets = 389111906
	  num_total_integers = 1479253359128
	  total bits for ints  = 3835033939822
	  total bits for offs  = 6134587136
	[emit color_sets] 34.083 s
	  [emit-colors peak RSS] 41.03 GiB (+0) [live in 7.07 GiB -> out 7.07 GiB]
	[removing tmp files] 6.487 s
	[peak resident memory] 41.03 GiB  (within budget of 64.00 GiB)
	done. wrote /mnt/hd2/pibiri/DNA/bw-661k.fa, /mnt/hd2/pibiri/DNA/bw-661k.u2c, and /mnt/hd2/pibiri/DNA/bw-661k.color_sets
	[total construction time] 54067.3 s