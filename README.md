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

By default k-mers are stored in a 64-bit integer, supporting `k ≤ 31` (cheaper
arithmetic and half the memory of k-mer-keyed structures). To support `k ≤ 63`,
configure with 128-bit k-mers:

```bash
cmake .. -DCDBG_LARGE_K=ON
```

The build produces a single executable, `cdbg-build`, in the `build/` directory.

## Usage

```
cdbg-build -i <filenames_list> -k <k> -o <out_basename> [-t <num_threads>] [-g <GiB>]
```

Options:

| Flag         | Description                                                       | Default |
|--------------|-------------------------------------------------------------------|---------|
| `-i PATH`    | Text file with one input path per line (one color each)           | —       |
| `-k INT`     | k-mer length (≤ 31 by default; build with `-DCDBG_LARGE_K=ON` for ≤ 63) | —  |
| `-o NAME`    | Output basename; writes `NAME.fa`, `NAME.u2c`, and `NAME.color_sets` | —    |
| `-t INT`     | Number of worker threads                                          | 1       |
| `-m INT`     | Minimizer length used for bucketing                               | auto    |
| `-b INT`     | log2 of the bucket count (`2^N` bucket files on disk)             | auto (derived from `-g` if set, else 10) |
| `-d PATH`    | Scratch directory for the bucket files (created if missing; if it already exists it must be empty; removed on success) | mkdtemp |
| `-g FLOAT`   | RAM budget in GiB. Sizes bucket count + spill thresholds and drives the spill machinery so peak RSS is held to `-g` (reported at the end). Honored as a hard cap for all phases **except** two non-spillable structures in bucket-process (the per-bucket walk set and the color-sets-dedup-map) — see `algorithm.md` §9.1; at large scale (e.g. 661k) the dedup-map can still require a larger `-g`. | unset |
| `--alpha FLOAT` | RAM-model per-thread buffer overhead multiplier (re-calibrates the bucket-count model) | 2.0 |
| `--beta FLOAT`  | RAM-model per-bucket compactor overhead multiplier (re-calibrates the bucket-count model) | 7.0 |
| `--flush INT`   | Per-thread→compactor handoff size in bytes (smaller → more, smaller buckets) | 4096 |
| `--spill INT`   | Compactor dedup window in bytes before a disk frame (smaller → more buckets, weaker dedup) | 65536 |
| `--verbose`  | Verbose output                                                    | off     |

The build runs as four sequential, GGCAT-style phases (see `algorithm.md`
for details), each communicating with the next only through on-disk
artifacts:

1. **bucket-write** — stream input, write super-k-mers into per-minimizer
   bucket files on disk (minimizer over (k−1)-mers, so every dBG edge/branch
   co-locates in one bucket).
2. **bucket-process** — walk each bucket independently into open-ended
   fragments and intern the global color sets.
3. **stitch** — join fragments across buckets by matching their shared
   **full boundary k-mer** (id-only doubling, then a single base/color
   assembly pass).
4. **emit** — write `.fa` + `.u2c` and finalize `.color_sets`.

When `-g` is set, every phase sizes itself against that budget (and the
phases that can, spill to disk), so peak RAM stays within the budget
rather than scaling with the full k-mer set.

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

	/usr/bin/time -v ./cdbg-build -i /mnt/hd2/pibiri/DNA/filenames/blackwell-661k_filenames.txt -o /mnt/hd2/pibiri/DNA/bw-661k -k 31 -t 48 --verbose -d /mnt/hd2/pibiri/DNA/tmp_dir/ -g 64
	
	k = 31, m = 12, num_colors = 661405, num_threads = 48, num_buckets = 40329
	  bucket-write model: flush_bases=4.00 KiB, spill_bytes=64.00 KiB, alpha=2, beta=7
	  tmp_dir = /mnt/hd2/pibiri/DNA/tmp_dir/
	[bucket-write] 699.17 GiB/699.17 GiB (100.0%) | files 661405/661405 (100.0%) 15486.6s
	[bucket-write] 15493.7 s
	  [bucket-write peak RSS] 34.87 GiB (+34.76 GiB) [live in 119.92 MiB -> out 34.85 GiB]
	  [bucket-write ctx-switch] voluntary 17924398  involuntary 3059151
	  bucket bytes written: 433927435329 (compressed; 508906060943 uncompressed)
	[bucket size dist] buckets=40329 nonempty=40329  uncompressed bytes:
	  mean=12.03 MiB  p50=11.57  p90=16.62  p99=22.62  p999=30.65  max=88.04 MiB
	  largest-48-sum=1.68 GiB (worst-case co-resident bucket bytes; kmer_info RAM is a ~constant multiple of this)
	  max-bucket=92319292 bytes  total=473.96 GiB
	[bucket-write profile] (per-thread time, ns/threads -> wall-equiv):
	  seq_read      30.84s   (libdeflate gzip + kseq parsing)
	  compute      1658.80s   (ACGT scan + 2-bit + ntHash + minimizer + per-thread append)
	    scan2bit    25.31s   (ACGT-run scan + 2-bit only; helicase's domain)
	  flush        3575.46s   (writer.flush total = lock + hashmap + spill)
	    lock_wait  198.79s
	    hashmap    496.74s
	    spill      2859.65s   (sort+unique + write_super_kmer_packed + gzwrite)
	  counts: files=661405 records=275184679404 flushes=2700094595 spills=19299630 inserts=14734927189
	[bucket-process] 40329/40329 (100.0%) 18028.9s
	  bucket fragments: 6079027539
	  distinct color classes: 389111906
	  color-sets-dedup-map: ~13.92 GiB resident / 32.00 GiB budget (within -> no spill)
	[bucket-process] 18028.9 s
	  [bucket-process peak RSS] 34.87 GiB (+0) [live in 867.04 MiB -> out 25.45 GiB]
	  [bucket-process ctx-switch] voluntary 155193662  involuntary 1243017
	[bucket-process profile] (per-thread time, ns/threads -> wall-equiv):
	  load        5561.28s   (bucket_reader + LZ4 + record intern + kmer hashmap roll)
	  resolve     1157.15s   (rsid -> local cid; record_sets.at + local_dict.intern)
	  walk         399.80s   (classify_left_end + extend_and_emit)
	  pre_decode   137.52s   (local_dict.at, lock-free)
	  pre_hash     114.60s   (wyhash + fnv1a on decoded class, lock-free)
	  merge_wait     0.00s   (waiting on global_mu)
	  merge       2822.03s   (global_dict.intern_with_hashes under global_mu)
	  counts: buckets=40329 records=14734927189 kmers=157408583654 local_classes=2282712132 unitigs=6079027539
	  stitch buckets: 2328, threads: 48
	  compact stitch: scalable, RAM-first cap 17.28 GiB
	  [stitch] round 0: joined=1655257190, 2364.9s
	  [stitch] round 1: joined=1102582898, 601.178s
	  [stitch] round 2: joined=723853801, 414.466s
	  [stitch] round 3: joined=450132365, 305.792s
	  [stitch] round 4: joined=281920720, 240.161s
	  [stitch] round 5: joined=183979910, 194.965s
	  ... (rounds 6-29 omitted; joined halves each round) ...
	  [stitch] round 30: joined=12, 0.307654s
	  [stitch] round 31: joined=1, 1.26333s
	  [stitch] round 32: joined=1, 1.23173s
	  [stitch] round 33: joined=1, 0.170981s
	  [stitch] round 34: joined=0, 0.156067s
	  [stitch] 35 rounds, 5008.88s total (3926.5s in rounds 0-4)
	[stitch] 6079027539/6079027539 (100.0%) 10717.8s
	  unitigs after stitching: 1313598628
	[stitch] 10717.8 s
	  [stitch peak RSS] 37.76 GiB (+2.89 GiB) [live in 1.75 GiB -> out 8.14 GiB]
	  [stitch ctx-switch] voluntary 14093458  involuntary 459244
	  [emit-fasta] 68 buckets (68 spilled to disk, rest read from RAM), cid-sort RAM cap 6.40 GiB
	[emit fasta] 1366.94 s
	  [emit-fasta peak RSS] 37.76 GiB (+0) [live in 308.95 MiB -> out 6.15 GiB]
	  [emit-fasta ctx-switch] voluntary 400024  involuntary 8238
	  num_color_sets = 389111906
	  num_total_integers = 1479253359128
	  total bits for ints  = 3835033939822
	  total bits for offs  = 6134587136
	[emit color_sets] 33.5212 s
	  [emit-colors peak RSS] 37.76 GiB (+0) [live in 6.15 GiB -> out 6.15 GiB]
	  [emit-colors ctx-switch] voluntary 24918  involuntary 105
	[removing tmp files] 0.0505843 s
	[peak resident memory] 37.76 GiB  (within budget of 64.00 GiB)
	done. wrote /mnt/hd2/pibiri/DNA/bw-661k.fa, /mnt/hd2/pibiri/DNA/bw-661k.u2c, and /mnt/hd2/pibiri/DNA/bw-661k.color_sets
	[total construction time] 45675.5 s
	
	        Command being timed: "/usr/bin/time -v ./cdbg-build -i /mnt/hd2/pibiri/DNA/filenames/blackwell-661k_filenames.txt -o /mnt/hd2/pibiri/DNA/bw-661k -k 31 -t 48 --verbose -d /mnt/hd2/pibiri/DNA/tmp_dir/ -g 64"
	        User time (seconds): 237533.49
	        System time (seconds): 6907.00
	        Percent of CPU this job got: 535%
	        Elapsed (wall clock) time (h:mm:ss or m:ss): 12:41:16
	        Maximum resident set size (kbytes): 39592348
	        Major (requiring I/O) page faults: 661406
	        Minor (reclaiming a frame) page faults: 2135645097
	        Voluntary context switches: 187637593
	        Involuntary context switches: 4771183
	        Swaps: 0
	        File system inputs: 3537004528
	        File system outputs: 4988045168
	        Page size (bytes): 4096
	        Exit status: 0