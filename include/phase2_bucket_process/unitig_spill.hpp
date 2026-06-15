#pragma once

// Disk-backed sinks/readers used to keep unitig and fragment seq
// bytes off the heap across the bucket-process / stitch / emit
// boundaries. Three classes live here:
//
// 1) unitig_bucket_writer  -- stitch -> emit
//    Sink for the stitch. Partitions finished merged unitigs into K
//    cid-range buckets: bucket b holds every unitig with cid in
//    [b * S, (b+1) * S) where S = ceil(num_color_classes / num_buckets).
//    emit_fasta iterates buckets in ascending order and sorts each
//    bucket's records by cid before writing FASTA, so the .fa output is
//    strictly cid-ascending (which the u2c bit_vector consumer requires).
//
//    Under a -g budget (the production case) it is WRITE-THROUGH: every
//    finished unitig is streamed straight to its bucket's disk file and no
//    unitig bytes are retained in RAM, so the stitch's only resident set is
//    the round store + color dict (the proven README path: 661k stitch peak
//    41 GiB within -g 64). emit then external-merge-sorts each bucket within
//    a small cid-sort cap. Holding unitigs in RAM instead (the old RAM-first
//    optimization, to skip the temp round-trip) blew past -g at 661k -- 1.31e9
//    sunk unitigs over a few cid-skewed buckets that the bucket-granular
//    governor spill could not hold -- so it is now used ONLY with NO -g
//    (ram_budget == SIZE_MAX), where there is no budget to honor and inputs
//    are small enough to keep everything in RAM and skip the disk round-trip.
//
//    Per-bucket on-disk format (spilled buckets only), a sequence of
//    records until EOF:
//      [u64 cid][u32 seq_len][seq_len bytes raw ACGT]
//
// 2) frag_unitig_writer  -- bucket-process -> stitch
//    Sink for process_buckets. Streams every finished fragment
//    (cid + open_flags + bases) to a single mutex-protected append-
//    only file. With this in place the in-RAM accumulator that used
//    to hold every fragment never grows past per-thread buffers.
//
//    File format: a sequence of records, one per fragment, until EOF:
//      [u32 cid][u8 open_flags][u32 seq_len][seq_len bytes raw ACGT]
//
// 3) frag_unitig_stream_reader  -- sequential, no-index view over (2)'s file
//    Used by the production stitch. Pulls fragments one at a time via
//    next(), holding only ONE fragment resident -- no per-fragment index
//    at all. The stitch seeds round 0 by reading every fragment once in
//    order and never revisits the spill, so a random-access index would be
//    pure overhead; at 5.7e9 fragments a 24 B/entry index was ~137 GB,
//    the last O(num_fragments) structure in the stitch. This removes it.

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include "phase2_bucket_process/bucket_walker.hpp"  // stitchable_unitig
#include "ram_governor.hpp"                          // ram_spillable

namespace cdbg {

// ram_spillable: the unitig writer is registered with the RAM governor during
// stitch, so the global measured-RSS backstop can reclaim its in-RAM buckets
// under pressure (this is the buffer the per-phase stitch watcher could not see,
// which busted -g at scale). The spill paths take m_spill_mu so the governor's
// poller thread can spill concurrently with the stitch's (externally-serialized)
// operator() calls.
class unitig_bucket_writer : public ram_spillable {
public:
    // num_color_classes is the final cid space (size of global_dict).
    // num_buckets controls how the cid range is partitioned. Picked
    // so per-bucket peak fits comfortably in the budget remainder.
    //
    // ram_budget bounds how many bytes of stitched unitigs we keep in MEMORY
    // (across all buckets) before spilling the largest buckets to disk. When the
    // whole unitig set fits, nothing is written here and emit reads it straight
    // from RAM -- skipping the temp-bucket write (here) AND the read-back (emit),
    // i.e. the unitig sequences never round-trip through disk. SIZE_MAX = keep
    // everything in RAM (the no-`-g` default); a finite budget spills the
    // overflow so the `-g` bound is honored regardless of unitig volume.
    unitig_bucket_writer(std::string dir, uint64_t num_color_classes, uint32_t num_buckets,
                         size_t ram_budget = SIZE_MAX)
        : m_dir(std::move(dir))
        , m_num_color_classes(num_color_classes)
        , m_num_buckets(num_buckets == 0 ? 1 : num_buckets)
        , m_ram_budget(ram_budget) {
        if (num_color_classes == 0) {
            // No unitigs ever produced -> a single empty bucket.
            m_num_buckets = 1;
        }
        // Files are opened lazily, only when a bucket is spilled (see operator()
        // / spill_bucket_). A run whose unitigs fit ram_budget creates none.
        m_files.assign(m_num_buckets, nullptr);
        m_ram.resize(m_num_buckets);
        m_ram_bytes.assign(m_num_buckets, 0);
    }

    ~unitig_bucket_writer() {
        for (auto*& f : m_files) {
            if (f) {
                std::fclose(f);
                f = nullptr;
            }
        }
    }

    unitig_bucket_writer(unitig_bucket_writer const&) = delete;
    unitig_bucket_writer& operator=(unitig_bucket_writer const&) = delete;

    // The sink callable used by the stitch. The stitch calls it serially under
    // its own sink mutex; m_spill_mu additionally guards against the RAM
    // governor's poller spilling buckets concurrently.
    void operator()(stitchable_unitig&& u) {
        std::lock_guard<std::mutex> lk(m_spill_mu);
        const uint64_t cid = u.mono_cid();  // post-split: monochromatic
        const uint32_t b = bucket_for_cid(cid);
        ++m_total_unitigs;
        // Write-through mode (a -g budget is in force): NEVER retain unitig bytes
        // in RAM -- stream every record straight to its bucket's disk file, exactly
        // like the proven README path (661k stitch peak 41 GiB). The RAM-first
        // optimization (hold unitigs in RAM, skip the temp round-trip) trades disk
        // I/O for RAM and at 661k scale -- 1.31e9 sunk unitigs over few, cid-skewed
        // buckets -- the bucket-granular governor spill cannot hold -g, so RSS ran
        // to 170 GiB. RAM-first is therefore used ONLY with no -g (ram_budget ==
        // SIZE_MAX), where there is no budget to honor and inputs are small.
        if (m_ram_budget != SIZE_MAX) {
            write_record_(ensure_bucket_file_(b), cid, u.seq);
            std::string().swap(u.seq);
            return;
        }
        if (m_files[b]) {
            // Bucket already spilled (one-way): write this record straight
            // through to disk and drop its bytes. A bucket is therefore EITHER
            // fully in RAM (m_files[b] == nullptr, records in m_ram[b]) OR fully
            // on disk -- never split -- so emit needs no RAM+disk merge.
            write_record_(m_files[b], cid, u.seq);
            std::string().swap(u.seq);
            return;
        }
        // RAM-resident bucket: keep the record in memory (this is the byte we are
        // trying NOT to round-trip through disk). est approximates its heap cost.
        const size_t est = sizeof(record) + u.seq.size();
        m_ram[b].push_back(record{cid, std::move(u.seq)});
        m_ram_bytes[b] += est;
        m_buffered_bytes += est;
        if (m_buffered_bytes > m_ram_budget) spill_largest_until_(m_ram_budget);
    }

    // ---- ram_spillable: governed during stitch (see builder) -----------------
    // Free about target_bytes by spilling the largest in-RAM buckets to disk;
    // returns the bytes actually freed. Thread-safe (m_spill_mu).
    uint64_t spill_under_pressure(uint64_t target_bytes) override {
        std::lock_guard<std::mutex> lk(m_spill_mu);
        const size_t before = m_buffered_bytes;
        const size_t floor = target_bytes >= before ? 0 : before - (size_t)target_bytes;
        spill_largest_until_(floor);
        return before - m_buffered_bytes;
    }
    uint64_t spillable_bytes() const override {
        std::lock_guard<std::mutex> lk(m_spill_mu);
        return m_buffered_bytes;
    }

    uint64_t total_unitigs() const { return m_total_unitigs; }
    uint32_t num_buckets() const { return m_num_buckets; }
    uint64_t num_color_classes() const { return m_num_color_classes; }

    std::string bucket_path(uint32_t b) const {
        return m_dir + "/unitig_bucket_" + std::to_string(b) + ".bin";
    }

    // One stitched unitig: its color-set id and its raw ACGT bases. Held in
    // m_ram[b] for RAM-resident buckets; serialized as [u64 cid][u32 len][bytes]
    // for spilled ones.
    struct record {
        uint64_t cid;
        std::string seq;
    };

    // Stream bucket b's records in cid-ASCENDING order, calling
    // emit(cid, seq_view) for each.
    //
    // RAM-resident bucket (the common case when the unitigs fit ram_budget):
    // sort m_ram[b] in place and emit -- ZERO disk I/O, no temp round-trip.
    //
    // Spilled bucket: external cid-sort from its disk file, using at most
    // ~mem_cap bytes regardless of bucket size (read in <=mem_cap chunks, sort
    // each, spill as a sorted run file, k-way merge the runs to the sink). This
    // bounds the emit peak even under cid skew (unitigs cluster in low cids;
    // 661k emit-fasta saw one equal-width range hit +36 GiB). Order WITHIN a cid
    // is irrelevant (verify compares unitig sets; u2c only needs each cid's
    // unitigs contiguous, which cid-order guarantees).
    template <typename Emit>
    void read_bucket_sorted(uint32_t b, uint64_t mem_cap, Emit&& emit) {
        if (!m_files[b]) {
            // RAM-resident: sort + emit directly, then free the bucket's memory.
            auto& recs = m_ram[b];
            std::sort(recs.begin(), recs.end(),
                      [](record const& x, record const& y) { return x.cid < y.cid; });
            for (auto& r : recs) emit(r.cid, std::string_view(r.seq));
            std::vector<record>().swap(recs);
            return;
        }
        std::FILE* f = m_files[b];
        std::fflush(f);
        std::rewind(f);
        if (mem_cap < (1u << 20)) mem_cap = 1u << 20;  // sane floor

        std::vector<record> chunk;
        uint64_t chunk_bytes = 0;
        std::vector<std::string> run_paths;  // spilled sorted-run files

        auto read_one = [&](std::FILE* fp, record& r) -> bool {
            uint32_t seq_len = 0;
            if (std::fread(&r.cid, sizeof(r.cid), 1, fp) != 1) {
                if (std::feof(fp)) return false;
                throw std::runtime_error("short read of cid in bucket " + std::to_string(b));
            }
            if (std::fread(&seq_len, sizeof(seq_len), 1, fp) != 1)
                throw std::runtime_error("short read of len in bucket " + std::to_string(b));
            r.seq.resize(seq_len);
            if (seq_len > 0 and std::fread(r.seq.data(), 1, seq_len, fp) != (size_t)seq_len)
                throw std::runtime_error("short read of seq in bucket " + std::to_string(b));
            return true;
        };
        auto write_one = [&](std::FILE* fp, record const& r) {
            uint32_t seq_len = (uint32_t)r.seq.size();
            if (std::fwrite(&r.cid, sizeof(r.cid), 1, fp) != 1 or
                std::fwrite(&seq_len, sizeof(seq_len), 1, fp) != 1 or
                (seq_len > 0 and std::fwrite(r.seq.data(), 1, seq_len, fp) != (size_t)seq_len))
                throw std::runtime_error("short write of run record for bucket " +
                                         std::to_string(b));
        };
        auto sort_chunk = [&]() {
            std::sort(chunk.begin(), chunk.end(),
                      [](record const& x, record const& y) { return x.cid < y.cid; });
        };
        auto spill_chunk = [&]() {
            std::string p = m_dir + "/unitig_run_" + std::to_string(b) + "_" +
                            std::to_string(run_paths.size()) + ".bin";
            std::FILE* rf = std::fopen(p.c_str(), "wb");
            if (!rf) throw std::runtime_error("cannot open unitig run file: " + p);
            for (auto const& r : chunk) write_one(rf, r);
            std::fclose(rf);
            run_paths.push_back(std::move(p));
            chunk.clear();
            chunk_bytes = 0;
        };

        // Pass 1: read in chunks; sort each. If it all fits one chunk and no
        // prior runs, emit directly (the no-spill fast path).
        record r;
        while (read_one(f, r)) {
            chunk_bytes += sizeof(record) + r.seq.size();
            chunk.push_back(std::move(r));
            r.seq.clear();
            if (chunk_bytes >= mem_cap) {
                sort_chunk();
                spill_chunk();
            }
        }
        if (run_paths.empty()) {
            // Whole bucket fit in RAM: sort + emit, no temp files.
            sort_chunk();
            for (auto const& rec : chunk) emit(rec.cid, std::string_view(rec.seq));
            return;
        }
        if (!chunk.empty()) {
            sort_chunk();
            spill_chunk();
        }  // final partial run

        // Pass 2: k-way merge the sorted runs. One record per run resident +
        // a heap -- bounded by (#runs * one record). #runs = bucket/mem_cap.
        struct cursor {
            std::FILE* fp;
            record cur;
            bool live;
        };
        std::vector<cursor> cur(run_paths.size());
        for (size_t i = 0; i < run_paths.size(); ++i) {
            cur[i].fp = std::fopen(run_paths[i].c_str(), "rb");
            if (!cur[i].fp) throw std::runtime_error("cannot reopen run: " + run_paths[i]);
            cur[i].live = read_one(cur[i].fp, cur[i].cur);
        }
        // Min-heap of run indices by current cid.
        auto worse = [&](size_t a, size_t c) { return cur[a].cur.cid > cur[c].cur.cid; };
        std::vector<size_t> heap;
        heap.reserve(cur.size());
        for (size_t i = 0; i < cur.size(); ++i)
            if (cur[i].live) heap.push_back(i);
        std::make_heap(heap.begin(), heap.end(), worse);
        while (!heap.empty()) {
            std::pop_heap(heap.begin(), heap.end(), worse);
            size_t i = heap.back();
            heap.pop_back();
            emit(cur[i].cur.cid, std::string_view(cur[i].cur.seq));
            if (read_one(cur[i].fp, cur[i].cur)) {
                heap.push_back(i);
                std::push_heap(heap.begin(), heap.end(), worse);
            }
        }
        for (size_t i = 0; i < cur.size(); ++i) {
            if (cur[i].fp) std::fclose(cur[i].fp);
            std::error_code ec;
            std::filesystem::remove(run_paths[i], ec);
        }
    }

    void close_and_unlink() {
        for (uint32_t b = 0; b < m_num_buckets; ++b) {
            if (m_files[b]) {
                std::fclose(m_files[b]);
                m_files[b] = nullptr;
                std::error_code ec;
                std::filesystem::remove(bucket_path(b), ec);  // only spilled buckets exist
            }
            std::vector<record>().swap(m_ram[b]);  // free any RAM-resident records
        }
        m_buffered_bytes = 0;
    }

    // How many buckets had to spill to disk (0 == the whole unitig set stayed in
    // RAM and the temp round-trip was skipped entirely). For the emit log line.
    uint32_t spilled_buckets() const { return m_spilled_buckets; }

private:
    // Write one record to an already-open spilled bucket file.
    void write_record_(std::FILE* f, uint64_t cid, std::string const& seq) {
        const uint32_t seq_len = (uint32_t)seq.size();
        if (std::fwrite(&cid, sizeof(cid), 1, f) != 1 or
            std::fwrite(&seq_len, sizeof(seq_len), 1, f) != 1 or
            (seq_len > 0 and std::fwrite(seq.data(), 1, seq_len, f) != (size_t)seq_len))
            throw std::runtime_error("short write of unitig record to a bucket file");
    }

    // Lazily open bucket b's disk file (one-way: once a bucket has a file it is
    // disk-resident). Shared by the write-through sink path and spill_bucket_.
    std::FILE* ensure_bucket_file_(uint32_t b) {
        if (!m_files[b]) {
            std::string path = bucket_path(b);
            m_files[b] = std::fopen(path.c_str(), "wb+");
            if (!m_files[b])
                throw std::runtime_error("cannot open unitig bucket file: " + path + ": " +
                                         std::strerror(errno));
            ++m_spilled_buckets;
        }
        return m_files[b];
    }

    // Spill bucket b's RAM-resident records to its (lazily-opened) disk file and
    // free the memory. One-way: after this, operator() routes b straight to disk.
    void spill_bucket_(uint32_t b) {
        ensure_bucket_file_(b);
        for (auto const& r : m_ram[b]) write_record_(m_files[b], r.cid, r.seq);
        m_buffered_bytes -= m_ram_bytes[b];
        m_ram_bytes[b] = 0;
        std::vector<record>().swap(m_ram[b]);
    }

    // Spill largest-first until m_buffered_bytes <= floor: the fattest in-RAM
    // buckets (cid-skew hot spots) go to disk, leaving the small ones resident.
    // Bucket granularity, but correct and bounded; degrades to all-disk if
    // nothing fits. Caller holds m_spill_mu. Used both for the model budget
    // (floor = m_ram_budget) and the governor's pressure spill.
    void spill_largest_until_(size_t floor) {
        while (m_buffered_bytes > floor) {
            uint32_t big = UINT32_MAX;
            size_t best = 0;
            for (uint32_t b = 0; b < m_num_buckets; ++b) {
                if (!m_files[b] and m_ram_bytes[b] > best) {
                    best = m_ram_bytes[b];
                    big = b;
                }
            }
            if (big == UINT32_MAX) break;  // nothing left in RAM to evict
            spill_bucket_(big);
        }
    }

    uint32_t bucket_for_cid(uint64_t cid) const {
        if (m_num_color_classes == 0 or m_num_buckets <= 1) return 0;
        // Cids in [b * S, (b+1) * S) live in bucket b, where
        // S = ceil(num_color_classes / num_buckets). Last bucket
        // mops up any rounding remainder.
        const uint64_t S = (m_num_color_classes + m_num_buckets - 1) / m_num_buckets;
        uint32_t b = (uint32_t)(cid / S);
        if (b >= m_num_buckets) b = m_num_buckets - 1;
        return b;
    }

    std::string m_dir;
    uint64_t m_num_color_classes;
    uint32_t m_num_buckets;
    std::vector<std::FILE*> m_files;  // null until the bucket spills
    uint64_t m_total_unitigs = 0;

    // In-RAM unitig buffer (skips the temp round-trip when it fits the budget).
    std::vector<std::vector<record>> m_ram;  // RAM-resident records, per bucket
    std::vector<size_t> m_ram_bytes;         // approx heap bytes held per bucket
    size_t m_ram_budget = SIZE_MAX;          // cap on total m_buffered_bytes
    size_t m_buffered_bytes = 0;             // sum of m_ram_bytes over RAM buckets
    uint32_t m_spilled_buckets = 0;          // count of buckets evicted to disk
    // Guards the mutating spill paths (operator(), spill_under_pressure) so the
    // RAM governor's poller can spill concurrently with the stitch's writes.
    // mutable: spillable_bytes() is const. Registered only during stitch, so
    // emit's single-threaded read_bucket_sorted runs lock-free.
    mutable std::mutex m_spill_mu;
};

// ----------------------------------------------------------------------------
// Disk-backed sink for the bucket-process -> stitch boundary.
//
// process_buckets feeds each finished fragment into a caller-provided
// sink. Used as that sink, frag_unitig_writer streams every fragment
// to a single append-only file in tmp_dir, so the accumulating
// std::vector<stitchable_unitig> never grows to its multi-GB peak
// during bucket-process. After bucket-process completes, the stitch
// reads it lazily via frag_unitig_stream_reader (one fragment at a time).
//
// File format: a sequence of records, one per finished fragment.
//   [u32 cid]
//   [u8  open_flags]
//   [u32 seq_len]
//   [seq_len bytes]    raw ACGT
//
// operator() is mutex-protected: process_buckets calls it from
// num_threads worker threads concurrently. Each call writes one
// record's worth of bytes (~50-150 B typical) under the lock; mutex
// contention on a workload of ~50 M fragments is negligible compared
// to the per-bucket walker work.

class frag_unitig_writer {
public:
    explicit frag_unitig_writer(std::string path) : m_path(std::move(path)) {
        m_file = std::fopen(m_path.c_str(), "wb+");
        if (!m_file) throw std::runtime_error("cannot open frag spill: " + m_path);
    }

    ~frag_unitig_writer() {
        if (m_file) std::fclose(m_file);
    }

    frag_unitig_writer(frag_unitig_writer const&) = delete;
    frag_unitig_writer& operator=(frag_unitig_writer const&) = delete;

    void operator()(stitchable_unitig&& u) {
        const uint8_t flags = u.open_flags;
        const uint32_t nruns = (uint32_t)u.runs.size();
        const uint32_t seq_len = (uint32_t)u.seq.size();

        std::lock_guard<std::mutex> lk(m_mu);
        if (std::fwrite(&flags, sizeof(flags), 1, m_file) != 1 or
            std::fwrite(&nruns, sizeof(nruns), 1, m_file) != 1)
            throw std::runtime_error("short write to " + m_path);
        for (auto const& r : u.runs) {
            if (std::fwrite(&r.cid, sizeof(r.cid), 1, m_file) != 1 or
                std::fwrite(&r.num_kmers, sizeof(r.num_kmers), 1, m_file) != 1)
                throw std::runtime_error("short write to " + m_path);
        }
        if (std::fwrite(&seq_len, sizeof(seq_len), 1, m_file) != 1)
            throw std::runtime_error("short write to " + m_path);
        if (seq_len > 0 and std::fwrite(u.seq.data(), 1, seq_len, m_file) != (size_t)seq_len)
            throw std::runtime_error("short write to " + m_path);
        ++m_count;
        m_total_seq_bytes += seq_len;
        // Free the merged seq's backing storage in place: the caller
        // already moved into us.
        std::string().swap(u.seq);
    }

    // Write a whole bucket's fragments under ONE lock acquisition instead of one
    // per fragment. process_buckets drains thousands of fragments per bucket from
    // 32 threads; calling operator() per fragment took the frag-sink mutex
    // O(num_fragments) times (~92.5 M on bw20k), and that lock churn was the
    // dominant source of bucket-process voluntary context switches. Here the frag
    // bytes are serialized into a thread-local buffer OUTSIDE the lock; the lock
    // then guards only one fwrite. frag_id order is the order batches arrive.
    void write_batch(std::vector<stitchable_unitig>& batch) {
        if (batch.empty()) return;
        auto put = [](std::vector<uint8_t>& b, void const* p, size_t n) {
            uint8_t const* q = (uint8_t const*)p;
            b.insert(b.end(), q, q + n);
        };
        std::vector<uint8_t> fbuf;
        uint64_t add_count = 0, add_seq = 0;
        for (auto const& u : batch) {
            const uint8_t flags = u.open_flags;
            const uint32_t nruns = (uint32_t)u.runs.size();
            const uint32_t seq_len = (uint32_t)u.seq.size();
            put(fbuf, &flags, sizeof(flags));
            put(fbuf, &nruns, sizeof(nruns));
            for (auto const& r : u.runs) {
                put(fbuf, &r.cid, sizeof(r.cid));
                put(fbuf, &r.num_kmers, sizeof(r.num_kmers));
            }
            put(fbuf, &seq_len, sizeof(seq_len));
            put(fbuf, u.seq.data(), seq_len);
            ++add_count;
            add_seq += seq_len;
        }
        {
            std::lock_guard<std::mutex> lk(m_mu);
            if (!fbuf.empty() and std::fwrite(fbuf.data(), 1, fbuf.size(), m_file) != fbuf.size())
                throw std::runtime_error("short write to " + m_path);
            m_count += add_count;
            m_total_seq_bytes += add_seq;
        }
        for (auto& u : batch) std::string().swap(u.seq);
    }

    uint64_t count() const { return m_count; }
    // Total raw seq bytes across all fragments written. Lets the stitch
    // size its bucket count without an index-walk over the spill (the old
    // O(num_fragments) frag_unitig_reader index).
    uint64_t total_seq_bytes() const { return m_total_seq_bytes; }
    std::string const& path() const { return m_path; }

    // Close the writer side. Call before reads.
    void close_for_writing() {
        std::lock_guard<std::mutex> lk(m_mu);
        if (m_file) {
            std::fflush(m_file);
            std::fclose(m_file);
            m_file = nullptr;
        }
    }

    void unlink() {
        std::error_code ec;
        if (!m_path.empty()) {
            std::filesystem::remove(m_path, ec);
            m_path.clear();
        }
    }

private:
    std::string m_path;
    std::FILE* m_file = nullptr;
    std::mutex m_mu;
    uint64_t m_count = 0;
    uint64_t m_total_seq_bytes = 0;
};

// ----------------------------------------------------------------------------
// Streaming (no-index) reader over a frag_unitig_writer's spill file.
//
// The stitch seeds round 0 by reading every fragment exactly once, in order,
// and nothing afterward revisits the original spill (later rounds live in the
// round_store). So the random-access frag_unitig_reader index -- a
// std::vector<entry> at 24 B/fragment, ~137 GB at 5.7e9 fragments, the last
// O(num_fragments) in-RAM structure in the stitch -- is pure overhead in
// production. This reader keeps only ONE fragment resident: a pull `next()`
// decodes the next record into reused buffers. Peak RAM is one fragment, not
// num_fragments.
//
// It deliberately does NOT satisfy the random-access Source interface (size /
// runs(i) / seq_view(i)); it exposes a pull `next(open_flags, runs, seq)` that
// the stitch's streaming seed path consumes. Record format matches
// frag_unitig_writer exactly:
//   [u8 flags][u32 nruns][nruns x (u64 cid + u32 num_kmers)][u32 seq_len][seq]
class frag_unitig_stream_reader {
public:
    explicit frag_unitig_stream_reader(std::string path) : m_path(std::move(path)) {
        m_file = std::fopen(m_path.c_str(), "rb");
        if (!m_file) throw std::runtime_error("cannot open frag spill: " + m_path);
    }
    ~frag_unitig_stream_reader() {
        if (m_file) std::fclose(m_file);
    }
    frag_unitig_stream_reader(frag_unitig_stream_reader const&) = delete;
    frag_unitig_stream_reader& operator=(frag_unitig_stream_reader const&) = delete;

    // Decode the next fragment into the caller's buffers. Returns false at EOF.
    // `runs` and `seq` are resized to this fragment's contents (capacity reused
    // across calls, so steady-state is alloc-free).
    bool next(uint8_t& open_flags, std::vector<color_run>& runs, std::string& seq) {
        uint8_t flags = 0;
        size_t got = std::fread(&flags, sizeof(flags), 1, m_file);
        if (got != 1) {
            if (std::feof(m_file)) return false;
            throw std::runtime_error("short read of flags from " + m_path);
        }
        uint32_t nruns = 0;
        if (std::fread(&nruns, sizeof(nruns), 1, m_file) != 1)
            throw std::runtime_error("short read of nruns from " + m_path);
        runs.resize(nruns);
        for (uint32_t r = 0; r < nruns; ++r) {
            if (std::fread(&runs[r].cid, sizeof(uint64_t), 1, m_file) != 1 or
                std::fread(&runs[r].num_kmers, sizeof(uint32_t), 1, m_file) != 1)
                throw std::runtime_error("short read of run from " + m_path);
        }
        uint32_t seq_len = 0;
        if (std::fread(&seq_len, sizeof(seq_len), 1, m_file) != 1)
            throw std::runtime_error("short read of seq_len from " + m_path);
        open_flags = flags;
        seq.resize(seq_len);
        if (seq_len > 0 and std::fread(seq.data(), 1, seq_len, m_file) != (size_t)seq_len)
            throw std::runtime_error("short read of seq from " + m_path);
        return true;
    }

private:
    std::string m_path;
    std::FILE* m_file = nullptr;
};

}  // namespace cdbg
