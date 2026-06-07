#pragma once

// Disk-backed sinks/readers used to keep unitig and fragment seq
// bytes off the heap across the bucket-process / stitch / emit
// boundaries. Three classes live here:
//
// 1) unitig_bucket_writer  -- stitch -> emit
//    Sink for stitch_unitigs_streaming. Partitions finished merged
//    unitigs into K cid-range bucket files: bucket b holds every
//    unitig with cid in [b * S, (b+1) * S) where
//    S = ceil(num_color_classes / num_buckets). emit_fasta iterates
//    buckets in ascending order and sorts each bucket's records by
//    cid in memory before writing FASTA, so the .fa output is
//    strictly cid-ascending (which the u2c bit_vector consumer
//    requires). Per-bucket peak at emit is one bucket's seq payload
//    plus its record-vector overhead.
//
//    Per-bucket file format: a sequence of records, one per finished
//    unitig, until EOF:
//      [u32 cid][u32 seq_len][seq_len bytes raw ACGT]
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

namespace cdbg {

class unitig_bucket_writer {
public:
    // num_color_classes is the final cid space (size of global_dict).
    // num_buckets controls how the cid range is partitioned. Picked
    // so per-bucket peak fits comfortably in the budget remainder.
    unitig_bucket_writer(std::string dir, uint64_t num_color_classes, uint32_t num_buckets)
        : m_dir(std::move(dir))
        , m_num_color_classes(num_color_classes)
        , m_num_buckets(num_buckets == 0 ? 1 : num_buckets) {
        if (num_color_classes == 0) {
            // No unitigs ever produced -> a single empty bucket.
            m_num_buckets = 1;
        }
        m_files.assign(m_num_buckets, nullptr);
        for (uint32_t b = 0; b < m_num_buckets; ++b) {
            std::string path = bucket_path(b);
            m_files[b] = std::fopen(path.c_str(), "wb+");
            if (!m_files[b])
                throw std::runtime_error("cannot open unitig bucket file: " + path + ": " +
                                         std::strerror(errno));
        }
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

    // The sink callable used by stitch_unitigs_streaming.
    void operator()(stitchable_unitig&& u) {
        const uint64_t cid = u.mono_cid();  // post-split: monochromatic
        const uint32_t b = bucket_for_cid(cid);
        std::FILE* f = m_files[b];
        const uint32_t seq_len = (uint32_t)u.seq.size();
        if (std::fwrite(&cid, sizeof(cid), 1, f) != 1)
            throw std::runtime_error("short write of unitig cid to bucket " + std::to_string(b));
        if (std::fwrite(&seq_len, sizeof(seq_len), 1, f) != 1)
            throw std::runtime_error("short write of unitig len to bucket " + std::to_string(b));
        if (seq_len > 0 and std::fwrite(u.seq.data(), 1, seq_len, f) != (size_t)seq_len) {
            throw std::runtime_error("short write of unitig seq to bucket " + std::to_string(b));
        }
        ++m_total_unitigs;
        // Free the seq buffer eagerly: stitch already moved-from frag,
        // and the merged seq's backing storage isn't needed any more.
        std::string().swap(u.seq);
    }

    uint64_t total_unitigs() const { return m_total_unitigs; }
    uint32_t num_buckets() const { return m_num_buckets; }
    uint64_t num_color_classes() const { return m_num_color_classes; }

    std::string bucket_path(uint32_t b) const {
        return m_dir + "/unitig_bucket_" + std::to_string(b) + ".bin";
    }

    // Read all records from bucket `b`. Caller-owned `out` is cleared
    // and refilled. Records arrive in *write order* (which is stitch
    // order, NOT cid order) -- caller sorts.
    struct record {
        uint64_t cid;
        std::string seq;
    };

    void read_bucket(uint32_t b, std::vector<record>& out) {
        out.clear();
        std::FILE* f = m_files[b];
        if (!f) throw std::runtime_error("read_bucket: file already closed");
        std::fflush(f);
        std::rewind(f);
        for (;;) {
            uint64_t cid = 0;
            uint32_t seq_len = 0;
            size_t got = std::fread(&cid, sizeof(cid), 1, f);
            if (got != 1) {
                if (std::feof(f)) break;
                throw std::runtime_error("short read of unitig cid in bucket " + std::to_string(b));
            }
            if (std::fread(&seq_len, sizeof(seq_len), 1, f) != 1)
                throw std::runtime_error("short read of unitig len in bucket " + std::to_string(b));
            std::string seq;
            seq.resize(seq_len);
            if (seq_len > 0 and std::fread(seq.data(), 1, seq_len, f) != (size_t)seq_len) {
                throw std::runtime_error("short read of unitig seq in bucket " + std::to_string(b));
            }
            out.push_back({cid, std::move(seq)});
        }
    }

    // Stream bucket b's records in cid-ASCENDING order, calling
    // emit(cid, seq_view) for each, using at most ~mem_cap bytes of RAM
    // regardless of bucket size. This replaces "load whole bucket + in-RAM
    // sort", whose peak was the FATTEST bucket -- and cid skew (unitigs
    // cluster in low cids; cid ranges are equal-width) made one bucket
    // dominate (661k emit-fasta: +36 GiB).
    //
    // Common case (bucket <= mem_cap): one in-RAM chunk, sort, emit -- no
    // spill, same speed as before. Large bucket: read in <=mem_cap chunks,
    // sort each, spill as a sorted run file, then k-way merge the runs to the
    // sink. Order WITHIN a cid is irrelevant (verify compares unitig sets; u2c
    // only needs each cid's unitigs contiguous, which cid-order guarantees).
    template <typename Emit>
    void read_bucket_sorted(uint32_t b, uint64_t mem_cap, Emit&& emit) {
        std::FILE* f = m_files[b];
        if (!f) throw std::runtime_error("read_bucket_sorted: file already closed");
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
            }
            std::error_code ec;
            std::filesystem::remove(bucket_path(b), ec);
        }
    }

private:
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
    std::vector<std::FILE*> m_files;
    uint64_t m_total_unitigs = 0;
};

// ----------------------------------------------------------------------------
// Disk-backed sink for the bucket-process -> stitch boundary.
//
// process_buckets feeds each finished fragment into a caller-provided
// sink. Used as that sink, frag_unitig_writer streams every fragment
// to a single append-only file in tmp_dir, so the accumulating
// std::vector<stitchable_unitig> never grows to its multi-GB peak
// during bucket-process. After bucket-process completes, callers
// either load the file back into a vector before stitch (simple path)
// or read it lazily during stitch (memory-frugal path -- not
// implemented here).
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
    // When k >= 2 a companion "links" spill is written alongside the frag spill,
    // one fixed-size record per fragment IN THE SAME frag_id order (both writes
    // happen under the same lock in operator()). Each record is
    //   [u8 open_flags][kbytes kl][kbytes kr]   (little-endian, kbytes=(2k+7)/8)
    // where kl/kr are the FORWARD boundary k-mers of the open ends (0 otherwise),
    // bit-identical to detail::id_fwd_kmer in compact_extmem. The id-only stitch's
    // seed reads this tiny stream instead of re-reading every base out of the 7 GB
    // frag spill just to recover two k-mers per fragment. k==0 disables it.
    explicit frag_unitig_writer(std::string path, uint32_t k = 0)
        : m_path(std::move(path)), m_k(k >= 2 ? k : 0) {
        m_file = std::fopen(m_path.c_str(), "wb+");
        if (!m_file) throw std::runtime_error("cannot open frag spill: " + m_path);
        if (m_k) {
            m_kbytes = (2u * m_k + 7u) / 8u;  // <= 16 (kmer_int_t is 128-bit)
            m_links_path = m_path + ".links";
            m_links_file = std::fopen(m_links_path.c_str(), "wb+");
            if (!m_links_file) throw std::runtime_error("cannot open links spill: " + m_links_path);
        }
    }

    ~frag_unitig_writer() {
        if (m_file) std::fclose(m_file);
        if (m_links_file) std::fclose(m_links_file);
    }

    frag_unitig_writer(frag_unitig_writer const&) = delete;
    frag_unitig_writer& operator=(frag_unitig_writer const&) = delete;

    void operator()(stitchable_unitig&& u) {
        const uint8_t flags = u.open_flags;
        const uint32_t nruns = (uint32_t)u.runs.size();
        const uint32_t seq_len = (uint32_t)u.seq.size();
        // Build the links record BEFORE taking the lock: the boundary k-mer
        // compute is the only real CPU here, so keep it off the critical section
        // (otherwise it serializes across all bucket-process threads). frag_id
        // order is still identical to the frag spill because the fwrites below
        // run under the lock in call order. Compute the FORWARD boundary k-mer of
        // each open end exactly as detail::id_fwd_kmer does
        // (v = (v<<2)|2bit over the first/last k bases); store low kbytes
        // little-endian. open_flags != 0 implies seq_len >= k.
        uint8_t lrec[1 + 32];
        size_t lrsz = 0;
        if (m_links_file) {
            kmer_int_t kl = 0, kr = 0;
            if (flags & UNITIG_OPEN_LEFT) {
                kmer_int_t v = 0;
                for (uint32_t i = 0; i < m_k; ++i) v = (v << 2) | (kmer_int_t)nuc_to_2bit(u.seq[i]);
                kl = v;
            }
            if (flags & UNITIG_OPEN_RIGHT) {
                kmer_int_t v = 0;
                char const* p = u.seq.data() + (size_t)(seq_len - m_k);
                for (uint32_t i = 0; i < m_k; ++i) v = (v << 2) | (kmer_int_t)nuc_to_2bit(p[i]);
                kr = v;
            }
            lrec[0] = flags;
            for (uint32_t b = 0; b < m_kbytes; ++b) lrec[1 + b] = (uint8_t)(kl >> (8 * b));
            for (uint32_t b = 0; b < m_kbytes; ++b)
                lrec[1 + m_kbytes + b] = (uint8_t)(kr >> (8 * b));
            lrsz = (size_t)1 + 2 * m_kbytes;
        }

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
        if (lrsz and std::fwrite(lrec, 1, lrsz, m_links_file) != lrsz)
            throw std::runtime_error("short write to " + m_links_path);
        // Free the merged seq's backing storage in place: the caller
        // already moved into us.
        std::string().swap(u.seq);
    }

    // Write a whole bucket's fragments under ONE lock acquisition instead of one
    // per fragment. process_buckets drains thousands of fragments per bucket from
    // 32 threads; calling operator() per fragment took the frag-sink mutex
    // O(num_fragments) times (~92.5 M on bw20k), and that lock churn was the
    // dominant source of bucket-process voluntary context switches. Here the frag
    // + link bytes are serialized into thread-local buffers OUTSIDE the lock; the
    // lock then guards only two fwrites. frag_id order matches the link spill
    // (both written together per batch) and is the order batches arrive.
    void write_batch(std::vector<stitchable_unitig>& batch) {
        if (batch.empty()) return;
        auto put = [](std::vector<uint8_t>& b, void const* p, size_t n) {
            uint8_t const* q = (uint8_t const*)p;
            b.insert(b.end(), q, q + n);
        };
        std::vector<uint8_t> fbuf, lbuf;
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
            if (m_links_file) {
                kmer_int_t kl = 0, kr = 0;
                if (flags & UNITIG_OPEN_LEFT) {
                    kmer_int_t v = 0;
                    for (uint32_t i = 0; i < m_k; ++i)
                        v = (v << 2) | (kmer_int_t)nuc_to_2bit(u.seq[i]);
                    kl = v;
                }
                if (flags & UNITIG_OPEN_RIGHT) {
                    kmer_int_t v = 0;
                    char const* p = u.seq.data() + (size_t)(seq_len - m_k);
                    for (uint32_t i = 0; i < m_k; ++i) v = (v << 2) | (kmer_int_t)nuc_to_2bit(p[i]);
                    kr = v;
                }
                uint8_t rec[1 + 32];
                rec[0] = flags;
                for (uint32_t b = 0; b < m_kbytes; ++b) rec[1 + b] = (uint8_t)(kl >> (8 * b));
                for (uint32_t b = 0; b < m_kbytes; ++b)
                    rec[1 + m_kbytes + b] = (uint8_t)(kr >> (8 * b));
                put(lbuf, rec, (size_t)1 + 2 * m_kbytes);
            }
            ++add_count;
            add_seq += seq_len;
        }
        {
            std::lock_guard<std::mutex> lk(m_mu);
            if (!fbuf.empty() and std::fwrite(fbuf.data(), 1, fbuf.size(), m_file) != fbuf.size())
                throw std::runtime_error("short write to " + m_path);
            if (m_links_file and !lbuf.empty() and
                std::fwrite(lbuf.data(), 1, lbuf.size(), m_links_file) != lbuf.size())
                throw std::runtime_error("short write to " + m_links_path);
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
    // Companion links spill path, or "" when k < 2 (links disabled).
    std::string const& links_path() const { return m_links_path; }

    // Close the writer side. Call before reads.
    void close_for_writing() {
        std::lock_guard<std::mutex> lk(m_mu);
        if (m_file) {
            std::fflush(m_file);
            std::fclose(m_file);
            m_file = nullptr;
        }
        if (m_links_file) {
            std::fflush(m_links_file);
            std::fclose(m_links_file);
            m_links_file = nullptr;
        }
    }

    // Read every record back into `out`. Reserves up-front from
    // m_count so no reallocation churn. The frag spill file is
    // sequential-access friendly so this is one big sequential read.
    void load_to_vector(std::vector<stitchable_unitig>& out) {
        out.clear();
        out.reserve(m_count);
        std::FILE* f = std::fopen(m_path.c_str(), "rb");
        if (!f) throw std::runtime_error("cannot reopen frag spill: " + m_path);
        for (;;) {
            uint8_t flags = 0;
            size_t got = std::fread(&flags, sizeof(flags), 1, f);
            if (got != 1) {
                if (std::feof(f)) break;
                std::fclose(f);
                throw std::runtime_error("short read of flags from " + m_path);
            }
            uint32_t nruns = 0;
            if (std::fread(&nruns, sizeof(nruns), 1, f) != 1) {
                std::fclose(f);
                throw std::runtime_error("short read of nruns from " + m_path);
            }
            stitchable_unitig u;
            u.runs.resize(nruns);
            for (uint32_t r = 0; r < nruns; ++r) {
                if (std::fread(&u.runs[r].cid, sizeof(uint64_t), 1, f) != 1 or
                    std::fread(&u.runs[r].num_kmers, sizeof(uint32_t), 1, f) != 1) {
                    std::fclose(f);
                    throw std::runtime_error("short read of run from " + m_path);
                }
            }
            uint32_t seq_len = 0;
            if (std::fread(&seq_len, sizeof(seq_len), 1, f) != 1) {
                std::fclose(f);
                throw std::runtime_error("short read of seq_len from " + m_path);
            }
            u.open_flags = flags;
            u.seq.resize(seq_len);
            if (seq_len > 0 and std::fread(u.seq.data(), 1, seq_len, f) != (size_t)seq_len) {
                std::fclose(f);
                throw std::runtime_error("short read of seq from " + m_path);
            }
            out.push_back(std::move(u));
        }
        std::fclose(f);
    }

    void unlink() {
        std::error_code ec;
        if (!m_path.empty()) {
            std::filesystem::remove(m_path, ec);
            m_path.clear();
        }
        if (!m_links_path.empty()) {
            std::filesystem::remove(m_links_path, ec);
            m_links_path.clear();
        }
    }

private:
    std::string m_path;
    std::FILE* m_file = nullptr;
    std::mutex m_mu;
    uint64_t m_count = 0;
    uint64_t m_total_seq_bytes = 0;
    // Companion links spill (boundary k-mers per fragment). Disabled when m_k==0.
    uint32_t m_k = 0;
    uint32_t m_kbytes = 0;
    std::string m_links_path;
    std::FILE* m_links_file = nullptr;
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
