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

#include "bucket_walker.hpp"  // stitchable_unitig

namespace cdbg {

class unitig_bucket_writer {
public:
    // num_color_classes is the final cid space (size of global_dict).
    // num_buckets controls how the cid range is partitioned. Picked
    // so per-bucket peak fits comfortably in the budget remainder.
    unitig_bucket_writer(std::string dir, uint64_t num_color_classes, uint32_t num_buckets)
        : m_dir(std::move(dir)),
          m_num_color_classes(num_color_classes),
          m_num_buckets(num_buckets == 0 ? 1 : num_buckets) {
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
        if (seq_len > 0 &&
            std::fwrite(u.seq.data(), 1, seq_len, f) != (size_t)seq_len) {
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
                throw std::runtime_error("short read of unitig cid in bucket " +
                                         std::to_string(b));
            }
            if (std::fread(&seq_len, sizeof(seq_len), 1, f) != 1)
                throw std::runtime_error("short read of unitig len in bucket " +
                                         std::to_string(b));
            std::string seq;
            seq.resize(seq_len);
            if (seq_len > 0 &&
                std::fread(seq.data(), 1, seq_len, f) != (size_t)seq_len) {
                throw std::runtime_error("short read of unitig seq in bucket " +
                                         std::to_string(b));
            }
            out.push_back({cid, std::move(seq)});
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
        if (m_num_color_classes == 0 || m_num_buckets <= 1) return 0;
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
    explicit frag_unitig_writer(std::string path)
        : m_path(std::move(path)) {
        m_file = std::fopen(m_path.c_str(), "wb+");
        if (!m_file)
            throw std::runtime_error("cannot open frag spill: " + m_path);
    }

    ~frag_unitig_writer() {
        if (m_file) std::fclose(m_file);
    }

    frag_unitig_writer(frag_unitig_writer const&) = delete;
    frag_unitig_writer& operator=(frag_unitig_writer const&) = delete;

    void operator()(stitchable_unitig&& u) {
        std::lock_guard<std::mutex> lk(m_mu);
        const uint8_t flags = u.open_flags;
        const uint32_t nruns = (uint32_t)u.runs.size();
        const uint32_t seq_len = (uint32_t)u.seq.size();
        if (std::fwrite(&flags, sizeof(flags), 1, m_file) != 1 ||
            std::fwrite(&nruns, sizeof(nruns), 1, m_file) != 1)
            throw std::runtime_error("short write to " + m_path);
        for (auto const& r : u.runs) {
            if (std::fwrite(&r.cid, sizeof(r.cid), 1, m_file) != 1 ||
                std::fwrite(&r.num_kmers, sizeof(r.num_kmers), 1, m_file) != 1)
                throw std::runtime_error("short write to " + m_path);
        }
        if (std::fwrite(&seq_len, sizeof(seq_len), 1, m_file) != 1)
            throw std::runtime_error("short write to " + m_path);
        if (seq_len > 0 &&
            std::fwrite(u.seq.data(), 1, seq_len, m_file) != (size_t)seq_len)
            throw std::runtime_error("short write to " + m_path);
        ++m_count;
        m_total_seq_bytes += seq_len;
        // Free the merged seq's backing storage in place: the caller
        // already moved into us.
        std::string().swap(u.seq);
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
                if (std::fread(&u.runs[r].cid, sizeof(uint64_t), 1, f) != 1 ||
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
            if (seq_len > 0 &&
                std::fread(u.seq.data(), 1, seq_len, f) != (size_t)seq_len) {
                std::fclose(f);
                throw std::runtime_error("short read of seq from " + m_path);
            }
            out.push_back(std::move(u));
        }
        std::fclose(f);
    }

    void unlink() {
        if (m_path.empty()) return;
        std::error_code ec;
        std::filesystem::remove(m_path, ec);
        m_path.clear();
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
            if (std::fread(&runs[r].cid, sizeof(uint64_t), 1, m_file) != 1 ||
                std::fread(&runs[r].num_kmers, sizeof(uint32_t), 1, m_file) != 1)
                throw std::runtime_error("short read of run from " + m_path);
        }
        uint32_t seq_len = 0;
        if (std::fread(&seq_len, sizeof(seq_len), 1, m_file) != 1)
            throw std::runtime_error("short read of seq_len from " + m_path);
        open_flags = flags;
        seq.resize(seq_len);
        if (seq_len > 0 &&
            std::fread(seq.data(), 1, seq_len, m_file) != (size_t)seq_len)
            throw std::runtime_error("short read of seq from " + m_path);
        return true;
    }

private:
    std::string m_path;
    std::FILE* m_file = nullptr;
};


}  // namespace cdbg
