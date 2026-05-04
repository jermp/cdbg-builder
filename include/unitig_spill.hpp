#pragma once

// Disk-backed unitig sink for the stitch -> emit boundary.
//
// stitch_unitigs_streaming feeds each finished merged unitig into a
// caller-provided sink. We use this writer as that sink so the seq
// strings never have to all coexist in RAM during emit_fasta.
//
// Bucketing by cid range:
//   We partition cid -> bucket so that bucket b holds every unitig
//   with cid in [b * S, (b+1) * S), where S = ceil(num_color_classes /
//   num_buckets). bucket 0 has the smallest cids; bucket num_buckets-1
//   has the largest. emit_fasta iterates buckets in ascending order
//   and within each bucket sorts the records by cid; the resulting
//   .fa output order is therefore exactly cid-ascending, which is
//   what the u2c bit_vector (downstream consumer's run-end marker)
//   requires.
//
// Per-bucket file format: a sequence of records, one per finished
// unitig, until EOF.
//   [u32 cid]
//   [u32 seq_len]
//   [seq_len bytes]    raw ACGT (the stitch sink writes the merged
//                      unitig's seq bytes verbatim)
//
// Per-bucket peak in-memory footprint at read time is one bucket's
// worth of seq bytes -- with K = 64 buckets and 1.88 M total unitigs
// totalling ~190 MB of seq, that's ~3 MB per bucket plus ~400 KB of
// per-record bookkeeping. Comfortably bounded.

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

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "bucket_walker.hpp"  // stitchable_unitig

namespace cdgb {

class unitig_bucket_writer {
public:
    // num_color_classes is the final cid space (size of global_dict).
    // num_buckets controls how the cid range is partitioned. Picked
    // so per-bucket peak fits comfortably in the budget remainder.
    unitig_bucket_writer(std::string dir, uint32_t num_color_classes, uint32_t num_buckets)
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
        const uint32_t b = bucket_for_cid(u.cid);
        std::FILE* f = m_files[b];
        const uint32_t cid = u.cid;
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
    uint32_t num_color_classes() const { return m_num_color_classes; }

    std::string bucket_path(uint32_t b) const {
        return m_dir + "/unitig_bucket_" + std::to_string(b) + ".bin";
    }

    // Read all records from bucket `b`. Caller-owned `out` is cleared
    // and refilled. Records arrive in *write order* (which is stitch
    // order, NOT cid order) -- caller sorts.
    struct record {
        uint32_t cid;
        std::string seq;
    };

    void read_bucket(uint32_t b, std::vector<record>& out) {
        out.clear();
        std::FILE* f = m_files[b];
        if (!f) throw std::runtime_error("read_bucket: file already closed");
        std::fflush(f);
        std::rewind(f);
        for (;;) {
            uint32_t cid = 0, seq_len = 0;
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
    uint32_t bucket_for_cid(uint32_t cid) const {
        if (m_num_color_classes == 0 || m_num_buckets <= 1) return 0;
        // Cids in [b * S, (b+1) * S) live in bucket b, where
        // S = ceil(num_color_classes / num_buckets). Last bucket
        // mops up any rounding remainder.
        const uint64_t S =
            ((uint64_t)m_num_color_classes + m_num_buckets - 1) / m_num_buckets;
        uint32_t b = (uint32_t)((uint64_t)cid / S);
        if (b >= m_num_buckets) b = m_num_buckets - 1;
        return b;
    }

    std::string m_dir;
    uint32_t m_num_color_classes;
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
        const uint32_t cid = u.cid;
        const uint8_t flags = u.open_flags;
        const uint32_t seq_len = (uint32_t)u.seq.size();
        if (std::fwrite(&cid, sizeof(cid), 1, m_file) != 1 ||
            std::fwrite(&flags, sizeof(flags), 1, m_file) != 1 ||
            std::fwrite(&seq_len, sizeof(seq_len), 1, m_file) != 1)
            throw std::runtime_error("short write to " + m_path);
        if (seq_len > 0 &&
            std::fwrite(u.seq.data(), 1, seq_len, m_file) != (size_t)seq_len)
            throw std::runtime_error("short write to " + m_path);
        ++m_count;
        // Free the merged seq's backing storage in place: the caller
        // already moved into us.
        std::string().swap(u.seq);
    }

    uint64_t count() const { return m_count; }
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
            uint32_t cid = 0;
            size_t got = std::fread(&cid, sizeof(cid), 1, f);
            if (got != 1) {
                if (std::feof(f)) break;
                std::fclose(f);
                throw std::runtime_error("short read of cid from " + m_path);
            }
            uint8_t flags = 0;
            uint32_t seq_len = 0;
            if (std::fread(&flags, sizeof(flags), 1, f) != 1 ||
                std::fread(&seq_len, sizeof(seq_len), 1, f) != 1) {
                std::fclose(f);
                throw std::runtime_error("short read of header from " + m_path);
            }
            stitchable_unitig u;
            u.cid = cid;
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
};

// ----------------------------------------------------------------------------
// mmap-backed read-only view over a frag_unitig_writer's spill file.
// Builds an in-memory index of (cid, open_flags, seq_offset, seq_len) per
// fragment and mmaps the file so seq bytes are accessible as
// std::string_view directly into the kernel page cache. With this, stitch
// can run without ever loading all fragment seqs into RAM as
// std::string's: the OS pages in seq bytes on demand and reclaims them
// under memory pressure.
//
// Memory cost: ~24 B per fragment for the index; for 26.6 M fragments
// that's ~640 MB (vs ~2.8 GB for an in-memory std::vector<stitchable_unitig>
// at the same scale). Mmap pages don't count toward RSS until accessed
// and dirty pages aren't ours -- the kernel can evict file-backed read-
// only pages without involving us.
//
// Used by the stitch phase via the templated stitch_unitigs_streaming.
// frag_unitig_reader satisfies the Source interface (size / cid /
// open_flags / seq_view).

class frag_unitig_reader {
public:
    explicit frag_unitig_reader(std::string path) : m_path(std::move(path)) {
        m_fd = ::open(m_path.c_str(), O_RDONLY);
        if (m_fd < 0)
            throw std::runtime_error("open failed for frag spill: " + m_path + ": " +
                                     std::strerror(errno));
        struct stat st;
        if (::fstat(m_fd, &st) < 0) {
            ::close(m_fd);
            throw std::runtime_error("fstat failed on " + m_path);
        }
        m_mmap_size = (size_t)st.st_size;
        if (m_mmap_size > 0) {
            void* p = ::mmap(nullptr, m_mmap_size, PROT_READ, MAP_PRIVATE, m_fd, 0);
            if (p == MAP_FAILED) {
                ::close(m_fd);
                throw std::runtime_error("mmap failed on " + m_path);
            }
            m_mmap_base = (uint8_t const*)p;
            // Walk_chain accesses fragments in chain order, which has
            // no spatial locality with the on-disk layout. MADV_RANDOM
            // tells the kernel to stop prefetching surrounding pages
            // and to reclaim resident pages more aggressively under
            // memory pressure -- both reduce stitch-phase peak RSS.
            ::madvise((void*)m_mmap_base, m_mmap_size, MADV_RANDOM);
        }
        // Walk the file once to populate the index. Records are
        // [u32 cid][u8 flags][u32 seq_len][seq_len bytes]. Sequential
        // access -- kernel readahead handles the I/O cost.
        size_t off = 0;
        while (off < m_mmap_size) {
            if (off + 4 + 1 + 4 > m_mmap_size)
                throw std::runtime_error("truncated frag spill header in " + m_path);
            uint32_t cid;
            std::memcpy(&cid, m_mmap_base + off, 4);
            off += 4;
            uint8_t flags = m_mmap_base[off];
            off += 1;
            uint32_t seq_len;
            std::memcpy(&seq_len, m_mmap_base + off, 4);
            off += 4;
            if (off + seq_len > m_mmap_size)
                throw std::runtime_error("truncated frag spill seq in " + m_path);
            if (seq_len > MAX_SEQ_LEN)
                throw std::runtime_error("seq_len exceeds 30-bit cap in " + m_path);
            entry e;
            e.cid = cid;
            e.flags_and_len = ((uint32_t)flags << 30) | (seq_len & MAX_SEQ_LEN);
            e.seq_offset = off;
            m_entries.push_back(e);
            off += seq_len;
        }
    }

    ~frag_unitig_reader() {
        if (m_mmap_base) {
            ::munmap((void*)m_mmap_base, m_mmap_size);
        }
        if (m_fd >= 0) {
            ::close(m_fd);
        }
    }

    frag_unitig_reader(frag_unitig_reader const&) = delete;
    frag_unitig_reader& operator=(frag_unitig_reader const&) = delete;

    // Source interface for stitch_unitigs_streaming.
    size_t size() const { return m_entries.size(); }
    uint32_t cid(size_t i) const { return m_entries[i].cid; }
    uint8_t open_flags(size_t i) const {
        return (uint8_t)(m_entries[i].flags_and_len >> 30);
    }
    std::string_view seq_view(size_t i) const {
        auto const& e = m_entries[i];
        return std::string_view((char const*)m_mmap_base + e.seq_offset,
                                e.flags_and_len & MAX_SEQ_LEN);
    }

private:
    // 16-byte packed entry: cid (4) + flags-and-len (4: 2 bits flags +
    // 30 bits seq_len) + seq_offset (8). Saves 8 B/entry vs the
    // natural-aligned 24-byte struct -- on 26.6 M fragments that's
    // ~213 MB of index RAM.
    static constexpr uint32_t MAX_SEQ_LEN = (1u << 30) - 1;
    struct entry {
        uint32_t cid;
        uint32_t flags_and_len;  // top 2 bits: open_flags; bottom 30 bits: seq_len
        uint64_t seq_offset;
    };
    static_assert(sizeof(entry) == 16, "frag_unitig_reader::entry must be 16 bytes");

    std::string m_path;
    int m_fd = -1;
    uint8_t const* m_mmap_base = nullptr;
    size_t m_mmap_size = 0;
    std::vector<entry> m_entries;
};

}  // namespace cdgb
