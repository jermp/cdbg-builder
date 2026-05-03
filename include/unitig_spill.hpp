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
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

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

}  // namespace cdgb
