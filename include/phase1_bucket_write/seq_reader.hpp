#pragma once

// FASTA / FASTQ reader with STREAMING zlib-based gzip decompression.
//
// Replaces the previous mmap + libdeflate whole-buffer backend. That
// design decompressed the ENTIRE file into one RAM buffer (libdeflate
// has no streaming API) and kept the whole compressed file mmap'd, so
// per-reader memory scaled with file size: fine on bacterial genomes
// (~5-10 MiB), catastrophic on human-scale pangenome inputs -- on
// HPRC-like collections (~0.8 GiB gzipped / ~3 GiB raw per genome,
// growth-doubling on top) it cost 4-7 GiB PER THREAD, i.e. hundreds of
// GiB of RSS across a 32-48 thread ingest. The mmap'd compressed pages
// additionally inflated RSS and misled the RSS-watcher / RAM governor.
//
// This backend streams instead: zlib's gzread inflates a bounded window
// at a time (it also transparently passes through plain uncompressed
// files), and kseq parses records incrementally from that stream.
// Per-reader memory is now O(zlib state + gz buffer + largest single
// record), independent of file size. The largest-record term is
// unavoidable with this API (next() hands out one whole sequence);
// for a human chromosome that's ~250 MiB, vs multiple GiB before.
//
// Cost: zlib inflate is ~2-3x slower than libdeflate. Ingest threads
// each decompress their own file in parallel and the pipeline is
// disk-bound at scale, so the wall-time impact stays small; bounded
// memory at pangenome scale is worth far more than that.

#include <stdexcept>
#include <string>

#include <zlib.h>

extern "C" {
#include "external/kseq.h"
}

// Instantiate kseq directly over zlib's streaming gzFile. gzread
// returns raw bytes for non-gzip inputs, so plain FASTA works through
// the same path with no sniffing here.
KSEQ_INIT(gzFile, gzread)

namespace cdbg {

struct seq_reader {
    explicit seq_reader(std::string const& path) {
        m_file = gzopen(path.c_str(), "rb");
        if (!m_file) {
            throw std::runtime_error("could not open input file: " + path);
        }
        // Default gz buffer is 8 KiB; a larger window cuts the number of
        // inflate calls and read syscalls on multi-GiB inputs.
        gzbuffer(m_file, 128 * 1024);
        m_kseq = kseq_init(m_file);
    }

    ~seq_reader() {
        if (m_kseq) kseq_destroy(m_kseq);
        if (m_file) gzclose(m_file);
    }

    seq_reader(seq_reader const&) = delete;
    seq_reader& operator=(seq_reader const&) = delete;

    // Returns true and fills out_seq + out_len on success; false on EOF.
    // Pointer is owned by the reader and invalidated by the next call.
    bool next(char const*& out_seq, size_t& out_len) {
        int r = kseq_read(m_kseq);
        if (r < 0) return false;
        out_seq = m_kseq->seq.s;
        out_len = m_kseq->seq.l;
        return true;
    }

private:
    gzFile m_file = nullptr;
    kseq_t* m_kseq = nullptr;
};

}  // namespace cdbg
