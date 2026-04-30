#pragma once

#include <stdexcept>
#include <string>
#include <zlib.h>

extern "C" {
#include "external/kseq.h"
}

// Instantiate kseq for gzFile readers (handles plain text and .gz transparently).
KSEQ_INIT(gzFile, gzread)

namespace cdgb {

// Iterates sequences in a FASTA / FASTQ / .gz file. Owns the gzFile handle.
struct SeqReader {
    SeqReader(const std::string& path) {
        m_fp = gzopen(path.c_str(), "r");
        if (!m_fp) throw std::runtime_error("could not open input file: " + path);
        m_seq = kseq_init(m_fp);
    }

    ~SeqReader() {
        if (m_seq) kseq_destroy(m_seq);
        if (m_fp) gzclose(m_fp);
    }

    SeqReader(const SeqReader&) = delete;
    SeqReader& operator=(const SeqReader&) = delete;

    // Returns true and fills out_seq + out_len on success; false on EOF.
    // Pointer is owned by the reader and invalidated by the next call.
    bool next(const char*& out_seq, size_t& out_len) {
        int r = kseq_read(m_seq);
        if (r < 0) return false;
        out_seq = m_seq->seq.s;
        out_len = m_seq->seq.l;
        return true;
    }

private:
    gzFile m_fp = nullptr;
    kseq_t* m_seq = nullptr;
};

}  // namespace cdgb
