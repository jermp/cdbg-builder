#pragma once

// FASTA / FASTQ reader with libdeflate-based gzip decompression.
//
// Replaces the previous zlib (gzopen/gzread) backend. libdeflate's
// decompressor is ~2-3x faster than zlib's on bacterial-genome inputs;
// since seq_read was 33% of bucket-write wall time on the 50K workload,
// this directly cuts ~half of that.
//
// libdeflate has no streaming API by design (whole-buffer in, whole-
// buffer out), so we do:
//   1. mmap the file (zero-copy, no RAM cost for the compressed bytes).
//   2. If gzip-magic, libdeflate_gzip_decompress_ex into a per-reader
//      output buffer; loop for multi-member streams.
//   3. Feed kseq from a memory_stream that just memcpys out of the
//      decompressed buffer (or directly out of the mmap if the input
//      wasn't gzipped).
//
// Memory: ~5-10 MiB per concurrent reader for a typical bacterial
// genome. With 32 ingest threads that's ~250 MiB extra, comfortably
// within the bucket-write budget (which is dominated by per-bucket
// hashmaps). The mmap costs nothing -- it's page-cache-backed.

#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include <libdeflate.h>

namespace cdbg {

// Memory stream type used by kseq's KSEQ_INIT below: a pointer + size +
// position. mem_read copies up to `len` bytes into the caller's buffer
// and advances pos. Same shape as POSIX read(2) so kseq can plug into
// it via KSEQ_INIT.
struct mem_stream {
    uint8_t const* data = nullptr;
    size_t size = 0;
    size_t pos = 0;
};

inline ssize_t mem_stream_read(mem_stream* s, void* buf, size_t len) {
    if (!s) return -1;
    if (s->pos >= s->size) return 0;
    size_t avail = s->size - s->pos;
    size_t n = (len < avail) ? len : avail;
    std::memcpy(buf, s->data + s->pos, n);
    s->pos += n;
    return (ssize_t)n;
}

}  // namespace cdbg

extern "C" {
#include "external/kseq.h"
}

// Instantiate kseq for our memory stream. kseq's __KS_GETC fast path
// inlines mem_stream_read calls; for repeated 4096-byte buffer fills
// this collapses to memcpy + pos bookkeeping, which the optimizer
// reduces to a tight loop over the decompressed buffer.
KSEQ_INIT(cdbg::mem_stream*, cdbg::mem_stream_read)

namespace cdbg {

struct seq_reader {
    explicit seq_reader(std::string const& path) {
        // mmap the input so we don't pay a read() copy for the
        // compressed bytes. Also lets the kernel keep the file in the
        // page cache across runs.
        int fd = ::open(path.c_str(), O_RDONLY);
        if (fd < 0) {
            throw std::runtime_error("could not open input file: " + path + ": " +
                                     std::strerror(errno));
        }
        struct stat st;
        if (::fstat(fd, &st) != 0) {
            int e = errno;
            ::close(fd);
            throw std::runtime_error("fstat failed on " + path + ": " + std::strerror(e));
        }
        m_in_size = (size_t)st.st_size;
        if (m_in_size > 0) {
            void* p = ::mmap(nullptr, m_in_size, PROT_READ, MAP_PRIVATE, fd, 0);
            if (p == MAP_FAILED) {
                int e = errno;
                ::close(fd);
                throw std::runtime_error("mmap failed on " + path + ": " + std::strerror(e));
            }
            m_in = (uint8_t const*)p;
            // Hint sequential access -- input is read once front-to-back.
            ::madvise((void*)m_in, m_in_size, MADV_SEQUENTIAL);
        }
        ::close(fd);

        // gzip magic 1f 8b? Empty files fall through as "not gzipped";
        // kseq will then immediately hit EOF.
        bool is_gz = m_in_size >= 2 and m_in[0] == 0x1f and m_in[1] == 0x8b;

        if (is_gz) {
            decompress_gzip_(path);
            m_stream.data = m_decompressed.data();
            m_stream.size = m_decompressed.size();
        } else {
            // Plain text / FASTA -- feed kseq directly from the mmap.
            m_stream.data = m_in;
            m_stream.size = m_in_size;
        }
        m_stream.pos = 0;
        m_kseq = kseq_init(&m_stream);
    }

    ~seq_reader() {
        if (m_kseq) kseq_destroy(m_kseq);
        if (m_in and m_in_size > 0) ::munmap((void*)m_in, m_in_size);
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
    // Decompress a gzip-wrapped buffer into m_decompressed using
    // libdeflate. Handles multi-member gzip via the _ex variant: each
    // call decompresses one member and reports input bytes consumed;
    // we loop until all input is consumed. Output buffer grows
    // geometrically on LIBDEFLATE_INSUFFICIENT_SPACE (libdeflate has
    // no streaming, so the whole decompressed result must fit).
    //
    // Initial guess: 4x compressed size. Bacterial-genome FASTA gzip
    // ratios are typically 3.0-3.5x, so 4x lands in one shot most of
    // the time.
    void decompress_gzip_(std::string const& path) {
        struct libdeflate_decompressor* dec = libdeflate_alloc_decompressor();
        if (!dec) throw std::runtime_error("libdeflate_alloc_decompressor failed");

        size_t out_cap = m_in_size * 4 + 16;
        m_decompressed.resize(out_cap);
        size_t out_off = 0;
        size_t in_off = 0;

        while (in_off < m_in_size) {
            size_t in_consumed = 0;
            size_t out_produced = 0;
            for (;;) {
                if (m_decompressed.size() - out_off < 4096) {
                    m_decompressed.resize(m_decompressed.size() * 2 + 4096);
                }
                enum libdeflate_result r = libdeflate_gzip_decompress_ex(
                    dec, m_in + in_off, m_in_size - in_off, m_decompressed.data() + out_off,
                    m_decompressed.size() - out_off, &in_consumed, &out_produced);
                if (r == LIBDEFLATE_SUCCESS) break;
                if (r == LIBDEFLATE_INSUFFICIENT_SPACE) {
                    m_decompressed.resize(m_decompressed.size() * 2 + 4096);
                    continue;
                }
                libdeflate_free_decompressor(dec);
                throw std::runtime_error("libdeflate_gzip_decompress failed on " + path +
                                         " (result=" + std::to_string((int)r) + ")");
            }
            in_off += in_consumed;
            out_off += out_produced;
            if (in_consumed == 0) break;  // safety: no progress (shouldn't happen)
        }

        libdeflate_free_decompressor(dec);
        m_decompressed.resize(out_off);
    }

    uint8_t const* m_in = nullptr;
    size_t m_in_size = 0;
    std::vector<uint8_t> m_decompressed;
    mem_stream m_stream;
    kseq_t* m_kseq = nullptr;
};

}  // namespace cdbg
