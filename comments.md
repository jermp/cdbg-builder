- file kmer.hpp:

are these functions

`reverse_complement`
`canonical`
`is_canonical`
`string_to_kmer`

used in hot loops?

Beware that they might be costly and reverse complements can be computed inplace and in O(1); see here:
https://github.com/jermp/sshash/blob/master/include/kmer.hpp#L159.

- file super_kmer.hpp:

is this function `inline void write_super_kmer(uint8_t flags, uint32_t const* colors, uint32_t num_colors, uint8_t const* bases, uint32_t len, std::vector<uint8_t>& out)` used anymore?

- file seq_reader.hpp:

why mmap is used here? careful that mmap can silently make the RSS grow if RAM is available. it would be good to remove it entirely from the codebase. Is it used somewhere else?

- file bucket_ingester.hpp:

in the function `ingest_file_bucketed`, we have these two hot loops:

```
....
while (end < l and nuc_to_2bit(s[end]) != 0xff) ++end;
size_t run_len = end - pos;
bool have_run = run_len >= k;
if (have_run) {
    bases_buf.resize(run_len);
    for (size_t i = 0; i < run_len; ++i) { bases_buf[i] = nuc_to_2bit(s[pos + i]); }
}
...
```

but clearly they can be merged into one, we can write directly to `bases_buf` the transformed values and stop when we hit a non-ACGT char.
If then the run is >= k, we go ahead, otherwise we reset the buffer. that is, we avoid having the while+for.

Also, we can avoid this call `bases_buf.resize(run_len);` but just keep track of the current begin into bases_buf and current run_len.

WELL, reading ahead the function `emit_super_kmers`: since super-kmers are computed over bytes, not a 2-bit encoded stream, why do we need to first copy the bytes into the local buffer bases_buf??? we can use directly the buffer held inside the seq_reader object, no?
So this buffer is directly passed to the `emit_super_kmers` function that directly produces super-kmers, having care to stop when a non-ACGT char is met. Otherwise we keep copying bytes twice and transforming ACGT into 0123 which is not necessary, we can use these maps
(like the traditional encoding map since you're using this one in the codebase):

#ifdef SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING
    /*
    char decimal  binary
        A     65     01000001 -> 00
        C     67     01000011 -> 01
        G     71     01000111 -> 10
        T     84     01010100 -> 11

        a     97     01100001 -> 00
        c     99     01100011 -> 01
        g    103     01100111 -> 10
        t    116     01110100 -> 11
    */
    static uint64_t char_to_uint(char c) { return (((c >> 1) ^ (c >> 2)) & 3); }
#else
    /*
    char decimal  binary
        A     65     01000.00.1 -> 00
        C     67     01000.01.1 -> 01
        G     71     01000.11.1 -> 11
        T     84     01010.10.0 -> 10

        a     97     01100.00.1 -> 00
        c     99     01100.01.1 -> 01
        g    103     01100.11.1 -> 11
        t    116     01110.10.0 -> 10
    */
    static uint64_t char_to_uint(char c) { return (c >> 1) & 3; }
#endif

in the function `emit_super_kmers`, we have this test `if (L < k) return;` but we already enforced this, so we should actually put an assert
like assert(L >= k).

file bucket_io.hpp:

struct pending_record {
    uint64_t hash;
    uint32_t color;
    uint32_t key_off;
    uint32_t key_len;
    uint8_t flags;
};
```

currently weights 21 bytes, so we are wasting 3 bytes since sizeof will be 24.


we could use pragma pack directive to pack the structure to 1 byte to avoid wasting space here

```
struct per_thread_bucket_buffers {
    std::vector<std::vector<bucket_compactor::pending_record>> recs;
    ...
};
```

for example:

#pragma pack(push, 1)
struct pending_record {
    uint64_t hash;
    uint32_t color;
    uint32_t key_off;
    uint32_t key_len;
    uint8_t flags;
};
#pragma pack(pop)

or can this layout be optimized?