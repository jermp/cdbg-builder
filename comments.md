- file seq_reader.hpp:

why mmap is used here? careful that mmap can silently make the RSS grow if RAM is available. it would be good to remove it entirely from the codebase. Is it used somewhere else?

- file bucket_ingester.hpp:

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
