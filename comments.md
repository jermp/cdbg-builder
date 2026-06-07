are kmer.hpp, minimizer.hpp, seq_reader.hpp, super_kmer.hpp and util.hpp used across all phases?
if so, they should be put in the same directoy as builder.hpp

in the kmer.hpp file:
`using kmer_int_t = __uint128_t;` -- would defaulting to `using kmer_int_t = uint64_t;` speed up things?
In SSHash, I keep k<=31 as default and use uint64_t, then only at compile time, I change that to __uint128_t to allow k<=63.

in the minimizer.hpp file:
where is `min_queue` used? didn't we say that a re-scan method is actually faster than a monotonic queue-based approach?
See here https://github.com/jermp/sshash/blob/master/include/minimizer_iterator.hpp.

