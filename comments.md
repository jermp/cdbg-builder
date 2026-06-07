in bucket_ingester.hpp

the bucket assignment is currently using a modulo:

```
auto bucket_of = [&](uint64_t h) -> uint32_t {
    // Skip the bottom bit to match GGCAT's "uniqueness-flag" reservation.
    // num_buckets need not be a power of two -- it's sized from the RAM
    // model (B = frac*g / (alpha*T*flush + beta*spill)), so map with a
    // modulo. This runs once per super-k-mer boundary (not per base), so
    // the divide is negligible against the per-base ntHash/minimizer scan.
    return (uint32_t)((h >> 1) % num_buckets);
};
```

but could use a much faster way using a remap, see here: https://github.com/jermp/pthash/blob/master/include/utils/util.hpp

