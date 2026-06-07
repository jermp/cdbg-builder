# Comments

file: builder.hpp

```
builder() = default;
explicit builder(build_config const& cfg) : m_cfg(cfg) {}
```

why "default" and "explicit"?


m_num_colors should hust be uint64_t throughout the codebase; we enforce it's always < 2^32 but keep it 64 bit wise so that we don't have to cast



instead of this:

```
auto const t_start = std::chrono::steady_clock::now();
auto print_total = [&] {
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::steady_clock::now() - t_start)
                  .count();
    std::cout << "[total construction time] " << (ms / 1000.0) << " s\n";
};
```

could use the clock object in external/bits/external/essentials

```
auto writer =
    std::make_unique<bucket_writer>(tmp_dir, num_buckets, m_flush_bases, m_spill_bytes);
```

why the std::make_unique? why a pointer is needed here? to automatically call the destructor of bucket_writer? we can also do it with a scope, no?

same comment for this `std::unique_ptr<progress> prog`. Why do we need the unique_ptr?

"When -g is set" --- but -g must ALWAYS be set, no?

`std::vector<uint64_t> file_sizes` -- why do we need all sizes, file by file and not just the total sum?

```
std::unique_ptr<progress> prog =
total_bytes > 0
    ? std::make_unique<progress>("bucket-write", done_bytes, total_bytes,
                                 /*render_bytes=*/true, &done,
                                 (uint64_t)files.size(), "files")
    : std::make_unique<progress>("bucket-write", done, (uint64_t)files.size());
```

but total_bytes is always > 0, no? how can i be 0?

pleae don't use size_t but just uint64_t

`std::unique_ptr<unitig_bucket_writer> uwriter_ptr;`, another unique_ptr...why?

in this comment:
```
// Fraction of -g the bucket-write phase is sized against. 0.50 matches the
// long-standing BUCKET_WRITE_SHARE that kept bucket-write lean. Phase 1
// doesn't need most of -g -- only enough that B = M/(alpha*T*flush+beta*
// spill) is large enough for a 48-way bucket-process.
```
you mean, for a "T-way bucket process", no? why 48?

don't use `||` but `or`, it's modern C++. Likewise for `&&` and `and`.

```
uint64_t count = b < 1.0 ? 1 : (uint64_t)b;
if (count < MIN_AUTO_BUCKETS) count = MIN_AUTO_BUCKETS;
if (count > UINT32_MAX) count = UINT32_MAX;
return (uint32_t)count;
```

don't use cast please... it's ugly.

if this function `static std::vector<std::string> read_filenames(std::string const& path)` is only used once, then subsitute the body directly where it is called.

same comment as above for this function `static void cleanup_tmp_dir(std::string const& tmp_dir)`.

```
static uint32_t pick_unitig_bucket_count_(uint64_t num_color_classes,
                                          uint64_t total_seq_bytes_estimate,
                                          double max_ram_gb) {
    if (num_color_classes == 0) return 1;
    constexpr uint32_t MIN_K = 16;
    constexpr uint32_t MAX_K = 1024;
    constexpr uint32_t DEFAULT_K = 64;
    constexpr double OVERHEAD = 2.0;
    constexpr double SHARE = 0.10;

    if (max_ram_gb <= 0 || total_seq_bytes_estimate == 0) {
        return (uint32_t)std::min<uint64_t>(num_color_classes, DEFAULT_K);
    }
    const uint64_t budget_bytes = (uint64_t)(max_ram_gb * 1024.0 * 1024.0 * 1024.0 * SHARE);
    if (budget_bytes == 0) {
        return (uint32_t)std::min<uint64_t>(num_color_classes, DEFAULT_K);
    }
    const uint64_t needed = (uint64_t)((double)total_seq_bytes_estimate * OVERHEAD);
    uint64_t k = (needed + budget_bytes - 1) / budget_bytes;
    if (k < MIN_K) k = MIN_K;
    if (k > MAX_K) k = MAX_K;
    if (k > num_color_classes) k = num_color_classes;
    return (uint32_t)k;
}
```

also here avoid the casting, it's ugly. directly use the proper type.
Also, in this function `pick_stitch_buckets_`.

