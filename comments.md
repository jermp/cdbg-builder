comments on util.hpp:

- still size_t around: replace them all with uint64_t or explain why size_t is better there

```
// RAII phase timer: prints "[label] X.XX s" on destruction. Used to bracket
// individual stages inside cdbg_builder::build().
class timer {
public:
    timer(char const* label) : m_label(label), m_t0(std::chrono::steady_clock::now()) {}
    ~timer() {
        auto t1 = std::chrono::steady_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - m_t0).count();
        std::cout << "[" << m_label << "] " << (ms / 1000.0) << " s\n";
    }

    timer(timer const&) = delete;
    timer& operator=(timer const&) = delete;

private:
    char const* m_label;
    std::chrono::steady_clock::time_point m_t0;
};
```

is this class used now that we use essentials::timer? clearly we can get by without...


- be sure this function

```
// Pretty-printer: 3.42 GiB / 728 MiB / 12 KiB, picking the largest
// unit at which the number is >= 1.
inline std::string format_bytes(uint64_t b) {
    static char const* units[] = {"B", "KiB", "MiB", "GiB", "TiB"};
    double v = b;
    int u = 0;
    while (v >= 1024.0 and u + 1 < (int)(sizeof(units) / sizeof(units[0]))) {
        v /= 1024.0;
        ++u;
    }
    char buf[64];
    if (u == 0)
        std::snprintf(buf, sizeof(buf), "%llu B", (unsigned long long)b);
    else
        std::snprintf(buf, sizeof(buf), "%.2f %s", v, units[u]);
    return buf;
}
```

is not called within a hot loop...

- these two seem very similar `format_bytes` and `format_bytes_in`: are they both necessary?

- aren't these classes:
`progress`
`bucket_write_prof`
`bucket_process_prof`
a bit intrusive? that is, I'm worried they can harm efficiency...
