#pragma once

// ram_governor: one always-on, measured-RSS controller shared by ALL phases.
//
// WHY (see ram-governor.md): today each phase predicts its structures' sizes,
// sums them, and trusts the model to fit -g. That is fragile -- every new
// dataset finds a size the model under-counted, a co-peak it didn't expect, or
// glibc holding freed pages, and -g is silently exceeded. The governor moves
// CORRECTNESS off the model: it reads the process's REAL resident set (which
// already includes every structure + glibc + fragmentation) and, when RSS
// crosses the budget, asks every registered spillable structure to spill to
// disk and returns the freed pages to the OS, until RSS drops. The per-phase
// size models then only decide HOW MUCH each structure keeps in RAM for speed --
// they no longer decide whether -g holds.
//
// PROPERTIES it is built to have:
//   1. Authoritative -- driven by measured RSS, not predicted sizes.
//   2. Comprehensive -- EVERY RAM-first buffer registers (bucket-write compactor
//      maps, the stitch stores, the unitig writer, emit sort buffers, the
//      dedup-map overflow). Nothing is invisible to it.
//   3. Always-on -- one controller runs for the whole build; no unwatched
//      windows between/within phases.
//   4. Reclaims glibc -- each spill is followed by malloc_trim, so a spill
//      actually LOWERS RSS (without this, freeing into glibc's arena leaves RSS
//      unchanged -- exactly why the per-phase stitch watcher couldn't hold -g).
//
// WHAT IT DOES NOT DO: it cannot spill a NON-spillable structure. So -g holds
// iff (the governor drives all spillable RAM down) AND (the non-spillable floor
// < -g). The non-spillable floor is the small, enumerable list in algorithm.md
// §9.1 (per-bucket walk set, color-sets-dedup-map, stitch transients), each
// handled separately by scaling or abort. The governor turns an open-ended
// surprise surface into that closed list.
//
// This is the generic foundation; phases are retrofitted to register in stages
// (ram-governor.md). Not yet wired into the pipeline -- defining the seam first.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <mutex>
#include <thread>
#include <vector>

#if defined(__GLIBC__)
#include <malloc.h>
#endif

#include "util.hpp"  // current_rss_bytes, process_peak_rss_bytes

namespace cdbg {

// A participant the governor can ask to free RAM under memory pressure.
// Implemented by every RAM-first structure (the stores, the unitig writer, the
// dedup-map overflow, ...). All methods MUST be thread-safe: they are called
// from the governor's poller thread while the participant's own worker threads
// run.
struct ram_spillable {
    virtual ~ram_spillable() = default;

    // Free about `target_bytes` of RAM by spilling to disk. Returns the bytes
    // actually freed (best-effort; 0 if nothing is left to spill). The governor
    // sums these across participants and re-measures RSS, so a participant may be
    // asked repeatedly as long as RSS stays high.
    virtual uint64_t spill_under_pressure(uint64_t target_bytes) = 0;

    // Optional hint: how much spillable RAM the participant currently holds, so
    // the governor can target the biggest first. 0 = unknown (governed in
    // registration order).
    virtual uint64_t spillable_bytes() const { return 0; }
};

// The controller. Construct once per build with the -g budget, start() it,
// register spillables (RAII), and it holds real RSS under budget for the whole
// run. budget == 0 disables it entirely (no thread, no-op) -- the "no -g" path.
struct ram_governor {
    // Trips (asks participants to spill) when RSS >= high_frac*budget; clears the
    // pull-style pressure signal when RSS <= low_frac*budget. Poll `interval`
    // matches the per-phase watchers (25 ms) so the overshoot above `high` --
    // allocation_rate * (poll + spill lag) -- stays small.
    explicit ram_governor(uint64_t budget, double high_frac = 0.85, double low_frac = 0.70,
                          std::chrono::milliseconds interval = std::chrono::milliseconds(25))
        : m_budget(budget)
        , m_high((uint64_t)(high_frac * (double)budget))
        , m_low((uint64_t)(low_frac * (double)budget))
        , m_interval(interval) {}

    ~ram_governor() { stop(); }
    ram_governor(ram_governor const&) = delete;
    ram_governor& operator=(ram_governor const&) = delete;

    // RAII registration: the participant is governed for the handle's lifetime
    // (register on phase entry, the handle's destructor unregisters on exit).
    struct registration {
        registration() = default;
        registration(ram_governor* g, ram_spillable* p) : m_g(g), m_p(p) {}
        registration(registration&& o) noexcept : m_g(o.m_g), m_p(o.m_p) {
            o.m_g = nullptr;
            o.m_p = nullptr;
        }
        registration& operator=(registration&& o) noexcept {
            if (this != &o) {
                reset();
                m_g = o.m_g;
                m_p = o.m_p;
                o.m_g = nullptr;
                o.m_p = nullptr;
            }
            return *this;
        }
        registration(registration const&) = delete;
        registration& operator=(registration const&) = delete;
        ~registration() { reset(); }
        void reset() {
            if (m_g) m_g->unregister_(m_p);
            m_g = nullptr;
            m_p = nullptr;
        }

    private:
        ram_governor* m_g = nullptr;
        ram_spillable* m_p = nullptr;
    };

    [[nodiscard]] registration add(ram_spillable* p) {
        if (m_budget == 0 or p == nullptr) return registration{};
        std::lock_guard<std::mutex> lk(m_mu);
        m_parts.push_back(p);
        return registration{this, p};
    }

    void start() {
        if (m_budget == 0) return;
        bool expected = false;
        if (!m_run.compare_exchange_strong(expected, true)) return;  // already running
        m_thread = std::thread([this] { loop_(); });
    }

    void stop() {
        if (!m_run.exchange(false)) return;
        if (m_thread.joinable()) m_thread.join();
    }

    // Pull-style signal for participants that poll a flag instead of registering
    // a spill callback (e.g. the existing stores during their transition): true
    // between the high and low crossings.
    bool under_pressure() const { return m_pressure.load(std::memory_order_relaxed); }

    uint64_t budget() const { return m_budget; }
    // Diagnostics (reported at phase boundaries to attribute spill activity).
    bool was_engaged() const { return m_events.load(std::memory_order_relaxed) > 0; }
    uint64_t pressure_events() const { return m_events.load(std::memory_order_relaxed); }
    uint64_t observed_rss_high() const { return m_rss_high.load(std::memory_order_relaxed); }

private:
    void unregister_(ram_spillable* p) {
        std::lock_guard<std::mutex> lk(m_mu);
        for (auto& e : m_parts) {
            if (e == p) {
                e = m_parts.back();
                m_parts.pop_back();
                break;
            }
        }
    }

    void loop_() {
        while (m_run.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(m_interval);
            if (!m_run.load(std::memory_order_relaxed)) break;
            uint64_t rss = current_rss_bytes();
            if (rss == 0) rss = process_peak_rss_bytes();  // no /proc: monotonic fallback
            if (rss == 0) continue;
            uint64_t prev = m_rss_high.load(std::memory_order_relaxed);
            while (rss > prev and
                   !m_rss_high.compare_exchange_weak(prev, rss, std::memory_order_relaxed)) {
            }
            if (rss >= m_high) {
                m_pressure.store(true, std::memory_order_relaxed);
                m_events.fetch_add(1, std::memory_order_relaxed);
                apply_pressure_(rss);
            } else if (rss <= m_low) {
                m_pressure.store(false, std::memory_order_relaxed);
            }
        }
    }

    void apply_pressure_(uint64_t rss) {
        const uint64_t want = rss > m_low ? rss - m_low : 0;
        if (want > 0) {
            std::vector<ram_spillable*> parts;
            {
                std::lock_guard<std::mutex> lk(m_mu);
                parts = m_parts;
            }
            // Target the biggest spillable holders first (best-effort hint).
            std::sort(parts.begin(), parts.end(), [](ram_spillable* a, ram_spillable* b) {
                return a->spillable_bytes() > b->spillable_bytes();
            });
            uint64_t freed = 0;
            for (auto* p : parts) {
                if (freed >= want) break;
                freed += p->spill_under_pressure(want - freed);
            }
        }
        // Return freed pages to the OS so RSS actually drops -- this is what
        // reclaims glibc-retained arena memory, INDEPENDENT of any participant
        // spill (freed pages strand in per-thread arenas otherwise). Done ONLY
        // under pressure and throttled, so it never degrades into the
        // per-interval full-heap trim that backfired in the per-phase watchers.
        maybe_trim_();
    }

    void maybe_trim_() {
#if defined(__GLIBC__)
        auto now = std::chrono::steady_clock::now();
        if (now - m_last_trim < std::chrono::milliseconds(250)) return;
        m_last_trim = now;
        ::malloc_trim(0);
#endif
    }

    uint64_t m_budget;
    uint64_t m_high;
    uint64_t m_low;
    std::chrono::milliseconds m_interval;

    std::atomic<bool> m_run{false};
    std::atomic<bool> m_pressure{false};
    std::atomic<uint64_t> m_events{0};
    std::atomic<uint64_t> m_rss_high{0};
    std::chrono::steady_clock::time_point m_last_trim{};  // throttle the glibc trim

    std::thread m_thread;
    std::mutex m_mu;
    std::vector<ram_spillable*> m_parts;
};

}  // namespace cdbg
