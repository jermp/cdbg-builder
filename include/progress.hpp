#pragma once

#include <atomic>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <string>
#include <thread>
#include <unistd.h>

namespace cdgb {

// Lightweight progress printer for long-running phases.
//
// Spawns a background thread that periodically samples a worker-owned atomic
// counter and prints a status line to stderr. When stderr is a TTY it updates
// a single line in place via '\r'; otherwise it emits one line per tick so
// progress is still visible in piped logs.
//
// Workers don't talk to Progress directly — they bump `counter`, and Progress
// only reads it. This keeps the hot path lock-free and Progress optional.
class Progress {
public:
    Progress(std::string label, std::atomic<uint64_t>& counter, uint64_t total,
             std::ostream& os = std::cerr,
             std::chrono::milliseconds interval = std::chrono::milliseconds(500))
        : m_label(std::move(label))
        , m_counter(counter)
        , m_total(total)
        , m_os(os)
        , m_interval(interval)
        , m_is_tty(::isatty(2) != 0)
        , m_start(std::chrono::steady_clock::now()) {
        m_thread = std::thread([this] { run(); });
    }

    ~Progress() { stop(); }

    Progress(const Progress&) = delete;
    Progress& operator=(const Progress&) = delete;

    void stop() {
        bool was = m_stopped.exchange(true);
        if (was) return;
        if (m_thread.joinable()) m_thread.join();
        print_line(true);
    }

private:
    void run() {
        // Print one line up front so the user sees the phase started.
        print_line(false);
        while (!m_stopped.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(m_interval);
            if (m_stopped.load(std::memory_order_relaxed)) break;
            print_line(false);
        }
    }

    void print_line(bool final) {
        uint64_t done = m_counter.load(std::memory_order_relaxed);
        if (final) done = m_total ? m_total : done;

        auto now = std::chrono::steady_clock::now();
        double sec =
            std::chrono::duration_cast<std::chrono::milliseconds>(now - m_start).count() / 1000.0;

        char buf[256];
        if (m_total) {
            double pct = 100.0 * (double)done / (double)m_total;
            std::snprintf(buf, sizeof(buf), "[%s] %llu/%llu (%.1f%%) %.1fs", m_label.c_str(),
                          (unsigned long long)done, (unsigned long long)m_total, pct, sec);
        } else {
            std::snprintf(buf, sizeof(buf), "[%s] %llu %.1fs", m_label.c_str(),
                          (unsigned long long)done, sec);
        }

        if (m_is_tty) {
            m_os << '\r' << buf << "      ";
            if (final) m_os << '\n';
        } else {
            // Skip mid-flight prints if nothing changed since last line.
            if (!final && done == m_last_done) return;
            m_last_done = done;
            m_os << buf << '\n';
        }
        m_os.flush();
    }

    std::string m_label;
    std::atomic<uint64_t>& m_counter;
    uint64_t m_total;
    std::ostream& m_os;
    std::chrono::milliseconds m_interval;
    bool m_is_tty;
    std::chrono::steady_clock::time_point m_start;
    std::atomic<bool> m_stopped{false};
    std::thread m_thread;
    uint64_t m_last_done = (uint64_t)-1;
};

}  // namespace cdgb
