# The RAM governor (design)

Status: **design + foundation struct** (branch `claude/ram-governor`, off
`claude/externalize-colorset-dedup`). The struct lives in `include/ram_governor.hpp`;
the phases are retrofitted to use it in stages (below). Nothing is wired into the
pipeline yet — this lands the seam first.

## The problem it ends

Every phase today *predicts* its structures' sizes, sums them, and trusts the
model to fit `-g`. That is structurally fragile: each new dataset finds a size
the model under-counted, a **co-peak** it didn't expect (e.g. the unitig writer
+ stores + transient all full at the stitch phase boundary — the 100k/`-g16`
overshoot), or **glibc** holding freed pages. We keep patching instances; it
doesn't converge, because the *model* is load-bearing for *correctness*.

## The fix: make the model non-load-bearing

One always-on controller (`ram_governor`) reads the process's **real RSS** —
which already includes every structure, glibc, and fragmentation — and, when RSS
crosses the budget, asks **every** registered spillable structure to spill to
disk and returns the freed pages to the OS, until RSS drops. The per-phase size
models then only decide *how much each structure keeps in RAM for speed*; they no
longer decide whether `-g` holds. Four properties:

1. **Authoritative** — driven by a measurement, not a prediction. A model
   under-count can't break `-g`; it just makes the governor act.
2. **Comprehensive** — every RAM-first buffer registers (bucket-write compactor
   maps, stitch stores, **the unitig writer**, emit sort buffers, the dedup-map
   overflow). The gap that bit us in stitch — the watcher couldn't see the unitig
   writer — cannot exist: if it holds RAM, it's registered.
3. **Always-on** — one controller for the whole build; no unwatched windows
   (emit has no watcher today; the inter-phase boundaries are uncovered).
4. **Reclaims glibc** — each spill is followed by `malloc_trim`, so a spill
   actually *lowers* RSS. Without this, `free()` returns memory to glibc's arena
   and RSS doesn't move — exactly why the per-phase stitch watcher couldn't hold
   16 GiB.

## Why this is "definitive"

It converts an **open-ended** surprise surface into a **closed, enumerable** one.
After it's in place, `-g` can be violated in only one way: a **non-spillable**
structure that alone exceeds `-g`. Those are the short, code-inspectable list in
`algorithm.md` §9.1 (per-bucket walk set, color-sets-dedup-map, stitch
transients), each already handled by *scaling* (frag_ranges, chain_buckets, the
planned bucket re-split, the dedup externalization) or *abort*. The failure modes
go from "unknown, dataset-dependent, keeps appearing" to "this known list, each
with a bound or an abort." The goal isn't "the model is finally right" — it's
"the model no longer has to be right."

## The interface (`include/ram_governor.hpp`)

```cpp
struct ram_spillable {                          // every RAM-first structure
    virtual uint64_t spill_under_pressure(uint64_t target_bytes) = 0;  // free ~target, return freed
    virtual uint64_t spillable_bytes() const { return 0; }            // size hint (target biggest first)
};

struct ram_governor {
    ram_governor(uint64_t budget, double high=0.85, double low=0.70, ms interval=25);
    registration add(ram_spillable*);   // RAII: governed for the handle's lifetime
    void start(); void stop();
    bool under_pressure() const;        // pull-style flag for poll-based participants
    // diagnostics: was_engaged(), pressure_events(), observed_rss_high()
};
```

Poller loop: every 25 ms read RSS; on `>= high*budget`, ask participants
(biggest `spillable_bytes()` first) to free down to the low watermark, then
`malloc_trim`; clear the pull-flag on `<= low*budget`. `budget == 0` (no `-g`)
disables it entirely.

`spill_under_pressure` and `spillable_bytes` MUST be thread-safe — they run on
the governor thread while the participant's workers run. (The existing stores
already spill themselves off a pressure flag; they can either keep that via
`under_pressure()` or implement the callback.)

## What it does NOT do

It can't spill a non-spillable structure, so it doesn't remove the need to
scale/abort those — it makes them the *only* remaining risk, and a closed one.
(This is why the dedup-map externalization still matters: it converts that big
non-spillable into a bounded/spillable one, shrinking the residual list.)

## Staged retrofit (each stage builds + passes `test_stitch` and 86630/171)

1. **Foundation** (this commit): the `ram_governor` + `ram_spillable` struct in
   `include/`, compiled but not wired.
2. **Own it in `builder`**: construct one `ram_governor(g)` at build start,
   `start()`/`stop()` around the run; replace the ad-hoc `cap_malloc_arenas_` /
   per-phase `release_free_heap_to_os_` calls with the governor's reclaim.
3. **Retrofit stitch**: the stores + the unitig writer implement `ram_spillable`
   (or subscribe to `under_pressure()`), register for the stitch's duration. This
   is the one that directly fixes the 100k overshoot at the *backstop* level
   (the unitig writer becomes reclaimable), letting `STITCH_BUDGET_FRAC` rise
   back toward speed.
4. **Retrofit bucket-write** (compactor maps) and **bucket-process** (admission
   gate reads the governor's pressure) onto the same controller, retiring their
   private watchers.
5. **Retrofit emit** (cid-sort buffers) and the **dedup-map overflow** once it
   exists.
6. **Report** governor activity per phase (events, observed RSS high) so spill
   pressure is visible, and confirm on 100k/661k that real RSS tracks `-g`.

## Open questions / knobs

- `high`/`low` watermarks (0.85/0.70 default) vs. the per-phase fractions — the
  governor is the *backstop*; the fractions keep it from firing in the common
  case. Once the governor is comprehensive, the fractions can be more generous
  (use more RAM) because the backstop is real.
- Spill ordering / fairness across participants under sustained pressure.
- Whether a participant that has spilled to empty while RSS is still over budget
  should escalate (it means the non-spillable floor exceeds `-g` — today the
  model-level checks abort; the governor could surface a precise diagnostic).
