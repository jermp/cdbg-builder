# The RAM governor (design + status)

Status: **MERGED, scoped to stitch.** `include/ram_governor.hpp` is in `main`,
constructed every build and `start()`ed before the stitch phase, with the
**unitig writer** as its one registered participant. The sections below
describe the full *vision* (one always-on controller for every phase); the
"Current status" and "Lessons" sections record how much of it is real today and
why the rest was deliberately *not* turned on. Read those two first if you only
want the truth on the ground.

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

## Current status (what's real today, vs the vision above)

The four properties above are the **goal**, not all delivered. As merged:

- **Scoped to stitch, NOT always-on.** The governor is constructed at build
  start but `start()`ed only before stitch. Bucket-write and bucket-process keep
  their own RAM controls (the bucket-write RSS watcher, the bucket-process
  admission gate) and run with the governor *off* — it never fires there
  anyway (their RSS stays below its high-watermark), and an always-on poller
  there is pure overhead. So "always-on / no unwatched windows" is **not** the
  current behavior; it's the eventual target.
- **One participant, not comprehensive.** Only the **unitig writer** implements
  `ram_spillable` and registers (during stitch). The stitch stores still spill
  via their own `ram_budget`; bucket-write/process/emit buffers are not yet
  registered. So "every RAM-first buffer registers" is also aspirational.
- **Watermarks** in use: `0.80·g` high / `0.62·g` low (a touch below `-g` for
  spinning-disk lag), not the `0.85/0.70` defaults in the struct.
- **Proven** at 100k/`-g 16`: stitch peaks at 13.11 GiB (held under 16), the
  governor force-spilling the unitig writer near the cap. Not yet validated at
  661k.

The remaining "always-on + retire the per-phase watchers" unification is a
deliberate **future** step (see Lessons), to be done as a *separately measured*
change rather than bundled in.

## Lessons from the integration (why it's only scoped to stitch)

Getting here cost a real detour; the design above survived but several
"obvious" pieces did **not**, and the reasons are load-bearing:

- **`malloc_trim` is the only way a spill lowers RSS** (freed glibc-arena pages
  stay resident otherwise) — but it is **expensive** (takes the malloc lock,
  walks the heap). It is safe only when fired under *genuine* pressure. Adding
  the same trim to the **bucket-write** watcher — which sits over its HIGH
  threshold for most of the phase — produced a ~40% slowdown (`131 trims × ~3 s`)
  and was reverted. So reclaim must be confined to where RSS truly approaches
  `-g`; that is *why* the governor is scoped to stitch, not always-on.
- **Global allocator knobs are dangerous.** An `M_ARENA_MAX` cap (one arena per
  thread) was tried to shrink glibc retention; it serialized allocation and
  ~halved bucket-write throughput. Dropped. The governor's *reactive* reclaim is
  the right tool; static global caps are not.
- **The governor's *polling* is free; its *action* is not.** A 25 ms
  `/proc/self/status` read costs nothing and never slowed any phase — every
  bucket-write/process slowdown traced to the trim or the arena cap, not to the
  governor watching. This is what makes "always-on" *eventually* safe: cheap
  polling everywhere, reclaim only near `-g`.
- **Measure, don't theorize.** The spill-count blowups were chased with several
  wrong stories (machine load, arena cap) before the data settled it. The fix
  process now leans on counters (`spills=`, sweep counts), not narratives.

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

1. ✅ **Foundation** (merged): the `ram_governor` + `ram_spillable` struct in
   `include/`.
2. ~ **Own it in `builder`** (merged, but **scoped to stitch**, not the whole
   build — see Current status). The `cap_malloc_arenas_` idea was tried here and
   **dropped** (Lessons); `release_free_heap_to_os_` at phase boundaries stays.
3. ✅ **Retrofit stitch** (merged): the **unitig writer** implements
   `ram_spillable` and registers for stitch's duration — the one piece that fixes
   the 100k overshoot at the backstop level (peak 13.11 GiB ≤ 16). The stitch
   *stores* still use their own `ram_budget` (not yet folded in).
4. ☐ **Retrofit bucket-write** (compactor maps) and **bucket-process** (admission
   gate reads the governor's pressure) onto the same controller, retiring their
   private watchers.
5. ☐ **Retrofit emit** (cid-sort buffers) and the **dedup-map overflow** once it
   exists.
6. ☐ **Report** governor activity per phase (events, observed RSS high) so spill
   pressure is visible, and confirm on 100k/661k that real RSS tracks `-g`.

Stages 4–6 are the "always-on + retire the per-phase watchers" unification, held
as a separate measured change (Lessons). Stages 1–3 (scoped to stitch) are what
honors `-g 16` at 100k today.

## Open questions / knobs

- `high`/`low` watermarks (0.85/0.70 default) vs. the per-phase fractions — the
  governor is the *backstop*; the fractions keep it from firing in the common
  case. Once the governor is comprehensive, the fractions can be more generous
  (use more RAM) because the backstop is real.
- Spill ordering / fairness across participants under sustained pressure.
- Whether a participant that has spilled to empty while RSS is still over budget
  should escalate (it means the non-spillable floor exceeds `-g` — today the
  model-level checks abort; the governor could surface a precise diagnostic).
