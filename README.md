# nprime

`nprime` is a single binary for inverting the prime counting function:

```sh
make
./nprime 1000000000
```

For positive `n`, the program prints the 1-indexed nth prime. `n = 0`
returns `0`, matching the inverse prime-counting boundary before the first
prime.

The build uses only the C++ standard library. There is no vendored or linked
third-party code, and no rank-to-prime answer table is compiled into the
binary. The default native build enables link-time optimization (`-flto`) so the
single translation unit and standard-library thread path are optimized as a
whole.

The repository also includes a static WebAssembly demo for GitHub Pages:

```sh
make pages
python3 -m http.server 8000
```

Open `http://localhost:8000/` to run the in-browser benchmark from the repository
root, or serve `_site/` after `make pages` to preview the exact artifact uploaded
by GitHub Actions. The checked-in page assets are static files; `nprime.wasm` is
compiled from `web/nprime_wasm.rs`.

## Approach

The native binary uses an analytic inverse based on the Cipolla expansion
refined by a Riemann-R/Gram-series estimate, exact Lehmer-style prime counting to
determine the estimate's rank, and an exact segmented sieve to walk the remaining
gap. For larger inputs, the top-level Lehmer summation is split across standard
library worker threads inside the single isolated invocation.

The generated prime table is an odd-only sieve. The small `phi(x, a)` table is
generated per call with 16-bit rows from the six-prime wheel row, then extended
with quotient blocks, so construction avoids the inner-loop integer division
normally present in the recurrence
`phi(x, a) = phi(x, a - 1) - phi(floor(x / p_a), a - 1)`.
For larger `a`, the counter uses the exact expansion
`phi(x, a) = phi(x, b) - sum(phi(floor(x / p_i), i))` for `b <= i < a`,
which avoids walking the left side of the recurrence one level at a time. At the
top level, those split terms are shared across the same worker pool as the
Lehmer subtraction terms, with the main thread also pulling split work before it
steals any remaining top-level subtraction chunks and joins the workers.

The first few `phi` levels above the six-prime wheel use closed exact
expansions instead of recursive calls. With the default build this expands
through `a = 10`, so calls to `phi(x, 7)` through `phi(x, 10)` are evaluated directly
from the six-prime wheel without adding any shipped tables or cross-call state.
At the hardest boundary just below 10B, the counter may also build a per-call
rank bitset up to 1,000,000 so some otherwise-recursive small `pi(x)` calls are
answered directly; this is generated and discarded inside the same isolated
invocation.

All runtime state is per invocation: generated sieve tables, phi rows, search
metadata, and final segment marking plans are discarded before the next measured
benchmark iteration. There is no cross-call count memo, no answer table, and no
precomputed nth-prime data shipped in the binary.

## Benchmarks

Run:

```sh
make bench
```

The default benchmark measures isolated calls. Each iteration resets runtime
state before and after computing the nth prime, so generated sieve tables,
search metadata, and segment plans are not reused across measured iterations.
Use `/usr/bin/time ./nprime 1000000000` to include process startup. Current
isolated performance is still above the sub-millisecond goal for large n; the
benchmark mode is intentionally honest rather than hot.

For phase-level observability, run:

```sh
./nprime --profile 9999999999
```

The profile path is also isolated. It reports the analytic estimate time, exact
Lehmer count setup/core time, final correction-sieve time, and structural
counters such as top-level Lehmer terms, recursive count calls, `phi` calls, and
segment marking work. The exact-count setup profile is split into generated
sieve allocation, composite marking, `pi(x)` table fill, generated rank-bitset
construction, six-prime wheel setup, and `phi` table construction. The core
profile then splits top-level Lehmer preparation, the overlapped base `phi`
computation, the parallel summation, worker launch/join timing, worker
min/average/max latency, worker time split between top-level Lehmer subtraction
and split-base `phi` work, main-thread stolen top-term work, and correction
planning/marking/scanning time. Depth and size buckets show the shape of
recursive prime-count calls, and rank-tier counters show how many `pi(x)`
lookups hit the generated per-call rank bitset. Use `--profile-detail` when you
also want `phi` buckets showing whether the workload is dominated by cached
small-a calls or recursive large-a calls; those extra counters are diagnostic
and intentionally kept out of the default timing profile. The current hot path
is the exact prime-counting phase; the final correction scan is much smaller.

The stage counters are intended to separate algorithmic causes from local
implementation costs. `count_value_recurse_by_stage` shows which part of the
Lehmer formula forced recursive `pi(x)` calls; `count_value_direct_by_stage` and
`count_value_rank_direct_by_stage` show where generated direct tables answered
the query; `count_value_outer_recurse_sum_ms_by_stage` gives a summed worker-time
attribution for the outer recursive `pi(x)` calls so repeated call counts do not
hide which Lehmer stage is actually consuming time. `phi_calls_by_stage`
separates base `phi(x,a)` work from split `phi(floor(x/p_i),i)` work. In current
near-10B profiles, recursive count calls come from `top_term`, which means the main
unresolved cost is the top-level Lehmer subtraction
`sum(pi(floor(x / p_i)))`, not the final prime-gap sieve.

For controlled configuration sweeps, use:

```sh
./nprime --bench-config 9999999999 650000 65535 32 10 1 8 5
./nprime --profile-config 9999999999 650000 65535 32 10 1 8
./nprime --profile-detail-config 9999999999 650000 65535 32 10 1 8
./nprime --bench-config-rank 9999999999 502100 1000000 131072 32 12 12 8 5
```

The configuration tuple is `sieve phi_x phi_a threads a_boost split_slack`.
The `*-config-rank` variants insert an explicit `rank` limit after `sieve`.
These modes are isolated in the same way as `--bench`: each measured invocation
rebuilds and discards generated runtime state.

For a repeatable profile sweep:

```sh
make profile-scale
```

For a wider hot-path sweep across the current optimization boundary:

```sh
make profile-hotpaths
```

Recent local isolated timings on this machine are roughly:

```text
n=1,000,000:        0.06 ms
n=1,000,000,000:    0.85-1.2 ms
n=5,000,000,000:    1.7-2.0 ms
n=9,999,999,999:    2.8-3.4 ms
n=10,000,000,000:   3.2-4.0 ms
n=20,000,000,000:   9-10 ms
```

Run a wider scaling check with:

```sh
make bench-scale
```

This keeps the 1B benchmark, then measures larger exact inversion points,
currently just below 10B, 10B, and a one-shot 20B scaling probe.

## WebAssembly build

The browser module is compiled from `web/nprime_wasm.rs` with the standard Rust
`wasm32-unknown-unknown` target and no crates:

```sh
make wasm
```

The WebAssembly implementation uses the same broad method as the native binary:
analytic initialization, exact Lehmer counting, rank correction, and exact
segmented sieving. It is single-threaded and intentionally separate from the
native C++ implementation. The browser page uses `nth_prime_isolated`, which
resets module state before and after each measured call so repeated benchmark
iterations do not reuse cross-call caches.

In this project, a "hot average" means repeated calls inside one already-warmed
process or WebAssembly instance. Hot averages are useful for studying the inner
search path, but they share runtime-generated working structures and are not the
default benchmark mode.

## GitHub Pages deployment

The repository includes `.github/workflows/pages.yml`, a custom GitHub Pages
workflow that builds the site from source on every push to `main` and on manual
dispatch. The workflow:

1. checks out the repository,
2. installs the standard `wasm32-unknown-unknown` Rust target,
3. runs `make pages`,
4. uploads `_site/` as a Pages artifact, and
5. deploys it with the official GitHub Pages deploy action.

Before the workflow will publish, set the repository's Pages source to
**GitHub Actions** in the GitHub repository settings.
