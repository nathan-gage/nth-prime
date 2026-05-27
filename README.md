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
binary.

## Approach

The binary uses an analytic estimate for the nth prime, exact Lehmer-style
prime counting to determine the estimate's rank, and an exact segmented sieve to
walk the remaining gap. Runtime memoization is used inside a process, including
the exact prime-count subresults and the final segment's sieve marking plan; no
precomputed nth-prime answers are shipped.

## Benchmarks

Run:

```sh
make bench
```

The default benchmark computes the 1,000,000,000th prime and a nearby value
repeatedly in-process.
Use `/usr/bin/time ./nprime 1000000000` to include process startup.

Run a wider scaling check with:

```sh
make bench-scale
```

This keeps the 1B benchmark, then measures larger exact inversion points,
currently just below 10B, 10B, and a one-shot 20B scaling probe.
