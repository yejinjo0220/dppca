# Private scalar scale proxy from paired differences

Internal helper used by the Huber scree estimator to obtain a private
scale proxy for a one-dimensional sample. The input vector `w` is first
paired, transformed into squared paired differences, aggregated by block
medians, and then fed into a noisy dyadic histogram.

## Usage

``` r
dp_m2(w, eps_m2, M = NULL, k_min_m2 = -20, k_max_m2 = 40)
```

## Arguments

- w:

  Numeric vector, typically the centered squared projected scores for a
  single principal component.

- eps_m2:

  Positive privacy epsilon used for the scale-proxy step.

- M:

  Optional integer number of blocks used in the block-median step. If
  `NULL`, a default based on `sqrt(n / 2)` is used.

- k_min_m2:

  Integer lower bound for the dyadic bin index.

- k_max_m2:

  Integer upper bound for the dyadic bin index.

## Value

A nonnegative numeric scalar private scale proxy.
