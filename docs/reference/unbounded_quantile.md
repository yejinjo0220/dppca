# Unbounded private quantile via log-binning

Internal wrapper around
[`unbounded_quantile_upper()`](https://yejinjo0220.github.io/dppca/reference/unbounded_quantile_upper.md)
that handles both upper and lower quantiles on a finite target interval
`[l, u]`.

If `q < 0.5`, the lower-tail quantile is computed by sign-flipping the
data and calling the upper-tail routine.

## Usage

``` r
unbounded_quantile(
  data,
  l,
  u,
  beta,
  q,
  eps_1,
  delta_1,
  eps_2,
  delta_2,
  max_extra_bins = 1000
)
```

## Arguments

- data:

  Numeric vector.

- l:

  Finite lower truncation bound.

- u:

  Finite upper truncation bound.

- beta:

  Log-binning base used in the geometric grid. Must satisfy `beta > 1`.

- q:

  Quantile level in `(0, 1)`.

- eps_1:

  Privacy epsilon for the noisy threshold step.

- delta_1:

  Privacy delta for the noisy threshold step.

- eps_2:

  Privacy epsilon for the noisy cumulative scan.

- delta_2:

  Privacy delta for the noisy cumulative scan.

- max_extra_bins:

  Nonnegative integer giving the extra number of bins searched past the
  largest occupied bin.

## Value

Numeric scalar private quantile estimate.
