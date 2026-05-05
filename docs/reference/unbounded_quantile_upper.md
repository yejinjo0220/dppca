# Unbounded private upper quantile via log-binning

Internal helper that computes a private upper-tail quantile estimate
using a logarithmic grid together with noisy thresholding and noisy
cumulative counts.

This function is used in the PMW scree routine to estimate an upper
clipping bound privately. The search is carried out over a geometric
grid determined by the lower anchor `l` and the log-binning base `beta`.

## Usage

``` r
unbounded_quantile_upper(
  data,
  l,
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

  Finite lower anchor for the search grid.

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

Numeric scalar private upper-quantile estimate.
