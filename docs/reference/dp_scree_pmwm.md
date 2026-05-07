# Differentially private scree estimation via PMWM

Internal helper that estimates the leading scree values using a
PMWM-style procedure based on:

1.  PCA preprocessing and direction estimation,

2.  centered squared scores for each component,

3.  private lower and upper quantile estimation,

4.  winsorization of the squared scores,

5.  Gaussian release of the winsorized mean,

6.  the \\n/(n-1)\\ scree correction.

If `split_mode = TRUE`, one subset of the sample is used for private
quantile estimation and the other subset is used for the winsorized mean
step. If `split_mode = FALSE`, the full sample is reused in both steps.

The clipping proportion is the practical choice
`max(trim_const / n_q, eta)`, truncated above at `0.49`.

## Usage

``` r
dp_scree_pmwm(
  X,
  k,
  eps_total,
  delta_total,
  g_dppca = FALSE,
  cpp.option = FALSE,
  split_mode = TRUE,
  center = TRUE,
  standardize = TRUE,
  beta = 1.01,
  a,
  b,
  trim_const = 10,
  eta = 0.01,
  mono = TRUE,
  max_extra_bins = 1000
)
```

## Arguments

- X:

  Numeric data matrix with observations in rows.

- k:

  Integer number of leading principal components.

- eps_total:

  Total privacy epsilon allocated to the scree routine.

- delta_total:

  Total privacy delta allocated to the scree routine.

- g_dppca:

  Logical; whether to privatize the PCA direction matrix.

- cpp.option:

  Logical passed to
  [`mech_tau_sph()`](https://yejinjo0220.github.io/dppca/reference/mech_tau_sph.md)
  when private directions are computed.

- split_mode:

  Logical; if `TRUE`, split the sample into quantile and mean subsets.

- center:

  Logical; whether to center columns before PCA.

- standardize:

  Logical; whether to standardize columns before PCA.

- beta:

  Log-binning base used in the private quantile estimator.

- a:

  Finite lower support bound supplied to the private quantile routine.

- b:

  Finite upper support bound supplied to the private quantile routine.

- trim_const:

  Positive constant controlling the practical clipping proportion
  `max(trim_const / n_q, eta)`.

- eta:

  Lower bound in the practical clipping proportion.

- mono:

  Logical; whether to post-process the final scree vector so that it is
  nonnegative and nonincreasing.

- max_extra_bins:

  Nonnegative integer controlling how far the unbounded quantile scan
  extends past the largest occupied bin.

## Value

A list with components:

- `scree`:

  Private scree estimates.

- `scree_np`:

  Non-private winsorized scree estimates.

- `pve`:

  Proportions of variance explained computed from the private scree
  estimates.
