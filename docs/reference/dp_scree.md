# Compute differentially private scree estimates

User-facing function for differentially private scree estimation. It
computes the leading non-private scree values together with one selected
DP scree estimator.

Supported methods are:

- `"clipped"`: clipped mean-based scree estimator,

- `"pmw"`: PMW-style scree estimator based on private quantiles and
  winsorized means,

- `"huber"`: Huber-type robust private mean-based scree estimator.

## Usage

``` r
dp_scree(
  X,
  k,
  method = c("clipped", "pmw", "huber"),
  eps_total,
  delta_total,
  center = TRUE,
  standardize = FALSE,
  g_dppca = FALSE,
  cpp.option = FALSE,
  mono = TRUE,
  C_clip = 3,
  beta = 1.01,
  a = NULL,
  b = NULL,
  trim_const = 10,
  eta = 0.01,
  split_mode = TRUE,
  mu0 = 0,
  eta0 = 1,
  T = NULL,
  M = NULL,
  k_min_m2 = -40,
  k_max_m2 = 40,
  m2_frac = 1/4
)
```

## Arguments

- X:

  Numeric data matrix with observations in rows.

- k:

  Integer number of leading principal components.

- method:

  One of `"clipped"`, `"pmw"`, or `"huber"`.

- eps_total:

  Total privacy epsilon allocated to the full scree routine.

- delta_total:

  Total privacy delta allocated to the full scree routine.

- center:

  Logical; whether to center columns before PCA.

- standardize:

  Logical; whether to standardize columns before PCA.

- g_dppca:

  Logical; whether to privatize the PCA direction matrix.

- cpp.option:

  Logical passed to
  [`mech_tau_sph()`](https://yejinjo0220.github.io/dppca/reference/mech_tau_sph.md)
  when private directions are computed.

- mono:

  Logical; whether to apply monotone post-processing to the final DP
  scree vector.

- C_clip:

  Positive clipping threshold used by the clipped estimator.

- beta:

  Log-binning base used by the PMW quantile estimator.

- a, b:

  Finite support bounds supplied to the PMW quantile routine.

- trim_const, eta:

  PMW practical clipping parameters.

- split_mode:

  Logical; whether the PMW estimator splits the sample.

- mu0:

  Initial value used by the Huber noisy gradient descent.

- eta0:

  Fixed step size used by the Huber noisy gradient descent.

- T:

  Optional integer number of Huber gradient descent iterations.

- M:

  Optional integer number of blocks used in
  [`dp_m2()`](https://yejinjo0220.github.io/dppca/reference/dp_m2.md).

- k_min_m2:

  Integer lower bound for dyadic histogram bins used in the Huber
  scale-proxy step.

- k_max_m2:

  Integer upper bound for dyadic histogram bins used in the Huber
  scale-proxy step.

- m2_frac:

  Fraction of the Huber scree privacy budget allocated to the private
  scale-proxy step.

## Value

A list containing `method`, `scree_np`, `evr_np`, `scree`, and `evr`.
