# Compute differentially private scree estimates

User-facing function for differentially private scree estimation. It
computes the leading non-private scree values together with one selected
DP scree estimator.

Supported methods are:

- `"clipped"`: clipped mean-based scree estimator,

- `"pmwm"`: PMWM-style scree estimator based on private quantiles and
  winsorized means,

- `"huber"`: Huber-type robust private mean-based scree estimator.

Method-specific tuning parameters are supplied through `control`:

- `control = clipped_control(C_clip = ...)` for `method = "clipped"`,

- `control = pmwm_control(beta = ..., a = ..., b = ..., trim_const = ..., eta = ..., split_mode = ...)`
  for `method = "pmwm"`,

- `control = huber_control(mu0 = ..., eta0 = ..., T = ..., M = ..., k_min_m2 = ..., k_max_m2 = ..., m2_frac = ...)`
  for `method = "huber"`.

## Usage

``` r
dp_scree(
  X,
  k,
  method = c("clipped", "pmwm", "huber"),
  eps_total,
  delta_total,
  center = TRUE,
  standardize = FALSE,
  g_dppca = FALSE,
  cpp.option = FALSE,
  mono = TRUE,
  control = NULL
)
```

## Arguments

- X:

  Numeric data matrix with observations in rows.

- k:

  Integer number of leading principal components.

- method:

  One of `"clipped"`, `"pmwm"`, or `"huber"`.

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

- control:

  Method-specific control list created by
  [`clipped_control()`](https://yejinjo0220.github.io/dppca/reference/clipped_control.md),
  [`pmwm_control()`](https://yejinjo0220.github.io/dppca/reference/pmwm_control.md),
  or
  [`huber_control()`](https://yejinjo0220.github.io/dppca/reference/huber_control.md).

## Value

A list containing `method`, `scree_np`, `evr_np`, `scree`, and `evr`.

## Examples

``` r
if (FALSE) { # \dontrun{
dp_scree(X, k = 3, method = "clipped", eps_total = 1, delta_total = 1e-6,
         control = clipped_control(C_clip = 3))

dp_scree(X, k = 3, method = "pmwm", eps_total = 1, delta_total = 1e-6,
         control = pmwm_control(a = 0, b = 10))

dp_scree(X, k = 3, method = "huber", eps_total = 1, delta_total = 1e-6,
         control = huber_control(T = 50, M = 20))
} # }
```
