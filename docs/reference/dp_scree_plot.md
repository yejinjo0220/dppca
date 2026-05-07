# Plot differentially private scree curves

User-facing function that computes and plots non-private and/or
differentially private scree curves using one or more DP scree
estimators.

By default, proportions of variance explained (PVE) are plotted. When
`type = "scree"`, raw scree values are plotted instead.

Method-specific tuning parameters are supplied through `control`. For a
single method, use a single control object, for example
`control = pmwm_control(a = ..., b = ...)`. For
`dp_scree_method = "all"`, use a named list, for example
`control = list(clipped = clipped_control(...), pmwm = pmwm_control(...), huber = huber_control(...))`.

## Usage

``` r
dp_scree_plot(
  X,
  k,
  dp_scree_method = c("clipped", "pmwm", "huber", "all"),
  eps_total,
  delta_total,
  center = TRUE,
  standardize = FALSE,
  control = NULL,
  g_dppca = FALSE,
  cpp.option = FALSE,
  mono = TRUE,
  type = c("pve", "scree")
)
```

## Arguments

- X:

  Numeric data matrix with observations in rows.

- k:

  Integer number of leading principal components.

- dp_scree_method:

  Which DP estimator(s) to plot. One of `"all"`, `"clipped"`, `"pmwm"`,
  or `"huber"`.

- eps_total:

  Total privacy epsilon allocated to the full scree routine.

- delta_total:

  Total privacy delta allocated to the full scree routine.

- center:

  Logical; whether to center columns before PCA.

- standardize:

  Logical; whether to standardize columns before PCA.

- control:

  Method-specific control list, or a named list of control lists when
  `dp_scree_method = "all"`.

- g_dppca:

  Logical; whether to privatize the PCA direction matrix.

- cpp.option:

  Logical passed to
  [`mech_tau_sph()`](https://yejinjo0220.github.io/dppca/reference/mech_tau_sph.md)
  when private directions are computed.

- mono:

  Logical; whether to apply monotone post-processing.

- type:

  Either `"pve"` or `"scree"`.

## Value

Invisibly returns a list containing non-private results and the
method-specific
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
outputs that were plotted.

## Details

Plot appearance is handled internally. The non-private curve is
overlaid, points are shown on each curve, and colors, line types, point
shapes, labels, limits, and legend position are set automatically.

## Examples

``` r
if (FALSE) { # \dontrun{
dp_scree_plot(
  X, k = 3, dp_scree_method = "all",
  eps_total = 1, delta_total = 1e-6,
  control = list(
    clipped = clipped_control(C_clip = 3),
    pmwm = pmwm_control(a = 0, b = 10),
    huber = huber_control(T = 50, M = 20)
  )
)
} # }
```
