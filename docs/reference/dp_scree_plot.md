# Plot differentially private scree curves

User-facing function that computes and plots non-private and/or
differentially private scree curves using one or more DP scree
estimators.

By default, explained variance ratios (EVR) are plotted. When
`type = "scree"`, raw scree values are plotted instead.

## Usage

``` r
dp_scree_plot(
  X,
  k,
  dp_scree_method = c("all", "clipped", "pmw", "huber"),
  eps_total,
  delta_total,
  center = TRUE,
  standardize = TRUE,
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
  m2_frac = 1/4,
  g_dppca = FALSE,
  cpp.option = FALSE,
  mono = TRUE,
  type = c("evr", "scree"),
  show_nonprivate = TRUE,
  show_points = TRUE,
  col_map = NULL,
  lty_map = NULL,
  pch_map = NULL,
  main = NULL,
  xlab = "Component",
  ylab = NULL,
  ylim = NULL,
  legend_pos = "topright",
  ...
)
```

## Arguments

- X:

  Numeric data matrix with observations in rows.

- k:

  Integer number of leading principal components.

- dp_scree_method:

  Which DP estimator(s) to plot. One of `"all"`, `"clipped"`, `"pmw"`,
  or `"huber"`.

- eps_total:

  Total privacy epsilon allocated to the full scree routine.

- delta_total:

  Total privacy delta allocated to the full scree routine.

- center:

  Logical; whether to center columns before PCA.

- standardize:

  Logical; whether to standardize columns before PCA.

- C_clip:

  Positive clipping threshold used by the clipped estimator.

- beta:

  Log-binning base used by the PMW quantile estimator.

- a, b:

  Finite support bounds supplied to the PMW quantile routine.

- trim_const, eta:

  PMW practical clipping parameters.

- split_mode:

  Logical; whether the PMW estimator splits the sample into quantile and
  mean subsets.

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

- g_dppca:

  Logical; whether to privatize the PCA direction matrix.

- cpp.option:

  Logical passed to
  [`mech_tau_sph()`](https://yejinjo0220.github.io/dppca/reference/mech_tau_sph.md)
  when private directions are computed.

- mono:

  Logical; whether to apply monotone post-processing.

- type:

  Either `"evr"` or `"scree"`.

- show_nonprivate:

  Logical; whether to overlay the non-private curve.

- show_points:

  Logical; whether to draw points on the plotted curves.

- col_map:

  Optional named color vector.

- lty_map:

  Optional named line-type vector.

- pch_map:

  Optional named point-shape vector.

- main, xlab, ylab, ylim:

  Plot controls.

- legend_pos:

  Legend position.

- ...:

  Additional graphical arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns a list containing non-private results and the
method-specific
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
outputs that were plotted.
