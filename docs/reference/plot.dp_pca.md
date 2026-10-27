# Plot an integrated PCA fit

Uses the stored estimates without estimating loading, scree or score
again. Use [`set.seed()`](https://rdrr.io/r/base/Random.html) before
plotting to reproduce synthetic samples.

## Usage

``` r
# S3 method for class 'dp_pca'
plot(
  x,
  display = "all",
  scree = list(),
  score = list(),
  loading = list(),
  variable = list(),
  biplot = list(),
  sampling = NULL,
  center = TRUE,
  variables = 10,
  base_size = 12,
  ...
)
```

## Arguments

- x:

  A fitted
  [`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
  object.

- display:

  `"all"` (default) or a vector selecting `"scree"`, `"score"`,
  `"loading"`, `"variable"`, and `"biplot"`.

- scree:

  Plot options: `type = "pve"` or `"scree"`, and `plot_control = NULL`
  ([`scree_plot_control()`](https://yejinjo0220.github.io/dppca/reference/scree_plot_control.md)).

- score:

  Plot options: `type = "histogram"` or `"sample"`, and
  `plot_control = NULL`
  ([`score_plot_control()`](https://yejinjo0220.github.io/dppca/reference/score_plot_control.md)).
  For grouped plots, the color map and legend settings also apply to the
  biplot.

- loading:

  Plot options: `type = "vector"`, `"point"`, or `"heatmap"`,
  `heatmap_components = NULL`, and `plot_control = NULL`
  ([`loading_plot_control()`](https://yejinjo0220.github.io/dppca/reference/loading_plot_control.md)).
  Multiple types require `display = "loading"`.

- variable:

  Plot options: `plot_control = NULL`. Shows labeled arrows and the
  correlation-style unit circle.

- biplot:

  Plot options: `alpha = 0` (between 0 and 1), `labels = TRUE`, and
  `plot_control = NULL`.

- sampling:

  A
  [`sampling_control()`](https://yejinjo0220.github.io/dppca/reference/sampling_control.md)
  object with one method, shared by the score sample plot and biplot. By
  default, DP sample size is the input sample size, and non-private
  plots use observed scores. An explicit `sample_size` samples both
  versions from their histograms. For groups, sizes are allocated across
  groups; zero-mass groups yield no points.

- center:

  Whether to shift displayed private scores by the midpoint of the
  private frame (`TRUE` by default). Shared by score and biplot; the
  biplot arrows start at the origin when `TRUE`. Ordinary reference
  scores remain centered at zero. This is display-only post-processing.

- variables:

  Variables to display: positive top-N count (`10` by default), `"all"`,
  `NULL`, names, or indices. Top-N selection uses norms of private
  loading coefficients on the two selected axes. Selection is shared by
  loading, variable and biplot panels, including non-private
  comparisons.

- base_size:

  Base font size (`12` by default); local plot controls override it.

- ...:

  Unused arguments are rejected.

## Value

Invisibly returns `private` and, when requested, `nonprivate` lists.
Each contains `plots`, `page`, `sampled_points` (unshifted), `biplot`,
`frames` and `direction_sign`. Biplot metadata includes the
alpha-transformed and displayed coordinates, `scale`, `arrow_scale`,
`score_center`, `score_radius`, `estimated_center`, `centered`,
`arrow_start`, `arrow_endpoints` and display frame. The plots are always
drawn; `out <- plot(fit)` also keeps these results for reuse. Multiple
loading displays with both versions additionally return a `comparison`
patchwork.

## Details

Plotting reuses the stored estimates; non-private panels appear
automatically if the fit contains them. Synthetic sampling is
post-processing, and [`set.seed()`](https://rdrr.io/r/base/Random.html)
can reproduce sampled displays.

Variable-plot coordinates use the covariance reconstructed from the `k`
private scree estimates and directions; they need not equal observed
variable–score correlations. Variables with zero implied variance are
omitted.

Biplot coordinates divide sampled scores by `d^alpha` and multiply
private direction coefficients by `d^alpha`, where
`d = sqrt((n - 1) * scree)` on the selected axes. A zero scree value
with positive alpha suppresses that axis. Arrows use one display
multiplier computed from the private sample and shared with the ordinary
comparison; empty or zero-radius samples omit arrows. The comparison
uses the same square biplot frame.

Centering translates only displayed coordinates. Combined grouped pages
show one group legend, while individual ggplot objects retain their own.
The returned object stores both the combined pages and individual plots.
