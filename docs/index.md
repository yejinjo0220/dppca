# dppca

`dppca` provides tools for **differentially private principal component
analysis (PCA) visualization**. The package focuses on
privacy-preserving versions of common PCA exploratory tools, including
PC direction estimation, scree/PVE plots, score plots, and
histogram-based summaries of PCA scores.

The main goal of `dppca` is to make it easy to compare non-private PCA
summaries with differentially private alternatives while keeping the
interface close to the standard PCA workflow in R.

## Installation

You can install the development version from GitHub with:

``` r
devtools::install_github("yejinjo0220/dppca")
```

After the package is available on CRAN, you will be able to install it
with:

``` r
install.packages("dppca")
```

## Overview

The package contains four main groups of functions.

### 1. Differentially private PC directions

[`dp_pc_dir()`](https://yejinjo0220.github.io/dppca/reference/dp_pc_dir.md)
estimates leading principal component directions under differential
privacy.

``` r
V_dp <- dp_pc_dir(
  X,
  k = 2,
  eps_total = 1,
  delta_total = 1e-6,
  center = TRUE,
  standardize = FALSE
)
```

The input `X` is the original data matrix with observations in rows.
Centering and standardization are handled internally through `center`
and `standardize`.

### 2. Differentially private scree and PVE plots

[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
computes private scree values or proportions of variance explained using
one of several private estimators.

``` r
res <- dp_scree(
  X,
  k = 5,
  method = "pmwm",
  eps_total = 1,
  delta_total = 1e-6,
  center = TRUE,
  standardize = FALSE,
  control = pmwm_control(a = 0, b = 10)
)
```

Available methods are:

- `"clipped"`: clipped mean based estimator
- `"pmwm"`: private modified winsorized mean based estimator
- `"huber"`: Huber-type robust estimator

Method-specific tuning parameters are supplied through control helper
functions:

``` r
clipped_control(C_clip = 3)

pmwm_control(
  beta = 1.01,
  a = 0,
  b = 10,
  trim_const = 10,
  eta = 0.01,
  split_mode = TRUE
)

huber_control(
  mu0 = 0,
  eta0 = 1,
  T = 50,
  M = 20,
  k_min_m2 = -40,
  k_max_m2 = 40,
  m2_frac = 1 / 4
)
```

To draw a private scree/PVE plot, use
[`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md):

``` r
dp_scree_plot(
  X,
  k = 5,
  dp_scree_method = "all",
  eps_total = 1,
  delta_total = 1e-6,
  center = TRUE,
  standardize = FALSE,
  type = "pve",
  control = list(
    clipped = clipped_control(C_clip = 3),
    pmwm = pmwm_control(a = 0, b = 10),
    huber = huber_control(T = 50, M = 20)
  )
)
```

Use `type = "pve"` to plot proportions of variance explained, or
`type = "scree"` to plot raw scree values. Plot styling is handled
internally so that the user-facing interface stays simple.

### 3. Differentially private score plots

[`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
computes differentially private PCA scores, and
[`dp_score_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot.md)
visualizes the first two private score coordinates.

``` r
score_res <- dp_score(
  X,
  k = 2,
  eps_total = 1,
  delta_total = 1e-6,
  center = TRUE,
  standardize = FALSE,
  bin_method = "WZ"
)

p <- dp_score_plot(
  X,
  k = 2,
  eps_total = 1,
  delta_total = 1e-6,
  center = TRUE,
  standardize = FALSE,
  bin_method = "WZ"
)

p
```

For bin recommendation in score-based visualizations, the available
options are:

- `"WZ"`
- `"Lei"`
- `"none"`

Grouped score plots are also supported:

``` r
dp_score_plot_group(
  X,
  group = group_labels,
  k = 2,
  eps_total = 1,
  delta_total = 1e-6,
  center = TRUE,
  standardize = FALSE,
  bin_method = "Lei"
)
```

### 4. Differentially private histogram summaries

The package also provides helper functions for histogram-based
visualization of private PCA scores. These functions are useful when the
goal is to summarize the distribution of projected observations rather
than release individual private score points.

The histogram routines support both single-sample and grouped
visualizations and are designed to work together with the score plot
workflow.

## Example datasets

`dppca` includes several example datasets for demonstrations and
vignettes:

- `adult`
- `eur_map`
- `eur_map_g`
- `gau`
- `gau_g`
- `gaussian_groups`

For example:

``` r
data(gau)
data(gau_g)

head(gau)
head(gau_g)
```

## Typical workflow

A typical analysis consists of the following steps.

``` r
library(dppca)

# 1. Load or prepare a numeric data matrix
X <- as.matrix(gau)

# 2. Estimate private PC directions
V_dp <- dp_pc_dir(
  X,
  k = 2,
  eps_total = 1,
  delta_total = 1e-6,
  center = TRUE,
  standardize = FALSE
)

# 3. Plot private PVE curves
pve_res <- dp_scree_plot(
  X,
  k = 5,
  dp_scree_method = "all",
  eps_total = 1,
  delta_total = 1e-6,
  type = "pve",
  control = list(
    clipped = clipped_control(C_clip = 3),
    pmwm = pmwm_control(a = 0, b = 10),
    huber = huber_control(T = 50, M = 20)
  )
)

# 4. Draw a private score plot
score_plot <- dp_score_plot(
  X,
  k = 2,
  eps_total = 1,
  delta_total = 1e-6,
  bin_method = "WZ"
)

score_plot

dp_score_plot_group(
  X = gau_g,
  G = "color",
  eps_total = 5,
  delta_total = 1e-6,
  center = TRUE,
  standardize = FALSE,
  bin_method = "WZ"
)
```

## Notes on preprocessing

Most user-facing PCA functions accept:

``` r
center = TRUE
standardize = FALSE
```

By default, columns are centered but not standardized. This corresponds
to a covariance-based PCA convention. Set `standardize = TRUE` when
correlation-based PCA is desired.

## Citation

If you use `dppca` in your work, please cite the package and the
associated methodological references described in the vignettes.

## License

This package is released under the license specified in the
`DESCRIPTION` file.
