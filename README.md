# dppca

`dppca` is an R package for differentially private PCA visualization.

The package focuses on two main tasks:

- **DP scree estimation** for privatizing leading scree values and explained variance ratios
- **DP score visualization** for constructing private 2D score histograms instead of releasing individual PCA score points

## Main features

### 1. DP scree estimation

The package provides several approaches for differentially private scree estimation:

- **Clipped mean**
- **Huber-type private mean**
- **PMW-style private winsorized mean**

User-facing functions:

- `dp_scree()`
- `dp_scree_plot()`

### 2. DP score histograms

Instead of plotting raw PCA score points, the package constructs a private 2D histogram on the score space.

It includes:

- empirical histogram
- Gaussian additive DP histogram
- sparse DP histogram
- group-wise DP score histograms

User-facing functions:

- `dp_score()`
- `dp_score_plot()`
- `dp_score_group()`
- `dp_score_plot_group()`

## Installation

You can install the package from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("yejinjo0220/dppca")
```

## Package structure

Main R files:

- `R/helper_pca.R`
- `R/helper_scree.R`
- `R/dp_scree.R`
- `R/helper_score.R`
- `R/dp_score.R`

## Example: DP scree estimation

```r
library(dppca)

set.seed(1)
n <- 200
d <- 8
X <- matrix(rnorm(n * d), nrow = n, ncol = d)

res_huber <- dp_scree(
  X = X,
  k = 3,
  method = "huber",
  eps_total = 1,
  delta_total = 1e-5,
  mu0 = 0,
  eta0 = 0.5,
  T = 5,
  M = 5,
  k_min_m2 = -10,
  k_max_m2 = 10,
  m2_frac = 0.25
)

res_huber
```

## Example: DP scree plot

```r
dp_scree_plot(
  X = X,
  k = 3,
  dp_scree_method = "all",
  eps_total = 1,
  delta_total = 1e-5,
  C_clip = 3,
  beta = 1.05,
  a = 0,
  b = 10,
  trim_const = 10,
  eta = 0.01,
  split_mode = TRUE,
  mu0 = 0,
  eta0 = 0.5,
  T = 5,
  M = 5,
  k_min_m2 = -10,
  k_max_m2 = 10,
  m2_frac = 0.25
)
```

## Example: DP score histogram

```r
library(dppca)

set.seed(1)
n <- 200
d <- 6
X <- matrix(rnorm(n * d), nrow = n, ncol = d)

res_score <- dp_score(
  X = X,
  center = TRUE,
  standardize = FALSE,
  g_dppca = FALSE,
  axes = c(1, 2),
  eps_total = 1,
  delta_total = 1e-5,
  inflate = 0.2,
  bin_method = "J"
)

str(res_score)
```

## Example: DP score plot

```r
res_plot <- dp_score_plot(
  X = X,
  center = TRUE,
  standardize = FALSE,
  g_dppca = FALSE,
  axes = c(1, 2),
  eps_total = 1,
  delta_total = 1e-5,
  inflate = 0.2,
  bin_method = "J"
)

res_plot$plot$all
```

## Example: group-wise DP score histogram

```r
G <- sample(c("A", "B"), n, replace = TRUE)

res_group <- dp_score_group(
  X = X,
  G = G,
  center = TRUE,
  standardize = FALSE,
  g_dppca = FALSE,
  axes = c(1, 2),
  eps_total = 1,
  delta_total = 1e-5,
  inflate = 0.2,
  bin_method = "J"
)

str(res_group)
```

## Notes

- `g_dppca = TRUE` allows privatization of the PCA direction matrix.
- `dp_scree_plot()` draws the scree figure directly and invisibly returns the underlying results.
- `dp_score_plot()` and `dp_score_plot_group()` return both calculation outputs and plot objects.
- For PMW-based scree estimation, you must provide support bounds such as `a` and `b`.
- For the clipped scree estimator, you must provide `C_clip`.

## License

MIT
