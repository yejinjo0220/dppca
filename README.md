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
make_spiked_data <- function(n, d, spikes = c(9, 6, 4, 2), seed = 1) {
  set.seed(seed)
  lambda_pop <- c(spikes, rep(1, d - length(spikes)))
  A <- matrix(rnorm(d * d), d, d)
  Q <- qr.Q(qr(A))
  Sigma_pop <- Q %*% diag(lambda_pop) %*% t(Q)
  Z <- matrix(rnorm(n * d), n, d)
  X <- Z %*% chol(Sigma_pop)
  X
}

n <- 5000
d <- 20
k <- 8

X <- make_spiked_data(n, d, seed = 1)


dp_scree_plot(
  X = X,
  k =4,
  dp_scree_method = "all",
  eps_total = 5,
  delta_total = 1e-5,
  a = 0,
  b = 100
)
```

## Example: DP score plot

```r
library(dppca)

data("eur_map")

res <- dp_score_plot(
  X = eur_map,
  eps_total = 4, delta_total = 1e-4,
)

res
res$score$s
res$plot$sparse
```

## Example: group-wise DP score histogram

```r
data("eur_map_g")

res_eur_g <- dp_score_plot_group(
  X = eur_map_g,
  G = "color",
  eps_total = 4, delta_total = 1e-4,
)

res_eur_g$plot$all
```

## Notes

- `g_dppca = TRUE` allows privatization of the PCA direction matrix.
- `dp_scree_plot()` draws the scree figure directly and invisibly returns the underlying results.
- `dp_score_plot()` and `dp_score_plot_group()` return both calculation outputs and plot objects.
- For PMW-based scree estimation, you must provide support bounds such as `a` and `b`.
- For the clipped scree estimator, you must provide `C_clip`.

## License

MIT
