# dppca: Differentially Private PCA and Histogram-Based Visualization

`dppca` is an R package that enables differentially private visualization of high-dimensional datasets through PCA projection and DP histograms.  
The package provides DP-PCA, DP quantile/frame estimation, DP histograms, sampling from histograms, and visualization tools for both single and multi-group datasets.

---

## Key Features

### Differential Privacy
- Differentially Private PCA using the spherical Gaussian mechanism.
- DP min/max estimation using smooth sensitivity (Nissim–Raskhodnikova–Smith).
- DP histograms using additive Gaussian noise or sparse Laplace mechanisms.

### Visualization Tools
- `dp_score_plot()` — DP scatter + DP histogram visualization.
- `dp_score_plot_group()` — Same as above, but with group labels.
- Both functions support DP-PCA or non-DP PCA.

### Histogram Utilities
- `dp_hist()` — Core DP histogram estimator.
- `dp_hist_group()` — Histogram per group.
- `sample_from_hist()` — Sample synthetic data from DP histograms.

### Helper Functions
- `dp_quantile_ss()` — Smooth-sensitivity based DP quantile estimation.
- `dp_frame()` — Generates DP plotting frame.
- `number_bins()` — Automatically chooses 2D bin counts.
- Internal helper functions for plotting (not exported).

---

## Installation

Since the package is under development, install from GitHub:

```r
# install.packages("devtools")
devtools::install_github("yejinjo0220/dppca")
library(dppca)


## Example Usage

```r
library(dppca)

data("eur_map")
data("eur_map_g")
data("gau")
data("gau_g")
data("adult")

# 1. European Map Example
res_eur <- dp_score_plot(
  X = eur_map,
  eps_total = 4, delta_total = 1e-4,
  sampling  = TRUE
)
res_eur$plot$all

res_eur_g <- dp_score_plot_group(
  X = eur_map_g,
  G = "color",
  eps_total = 4, delta_total = 1e-4,
  sampling  = TRUE
)
res_eur_g$plot$all

# 2. Gaussian Mixture Example
res_gau <- dp_score_plot(
  X = gau,
  eps_total = 4, delta_total = 1e-5,
  sampling  = TRUE
)
res_gau$plot$all

res_gau_g <- dp_score_plot_group(
  X = gau_g,
  G = "color",
  eps_total = 4, delta_total = 1e-5,
  sampling  = TRUE
)
res_gau_g$plot$all

# 3. Adult Census Example
res_adult <- dp_score_plot(
  X = adult,
  scale. = TRUE,
  eps_total = 4,
  delta_total = 1e-5,
  sampling = TRUE
)
res_adult$plot$all

# 4. Run shiny with Built-in data
run_dp_score()

# 5. Run shiny with your data
X <- matrix(rnorm(500*10), 500, 10)
G <- sample(c("A","B","C"), 500, TRUE)

dp_score_app(X, G)
