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
