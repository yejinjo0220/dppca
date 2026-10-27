# dppca

`dppca` is an R package for **differentially private principal component
analysis (PCA) visualization**. It estimates private PC directions,
scree values, and two-dimensional score histograms, and uses these
quantities to create loading plots, scree plots, score plots, variable
plots, and biplots.

The package supports **standalone analyses** with
[`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md),
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md),
and
[`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md),
as well as a **shared PCA workflow** with
[`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md).
The fitting functions return S3 objects that can be inspected with
[`summary()`](https://rdrr.io/r/base/summary.html) and visualized with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Installation

Install the released package from CRAN:

``` r
install.packages("dppca")
```

To use the latest development code shown in the tutorials, install from
GitHub:

``` r
# install.packages("remotes")
remotes::install_github("yejinjo0220/dppca")
```

``` r
library(dppca)
```

## Quick start: one fit, five PCA visualizations

The integrated
[`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
function estimates a private PC direction matrix **once** and reuses it
for private scree and score estimation. It requires **one scree
estimator** (`"clipped"`, `"pmwm"`, or `"huber"`) and **one score
histogram estimator** (`"add"` or `"sparse"`).

``` r
library(dppca)
data("gau", package = "dppca")
X <- gau  # 1,000 observations, 20 numeric variables

set.seed(123)
pca_fit <- dp_pca(
  X,
  privacy = privacy_control(eps = 6, delta = 1e-4),
  scree = list(
    k = 5,
    method = "clipped",
    control = clipped_control(C_clip = 3)
  ),
  score = list(method = "add", bins = c(10, 10)),
  axes = c(1, 2)
)

pca_fit
summary(pca_fit)

set.seed(456)
pca_plots <- plot(pca_fit)  # Five DP panels
```

By default, `plot(pca_fit)` displays **PVE**, a **score histogram**, a
**loading vector plot**, a **variable plot**, and a **biplot**.

### Inspect the returned results

| Component | Description |
|:---|:---|
| `pca_fit$directions` | Private PC directions (columns of an orthonormal matrix). |
| `pca_fit$scree` | Estimated variances along the retained private directions. |
| `pca_fit$pve` | Proportions of estimated variation over the retained components. |
| `pca_fit$frame` | Privately estimated two-dimensional score frame. |
| `pca_fit$hist` | Private pooled score histogram. |
| `pca_fit$score_groups` | Per-group histograms, when a `group` is provided. |
| `pca_fit$metadata` | Dimensions, selected axes, methods, and privacy allocations. |
| `pca_fit$nonprivate` | Optional ordinary PCA reference; `NULL` by default. |

You can choose individual panels and their display settings without
refitting:

``` r
set.seed(456)
selected_plots <- plot(
  pca_fit,
  display = c("score", "biplot"),
  score = list(type = "sample"),
  biplot = list(alpha = 0, labels = TRUE),
  sampling = sampling_control(method = "uniform", sample_size = 300),
  center = TRUE,  # Shared by score and biplot (default)
  variables = 10,
  base_size = 12
)

# The individual ggplot objects are also stored:
selected_plots$private$plots$score
selected_plots$private$plots$biplot
```

The `variables` argument defaults to the **10 variables with the largest
Euclidean norms of their private loading coefficients on the displayed
two PC axes**. T he same selected variables are reused in loading,
variable, and biplot panels.

### Optional non-private comparison

To compare against ordinary PCA, add `non_private = TRUE` to the initial
[`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
call and fit again. Then, `plot(pca_fit)` automatically draws the
**ordinary reference page first**, followed by the DP page.

## Standalone estimators and plots

The standalone workflows offer additional method and display choices.

| Function | What it estimates | Plot interface |
|:---|:---|:---|
| [`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md) | Private PC directions | `plot(fit, type = "vector")`, `"point"`, or `"heatmap"` |
| [`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md) | Private scree and PVE with one or more of clipped, PMWM, and Huber methods | `plot(fit, type = "scree")` or `"pve"` |
| [`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md) | Two-dimensional additive and/or sparse score histograms | `plot(fit, type = "histogram")` or `"sample"` |
| [`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md) | Shared directions, scree, and histogram for five plots | `plot(fit, display = "all")` or selected panel names |

### Private loading directions

``` r
loading_fit <- dp_loading(
  X,
  privacy = privacy_control(eps = 3, delta = 1e-4)
)

loading_fit$directions
loading_plot <- plot(loading_fit, type = "vector", variables = 10)
```

### Scree / PVE

``` r
scree_fit <- dp_scree(
  X,
  k = 5,
  privacy = privacy_control(eps = 3, delta = 1e-4),
  method = "clipped",
  control = clipped_control(C_clip = 3)
)

summary(scree_fit)
plot(scree_fit, type = "pve")
```

For multiple methods,
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
also accepts `method = c("clipped", "pmwm", "huber")` with suitable
method-specific controls.

### Score histograms and synthetic samples

``` r
score_fit <- dp_score(
  X,
  privacy = privacy_control(eps = 3, delta = 1e-4),
  bins = c(10, 10),
  method = c("add", "sparse")
)

summary(score_fit)
set.seed(456)
score_plots <- plot(
  score_fit,
  type = c("histogram", "sample"),
  center = FALSE,  # Optional: disable the default display centering
  sampling_control = sampling_control(method = "uniform", sample_size = 300)
)
```

By default, `center = TRUE` translates the **displayed private score
coordinates** by the midpoint of the private frame. Set `center = FALSE`
to leave those coordinates unshifted. This affects the plot, not the
stored frame or histogram probabilities. Grouped score plots use the
same `center` option.

## Grouped PCA visualizations

A grouping column is excluded from PCA features and used to color
synthetic score points and form group histograms. The `gau_g` dataset
contains the same 20 numeric features as `gau`, plus a `group` label:

``` r
data("gau_g", package = "dppca")

set.seed(123)
group_fit <- dp_pca(
  gau_g,
  privacy = privacy_control(eps = 6, delta = 1e-4),
  scree = list(k = 5, method = "clipped",
               control = clipped_control(C_clip = 3)),
  score = list(method = "add", bins = c(10, 10)),
  group = "group"
)

set.seed(456)
plot(group_fit, display = c("score", "biplot"),
     score = list(type = "sample"), center = FALSE)
```

Standalone grouped score analysis is available with
`dp_score(..., group = "group")`.

**Privacy note:** For the grouped histogram privacy argument, group
labels and the group partition are assumed public and fixed. Private
sampled points are not original score records.

## Privacy-budget allocation

[`privacy_control()`](https://yejinjo0220.github.io/dppca/reference/privacy_control.md)
specifies the total `(eps, delta)` and optional named `split`. If no
split is given, each fitting function uses an **equal allocation among
its required components**:

| Fitting function | Budget components | Implicit allocation |
|:---|:---|:---|
| [`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md) | `loading` | 100% |
| [`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md) | `loading`, `scree` | 50% / 50% |
| [`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md) | `loading`, `score` | 50% / 50% |
| [`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md) | `loading`, `scree`, `score` | 1/3 each |

For example, customize the integrated budget as follows:

``` r
privacy <- privacy_control(
  eps = 6,
  delta = 1e-4,
  split = c(loading = 0.2, scree = 0.4, score = 0.4)
)

privacy$total
privacy$budgets
```

The selected `score` epsilon allocation is further divided **35% for the
private frame** and **65% for the histogram**; its delta allocation is
assigned to the histogram.

## Example datasets

| Dataset | Shape | Description |
|:---|:---|:---|
| `gau` | 1,000 × 20 | Five synthetic Gaussian clusters with numeric features `V1`–`V20`. |
| `gau_g` | 1,000 × 21 | The same observations and features, plus the `group` column. |
| `vowel` | 990 × 14 | British English vowel data: 10 acoustic measurements and four categorical descriptors. |

The `vowel` dataset originates from the UCI Machine Learning Repository
and is distributed under the terms described in its package
documentation. See
[`?vowel`](https://yejinjo0220.github.io/dppca/reference/vowel.md) for
the data source and license. Dataset details and preprocessing examples
are provided in the Getting Started tutorial.

## Tutorials and theory

The **practical tutorials** are:

- [Getting
  Started](https://yejinjo0220.github.io/dppca/articles/tutorials-getting-started.html):
  installation, datasets, preprocessing, and privacy budgets.
- [DP Loading
  Plot](https://yejinjo0220.github.io/dppca/articles/tutorials-loading-plots.html):
  private PC directions and loading displays.
- [DP Scree
  Plot](https://yejinjo0220.github.io/dppca/articles/tutorials-scree-plots.html):
  clipped, PMWM, and Huber estimation; scree and PVE.
- [DP Score
  Plot](https://yejinjo0220.github.io/dppca/articles/tutorials-score-plots.html):
  histogram methods, synthetic points, and grouped scores.
- [DP PCA
  Visualization](https://yejinjo0220.github.io/dppca/articles/tutorials-dp-visualization.html):
  five plots built from a shared fit.

For the mathematics, see [PCA
Background](https://yejinjo0220.github.io/dppca/articles/theory-pca-background.html),
[Theory of DP
Loading](https://yejinjo0220.github.io/dppca/articles/theory-dp-loading.html),
[Theory of DP
Scree](https://yejinjo0220.github.io/dppca/articles/theory-dp-scree.html),
[Theory of DP
Score](https://yejinjo0220.github.io/dppca/articles/theory-dp-score.html),
and [Theory of DP PCA
Visualizations](https://yejinjo0220.github.io/dppca/articles/theory-dp-visualization.html).

For complete argument descriptions, use the R help pages:
[`?dp_loading`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md),
[`?dp_scree`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md),
[`?dp_score`](https://yejinjo0220.github.io/dppca/reference/dp_score.md),
[`?dp_pca`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md),
and the corresponding `plot.*` methods. The development implementation
is available on [GitHub](https://github.com/yejinjo0220/dppca).

## Online Shiny app

Explore the `dppca` visualizations interactively in your browser,
without installing R:

**[Launch the online DP-PCA Shiny
app](https://yejinjo0220-dppca.share.connect.posit.cloud/)**

See the [Shiny App
guide](https://yejinjo0220.github.io/dppca/articles/shiny-app.html) for
the available analyses, data options, and privacy considerations.

## References

- Kim, M. and Jung, S. (2025). Robust and Differentially Private
  Principal Component Analysis. *Statistical Analysis and Data Mining:
  An ASA Data Science Journal*, 18(6), e70053.
  <https://doi.org/10.1002/sam.70053>
- Dwork, C. and Roth, A. (2014). The Algorithmic Foundations of
  Differential Privacy. *Foundations and Trends in Theoretical Computer
  Science*, 9(3–4), 211–407. <https://doi.org/10.1561/0400000042>
- Ramsay, K. and Spicker, D. (2025). Improved subsample-and-aggregate
  via the private modified winsorized mean. arXiv:2501.14095.
- Yu, M., Ren, Z., and Zhou, W.-X. (2024). Gaussian differentially
  private robust mean estimation and inference. *Bernoulli*, 30(4),
  3059–3088.
- Nissim, K., Raskhodnikova, S., and Smith, A. (2007). Smooth
  Sensitivity and Sampling in Private Data Analysis. In *STOC ’07*,
  75–84. <https://doi.org/10.1145/1250790.1250803>
- Wasserman, L. and Zhou, S. (2010). A Statistical Framework for
  Differential Privacy. *Journal of the American Statistical
  Association*, 105(489), 375–389.
  <https://doi.org/10.1198/jasa.2009.tm08651>
- Karwa, V. and Vadhan, S. P. (2017). Finite Sample Differentially
  Private Confidence Intervals. arXiv:1711.03908.
