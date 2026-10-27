# Differentially private score histograms

This function projects observations onto supplied private directions or
directions estimated by
[`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
and computes two-dimensional score histograms. Observed coordinates and
empirical reference histograms are returned only when
`non_private = TRUE`.

## Usage

``` r
dp_score(
  X,
  privacy,
  bins,
  method = c("add", "sparse"),
  cpp.option = TRUE,
  axes = c(1, 2),
  V_dp = NULL,
  non_private = FALSE,
  group = NULL
)

# S3 method for class 'dp_score'
summary(object, ...)

# S3 method for class 'dp_score_group'
summary(object, ...)

# S3 method for class 'summary.dp_score'
print(x, digits = 4, ...)

# S3 method for class 'dp_score'
print(x, ...)

# S3 method for class 'dp_score_group'
print(x, ...)
```

## Arguments

- X:

  A numeric matrix or data frame with observations in rows and variables
  in columns. Non-numeric data-frame columns are excluded with a
  warning; numeric values must be finite. At least two observations and
  two numeric features are required. If directions must be estimated,
  the pairwise-separation condition of
  [`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
  also applies.

- privacy:

  A
  [`privacy_control()`](https://yejinjo0220.github.io/dppca/reference/privacy_control.md)
  object allocating budget to `loading` and `score` (equal shares by
  default). Each requested histogram method receives the full histogram
  allocation.

- bins:

  Positive integer vector of length two giving the number of histogram
  bins along the selected score axes.

- method:

  One or both histogram methods: `"add"` (Gaussian additive noise) and
  `"sparse"` (thresholded sparse histogram). Default:
  `c("add", "sparse")`.

- cpp.option:

  Whether to use the compiled spherical-Kendall implementation (`TRUE`
  by default) when estimating private directions; ignored if `V_dp` is
  supplied.

- axes:

  Two distinct PC indices for the horizontal and vertical score axes,
  respectively (`c(1, 2)` by default). Each must be between `1` and the
  number of numeric PCA features.

- V_dp:

  Optional `p` by `p` private PC direction matrix (`NULL` by default),
  where `p` is the number of numeric PCA features. Columns must be
  orthonormal and rows must follow the feature order in `X`; matching
  row names are required when both sets of names exist.

- non_private:

  Whether to compute and store ordinary PCA reference results (`FALSE`
  by default). For score analysis, these include mean-centered ordinary
  PCA scores and empirical reference histograms.

- group:

  Optional group column name in `X` or vector of group labels (`NULL` by
  default). A selected group column is excluded from PCA features, even
  when numeric; group labels and membership are assumed public. All
  groups share the same private directions, frame, and histogram grid.

- object:

  A result from `dp_score()`, with or without `group`.

- ...:

  Reserved for S3 generic compatibility.

- x:

  A result from `dp_score()`, with or without `group`.

- digits:

  Number of significant digits shown in tables.

## Value

A `dp_score` object containing the private `frame` and requested
histogram methods (`add` and/or `sparse`), or `groups` for grouped fits.
Each histogram has `xmin`, `xmax`, `ymin`, `ymax` and `prob` columns.
`nonprivate` is `NULL` by default; when requested it contains `scores`,
`frame` and the empirical `histogram` (or group histograms under
`groups`). `metadata` contains data dimensions, axes, bins, methods and
privacy settings. Use [`summary()`](https://rdrr.io/r/base/summary.html)
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) to inspect
or visualize stored results.

A summary containing the public sample size, axes, grid, private frame
(including its center), allocated privacy budget and a table of
positive-bin counts and probability sums. No data-dependent private
mechanisms are run.

Invisibly returns `x`.

Invisibly returns `x`.

## Details

Scores are projected onto the selected private PC directions. Rather
than releasing individual projected coordinates, the function estimates
a private two-dimensional plotting frame followed by private histogram
probabilities. The frame uses two private coordinate medians and an
inflated private 0.99 quantile of distances from their center (inflation
factor `1.2`), using unbounded quantile mechanisms.

For a fixed frame and directions under replacement adjacency, the
histogram count vector has `l1` sensitivity at most `2` and `l2`
sensitivity at most `sqrt(2)`. Out-of-frame scores are assigned to
boundary bins. Adaptive private directions and frame selection require
privacy composition. Input is not automatically centered or
standardized; data-dependent mean centering or standardization is not
covered by the count sensitivity argument.

Two histogram mechanisms are available:

- `"add"`: Add Gaussian noise to all counts, truncate negatives to zero,
  and normalize. The classical Gaussian calibration has a standard
  sufficient guarantee for `0 < eps_hist < 1`; larger values require a
  separate calibration check. All-zero truncated counts cause an error.

- `"sparse"`: Add Laplace noise to nonzero bin proportions, retain bins
  above a privacy-dependent threshold, and normalize positive mass. If
  none survive, return an all-zero histogram.

The `score` epsilon allocation is split 35% for the private frame and
65% for the histogram; all score delta goes to the histogram. The frame
epsilon is divided equally among two private medians and one radial
quantile. Directions and frame are shared across methods. Each requested
histogram method receives the full histogram allocation, so jointly
releasing `"add"` and `"sparse"` requires composition across methods.

For grouped data, the partition and labels are assumed public and fixed;
groups share directions, frame and grid. Parallel composition applies to
disjoint group histograms under this assumption, while different methods
still require composition. `non_private = FALSE` does not privatize
group membership. Input dimensions, feature schema, bins, and tuning
choices are assumed public. Repeated releases also require composition.

The optional ordinary PCA reference uses sample-mean-centered scores and
a zero-centered comparison frame. It is not differentially private. See
[`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
for direction-estimation and Gaussian calibration conditions.

## References

Dwork C, Roth A (2014). “The Algorithmic Foundations of Differential
Privacy.” *Found. Trends Theor. Comput. Sci.*, **9**(3–4), 211–407. ISSN
1551-305X, [doi:10.1561/0400000042](https://doi.org/10.1561/0400000042)
.

Durfee D (2023). “Unbounded Differentially Private Quantile and Maximum
Estimation.” In *Advances in Neural Information Processing Systems*,
volume 36, 77691–77712.
[doi:10.52202/075280-3394](https://doi.org/10.52202/075280-3394) ,
<https://proceedings.neurips.cc/paper_files/paper/2023/hash/f4b6ef2a78684dca2fb3f1c09372e041-Abstract-Conference.html>.

Wasserman L, Zhou S (2010). “A Statistical Framework for Differential
Privacy.” *Journal of the American Statistical Association*,
**105**(489), 375–389.
[doi:10.1198/jasa.2009.tm08651](https://doi.org/10.1198/jasa.2009.tm08651)
.

Karwa V, Vadhan S (2018). “Finite Sample Differentially Private
Confidence Intervals.” In *Proceedings of the 9th Innovations in
Theoretical Computer Science Conference*, volume 94 of *Leibniz
International Proceedings in Informatics*, 44:1–44:9.
[doi:10.4230/LIPIcs.ITCS.2018.44](https://doi.org/10.4230/LIPIcs.ITCS.2018.44)
.

Kim M, Jung S (2025). “Robust and Differentially Private Principal
Component Analysis.” *Statistical Analysis and Data Mining: An ASA Data
Science Journal*, **18**(6), e70053.
[doi:10.1002/sam.70053](https://doi.org/10.1002/sam.70053) .

## See also

[`plot.dp_score()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score.md)
for plotting stored results without new private estimates.
[`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
and
[`plot.dp_pca()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_pca.md)
for estimation followed by plots of a stored fit.
[`plot.dp_score_group()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score_group.md)
for displaying a grouped result.
[`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
for private principal component direction estimation.

## Examples

``` r
data(gau, package = "dppca")

# Use a small subset to keep the example fast.
X <- gau[1:300, ]

# Release a private histogram with fixed preprocessing and no raw scores.
set.seed(123)
score_gau <- dp_score(
  X,
  privacy = privacy_control(eps = 1.8, delta = 1e-3),
  method = "add",
  bins = c(10, 10),
  non_private = FALSE
)

score_gau$frame
#> $xlim
#> [1] -12.56812  12.50752
#> 
#> $ylim
#> [1] -12.51772  12.55792
#> 
head(score_gau$add)
#>         xmin       xmax      ymin      ymax        prob
#> 1 -12.568121 -10.060557 -12.51772 -10.01016 0.000000000
#> 2 -10.060557  -7.552993 -12.51772 -10.01016 0.000000000
#> 3  -7.552993  -5.045429 -12.51772 -10.01016 0.005782674
#> 4  -5.045429  -2.537865 -12.51772 -10.01016 0.007195322
#> 5  -2.537865  -0.030301 -12.51772 -10.01016 0.000000000
#> 6  -0.030301   2.477263 -12.51772 -10.01016 0.023915061
summary(score_gau)
#> DP Score Histograms Summary
#> Data: 300 observations x 20 variables
#> Score axes: PC1, PC2 
#> Grid: 10 x 10 
#> Methods: score = add 
#> 
#> Private frame:
#>  Axis  Center  Lower Upper
#>   PC1 -0.0303 -12.57 12.51
#>   PC2  0.0201 -12.52 12.56
#> 
#> Histogram summary:
#>  Method Bins Positive bins Probability sum
#>     add  100            59               1
#> 
#> Allocated privacy budget:
#>                   eps delta
#> Loading         0.900 5e-04
#> Score frame     0.315 0e+00
#> Score histogram 0.585 5e-04
#> Total           1.800 1e-03
#> 
#> Non-private reference: not stored 
# plot(score_gau, type = "histogram")
# To use a group column: dp_score(X, privacy, bins, group = "region")
```
