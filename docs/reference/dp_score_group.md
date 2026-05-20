# Group-wise differentially private score histograms

This function computes two-dimensional principal component scores and
releases group-wise differentially private histograms on a common score
frame and grid. It is useful when observations have group labels and the
low-dimensional score distribution should be compared across groups.

## Usage

``` r
dp_score_group(
  X,
  group,
  eps,
  delta,
  bins,
  method = c("add", "sparse"),
  center = TRUE,
  standardize = FALSE,
  g_dppca = FALSE,
  cpp.option = FALSE,
  axes = c(1, 2),
  fixed_frame = NULL
)
```

## Arguments

- X:

  A matrix or data frame where rows correspond to observations and
  columns correspond to variables. `X` can additionally include a named
  column representing the group label for each observation.

- group:

  Group labels. This can be a vector of length `nrow(X)` or a single
  column name in `X`. If a column name is supplied, that column is used
  as the group label and removed from the feature matrix.

- eps:

  Positive number defining the total `epsilon` privacy parameter.

- delta:

  Number in `(0, 1)` defining the total `delta` privacy parameter.

- bins:

  Integer vector of length 2 defining the number of histogram bins along
  the first and second score axes, respectively.

- method:

  Character vector specifying which private histogram methods to
  compute. Use `"add"` for the additive Gaussian histogram and
  `"sparse"` for the sparse thresholded histogram. The default is
  `c("add", "sparse")`.

- center:

  A logical value indicating whether to center the columns of `X` before
  computing principal component directions. The default is `TRUE`.

- standardize:

  A logical value indicating whether to scale the columns of `X` by
  their sample standard deviations after optional centering. The default
  is `FALSE`.

- g_dppca:

  A logical value indicating whether to use private principal component
  directions. The default is `FALSE`. See
  [`dp_pc_dir()`](https://yejinjo0220.github.io/dppca/reference/dp_pc_dir.md)
  for details.

- cpp.option:

  A logical value passed to
  [`dp_pc_dir()`](https://yejinjo0220.github.io/dppca/reference/dp_pc_dir.md)
  when `g_dppca = TRUE`. The default is `FALSE`.

- axes:

  Integer vector of length 2 specifying the principal components used to
  construct the score coordinates. The default is `c(1, 2)`.

- fixed_frame:

  Optional fixed plotting frame. If supplied, it can be a numeric vector
  `c(lower, upper)` used for both axes, or a list with numeric
  components `xlim` and `ylim`. If `NULL`, a private square frame is
  estimated from the score coordinates.

## Value

A list with components:

- score:

  An \\n \times 2\\ matrix containing the PC scores for the two selected
  axes.

- frame:

  A list with components `xlim` and `ylim`.

- groups:

  A named list of group-specific histogram outputs.

- method:

  Character vector of private histogram methods used.

## Details

The score directions, plotting frame, and histogram grid are shared
across all groups. For each group \\g\\, the group-specific count in bin
\\B_k\\ is \\c_k^{(g)} = \sum_i 1\\s_i \in B_k, g_i = g\\\\. Private
histograms are then computed separately for each group on the common
grid. Because the groups form a partition of the rows, the group-wise
histograms use the same histogram privacy parameters for each group by
parallel composition.

## See also

[`dp_score_plot_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot_group.md)
for plotting group-wise score histograms.
[`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
for pooled score histograms.

## Examples

``` r
data(gau_g, package = "dppca")

# Compute private grouped PCA scores.
set.seed(123)
score_gau_g <- dp_score_group(
  gau_g,
  group = "group",
  eps = 3,
  delta = 1e-3,
  bins = c(8, 8),
)

head(score_gau_g$score)
#>            PC1        PC2
#> [1,] 0.5752441  0.9711984
#> [2,] 3.1987535  1.5765379
#> [3,] 2.0936989  1.0882390
#> [4,] 0.9725891  1.5285036
#> [5,] 1.0839100  0.4533152
#> [6,] 1.5430497 -1.0954815
head(score_gau_g$groups$group1$add)
#>         xmin       xmax      ymin     ymax         prob
#> 1 -24.964534 -18.961402 -24.96453 -18.9614 0.0000000000
#> 2 -18.961402 -12.958271 -24.96453 -18.9614 0.0102258822
#> 3 -12.958271  -6.955139 -24.96453 -18.9614 0.0004625692
#> 4  -6.955139  -0.952008 -24.96453 -18.9614 0.0008481902
#> 5  -0.952008   5.051123 -24.96453 -18.9614 0.0112516577
#> 6   5.051123  11.054255 -24.96453 -18.9614 0.0030238337
```
