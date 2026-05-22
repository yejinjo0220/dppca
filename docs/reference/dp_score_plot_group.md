# Plot group-wise differentially private score histograms

This function computes and visualizes group-wise differentially private
score histograms. It is a plotting wrapper around
[`dp_score_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_group.md)
and returns both the computed group-wise score output and `ggplot`
objects.

## Usage

``` r
dp_score_plot_group(
  X,
  group,
  eps,
  delta,
  bins,
  center = TRUE,
  standardize = FALSE,
  g_dppca = FALSE,
  cpp.option = FALSE,
  axes = c(1, 2),
  method = c("add", "sparse")
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

- method:

  Character vector specifying which private histogram methods to
  compute. Use `"add"` for the additive Gaussian histogram and
  `"sparse"` for the sparse thresholded histogram. The default is
  `c("add", "sparse")`.

## Value

A list with components:

- score:

  The output of
  [`dp_score_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_group.md).

- plot:

  A list containing group-wise histogram plots.

- group_colors:

  Named vector of colors used for the groups.

## See also

[`dp_score_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_group.md)
for computing group-wise score histograms without plotting.
[`dp_score_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot.md)
for pooled score histogram plots.

## Examples

``` r
data(gau_g, package = "dppca")

# Draw a private grouped score plot.
set.seed(123)
score_plot_gau_g <- dp_score_plot_group(
  gau_g,
  group = "group",
  eps = 3,
  delta = 1e-3,
  bins = c(8, 8)
)

score_plot_gau_g$plot$all

```
