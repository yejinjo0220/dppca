# Plot differentially private score histograms

This function computes and visualizes two-dimensional principal
component score histograms, including the original scatter plot, the
non-private empirical histogram, and one or more differentially private
histogram estimates. It is a plotting wrapper around
[`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
and returns both the computed score output and `ggplot` objects.

## Usage

``` r
dp_score_plot(
  X,
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

  A numeric matrix or data frame. Rows correspond to observations and
  columns correspond to variables.

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

  The output of
  [`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md).

- plot:

  A list containing the scatter plot, histogram panels, and the combined
  patchwork plot.

## See also

[`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
for computing score histograms without plotting.
[`dp_score_plot_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot_group.md)
for group-wise score histogram plots.

## Examples

``` r
data(gau, package = "dppca")

# Use a small subset to keep the example fast.
X <- gau[1:300, ]

# Draw a private score plot using the additive histogram method.
set.seed(123)
score_plot <- dp_score_plot(
  X,
  eps = 3,
  delta = 1e-3,
  bins = c(8, 8),
  method = "add"
)
score_plot$plot$add


# Draw score plots for all available histogram methods.
set.seed(123)
score_plot <- dp_score_plot(
  X,
  eps = 3,
  delta = 1e-3,
  bins = c(8, 8)
)
score_plot$plot$all

```
