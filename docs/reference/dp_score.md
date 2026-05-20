# Differentially private score histograms

This function computes two-dimensional principal component scores and
returns differentially private histogram estimates on the score space.
It returns the score coordinates, the plotting frame, the non-private
histogram, and the requested private histogram estimates.

## Usage

``` r
dp_score(
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

  An \\n \times 2\\ matrix containing the PC scores for the two selected
  axes.

- frame:

  A list with components `xlim` and `ylim`.

- none:

  Data frame for the non-private empirical histogram.

- add:

  Data frame for the additive Gaussian private histogram, or `NULL` if
  not requested.

- sparse:

  Data frame for the sparse private histogram, or `NULL` if not
  requested.

- method:

  Character vector of private histogram methods used.

## Details

Let \\v_a\\ and \\v_b\\ be the principal component directions selected
by `axes = c(a, b)` for some \\1 \le a \< b \le ncol(X)\\. After
preprocessing, the score point for \\i\\th observation is \\s_i =
(x_i^\top v_a, x_i^\top v_b)\\. A non-private score plot would display
the points \\s_1, \ldots, s_n\\ directly. This function instead
summarizes their empirical distribution by a two-dimensional histogram
and releases private versions of the histogram for the visualization.

If `fixed_frame = NULL`, the plotting frame is constructed privately.
The two score coordinates are stacked into one vector, private lower and
upper quantiles are estimated using a smooth-sensitivity based quantile
mechanism (Nissim et al. 2007) , and the resulting interval is used as a
common square frame for both axes. If `fixed_frame` is supplied, it is
treated as public and no privacy budget is spent on frame construction.

The private histogram is computed on the rectangular grid defined by
`fixed_frame` or by the private frame and the bin counts in `bins`.
Under row-level adjacency, changing one observation can increase one bin
count by one and decrease another by one, giving \\\ell_1\\ sensitivity
at most \\2\\ and \\\ell_2\\ sensitivity at most \\\sqrt{2}\\ for the
count vector.

Two private histogram mechanisms are supported:

- `"add"` constructs an additive differentially private histogram by
  adding Gaussian noise to all bin counts, clipping negative noisy
  counts to zero, and normalizing the result. This additive-noise
  approach is commonly used for private histograms; see Wasserman and
  Zhou (2010) .

- `"sparse"` constructs a sparse differentially private histogram for
  settings where many bins are empty. It perturbs only nonzero empirical
  bin proportions and keeps bins whose noisy values exceed a stability
  threshold, following the stability-based private histogram idea of
  Karwa and Vadhan (2017) .

The privacy parameters are allocated across the privacy-consuming steps.
If `g_dppca = FALSE` and `fixed_frame = NULL`, half of `eps` and `delta`
is used for private frame construction and half for the private
histogram. If `g_dppca = TRUE` and `fixed_frame = NULL`, the parameters
are split equally among private direction estimation, private frame
construction, and private histogram release. If `fixed_frame` is
supplied, the frame step is skipped and the remaining steps split the
privacy parameters equally.

For a detailed procedure and mathematical formulations, refer
<https://yejinjo0220.github.io/dppca/articles/dp_score>.

## References

Dwork C, Roth A (2014). “The Algorithmic Foundations of Differential
Privacy.” *Found. Trends Theor. Comput. Sci.*, **9**(3–4), 211–407. ISSN
1551-305X, [doi:10.1561/0400000042](https://doi.org/10.1561/0400000042)
.

Nissim K, Raskhodnikova S, Smith A (2007). “Smooth Sensitivity and
Sampling in Private Data Analysis.” In *STOC'07: Proceedings of the 39th
Annual ACM Symposium on Theory of Computing*, 75–84. ISBN 9781595936318,
[doi:10.1145/1250790.1250803](https://doi.org/10.1145/1250790.1250803) .

Wasserman L, Zhou S (2010). “A Statistical Framework for Differential
Privacy.” *Journal of the American Statistical Association*,
**105**(489), 375–389.
[doi:10.1198/jasa.2009.tm08651](https://doi.org/10.1198/jasa.2009.tm08651)
.

Karwa V, Vadhan SP (2017). “Finite Sample Differentially Private
Confidence Intervals.” *CoRR*, **abs/1711.03908**. 1711.03908.

Kim M, Jung S (2025). “Robust and Differentially Private Principal
Component Analysis.” *Statistical Analysis and Data Mining: An ASA Data
Science Journal*, **18**(6), e70053.
[doi:10.1002/sam.70053](https://doi.org/10.1002/sam.70053) .

## See also

[`dp_score_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot.md)
for plotting the output of this function.
[`dp_score_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_group.md)
and
[`dp_score_plot_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot_group.md)
for group-wise score histograms.
[`dp_pc_dir()`](https://yejinjo0220.github.io/dppca/reference/dp_pc_dir.md)
for private principal component direction estimation.

## Examples

``` r
data(gau, package = "dppca")

# Use a small subset to keep the example fast.
X <- gau[1:300, ]

# Compute private two-dimensional PCA scores using the additive histogram method.
set.seed(123)
score_gau <- dp_score(
  X,
  eps = 2,
  delta = 1e-3,
  method = "add",
  bins = c(10, 10)
)

head(score_gau$score)
#>          PC1        PC2
#> 1 -1.6418971 -2.9417503
#> 2  2.4192805 -1.9747774
#> 3 -1.5647289  2.4500389
#> 4  1.1818664  0.6632302
#> 5 -0.7668155 -2.7729387
#> 6  1.4701354  3.1142919
head(score_gau$add)
#>         xmin       xmax      ymin      ymax         prob
#> 1 -6.7554436 -5.3600789 -6.755444 -5.360079 0.0000000000
#> 2 -5.3600789 -3.9647142 -6.755444 -5.360079 0.0176699895
#> 3 -3.9647142 -2.5693495 -6.755444 -5.360079 0.0007993045
#> 4 -2.5693495 -1.1739848 -6.755444 -5.360079 0.0014656449
#> 5 -1.1739848  0.2213799 -6.755444 -5.360079 0.0194424961
#> 6  0.2213799  1.6167446 -6.755444 -5.360079 0.0052250857
```
