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
X <- head(gau, 50)

set.seed(123)
score_out <- dp_score(
  X,
  eps = 1,
  delta = 1e-2,
  bins = c(3, 3),
  method = "add"
)
head(score_out$none)
#>         xmin      xmax       ymin      ymax prob
#> 1 -13.289336 -2.395967 -13.289336 -2.395967 0.00
#> 2  -2.395967  8.497402 -13.289336 -2.395967 0.06
#> 3   8.497402 19.390771 -13.289336 -2.395967 0.00
#> 4 -13.289336 -2.395967  -2.395967  8.497402 0.08
#> 5  -2.395967  8.497402  -2.395967  8.497402 0.86
#> 6   8.497402 19.390771  -2.395967  8.497402 0.00
head(score_out$add)
#>         xmin      xmax       ymin      ymax        prob
#> 1 -13.289336 -2.395967 -13.289336 -2.395967 0.000000000
#> 2  -2.395967  8.497402 -13.289336 -2.395967 0.202924127
#> 3   8.497402 19.390771 -13.289336 -2.395967 0.007619123
#> 4 -13.289336 -2.395967  -2.395967  8.497402 0.059958028
#> 5  -2.395967  8.497402  -2.395967  8.497402 0.679692208
#> 6   8.497402 19.390771  -2.395967  8.497402 0.049806514

# Simulated low-rank data.
set.seed(123)
n <- 50
z1 <- rnorm(n)
z2 <- rnorm(n)
X_sim <- cbind(
  x1 = z1 + 0.2 * rnorm(n),
  x2 = 0.8 * z1 + 0.2 * rnorm(n),
  x3 = z2 + 0.2 * rnorm(n),
  x4 = 0.5 * z1 - 0.4 * z2 + 0.2 * rnorm(n)
)

set.seed(123)
dp_score(X_sim, eps = 1, delta = 1e-2, bins = c(3, 3))
#> $score
#>               PC1          PC2
#>  [1,]  0.99162675 -0.288962607
#>  [2,]  0.26756662 -0.059616641
#>  [3,] -2.08102800 -0.431050283
#>  [4,]  0.67608389 -1.255480719
#>  [5,] -0.32440649  0.506916368
#>  [6,] -1.57356213 -1.951067190
#>  [7,] -1.43161032  1.663396014
#>  [8,]  2.07318063  0.333144915
#>  [9,]  0.98049866 -0.032670737
#> [10,]  0.70034409  0.004699994
#> [11,] -1.44297007 -0.797121334
#> [12,] -0.64664961  0.421804405
#> [13,] -0.29225529  0.239643959
#> [14,] -1.11234956  1.072458305
#> [15,] -0.07158962  1.664922280
#> [16,] -2.35070173 -1.160022906
#> [17,] -0.57249730 -0.474910038
#> [18,]  2.82708564  0.998227532
#> [19,] -0.50413370 -0.822948711
#> [20,]  1.54523530 -1.397568542
#> [21,]  1.04887684  1.186149747
#> [22,] -0.72377795  2.461344053
#> [23,]  1.84746380 -0.568229874
#> [24,]  0.57899838  0.846442016
#> [25,]  0.21079536  1.075908809
#> [26,]  2.89686968 -0.117704009
#> [27,] -1.34946217  0.182844498
#> [28,] -0.87644481  1.369925298
#> [29,]  1.64361192  0.318372256
#> [30,] -1.74013894 -0.036355284
#> [31,] -0.41510932 -0.509207761
#> [32,]  0.27087800 -0.168180280
#> [33,] -1.36426814  0.119990175
#> [34,] -0.75924397 -0.667098498
#> [35,] -0.94866365  0.260810842
#> [36,] -1.00014259 -0.296616448
#> [37,] -0.20491448 -1.084477894
#> [38,]  0.02061806 -0.324861131
#> [39,] -0.07561976  0.347182899
#> [40,]  1.18208081 -0.572686094
#> [41,]  1.07110005 -0.420927467
#> [42,]  0.44196389 -0.152921070
#> [43,]  1.92627855  0.291747916
#> [44,] -2.87743776  0.177353526
#> [45,] -0.82788099 -1.405900189
#> [46,]  1.13550253  0.774676346
#> [47,]  1.46119153 -1.677920018
#> [48,]  1.36127691 -0.973896445
#> [49,] -1.41984828  0.057224396
#> [50,] -0.17242125  1.273215618
#> 
#> $frame
#> $frame$xlim
#> [1] -10.80660  14.12147
#> 
#> $frame$ylim
#> [1] -10.80660  14.12147
#> 
#> 
#> $none
#>         xmin      xmax       ymin      ymax prob
#> 1 -10.806602 -2.497246 -10.806602 -2.497246 0.00
#> 2  -2.497246  5.812110 -10.806602 -2.497246 0.00
#> 3   5.812110 14.121466 -10.806602 -2.497246 0.00
#> 4 -10.806602 -2.497246  -2.497246  5.812110 0.02
#> 5  -2.497246  5.812110  -2.497246  5.812110 0.98
#> 6   5.812110 14.121466  -2.497246  5.812110 0.00
#> 7 -10.806602 -2.497246   5.812110 14.121466 0.00
#> 8  -2.497246  5.812110   5.812110 14.121466 0.00
#> 9   5.812110 14.121466   5.812110 14.121466 0.00
#> 
#> $add
#>         xmin      xmax       ymin      ymax       prob
#> 1 -10.806602 -2.497246 -10.806602 -2.497246 0.00000000
#> 2  -2.497246  5.812110 -10.806602 -2.497246 0.24197555
#> 3   5.812110 14.121466 -10.806602 -2.497246 0.01094580
#> 4 -10.806602 -2.497246  -2.497246  5.812110 0.02785489
#> 5  -2.497246  5.812110  -2.497246  5.812110 0.64767063
#> 6   5.812110 14.121466  -2.497246  5.812110 0.07155313
#> 7 -10.806602 -2.497246   5.812110 14.121466 0.00000000
#> 8  -2.497246  5.812110   5.812110 14.121466 0.00000000
#> 9   5.812110 14.121466   5.812110 14.121466 0.00000000
#> 
#> $sparse
#>         xmin      xmax       ymin      ymax prob
#> 1 -10.806602 -2.497246 -10.806602 -2.497246    0
#> 2  -2.497246  5.812110 -10.806602 -2.497246    0
#> 3   5.812110 14.121466 -10.806602 -2.497246    0
#> 4 -10.806602 -2.497246  -2.497246  5.812110    0
#> 5  -2.497246  5.812110  -2.497246  5.812110    0
#> 6   5.812110 14.121466  -2.497246  5.812110    0
#> 7 -10.806602 -2.497246   5.812110 14.121466    0
#> 8  -2.497246  5.812110   5.812110 14.121466    0
#> 9   5.812110 14.121466   5.812110 14.121466    0
#> 
#> $method
#> [1] "add"    "sparse"
#> 
```
