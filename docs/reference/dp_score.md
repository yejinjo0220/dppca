# Differentially private score histograms

Computes two-dimensional PCA scores from the input data and constructs
three histogram versions on the resulting score space:

- the empirical histogram,

- the Gaussian additive DP histogram,

- the sparse Laplace-thresholded DP histogram.

## Usage

``` r
dp_score(
  X,
  center = TRUE,
  standardize = FALSE,
  g_dppca = FALSE,
  cpp.option = FALSE,
  axes = c(1, 2),
  eps_total,
  delta_total,
  eps_ratio = NULL,
  delta_ratio = NULL,
  inflate = 0.2,
  q_frame = NULL,
  frame = NULL,
  m_x = NULL,
  m_y = NULL,
  bin_method = c("J", "W", "none")
)
```

## Arguments

- X:

  Numeric matrix or an object coercible to a matrix with rows
  corresponding to observations and columns corresponding to variables.

- center:

  Logical; whether the variables should be centered before PCA.

- standardize:

  Logical; whether the variables should be scaled to unit variance
  before PCA.

- g_dppca:

  Logical; whether to use a differentially private PCA direction matrix.
  If `FALSE`, the usual sample PCA directions are used.

- cpp.option:

  Logical passed to
  [`mech_tau_sph()`](https://yejinjo0220.github.io/dppca/reference/mech_tau_sph.md)
  when private directions are computed.

- axes:

  Integer vector of length 2 indicating which principal components to
  use. Default is `c(1, 2)`.

- eps_total:

  Total privacy budget.

- delta_total:

  Total privacy parameter.

- eps_ratio:

  Optional privacy-budget split. If `g_dppca = TRUE`, this must have
  length 3 and corresponds to `(PCA, q, hist)`. If `g_dppca = FALSE`,
  this must have length 2 and corresponds to `(q, hist)`.

- delta_ratio:

  Optional delta split, with the same convention as `eps_ratio`.

- inflate:

  Non-negative numeric value controlling frame expansion.

- q_frame:

  Optional quantile level passed to
  [`dp_frame()`](https://yejinjo0220.github.io/dppca/reference/dp_frame.md).

- frame:

  Optional user-specified frame.

- m_x:

  Optional number of bins along the x-axis.

- m_y:

  Optional number of bins along the y-axis.

- bin_method:

  Character string specifying the heuristic used by
  [`recommend_bins()`](https://yejinjo0220.github.io/dppca/reference/recommend_bins.md)
  when `m_x` and/or `m_y` are not supplied. One of `"J"`, `"W"`, or
  `"none"`.

## Value

A list with components:

- `score`: \\n \times 2\\ score matrix.

- `frame`: a list with `xlim` and `ylim`.

- `none`: data frame for the empirical histogram.

- `add`: data frame for the Gaussian additive DP histogram.

- `sparse`: data frame for the sparse DP histogram.
