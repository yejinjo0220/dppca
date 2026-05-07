# Plot differentially private score histograms

Computes PCA score histograms through
[`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
and returns both the calculation results and the corresponding plots.

## Usage

``` r
dp_score_plot(
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
  bin_method = c("WZ", "Lei", "none"),
  color = "#6A5ACD"
)
```

## Arguments

- X:

  Numeric matrix or an object coercible to a matrix.

- center:

  Logical; whether the variables should be centered before PCA.

- standardize:

  Logical; whether the variables should be scaled before PCA.

- g_dppca:

  Logical; whether to use a DP PCA direction matrix.

- cpp.option:

  Logical passed to
  [`mech_tau_sph()`](https://yejinjo0220.github.io/dppca/reference/mech_tau_sph.md)
  when private directions are computed.

- axes:

  Integer vector of length 2 specifying the principal components used
  for score construction.

- eps_total:

  Total privacy budget.

- delta_total:

  Total privacy parameter.

- eps_ratio:

  Optional privacy-budget split.

- delta_ratio:

  Optional delta split.

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

  Character string specifying the bin recommendation rule.

- color:

  Base plotting color.

## Value

A list with components:

- `score`: the output of
  [`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md).

- `plot`: a list containing the combined plot and each individual panel.
