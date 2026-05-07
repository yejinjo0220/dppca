# Control parameters for the PMWM DP scree estimator

Creates a control list for `method = "pmwm"` in
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
and
[`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md).

## Usage

``` r
pmwm_control(
  beta = 1.01,
  a = NULL,
  b = NULL,
  trim_const = 10,
  eta = 0.01,
  split_mode = TRUE
)
```

## Arguments

- beta:

  Log-binning base used by the PMWM private quantile estimator.

- a, b:

  Finite support bounds supplied to the PMWM private quantile routine.

- trim_const, eta:

  Practical clipping parameters used by the PMWM estimator.

- split_mode:

  Logical; whether the PMWM estimator splits the sample into
  quantile-estimation and mean-estimation subsets.

## Value

A list of PMWM-estimator tuning parameters.
