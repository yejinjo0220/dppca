# Control parameters for the clipped DP scree estimator

Creates a control list for `method = "clipped"` in
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
and
[`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md).

## Usage

``` r
clipped_control(C_clip = 3)
```

## Arguments

- C_clip:

  Positive clipping threshold used by the clipped estimator.

## Value

A list of clipped-estimator tuning parameters.
