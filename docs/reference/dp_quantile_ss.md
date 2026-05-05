# Differentially private quantile via smooth sensitivity

Internal helper that computes a differentially private estimate of a
quantile using the smooth sensitivity framework. This function is
primarily used by
[`dp_frame()`](https://yejinjo0220.github.io/dppca/reference/dp_frame.md)
to construct a differentially private plotting frame for score-based
visualizations.

## Usage

``` r
dp_quantile_ss(x, q, epsilon, delta)
```

## Arguments

- x:

  Numeric vector of observations.

- q:

  Target quantile level in \\(0, 1)\\.

- epsilon:

  Privacy budget \\\varepsilon\\.

- delta:

  Privacy parameter \\\delta\\.

## Value

A numeric scalar giving the noisy quantile estimate.
