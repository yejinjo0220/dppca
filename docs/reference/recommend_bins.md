# Recommend the number of histogram bins

Internal helper that recommends the number of bins per axis for a
two-dimensional histogram based on the sample size.

## Usage

``` r
recommend_bins(X, method = c("WZ", "Lei"))
```

## Arguments

- X:

  Matrix or data frame. Only the number of rows is used.

- method:

  Character string specifying the binning rule. Supported options are
  `"WZ"` and `"Lei"`.

## Value

A positive integer giving the recommended number of bins per axis.
