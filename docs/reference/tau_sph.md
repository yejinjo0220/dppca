# Spherical Kendall matrix

Internal helper that computes the spherical Kendall matrix from pairwise
normalized differences of observations.

## Usage

``` r
tau_sph(X, cpp.option = FALSE)
```

## Arguments

- X:

  Numeric matrix.

- cpp.option:

  Logical; currently only `FALSE` is supported.

## Value

A symmetric matrix.
