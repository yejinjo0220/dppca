# Gaussian mechanism applied to the spherical Kendall matrix

Internal helper that adds Gaussian noise to the spherical Kendall matrix
to construct a noisy symmetric matrix for DP principal component
direction estimation.

## Usage

``` r
mech_tau_sph(X, sig, cpp.option = FALSE)
```

## Arguments

- X:

  Numeric matrix.

- sig:

  Noise standard deviation.

- cpp.option:

  Logical; currently only `FALSE` is supported.

## Value

A noisy symmetric matrix.
