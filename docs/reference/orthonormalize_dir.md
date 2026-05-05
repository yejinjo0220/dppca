# Re-orthonormalize principal component directions

Internal helper that re-orthonormalizes candidate principal component
directions using a QR decomposition and keeps the first `k` columns.

## Usage

``` r
orthonormalize_dir(V, k)
```

## Arguments

- V:

  Numeric matrix.

- k:

  Number of columns to retain.

## Value

A direction matrix with orthonormal columns.
