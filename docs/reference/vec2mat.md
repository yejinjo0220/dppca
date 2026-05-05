# Convert a packed symmetric vector into a matrix

Internal helper that converts a packed vector representation of a
symmetric matrix into the full matrix form.

## Usage

``` r
vec2mat(v, p)
```

## Arguments

- v:

  Numeric vector of length \\p(p+1)/2\\.

- p:

  Matrix dimension.

## Value

A symmetric \\p \times p\\ matrix.
