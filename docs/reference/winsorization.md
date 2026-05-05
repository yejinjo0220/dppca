# Clamp numeric values to a fixed interval

Internal helper that truncates numeric values to the closed interval
`[lo, hi]`. This is used in several places where intermediate values
must be kept inside valid numerical or algorithmic bounds.

## Usage

``` r
winsorization(x, lo, hi)
```

## Arguments

- x:

  Numeric vector.

- lo:

  Finite lower bound.

- hi:

  Finite upper bound satisfying `lo <= hi`.

## Value

Numeric vector of the same length as `x`, with all entries truncated to
`[lo, hi]`.
