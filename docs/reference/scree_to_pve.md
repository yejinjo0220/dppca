# Convert scree values to proportions of variance explained safely

Internal helper that converts a vector of scree estimates to proportions
of variance explained (PVE). If the total scree is positive and finite,
this returns the normalized vector `scree / sum(scree)`. If the total
scree is nonpositive or nonfinite, this helper returns a zero vector
instead of propagating problematic values.

## Usage

``` r
scree_to_pve(scree)
```

## Arguments

- scree:

  Numeric vector of scree estimates.

## Value

Numeric vector of the same length as `scree`, interpreted as proportions
of variance explained.
