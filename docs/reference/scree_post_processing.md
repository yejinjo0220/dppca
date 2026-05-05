# Enforce a nonnegative and nonincreasing scree sequence

Internal post-processing helper for scree estimates. Since scree values
are variances, they should be nonnegative, and in many workflows it is
desirable to enforce a monotone decreasing profile across components.
This function first truncates the input at zero and then applies
isotonic regression to the negated sequence so that the final output is
nonincreasing.

## Usage

``` r
scree_post_processing(x)
```

## Arguments

- x:

  Numeric vector of raw scree estimates.

## Value

Numeric vector of the same length as `x`, containing a nonnegative and
nonincreasing version of the input.
