# Plot stored private scree estimates

Plot stored private scree estimates

## Usage

``` r
# S3 method for class 'dp_scree'
plot(x, type = c("pve", "scree"), plot_control = NULL, ...)
```

## Arguments

- x:

  A `dp_scree` object.

- type:

  One or both quantities to plot: `"pve"` for proportions of variance
  explained and `"scree"` for scree values. The default draws both side
  by side. Supply one value to draw only that quantity.

- plot_control:

  Optional list from
  [`scree_plot_control()`](https://yejinjo0220.github.io/dppca/reference/scree_plot_control.md).
  If `NULL`, default plotting settings are used.

- ...:

  Further arguments, reserved for compatibility with the generic.

## Value

Invisibly returns `x` after drawing the requested plot(s).

## Details

Uses stored estimates only, so repeated plotting consumes no new privacy
budget. The ordinary reference is non-private comparison output and is
displayed automatically when present in the fitted result.
