# Plot a single-group histogram panel

Internal plotting helper used by
[`dp_score_plot_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot_group.md)
for one group's histogram display.

## Usage

``` r
make_hist_single_dp(df, xlim, ylim, col, title = NULL)
```

## Arguments

- df:

  Histogram data frame with bin coordinates and probabilities.

- xlim, ylim:

  Plotting limits.

- col:

  Fill color.

- title:

  Optional plot title.

## Value

A ggplot object or a patchwork spacer.
