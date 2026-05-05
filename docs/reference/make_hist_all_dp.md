# Plot all-group histogram panel

Internal plotting helper used by
[`dp_score_plot_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot_group.md)
for pooled group-wise histogram displays.

## Usage

``` r
make_hist_all_dp(df, xlim, ylim, col_map, title = NULL)
```

## Arguments

- df:

  Histogram data frame containing a `group` column.

- xlim, ylim:

  Plotting limits.

- col_map:

  Named color vector.

- title:

  Optional plot title.

## Value

A ggplot object or a patchwork spacer.
