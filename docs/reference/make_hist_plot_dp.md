# Plot a single histogram panel

Internal plotting helper used by
[`dp_score_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot.md)
for a single histogram panel.

## Usage

``` r
make_hist_plot_dp(hist_df, xlim, ylim, color, title = NULL)
```

## Arguments

- hist_df:

  Histogram data frame with bin coordinates and probabilities.

- xlim, ylim:

  Plotting limits.

- color:

  Fill color.

- title:

  Optional plot title.

## Value

A ggplot object.
