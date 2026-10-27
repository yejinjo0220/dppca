# Control loading plot appearance

Specify colors, text sizes, labels, and heatmap settings for loading
plots. Pass the resulting list through the `plot_control` argument of
[`plot.dp_loading()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_loading.md).

## Usage

``` r
loading_plot_control(
  nonprivate_color = "black",
  private_color = "#1F4E79",
  arrow_size = 0.7,
  point_size = 2.8,
  label_size = 4,
  heatmap_values = TRUE,
  heatmap_digits = 2L,
  heatmap_text_size = 3.5,
  base_size = 12,
  title_size = 12,
  nonprivate_title = "Non-private",
  private_title = "DP",
  xlab = NULL,
  ylab = NULL,
  heatmap_colors = c(low = "#2C5D8A", high = "#B44A4A")
)
```

## Arguments

- nonprivate_color, private_color:

  Colors of the non-private and private vectors or points, respectively.

- arrow_size:

  Width of arrows in vector plots.

- point_size:

  Size of points in point plots.

- label_size:

  Size of variable labels in vector and point plots.

- heatmap_values:

  Whether to display loading coefficients in heatmap cells. Default:
  `TRUE`.

- heatmap_digits:

  Number of decimal places for heatmap labels.

- heatmap_text_size:

  Size of coefficient labels in heatmaps.

- base_size:

  Base font size for all loading plots.

- title_size:

  Font size for panel titles.

- nonprivate_title, private_title:

  Text used for the corresponding panel titles. The display name is
  appended in parentheses; `NULL` hides a title.

- xlab, ylab:

  Axis labels for vector and point plots. `NULL` uses the selected PC
  names.

- heatmap_colors:

  Named character vector with `low` and `high` colors for negative and
  positive loading coefficients. White represents zero. Private and
  non-private heatmaps use separate, symmetric color ranges based on
  their own displayed coefficients.

## Value

A named list of validated loading plot settings.
