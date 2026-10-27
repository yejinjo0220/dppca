# Plot stored private PC directions

Create loading plots from an existing
[`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
fit. Plotting uses the stored direction matrices and does not
re-estimate PCs or spend additional privacy budget.

Draw the combined plot stored in a `dp_loading_plot` object. The
original directions remain available in `x$loading`.

Summarize the fitted PC directions stored in a `dp_loading_plot` object.

Display a previously constructed plot or redraw it using new plotting
options without re-estimating the PC directions.

## Usage

``` r
# S3 method for class 'dp_loading'
plot(
  x,
  axes = c(1, 2),
  type = c("vector", "point", "heatmap"),
  heatmap_components = NULL,
  variables = 10,
  plot_control = NULL,
  ...
)

# S3 method for class 'dp_loading_plot'
print(x, ...)

# S3 method for class 'dp_loading_plot'
summary(object, ...)

# S3 method for class 'dp_loading_plot'
plot(x, ...)
```

## Arguments

- x:

  A `dp_loading_plot` object.

- axes:

  Two distinct PC indices for vector and point plots, displayed on fixed
  axes from -1.2 to 1.2. Default: `c(1, 2)`.

- type:

  One or more of `"vector"`, `"point"`, and `"heatmap"`. The default
  shows all three types. With only one stored component, the default
  shows a heatmap because vector and point plots need two axes.

- heatmap_components:

  PC indices to include in a heatmap; `NULL` (default) uses `axes`.

- variables:

  Variables to display. A single positive integer selects the variables
  with the largest Euclidean norms across the displayed PC coefficients
  (top 10 by default). Use `"all"` or `NULL` for every variable, or
  provide variable names or a vector of column indices. Selection is
  performed separately for each version and display.

- plot_control:

  Optional appearance settings from
  [`loading_plot_control()`](https://yejinjo0220.github.io/dppca/reference/loading_plot_control.md)
  or a named list of its options.

- ...:

  Options passed to `plot.dp_loading()`. With no additional options,
  redraw the stored plot using its original settings.

- object:

  A `dp_loading_plot` object.

## Value

Draws the combined plot and invisibly returns a `dp_loading_plot` object
with `loading` (the original fit) and `plot` (individual plots and the
combined `plot$all`). Printing the saved object displays the plot.

Invisibly returns `x`.

A `summary.dp_loading` object.

Invisibly returns a `dp_loading_plot` object.

## Details

Vector and point plots use fixed axes from -1.2 to 1.2 and show variable
labels. Each heatmap scales colors over its own displayed coefficients,
so comparable coefficients can have different colors between versions.
An ordinary PCA reference is shown automatically only when it was stored
at fitting; it is not a private release.
