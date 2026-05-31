# Launch the dppca Shiny app

Launches an interactive Shiny application for exploring differentially
private PCA visualizations. The app provides a graphical interface for
[`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md)
and
[`dp_score_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot.md),
including private scree plot methods and histogram-based private score
plot methods.

## Usage

``` r
dppca_app(X = NULL, group = NULL)
```

## Arguments

- X:

  Optional numeric matrix or data frame. If supplied, the app opens with
  this data as the initial dataset.

- group:

  Optional group labels. This can be either a vector of length `nrow(X)`
  or a single column name in `X`.

## Value

No return value. This function opens a Shiny application.

## Details

The app can be opened with built-in example datasets or with a
user-supplied dataset. If `X` is supplied, the app starts with `X` as
the initial dataset. If `group` is supplied, the score plot can use the
group labels for coloring. The `group` argument can be either a vector
of length `nrow(X)` or the name of a column in `X`. When `group` is a
column name, that column is used as group labels and is removed from the
PCA feature matrix.

## Examples

``` r
if (interactive()) {
# Launch the app with built-in example datasets.
dppca_app()

# Launch the app with a user-supplied numeric dataset.
data(gau, package = "dppca")
dppca_app(gau)

# Launch the app with group labels stored in a column.
data(gau_g, package = "dppca")
dppca_app(gau_g, group = "group")
}
```
