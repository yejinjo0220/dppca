# Launch the dppca Shiny app

Launches an interactive Shiny application for exploring differentially
private PCA visualizations. The app provides interfaces for
[`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md)
and
[`dp_score_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot.md),
including private scree plot methods and histogram-based private score
plot methods.

If `X` is supplied, the app starts with `X` as the initial dataset. If
`group` is supplied, the score plot can use the group labels for
coloring. The argument `group` can be either a vector of length
`nrow(X)` or the name of a column in `X`. When `group` is a column name,
that column is used as group labels and removed from the PCA feature
matrix.

## Usage

``` r
dppca_app(X = NULL, group = NULL)
```

## Arguments

- X:

  Optional numeric matrix or data frame. If supplied, the app opens with
  this data as the initial data source.

- group:

  Optional group labels. This can be either a vector of length `nrow(X)`
  or a single column name in `X`.

## Value

No return value. This function is called for its side effect of
launching a Shiny application.

## Examples

``` r
if (interactive()) {
  dppca_app()

  data(gau, package = "dppca")
  dppca_app(gau)

  data(gau_g, package = "dppca")
  dppca_app(gau_g, group = "group")
}
```
