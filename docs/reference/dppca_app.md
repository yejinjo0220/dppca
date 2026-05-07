# Launch the dppca Shiny app

Opens an interactive Shiny app for exploring
[`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md)
and
[`dp_score_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot.md).
If `X` is supplied, the app opens using that data. If `group` is
supplied, the score plot uses the grouped version when the group labels
are available. If `group = "color"`, the column named `color` is used as
group labels and removed from the PCA feature matrix.

## Usage

``` r
dppca_app(X = NULL, group = NULL)
```

## Arguments

- X:

  Optional numeric matrix or data frame. If supplied, the app opens with
  this data selected as the default data source.

- group:

  Optional group labels. This can be either a vector of length `nrow(X)`
  or a single column name in `X`.

## Value

Invisibly launches a Shiny application.

## Examples

``` r
if (FALSE) { # \dontrun{
dppca_app()
dppca_app(gau)
dppca_app(gau_g, group = "color")
} # }
```
