# Launch the dppca Shiny app

Opens an interactive Shiny app for exploring
[`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md)
and
[`dp_score_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot.md)
with the updated scree and score method options. If `X` is supplied, the
app opens using that data. If `group` is supplied, the score plot uses
the grouped version when the group labels are available. If
`group = "group"`, the column named `group` is used as group labels and
removed from the PCA feature matrix. In the app, the scree menu uses
checkboxes and can overlay `method = "clipped"`, `method = "pmwm"`,
`method = "huber"`, or a vector such as
`method = c("clipped", "pmwm", "huber")`. The score menu can pass
`method = "add"`, `method = "sparse"`, or `method = c("add", "sparse")`
and bin counts through `bins = c(m_x, m_y)` to the score plotting
functions. The app uses a taller, centered scree plot panel so that the
scree curves are not visually flattened on wide screens.

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
if (interactive()) {
  dppca_app()

  data(gau, package = "dppca")
  dppca_app(gau)

  data(gau_g, package = "dppca")
  dppca_app(gau_g, group = "group")
}
```
