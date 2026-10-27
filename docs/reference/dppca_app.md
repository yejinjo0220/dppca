# Launch the dppca Shiny app

Launches an interactive 'shiny' application for exploring differentially
private PCA visualizations. The PCA tab uses
[`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
for integrated estimation, summaries, and scree, score, loading,
variable, and biplot displays, including grouped score plots when
grouping is enabled. The Loading, Scree, and Score tabs provide
standalone analyses, including method comparisons and grouped score
plots.

## Usage

``` r
dppca_app(X = NULL, group = NULL)
```

## Arguments

- X:

  A numeric matrix or data frame with observations in rows and variables
  in columns. Non-numeric data-frame columns are excluded with a
  warning; numeric values must be finite. Input is not automatically
  centered or standardized. If provided, the app uses this dataset
  initially. At least two observations and two numeric features are
  required.

- group:

  Optional group column name in `X` or vector of group labels (`NULL` by
  default). A selected group column is excluded from PCA features, even
  when numeric; group labels and membership are assumed public. Only
  used with supplied `X`; activate grouping in the PCA or Score tab to
  use these labels.

## Value

Invisibly returns the value produced by
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html) after the
application exits.

## Details

The app accepts built-in data, supplied `X`, or uploaded CSV files.
Grouping is available for PCA and Score analyses; a selected group
column is excluded from numeric PCA features. Loading and Scree analyses
use all numeric feature columns.

Running an analysis stores its estimates for subsequent plotting.
Display settings can redraw an existing fit without consuming another
privacy budget. Non-private reference output is optional in the app and
is not a private release. Budget tables display allocations but do not
track cumulative privacy cost across runs.

Input is not automatically centered or standardized. Any data-dependent
preprocessing requires separate privacy accounting; see
[`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
and
[`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
for applicable conditions.

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
