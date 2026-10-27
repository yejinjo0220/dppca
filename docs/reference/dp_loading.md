# Estimate differentially private PC directions

Estimate principal component (PC) directions using a noisy spherical
Kendall matrix, with optional ordinary PCA directions for comparison.

Display a concise overview of a `dp_loading()` fit, including its data
dimensions, privacy parameters, and whether a non-private reference is
stored.

Summarize a fitted `dp_loading()` object without recomputing its
estimates or consuming additional privacy budget.

## Usage

``` r
dp_loading(X, privacy, cpp.option = TRUE, non_private = FALSE)

# S3 method for class 'dp_loading'
print(x, ...)

# S3 method for class 'dp_loading'
summary(object, ...)

# S3 method for class 'summary.dp_loading'
print(x, digits = 3, ...)
```

## Arguments

- X:

  A numeric matrix or data frame with observations in rows and variables
  in columns. Non-numeric data-frame columns are excluded with a
  warning; numeric values must be finite. At least two observations and
  one numeric feature are required. Observations must be sufficiently
  separated to compute normalized pairwise differences.

- privacy:

  A
  [`privacy_control()`](https://yejinjo0220.github.io/dppca/reference/privacy_control.md)
  object allocating privacy budget to `loading` alone. Either specify
  total `eps` and `delta` or an absolute `loading` budget.

- cpp.option:

  Whether to use the compiled spherical-Kendall implementation (`TRUE`
  by default); `FALSE` uses the R implementation.

- non_private:

  Whether to compute and store ordinary PCA reference results (`FALSE`
  by default).

- x:

  A `summary.dp_loading` object.

- ...:

  Additional arguments.

- object:

  A `dp_loading` object.

- digits:

  Number of significant digits used for printing.

## Value

A `dp_loading` object with the following components:

- `directions`:

  A `p` by `p` matrix of private PC directions.

- `nonprivate`:

  A `p` by `p` matrix of ordinary PC directions, or `NULL` when
  `non_private = FALSE`.

- `metadata`:

  Sample size, variable count and privacy parameters.

Rows represent variables, and columns represent `PC1`, ..., `PCp`. The
columns are orthonormal eigenvectors, without eigenvalue or correlation
scaling.

Invisibly returns `x`.

A `summary.dp_loading` object containing the stored `private` and
`nonprivate` direction matrices, with sample size and privacy budget
metadata. The print method displays up to six variables and four PCs;
the full matrices remain stored. The non-private matrix remains `NULL`
if it was not stored.

Invisibly returns `x`.

## Details

The private PC directions are eigenvectors of a spherical Kendall matrix
perturbed with symmetric Gaussian noise. Directions are ordered by the
noisy eigenvalues; they need not be eigenvectors of the sample
covariance matrix. The spherical Kendall construction uses normalized
pairwise differences and is invariant to common translations.

With `non_private = TRUE`, ordinary covariance-based directions are
stored after sign alignment with the private directions. This comparison
is not a private release.

Under fixed-size replacement adjacency, the matrix has Frobenius
sensitivity at most `4 / n`. The diagonal Gaussian noise standard
deviation is `4 * sqrt(2 * log(1.25 / delta)) / (n * eps)`; off-diagonal
entries use this scale divided by `sqrt(2)`. This classical calibration
has a standard sufficient guarantee for `0 < eps < 1`; larger `eps`
requires an additional calibration check, which is not performed by the
function.

Duplicate or numerically indistinguishable observations are rejected,
but these checks do not privatize input validation. Data schema and
sample size are assumed public. Data-dependent preprocessing requires
separate privacy analysis and accounting.

## Examples

``` r
set.seed(123)
X <- matrix(rnorm(150), nrow = 50, ncol = 3)
fit <- dp_loading(
  X, privacy = privacy_control(eps = 0.5, delta = 1e-3),
  cpp.option = FALSE
)
print(fit)
#> Differentially Private PC Directions
#> Data: 50 observations x 3 variables
#> Components: 3 PC directions
#> Total privacy budget: eps = 0.5 ; delta = 0.001 
#> Non-private reference: not stored 
summary(fit)
#> DP PC Directions Summary
#> Data: 50 observations x 3 variables
#> Components: 3 PC directions
#> 
#> Private PC directions  (first 3 variables x first 3 PCs):
#>       PC1    PC2    PC3
#> X1  0.701 -0.409 -0.584
#> X2 -0.703 -0.259 -0.663
#> X3  0.120  0.875 -0.469
#> 
#> Total privacy budget: eps = 0.5 ; delta = 0.001 
#> Non-private reference: not stored 
#> The full direction matrix is available in $directions.
fit$directions
#>           PC1        PC2        PC3
#> X1  0.7012424 -0.4092929 -0.5837280
#> X2 -0.7027776 -0.2592186 -0.6625024
#> X3  0.1198444  0.8748057 -0.4694170
fit_with_reference <- dp_loading(
  X, privacy = privacy_control(loading = c(eps = 0.5, delta = 1e-3)),
  cpp.option = FALSE, non_private = TRUE
)
fit_with_reference$nonprivate
#>           PC1          PC2       PC3
#> X1  0.1793378  0.977818990 0.1082034
#> X2 -0.5006692 -0.003963068 0.8656297
#> X3  0.8468579 -0.209414222 0.4888531
if (FALSE) { # \dontrun{
loading_plot <- plot(fit, type = c("vector", "point", "heatmap"),
                     heatmap_components = 1:3, variables = "all")
loading_plot                 # Display the saved plot
loading_plot$loading$directions # Access the stored directions
loading_plot$plot$all       # Access the combined plot
} # }
```
