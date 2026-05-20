# Estimate principal component directions

This function returns the principal component directions used by the
plotting and estimation functions in this package. By default, it
computes the usual non-private directions from the sample covariance
matrix. If `g_dppca = TRUE`, it computes private directions using the
spherical-transformation version of g-DPPCA proposed by Kim and Jung
(2025) : it adds Gaussian noise to the spherical Kendall matrix and then
extracts its leading eigenvectors.

## Usage

``` r
dp_pc_dir(
  X,
  k,
  center = TRUE,
  standardize = FALSE,
  g_dppca = FALSE,
  eps_dir = NULL,
  delta_dir = NULL,
  cpp.option = FALSE
)
```

## Arguments

- X:

  A numeric matrix or data frame. Rows correspond to observations and
  columns correspond to variables.

- k:

  Number of principal component directions to return. Must be an integer
  between `1` and the number of columns in `X`.

- center:

  A logical value indicating whether to center the columns of `X` before
  computing principal component directions. The default is `TRUE`.

- standardize:

  A logical value indicating whether to scale the columns of `X` by
  their sample standard deviations after optional centering. The default
  is `FALSE`.

- g_dppca:

  Whether to compute private principal component directions using the
  spherical Kendall mechanism based on the g-DPPCA method. If `FALSE`,
  the usual non-private directions are computed from the sample
  covariance matrix. The default is `FALSE`.

- eps_dir:

  Positive number defining the `epsilon` privacy parameter for private
  principal component directions. Required when `g_dppca = TRUE`.

- delta_dir:

  Number in `(0, 1)` defining the `delta` privacy parameter for private
  principal component directions. Required when `g_dppca = TRUE`.

- cpp.option:

  A logical value reserved for a future C++ implementation of the
  spherical Kendall matrix. Currently only `FALSE` is supported.

## Value

A numeric matrix with `ncol(X)` rows and `k` columns. The columns are
orthonormal principal component directions.

## Details

The non-private option computes leading eigenvectors of the sample
covariance matrix of the preprocessed data. The private option is based
on the spherical Kendall mechanism of Kim and Jung (2025) : it first
forms the spherical Kendall matrix from pairwise normalized differences,
adds symmetric Gaussian noise, and then computes leading eigenvectors.
The final eigenvector matrix is re-orthonormalized by QR decomposition.
For a detailed procedure and mathematical formulations, refer
<https://yejinjo0220.github.io/dppca/articles/pc_direction>.

## References

Kim M, Jung S (2025). “Robust and Differentially Private Principal
Component Analysis.” *Statistical Analysis and Data Mining: An ASA Data
Science Journal*, **18**(6), e70053.
[doi:10.1002/sam.70053](https://doi.org/10.1002/sam.70053) .

## Examples

``` r
data(gau, package = "dppca")

# Use a small subset to keep the example fast.
X <- head(gau, 50)

# Non-private principal component directions
V <- dp_pc_dir(X, k = 2)
head(V)
#>             [,1]        [,2]
#> [1,] -0.04609090  0.36535470
#> [2,]  0.06004680 -0.03050088
#> [3,] -0.13446543  0.25881051
#> [4,]  0.34415959  0.24385316
#> [5,]  0.08629184 -0.16559258
#> [6,]  0.01282587  0.12106961

# Private principal component directions
set.seed(123)
V_private <- dp_pc_dir(
  X,
  k = 2,
  g_dppca = TRUE,
  eps_dir = 1,
  delta_dir = 1e-2
)
head(V_private)
#>             [,1]       [,2]
#> [1,] -0.25576159  0.2022550
#> [2,] -0.03711666 -0.3422293
#> [3,] -0.06825110 -0.1075148
#> [4,]  0.05285452  0.2211351
#> [5,]  0.15364390 -0.1798014
#> [6,]  0.25121749  0.2434837

# Generate a small low-rank dataset.
n <- 50
z1 <- rnorm(n)
z2 <- rnorm(n)
X <- cbind(
  x1 = z1 + 0.2 * rnorm(n),
  x2 = 0.8 * z1 + 0.2 * rnorm(n),
  x3 = z2 + 0.2 * rnorm(n),
  x4 = 0.5 * z1 - 0.4 * z2 + 0.2 * rnorm(n)
)

# Non-private principal component directions
V <- dp_pc_dir(X, k = 2)
head(V)
#>            [,1]       [,2]
#> [1,] -0.6367737 -0.3598352
#> [2,] -0.5229317 -0.2991188
#> [3,]  0.3449065 -0.8616297
#> [4,] -0.4495567  0.1965728

# Private principal component directions
V_private <- dp_pc_dir(
  X,
  k = 2,
  g_dppca = TRUE,
  eps_dir = 1,
  delta_dir = 1e-2
)
head(V_private)
#>              [,1]       [,2]
#> [1,] -0.273649625 -0.6916640
#> [2,]  0.003308571 -0.6553634
#> [3,]  0.876794547 -0.2730895
#> [4,] -0.395393803 -0.1323696
```
