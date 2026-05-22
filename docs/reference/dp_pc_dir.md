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
  eps = NULL,
  delta = NULL,
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

- eps:

  Positive number defining the `epsilon` privacy parameter for private
  principal component directions. Required when `g_dppca = TRUE`.

- delta:

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
X <- gau[1:200, ]

# Non-private principal component directions
V <- dp_pc_dir(X, k = 2)
head(V)
#>              [,1]        [,2]
#> [1,] -0.276390432 -0.05512405
#> [2,] -0.008155916 -0.03392589
#> [3,] -0.433823614 -0.26807870
#> [4,] -0.007115806  0.43519507
#> [5,]  0.406349436  0.12107397
#> [6,]  0.108264920  0.07128666


# Private principal component directions
set.seed(123)
V_private <- dp_pc_dir(
  X,
  k = 2,
  g_dppca = TRUE,
  eps = 2,
  delta = 1e-3
)
head(V_private)
#>              [,1]       [,2]
#> [1,] -0.223077596  0.3476042
#> [2,] -0.050923111 -0.2384854
#> [3,] -0.018308885  0.1395241
#> [4,]  0.091084919  0.1563828
#> [5,]  0.009964607 -0.2978056
#> [6,]  0.161973207  0.1422014
```
