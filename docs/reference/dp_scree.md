# Differentially private scree values

This function computes estimates of scree values, eigenvalues of
covariance matrix, for principal component analysis, including both the
usual non-private estimates and differentially private estimates. The
private estimates are computed as the private mean of the squared
principal component scores. See Details for the estimating equations and
method-specific construction.

## Usage

``` r
dp_scree(
  X,
  k,
  method = c("clipped", "pmwm", "huber"),
  control = NULL,
  eps,
  delta,
  center = TRUE,
  standardize = FALSE,
  g_dppca = FALSE,
  cpp.option = FALSE,
  mono = TRUE
)
```

## Arguments

- X:

  A numeric matrix or data frame. Rows correspond to observations and
  columns correspond to variables.

- k:

  Positive integer defining the number of leading principal components
  to estimate. Must be an integer between `1` and the number of columns
  in `X`.

- method:

  Scree value estimation method. One of `"clipped"`, `"pmwm"`, or
  `"huber"`.

- control:

  Optional method-specific control list created by
  [`clipped_control()`](https://yejinjo0220.github.io/dppca/reference/clipped_control.md),
  [`pmwm_control()`](https://yejinjo0220.github.io/dppca/reference/pmwm_control.md),
  or
  [`huber_control()`](https://yejinjo0220.github.io/dppca/reference/huber_control.md).

- eps:

  Positive number defining the total `epsilon` privacy parameter. If
  `g_dppca = TRUE`, it is split between private direction estimation and
  private scree estimation.

- delta:

  Number in `(0, 1)` defining the total `delta` privacy parameter. If
  `g_dppca = TRUE`, it is split between private direction estimation and
  private scree estimation.

- center:

  A logical value indicating whether to center the columns of `X` before
  computing principal component directions. The default is `TRUE`.

- standardize:

  A logical value indicating whether to scale the columns of `X` by
  their sample standard deviations after optional centering. The default
  is `FALSE`.

- g_dppca:

  A logical value indicating whether to use private principal component
  directions for scree estimation. The default is `FALSE`. See
  [`dp_pc_dir()`](https://yejinjo0220.github.io/dppca/reference/dp_pc_dir.md)
  for details.

- cpp.option:

  A logical value passed to
  [`dp_pc_dir()`](https://yejinjo0220.github.io/dppca/reference/dp_pc_dir.md)
  when `g_dppca = TRUE`. The default is `FALSE`.

- mono:

  A logical value indicating whether to apply monotone post-processing
  to the vector of private scree values. The default is `TRUE`.

## Value

A list with components:

- `method`: scree value estimation method.

- `scree_np`: non-private scree estimates.

- `pve_np`: non-private proportions of variance explained.

- `scree`: differentially private scree value estimates.

- `pve`: differentially private proportions of variance explained.

## Details

Let \\X\\ denote the preprocessed data matrix and let \\v_l\\ be the
\\l\\th principal component direction. The \\l\\th score vector is \\z_l
= X v_l\\. The corresponding sample scree value can be written as \$\$
\hat{\lambda}\_l = v_l^\top \widehat{\Sigma} v_l = \frac{1}{n -
1}\sum\_{i = 1}^n z\_{il}^2 = \frac{n}{n - 1}\left(\frac{1}{n}\sum\_{i =
1}^n w\_{il}\right), \qquad w\_{il} = z\_{il}^2. \$\$ Therefore, each
scree value is estimated by privately estimating the mean of \\w\_{1l},
\ldots, w\_{nl}\\ and multiplying by \\n/(n - 1)\\.

The supported methods differ in how this private mean is estimated:

- `"clipped"` clips the squared scores \\w\_{i\ell}\\ at `C_clip` and
  then applies the Gaussian mechanism (Dwork and Roth 2014) . This is
  the simplest option but depends directly on the clipping threshold.

- `"pmwm"` uses the private modified winsorized mean approach of Ramsay
  and Spicker (2025) , adapted from the accompanying Python
  implementation into R. It privately estimates tail cutoffs, winsorizes
  the squared scores \\w\_{i\ell}\\, and releases a noisy winsorized
  mean.

- `"huber"` uses a Huber-type private robust mean estimator based on
  noisy gradient descent, following Yu et al. (2024) .

The argument `g_dppca` controls how the principal component directions
are obtained. If `g_dppca = FALSE`, the directions are computed
non-privately as an eigenvector of sample covariance, and the full
privacy parameters `eps` and `delta` are used for private scree value
estimation. If `g_dppca = TRUE`, the directions are computed privately
using
[`dp_pc_dir()`](https://yejinjo0220.github.io/dppca/reference/dp_pc_dir.md).
In that case, the privacy parameters are split equally: `eps / 2` and
`delta / 2` are used for private direction estimation, and the remaining
`eps / 2` and `delta / 2` are used for private scree value estimation.
When `mono = TRUE`, the final monotone adjustment is a post-processing
step and does not change the privacy guarantee.

For a detailed procedure and mathematical formulations, refer
<https://yejinjo0220.github.io/dppca/articles/dp_scree>.

## References

Dwork C, Roth A (2014). “The Algorithmic Foundations of Differential
Privacy.” *Found. Trends Theor. Comput. Sci.*, **9**(3–4), 211–407. ISSN
1551-305X, [doi:10.1561/0400000042](https://doi.org/10.1561/0400000042)
.

Ramsay K, Spicker D (2025). “Improved subsample-and-aggregate via the
private modified winsorized mean.” Code available at
<https://github.com/12ramsake/PMWM>, 2501.14095,
<https://arxiv.org/abs/2501.14095>.

Yu M, Ren Z, Zhou W (2024). “Gaussian differentially private robust mean
estimation and inference.” *Bernoulli*, **30**(4), 3059–3088.

Kim M, Jung S (2025). “Robust and Differentially Private Principal
Component Analysis.” *Statistical Analysis and Data Mining: An ASA Data
Science Journal*, **18**(6), e70053.
[doi:10.1002/sam.70053](https://doi.org/10.1002/sam.70053) .

## See also

[`dp_pc_dir()`](https://yejinjo0220.github.io/dppca/reference/dp_pc_dir.md)
for principal component direction estimation.
[`clipped_control()`](https://yejinjo0220.github.io/dppca/reference/clipped_control.md),
[`pmwm_control()`](https://yejinjo0220.github.io/dppca/reference/pmwm_control.md),
and
[`huber_control()`](https://yejinjo0220.github.io/dppca/reference/huber_control.md)
for method-specific tuning parameters.

## Examples

``` r
data(gau, package = "dppca")

# Use a small subset to keep the example fast.
X <- gau[1:100, ]

# Estimate the private scree values using the clipped mean method.
set.seed(123)
dp_scree(
  X,
  k = 2,
  method = "clipped",
  control = clipped_control(C_clip = 3),
  eps = 2,
  delta = 1e-3
)
#> $method
#> [1] "clipped"
#> 
#> $scree_np
#> [1] 2.061105 1.771498
#> 
#> $pve_np
#> [1] 0.5377819 0.4622181
#> 
#> $scree
#> [1] 1.196163 1.196163
#> 
#> $pve
#> [1] 0.5 0.5
#> 

# Other scree methods can be used by changing `method` and `control`, e.g.,
# method = "pmwm",
# control = pmwm_control(a = 0, b = 50, trim_const = 10, eta = 0.01)
#
# method = "huber",
# control = huber_control(k_min_m2 = -10, k_max_m2 = 10, m2_frac = 1 / 4)
```
