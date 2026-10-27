# Differentially private scree values

Estimates scree values along private principal component directions
using means of squared differences from disjoint score pairs. An
ordinary PCA reference is computed only when `non_private = TRUE`.

## Usage

``` r
dp_scree(
  X,
  k,
  privacy,
  method = c("clipped", "pmwm", "huber"),
  control = NULL,
  cpp.option = TRUE,
  mono = TRUE,
  V_dp = NULL,
  non_private = FALSE
)

# S3 method for class 'dp_scree'
print(x, ...)

# S3 method for class 'dp_scree'
summary(object, ...)

# S3 method for class 'summary.dp_scree'
print(x, digits = 4, ...)
```

## Arguments

- X:

  A numeric matrix or data frame with observations in rows and variables
  in columns. Non-numeric data-frame columns are excluded with a
  warning; numeric values must be finite. At least two observations and
  one numeric feature are required. If directions must be estimated, the
  pairwise-separation condition of
  [`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
  also applies.

- k:

  Number of leading components to estimate; an integer from `1` to the
  number of retained numeric features.

- privacy:

  A
  [`privacy_control()`](https://yejinjo0220.github.io/dppca/reference/privacy_control.md)
  object allocating budget to `loading` and `scree` (equal shares by
  default). Each requested scree method receives the full scree
  allocation.

- method:

  One or more scree estimators: `"clipped"`, `"pmwm"`, or `"huber"`
  (`"clipped"` by default).

- control:

  Tuning parameters from
  [`clipped_control()`](https://yejinjo0220.github.io/dppca/reference/clipped_control.md),
  [`pmwm_control()`](https://yejinjo0220.github.io/dppca/reference/pmwm_control.md),
  or
  [`huber_control()`](https://yejinjo0220.github.io/dppca/reference/huber_control.md).
  For multiple methods, supply a named list of method-specific controls;
  required tuning values must be supplied.

- cpp.option:

  Whether to use the compiled spherical-Kendall implementation (`TRUE`
  by default) when estimating private directions; ignored if `V_dp` is
  supplied.

- mono:

  Whether to apply nonincreasing isotonic regression (`TRUE` by
  default). Negative scree estimates are truncated to zero in either
  case.

- V_dp:

  Optional `p` by `p` private PC direction matrix (`NULL` by default),
  where `p` is the number of numeric PCA features. Columns must be
  orthonormal and rows must follow the feature order in `X`; matching
  row names are required when both sets of names exist.

- non_private:

  Whether to compute and store ordinary PCA reference results (`FALSE`
  by default). The reference contains ordinary PCA scree and PVE.

- x:

  A `summary.dp_scree` object.

- ...:

  Further arguments, reserved for compatibility with the generic.

- object:

  A `dp_scree` object.

- digits:

  Number of significant digits used for printing.

## Value

A `dp_scree` list containing one entry per requested method (`clipped`,
`pmwm`, and/or `huber`). Each entry has `scree` and `pve` vectors named
`PC1`, ..., `PCk`. The `nonprivate` entry holds the optional ordinary
PCA reference (or `NULL` when omitted), and `metadata` stores input
dimensions, methods and the recorded privacy budget. PVE is normalized
over the `k` returned components in every entry. Use
[`summary()`](https://rdrr.io/r/base/summary.html) for scree, PVE and
cumulative PVE tables, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) to draw stored
estimates without a new randomized computation.

Invisibly returns `x`.

A `summary.dp_scree` object containing one table per stored method. Each
table gives scree values, PVE and cumulative PVE. PVE is normalized over
the fitted components, as in `dp_scree()`. Printing the summary shows
component importance by method and the privacy allocations.

## Details

Let `Y = X V` be scores along the private directions. For each method,
`m = floor(n / 2)` disjoint pairs of rows are formed and the
component-wise values are \$\$W\_{jl} = (Y\_{b_j,l} -
Y\_{a_j,l})^2/2.\$\$ An odd final observation is unused. The methods
privately estimate the column means of `W`, corresponding to variance
along private directions rather than necessarily to sample covariance
eigenvalues.

The estimators differ as follows:

- `"clipped"`: Clip the paired values at `C_clip`, then perturb their
  mean vector with Gaussian noise.

- `"pmwm"`: Estimate private lower and upper quantiles, winsorize the
  paired values, and perturb their mean vector with Gaussian noise.
  Quantile calls share half of the method's epsilon; the mean receives
  the remaining half and all delta.

- `"huber"`: Estimate private component scales and apply noisy Huber
  gradient updates. The scale step uses `m2_frac` of epsilon; the
  updates use the remainder and all delta. At least eight observations
  are required.

Directions are estimated once unless `V_dp` is supplied. The input is
not automatically centered or standardized; sensitivity arguments use
fixed-size replacement adjacency and condition on the directions.
Data-dependent scaling requires separate privacy analysis. Sample size,
schema, component selection, and tuning parameters are assumed public.

Each requested scree method uses the full scree allocation. Releasing
several methods together, or releasing results from repeated calls,
requires composition beyond the stored budget for a single method. The
optional ordinary PCA reference is not private.

Scree values are truncated at zero, optionally made nonincreasing, and
normalized over the returned `k` components for PVE. A zero total yields
zero PVE. These steps are post-processing and do not reorder the
directions.

## References

Dwork C, Roth A (2014). “The Algorithmic Foundations of Differential
Privacy.” *Found. Trends Theor. Comput. Sci.*, **9**(3–4), 211–407. ISSN
1551-305X, [doi:10.1561/0400000042](https://doi.org/10.1561/0400000042)
.

Ramsay K, Spicker D (2025). “Improved subsample-and-aggregate via the
private modified winsorized mean.” Code available at
<https://github.com/12ramsake/PMWM>, 2501.14095,
<https://arxiv.org/abs/2501.14095>.

Durfee D (2023). “Unbounded Differentially Private Quantile and Maximum
Estimation.” In *Advances in Neural Information Processing Systems*,
volume 36, 77691–77712.
[doi:10.52202/075280-3394](https://doi.org/10.52202/075280-3394) ,
<https://proceedings.neurips.cc/paper_files/paper/2023/hash/f4b6ef2a78684dca2fb3f1c09372e041-Abstract-Conference.html>.

Yu M, Ren Z, Zhou W (2024). “Gaussian differentially private robust mean
estimation and inference.” *Bernoulli*, **30**(4), 3059–3088.

Kim M, Jung S (2025). “Robust and Differentially Private Principal
Component Analysis.” *Statistical Analysis and Data Mining: An ASA Data
Science Journal*, **18**(6), e70053.
[doi:10.1002/sam.70053](https://doi.org/10.1002/sam.70053) .

## See also

[`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
for private principal component direction estimation.
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
privacy <- privacy_control(eps = 3, delta = 1e-4)

# Estimate clipped private scree values and an ordinary PCA comparison.
set.seed(123)
out <- dp_scree(
  X,
  k = 2,
  privacy = privacy,
  method = "clipped",
  control = clipped_control(C_clip = 3)
)

out$nonprivate
#> NULL
out$clipped
#> $scree
#>       PC1       PC2 
#> 0.7621344 0.7621344 
#> 
#> $pve
#> PC1 PC2 
#> 0.5 0.5 
#> 

# Multiple methods can be requested together by using a named control list.
# Loading is estimated once; each method receives the full scree budget.
# Change the split in privacy_control() to choose a different allocation.
```
