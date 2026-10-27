# Allocate privacy budgets for shared PCA estimation

Allocate privacy budgets for shared PCA estimation

## Usage

``` r
privacy_control(
  eps = NULL,
  delta = NULL,
  split = c(loading = 1/3, scree = 1/3, score = 1/3),
  loading = NULL,
  scree = NULL,
  score = NULL
)
```

## Arguments

- eps, delta:

  Total privacy parameters. `eps` must be positive and finite, and
  `delta` must lie in `(0, 1)`. Leave both `NULL` when supplying
  absolute component budgets instead.

- split:

  Named positive proportions that sum to one. The same proportions apply
  to `eps` and `delta`. Include `loading` alone, or `loading` with
  `scree`, `score`, or both. When `split` is omitted, the fitting
  function automatically divides the total budget equally among its
  required components: `loading` for
  [`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md),
  `loading` and `scree` for
  [`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md),
  `loading` and `score` for
  [`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md),
  or all three for
  [`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md).
  An explicitly supplied `split` is used exactly as specified, without
  automatic redistribution.

- loading, scree, score:

  Optional absolute component budgets, each a named
  `c(eps = ..., delta = ...)` vector. Supply `loading` alone or with
  `scree`, `score`, or both, without total budgets or a `split`
  argument. Each component must have positive finite `eps` and `delta`
  in `(0, 1)`; the component deltas must sum to less than one. Absolute
  budgets may allocate different proportions of `eps` and `delta`.

## Value

A `dppca_privacy_control` list with `total`, a named `c(eps, delta)`
vector, and `budgets`, a list of the supplied component budget vectors
in `loading`, `scree`, `score` order.

## Details

The estimation function resolves the required component budgets:
`loading` for
[`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md),
`loading` plus `scree` for
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md),
`loading` plus `score` for
[`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md),
and all three for
[`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md).
With total `eps` and `delta` and no explicit `split`, the constructor
initially shows equal thirds, but the fitting function resolves the
actual equal shares. An explicit split or absolute budget must match the
required components.

With supplied `V_dp`, the loading allocation records the prior direction
release; it is not spent again or redistributed. Score epsilon is
divided 35% to the private frame and 65% to the histogram; all score
delta goes to the histogram. When multiple scree or score estimators are
requested, each receives its respective full allocation. Releasing them
together requires privacy composition; the recorded totals are
allocations, not cumulative privacy expenditure. Preprocessing and
supplied directions require their own privacy justification.

## Examples

``` r
# With no split, the fitting function resolves equal shares:
# loading only (100%), scree/score (50%-50%), or PCA (one-third each).
privacy_control(eps = 6, delta = 1e-4)
#> $total
#>   eps delta 
#> 6e+00 1e-04 
#> 
#> $budgets
#> $budgets$loading
#>          eps        delta 
#> 2.000000e+00 3.333333e-05 
#> 
#> $budgets$scree
#>          eps        delta 
#> 2.000000e+00 3.333333e-05 
#> 
#> $budgets$score
#>          eps        delta 
#> 2.000000e+00 3.333333e-05 
#> 
#> 
#> attr(,"class")
#> [1] "dppca_privacy_control"
#> attr(,"implicit_split")
#> [1] TRUE

# Alternatively, specify the loading-only absolute budget.
privacy_control(loading = c(eps = 6, delta = 1e-4))
#> $total
#>   eps delta 
#> 6e+00 1e-04 
#> 
#> $budgets
#> $budgets$loading
#>   eps delta 
#> 6e+00 1e-04 
#> 
#> 
#> attr(,"class")
#> [1] "dppca_privacy_control"

# Explicit splits customize allocation for scree or score.
privacy_control(
  eps = 6,
  delta = 1e-4,
  split = c(loading = 0.5, scree = 0.5)
)
#> $total
#>   eps delta 
#> 6e+00 1e-04 
#> 
#> $budgets
#> $budgets$loading
#>   eps delta 
#> 3e+00 5e-05 
#> 
#> $budgets$scree
#>   eps delta 
#> 3e+00 5e-05 
#> 
#> 
#> attr(,"class")
#> [1] "dppca_privacy_control"

# Score estimation: a different allocation for its two components.
privacy_control(
  eps = 6,
  delta = 1e-4,
  split = c(loading = 0.4, score = 0.6)
)
#> $total
#>   eps delta 
#> 6e+00 1e-04 
#> 
#> $budgets
#> $budgets$loading
#>     eps   delta 
#> 2.40000 0.00004 
#> 
#> $budgets$score
#>     eps   delta 
#> 3.60000 0.00006 
#> 
#> 
#> attr(,"class")
#> [1] "dppca_privacy_control"

# Absolute budgets, without total budgets or split.
privacy_control(
  loading = c(eps = 2, delta = 2e-5),
  scree = c(eps = 4, delta = 8e-5)
)
#> $total
#>   eps delta 
#> 6e+00 1e-04 
#> 
#> $budgets
#> $budgets$loading
#>   eps delta 
#> 2e+00 2e-05 
#> 
#> $budgets$scree
#>   eps delta 
#> 4e+00 8e-05 
#> 
#> 
#> attr(,"class")
#> [1] "dppca_privacy_control"
```
