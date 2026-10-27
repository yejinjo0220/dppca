# Control options for Huber scree estimation

Creates a control list for `method = "huber"` in
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
or the scree configuration of
[`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md).

## Usage

``` r
huber_control(
  k_min_m2,
  k_max_m2,
  m2_frac,
  mu0 = 0,
  eta0 = 1,
  T = NULL,
  M = NULL
)
```

## Arguments

- k_min_m2, k_max_m2:

  Required public integer bounds for the dyadic scale grid, with
  `k_min_m2 < k_max_m2`.

- m2_frac:

  Required fraction in `(0, 1)` of epsilon used for the private scale
  step.

- mu0:

  Public initial mean estimate (`0` by default).

- eta0:

  Positive public gradient step size (`1` by default).

- T:

  Number of gradient iterations; `NULL` (default) uses
  `ceiling(log(m))`, where `m = floor(n / 2)`.

- M:

  Number of scale-estimation blocks; `NULL` (default) uses
  `floor(sqrt(m / 2))`, capped by available input pairs.

## Value

A control list for `method = "huber"`.

## Details

The Huber method privately estimates a robust scale for each component
using a noisy dyadic histogram, then performs projected noisy gradient
updates to estimate the mean of squared pair differences. The scale
proxy sets the Huber threshold; it is not a raw second moment. At least
eight observations are required.

The private scale step uses `m2_frac * epsilon / k` for each component
and no delta. Joint gradient updates use the remaining epsilon and all
delta, calibrated across components and iterations. Each update is
projected to nonnegative values; optional monotone adjustment follows
estimation. Tuning parameters must be public or have separate privacy
accounting.

## References

Yu M, Ren Z, Zhou W (2024). “Gaussian differentially private robust mean
estimation and inference.” *Bernoulli*, **30**(4), 3059–3088.

## See also

[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
for computing differentially private scree estimates using these control
options.
[`plot.dp_scree()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_scree.md)
for plotting scree estimates.

## Examples

``` r
huber_control(k_min_m2 = -10, k_max_m2 = 10, m2_frac = 1 / 4)
#> $mu0
#> [1] 0
#> 
#> $eta0
#> [1] 1
#> 
#> $T
#> NULL
#> 
#> $M
#> NULL
#> 
#> $k_min_m2
#> [1] -10
#> 
#> $k_max_m2
#> [1] 10
#> 
#> $m2_frac
#> [1] 0.25
#> 
```
