# Control options for Huber scree estimation

Creates a control list for the Huber-type private scree estimator used
by
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
and
[`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md)
when `method = "huber"`.

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

  Integers defining the lower and upper dyadic bin indices used in the
  private second-moment scale step. The histogram searches over scale
  levels \\2^k\\ for \\k\_{\min} \le k \le k\_{\max}\\. These values
  have no defaults because they should be chosen according to the scale
  of the data.

- m2_frac:

  Number in `(0, 1)` defining the fraction of the Huber scree privacy
  parameters allocated to the private second-moment scale step. This
  value has no default.

- mu0:

  Numeric initial value for Huber noisy gradient descent. The default is
  `0`.

- eta0:

  Positive number defining the fixed step size for Huber noisy gradient
  descent. The default is `1`.

- T:

  Optional positive integer defining the number of noisy gradient
  descent iterations. If `NULL`, the implementation uses \\\lceil \log n
  \rceil\\, where \\n\\ is the number of observations.

- M:

  Optional positive integer defining the number of blocks used in the
  private second-moment scale step. If `NULL`, the implementation uses
  \\\lfloor \sqrt{n} / 2 \rfloor\\, where \\n\\ is the number of
  observations.

## Value

A list of control options for `method = "huber"`.

## Details

The Huber method estimates the mean of squared principal component
scores by noisy gradient descent on the Huber loss. It follows the
Huber-type private robust mean approach of Yu et al. (2024) .

The method first privately estimates a scale proxy for the squared
scores, denoted by \\m_2\\. This scale proxy is then used to choose the
Huber robustification level for noisy gradient descent. The parameters
`k_min_m2`, `k_max_m2`, and `m2_frac` control this private scale-proxy
step, while `mu0`, `eta0`, `T`, and `M` control the subsequent noisy
gradient descent routine.

The dyadic indices `k_min_m2` and `k_max_m2` define the search range for
the private histogram used to estimate \\m_2\\. The histogram searches
over candidate scale levels \\2^k\\ satisfying \\k\_{\min} \le k \le
k\_{\max}\\. Because this range depends on the scale of the squared
scores, these arguments are intentionally not given defaults.

The argument `m2_frac` determines how the Huber scree privacy parameters
are split between the private scale-proxy step and the noisy gradient
descent step. If \\(\epsilon\_{\mathrm{scree}},
\delta\_{\mathrm{scree}})\\ denotes the privacy parameters available for
Huber scree estimation, then \\m2_frac \cdot
(\epsilon\_{\mathrm{scree}}, \delta\_{\mathrm{scree}})\\ is used to
privately estimate \\m_2\\, while \\(1 - m2_frac) \cdot
(\epsilon\_{\mathrm{scree}}, \delta\_{\mathrm{scree}})\\ is used for
Huber noisy gradient descent.

The remaining parameters have default values. The default `mu0 = 0` is
the initial value for noisy gradient descent, and the default `eta0 = 1`
is the fixed step size. If `T = NULL`, the number of noisy gradient
descent iterations is chosen as \\\lceil \log n \rceil\\. If `M = NULL`,
the number of blocks used in the private estimator of \\m_2\\ is chosen
as \\\lfloor \sqrt{n} / 2 \rfloor\\.

## References

Yu M, Ren Z, Zhou W (2024). “Gaussian differentially private robust mean
estimation and inference.” *Bernoulli*, **30**(4), 3059–3088.

## See also

[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
for computing differentially private scree estimates using these control
options.
[`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md)
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
