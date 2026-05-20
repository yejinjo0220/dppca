# Control options for private modified winsorized scree estimation

Creates a control list for the private modified winsorized mean scree
estimator used by
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
and
[`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md)
when `method = "pmwm"`.

## Usage

``` r
pmwm_control(a, b, trim_const, eta, beta = 1.001, split_mode = TRUE)
```

## Arguments

- a, b:

  Finite lower and upper search bounds supplied to the private quantile
  routine. The private lower and upper clipping cutoffs are searched
  within this range. These values have no defaults because they should
  be chosen on the scale of squared principal component scores.

- trim_const:

  Positive number controlling the baseline clipping level in the
  practical clipping proportion. This value has no default.

- eta:

  Nonnegative number controlling the expected contamination level in the
  practical clipping proportion. This value has no default.

- beta:

  Positive number greater than `1` defining the log-binning base used by
  the private quantile routine. The default is `1.001`.

- split_mode:

  A logical value indicating whether to split the sample into
  quantile-estimation and mean-estimation subsets. The default is
  `TRUE`.

## Value

A list of control options for `method = "pmwm"`.

## Details

The PMWM method privately estimates lower and upper tail cutoffs,
winsorizes the squared scores to those cutoffs, and then releases a
noisy winsorized mean. It is based on the private modified winsorized
mean of Ramsay and Spicker (2025) .

The implementation used here is an R adaptation of the publicly
available Python implementation accompanying Ramsay and Spicker (2025) .
The adaptation is used for scree estimation by applying the PMWM
estimator to squared principal component scores.

The PMWM scree estimator uses additional control parameters for private
quantile estimation and winsorization. The parameter `beta` determines
the spacing of the geometric search grid used by the private quantile
estimator and must satisfy \\\beta \> 1\\. Smaller values of `beta` give
a finer grid but may increase computation.

The bounds `a` and `b` define the lower and upper search range supplied
to the private quantile routine. The private lower and upper
winsorization cutoffs are searched within this range. These bounds
should be chosen on the scale of the squared principal component scores.

The parameters `trim_const` and `eta` determine the practical clipping
proportion used by the modified winsorized mean. If \\n_q\\ denotes the
number of observations used for private quantile estimation, the
clipping proportion is \$\$ p = \min\left\\
\max\left(\frac{\mathrm{trim\\const}}{n_q}, \eta\right), 0.49 \right\\.
\$\$ Here, `trim_const / n_q` controls the baseline clipping level,
while `eta` gives a lower bound reflecting the expected contamination
level.

If `split_mode = TRUE`, the sample is split into two parts: one part is
used for private quantile estimation and the other part is used for the
winsorized mean step. If `split_mode = FALSE`, all observations are used
in both steps.

The parameters `a`, `b`, `trim_const`, and `eta` are intentionally not
given defaults. They are data- and robustness-dependent choices and
should be set deliberately by the user.

## References

Ramsay K, Spicker D (2025). “Improved subsample-and-aggregate via the
private modified winsorized mean.” Code available at
<https://github.com/12ramsake/PMWM>, 2501.14095,
<https://arxiv.org/abs/2501.14095>.

## See also

[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
for computing differentially private scree estimates using these control
options.
[`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md)
for plotting scree estimates.

## Examples

``` r
pmwm_control(a = 0, b = 20, trim_const = 10, eta = 0.01)
#> $beta
#> [1] 1.001
#> 
#> $a
#> [1] 0
#> 
#> $b
#> [1] 20
#> 
#> $trim_const
#> [1] 10
#> 
#> $eta
#> [1] 0.01
#> 
#> $split_mode
#> [1] TRUE
#> 
pmwm_control(
  a = 0,
  b = 20,
  trim_const = 10,
  eta = 0.01,
  beta = 1.001,
  split_mode = FALSE
)
#> $beta
#> [1] 1.001
#> 
#> $a
#> [1] 0
#> 
#> $b
#> [1] 20
#> 
#> $trim_const
#> [1] 10
#> 
#> $eta
#> [1] 0.01
#> 
#> $split_mode
#> [1] FALSE
#> 
```
