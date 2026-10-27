# Control options for private modified winsorized scree estimation

Creates a control list for `method = "pmwm"` in
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
or the scree configuration of
[`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md).

## Usage

``` r
pmwm_control(a, b, trim_const, eta, beta = 1.001, split_mode = FALSE)
```

## Arguments

- a, b:

  Required finite public bounds for private quantile cutoffs, with
  `a < b`, on the scale of squared pair differences.

- trim_const:

  Required positive public trimming constant.

- eta:

  Required public tuning value in `[0, 0.5)`.

- beta:

  Geometric quantile-search base greater than `1` (`1.001` by default).

- split_mode:

  Whether to split score pairs between quantile and mean estimation
  (`FALSE` by default).

## Value

A control list for `method = "pmwm"`.

## Details

PMWM estimates lower and upper private quantiles of squared pair
differences, truncates the cutoffs to `[a, b]`, winsorizes the data, and
perturbs the mean vector. The cutoffs, not the original data, are
bounded. If the estimated upper cutoff is smaller than the lower, the
cutoffs are reset to `a` and `b`.

With `k` components, each of `2k` private quantiles receives
`epsilon / (4k)`; Gaussian mean estimation receives the remaining half
of epsilon and all delta. If `n_q` is the number of paired rows used for
quantile estimation, the trimming proportion is
\$\$p=\min\\\max(\mathrm{trim\\const}/n_q,\eta),0.49\\.\$\$ With
`split_mode = TRUE`, the pair rows are split between quantile and mean
estimation (at least four observations are required); otherwise all rows
are reused under the same privacy allocation. Tunings must be public or
selected through an accounted private procedure.

## References

Ramsay K, Spicker D (2025). “Improved subsample-and-aggregate via the
private modified winsorized mean.” Code available at
<https://github.com/12ramsake/PMWM>, 2501.14095,
<https://arxiv.org/abs/2501.14095>.

Durfee D (2023). “Unbounded Differentially Private Quantile and Maximum
Estimation.” In *Advances in Neural Information Processing Systems*,
volume 36, 77691–77712.
[doi:10.52202/075280-3394](https://doi.org/10.52202/075280-3394) ,
<https://proceedings.neurips.cc/paper_files/paper/2023/hash/f4b6ef2a78684dca2fb3f1c09372e041-Abstract-Conference.html>.

## See also

[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
for computing differentially private scree estimates using these control
options.
[`plot.dp_scree()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_scree.md)
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
#> [1] FALSE
#> 
```
