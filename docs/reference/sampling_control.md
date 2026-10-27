# Control options for synthetic score sampling

Creates a control list for generating synthetic score points from a
score histogram. Sampling from a released private histogram is
post-processing under the conditions below; it does not alter the
histogram mechanism.
[`plot.dp_pca()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_pca.md)
can also use these controls to sample an empirical reference histogram,
which does not make that reference private.

## Usage

``` r
sampling_control(
  method = "uniform",
  sample_size = NULL,
  bandwidth_scale = NULL
)
```

## Arguments

- method:

  Character vector specifying the sampling construction. `"center"`
  selects a histogram bin according to its histogram probability and
  samples around that bin center using Gaussian jitter. `"uniform"`
  selects a histogram bin according to its histogram probability,
  samples uniformly inside the selected bin, and then adds Gaussian
  jitter. The default is `"uniform"`. Both methods can be requested with
  `c("center", "uniform")` in
  [`plot.dp_score()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score.md);
  [`plot.dp_pca()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_pca.md)
  requires exactly one sampling method.

- sample_size:

  Optional positive integer giving the synthetic sample size for each
  requested sampling construction. For
  [`plot.dp_score()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score.md),
  this is the number of synthetic score points, or the total across all
  groups for a grouped fit. For
  [`plot.dp_pca()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_pca.md),
  it is the size shared by score and biplot samples within each version;
  an explicit value also samples non-private reference histograms when
  comparison panels are requested. With `NULL`, private samples use the
  number of input observations and integrated non-private panels use
  observed scores. The default private sample size assumes that the
  input sample size is public. See
  [`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
  and
  [`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
  for grouped sampling weights and the applicable privacy conditions.

- bandwidth_scale:

  Gaussian smoothing scale. If `NULL`, method-specific defaults are
  used: `0.5` for `"center"` and `0.25` for `"uniform"`. A single
  nonnegative number applies the same scale to all requested sampling
  methods. A named numeric vector can specify method-specific values,
  for example `c(center = 0.5, uniform = 0.2)`; it must include every
  requested method. A scale of zero disables Gaussian jitter.

## Value

A `dppca_sampling_control` list with components `method`, `sample_size`,
and `bandwidth_scale`. The bandwidth values are resolved to a named
numeric vector, while `sample_size = NULL` is resolved later by the
plotting function.

## Details

For released histogram probabilities `p_k` and cells `B_k`, first draw
\$\$J \sim \mathrm{Categorical}(p_1,\ldots,p_m).\$\$ With
`method = "center"`, sample `S^* = c_J + Z`, where `c_J` is the selected
cell center. With `method = "uniform"`, sample `U` uniformly in `B_J`
and use `S^* = U + Z`. In both cases Gaussian jitter has covariance
\$\$\mathrm{Cov}(Z)=c^2\\\mathrm{diag}(\Delta_x^2,\Delta_y^2),\$\$ where
`c` is `bandwidth_scale` and the deltas are bin widths.

Sampling uses only released histogram probabilities and bin geometry.
With publicly chosen settings it is post-processing and requires no new
privacy budget. Jitter can place points outside the frame; an all-zero
histogram yields no sampled points and a warning.

## Examples

``` r
sampling_control()
#> $method
#> [1] "uniform"
#> 
#> $sample_size
#> NULL
#> 
#> $bandwidth_scale
#> uniform 
#>    0.25 
#> 
#> attr(,"class")
#> [1] "dppca_sampling_control" "list"                  

sampling_control(
  method = c("center", "uniform"),
  sample_size = 1000,
  bandwidth_scale = c(center = 0.5, uniform = 0.2)
)
#> $method
#> [1] "center"  "uniform"
#> 
#> $sample_size
#> [1] 1000
#> 
#> $bandwidth_scale
#>  center uniform 
#>     0.5     0.2 
#> 
#> attr(,"class")
#> [1] "dppca_sampling_control" "list"                  
```
