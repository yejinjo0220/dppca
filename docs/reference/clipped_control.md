# Control options for clipped scree estimation

Creates a control list for the clipped-mean scree estimator used by
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
and
[`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md)
when `method = "clipped"`.

## Usage

``` r
clipped_control(C_clip)
```

## Arguments

- C_clip:

  Positive clipping threshold for squared principal component scores.
  This value has no default because it should be chosen according to the
  scale of the data.

## Value

A list of control options for `method = "clipped"`.

## Details

The clipped method estimates each scree value by clipping the squared
scores at `C_clip` and then applying a sensitivity-calibrated Gaussian
mechanism (Dwork and Roth 2014) . Larger values of `C_clip` reduce
clipping bias but increase the sensitivity, and therefore the scale of
the privacy noise.

## References

Dwork C, Roth A (2014). “The Algorithmic Foundations of Differential
Privacy.” *Found. Trends Theor. Comput. Sci.*, **9**(3–4), 211–407. ISSN
1551-305X, [doi:10.1561/0400000042](https://doi.org/10.1561/0400000042)
.

## See also

[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
for using clipped scree estimation.
[`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md)
for plotting private scree estimates.

## Examples

``` r
clipped_control(C_clip = 3)
#> $C_clip
#> [1] 3
#> 
```
