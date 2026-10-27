# Control options for clipped scree estimation

Creates a control list for `method = "clipped"` in
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
or the scree configuration of
[`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md).

## Usage

``` r
clipped_control(C_clip)
```

## Arguments

- C_clip:

  Required positive public clipping threshold for squared pair
  differences. Larger values reduce clipping bias but increase Gaussian
  noise; data-dependent selection requires separate privacy accounting.

## Value

A control list for `method = "clipped"`.

## Details

The method clips the paired quantities at `C_clip` and adds Gaussian
noise to their vector of column means. For \\m=\lfloor n/2\rfloor\\
pairs and \\k\\ components, the calibration uses vector sensitivity
\\\sqrt{k}\\C\_{clip}/m\\ under the conditions in
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md).
Larger thresholds reduce clipping bias but increase the noise scale.

## References

Dwork C, Roth A (2014). “The Algorithmic Foundations of Differential
Privacy.” *Found. Trends Theor. Comput. Sci.*, **9**(3–4), 211–407. ISSN
1551-305X, [doi:10.1561/0400000042](https://doi.org/10.1561/0400000042)
.

## See also

[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
for using clipped scree estimation.
[`plot.dp_scree()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_scree.md)
for plotting private scree estimates.

## Examples

``` r
clipped_control(C_clip = 3)
#> $C_clip
#> [1] 3
#> 
```
