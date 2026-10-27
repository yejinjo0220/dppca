# Control options for scree and PVE plots

Creates a control list for the appearance of plots produced by
[`plot.dp_scree()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_scree.md).
The plotting function uses base R graphics.

## Usage

``` r
scree_plot_control(
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  legend_position = "topright",
  legend_labels = NULL,
  col = NULL,
  lty = NULL,
  pch = NULL,
  lwd = 1,
  point_cex = 1,
  cex_main = 1.5,
  cex_lab = 1.2,
  cex_axis = 1.1,
  cex_legend = 1.1,
  legend_bty = "n"
)
```

## Arguments

- title:

  Optional plot title. If `NULL`,
  [`plot.dp_scree()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_scree.md)
  uses `"DP PVE Plot"` when `type = "pve"` and `"DP Scree Plot"` when
  `type = "scree"`.

- xlab:

  Optional x-axis label. If `NULL`, `"Component"` is used.

- ylab:

  Optional y-axis label. If `NULL`, `"PVE"` is used for PVE and
  `"Scree Value"` for scree values.

- legend_position:

  Legend position passed to
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html). The
  default is `"topright"`.

- legend_labels:

  Optional named character vector giving legend labels. Recognized names
  are `"nonprivate"`, `"clipped"`, `"pmwm"`, and `"huber"`.

- col:

  Optional named character vector of plotting colors, using the same
  series names as `legend_labels`.

- lty:

  Optional named vector of line types, using the same series names as
  `legend_labels`.

- pch:

  Optional named vector of point symbols, using the same series names as
  `legend_labels`.

- lwd:

  Positive line width. The default is `1`.

- point_cex:

  Positive point-size multiplier. The default is `1`.

- cex_main:

  Positive title-size multiplier. The default is `1.5`.

- cex_lab:

  Positive axis-label-size multiplier. The default is `1.2`.

- cex_axis:

  Positive axis-text-size multiplier. The default is `1.1`.

- cex_legend:

  Positive legend-text-size multiplier. The default is `1.1`.

- legend_bty:

  Box type passed to
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html). The
  default `"n"` suppresses the legend border.

## Value

A list of plotting options for
[`plot.dp_scree()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_scree.md).

## Examples

``` r
scree_plot_control()
#> $title
#> NULL
#> 
#> $xlab
#> NULL
#> 
#> $ylab
#> NULL
#> 
#> $legend_position
#> [1] "topright"
#> 
#> $legend_labels
#> NULL
#> 
#> $col
#> NULL
#> 
#> $lty
#> NULL
#> 
#> $pch
#> NULL
#> 
#> $lwd
#> [1] 1
#> 
#> $point_cex
#> [1] 1
#> 
#> $cex_main
#> [1] 1.5
#> 
#> $cex_lab
#> [1] 1.2
#> 
#> $cex_axis
#> [1] 1.1
#> 
#> $cex_legend
#> [1] 1.1
#> 
#> $legend_bty
#> [1] "n"
#> 

scree_plot_control(
  title = "PVE comparison",
  legend_position = "topright",
  cex_main = 1.4
)
#> $title
#> [1] "PVE comparison"
#> 
#> $xlab
#> NULL
#> 
#> $ylab
#> NULL
#> 
#> $legend_position
#> [1] "topright"
#> 
#> $legend_labels
#> NULL
#> 
#> $col
#> NULL
#> 
#> $lty
#> NULL
#> 
#> $pch
#> NULL
#> 
#> $lwd
#> [1] 1
#> 
#> $point_cex
#> [1] 1
#> 
#> $cex_main
#> [1] 1.4
#> 
#> $cex_lab
#> [1] 1.2
#> 
#> $cex_axis
#> [1] 1.1
#> 
#> $cex_legend
#> [1] 1.1
#> 
#> $legend_bty
#> [1] "n"
#> 
```
