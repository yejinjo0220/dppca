# Control options for score plots

Creates a control list for the appearance of plots produced by
[`plot.dp_score()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score.md),
including grouped score plots. Score plots are constructed with
`ggplot2` and combined with `patchwork`.

## Usage

``` r
score_plot_control(
  color = "#6A5ACD",
  scatter_alpha = 0.6,
  scatter_size = 1.8,
  hist_alpha_range = c(0, 1),
  base_size = 12,
  title_size = 14,
  scatter_title = "Non-Private Scatter",
  nonprivate_title = "Non-Private Histogram",
  add_title = "Add DP Histogram",
  sparse_title = "Sparse DP Histogram",
  xlab = NULL,
  ylab = NULL,
  group_colors = NULL,
  group_ncol = 3L,
  show_group_legend = TRUE
)
```

## Arguments

- color:

  Single color used for non-group score plots. The default is
  `"#6A5ACD"`.

- scatter_alpha:

  Number in `[0, 1]` controlling point transparency in non-private and
  sampled scatter panels. The default is `0.6`.

- scatter_size:

  Positive number controlling point size in non-private and sampled
  scatter panels. The default is `1.8`.

- hist_alpha_range:

  Numeric vector of length 2 giving the minimum and maximum alpha values
  used for histogram cells. Values must satisfy
  `0 <= hist_alpha_range[1] <= hist_alpha_range[2] <= 1`. The default is
  `c(0, 1)`.

- base_size:

  Positive base font size for score plots. The default is `12`.

- title_size:

  Positive plot-title font size. The default is `14`.

- scatter_title:

  Optional title for the non-private scatter panel. The default is
  `"Non-Private Scatter"`.

- nonprivate_title:

  Optional title for the non-private histogram panel. The default is
  `"Non-Private Histogram"`.

- add_title:

  Optional title for the additive private histogram panel. The default
  is `"Add DP Histogram"`.

- sparse_title:

  Optional title for the sparse private histogram panel. The default is
  `"Sparse DP Histogram"`.

- xlab:

  Optional x-axis label. If `NULL`, the selected PC name is used.

- ylab:

  Optional y-axis label. If `NULL`, the selected PC name is used.

- group_colors:

  Optional nonempty character vector of colors for grouped score plots.
  Named colors are matched to group labels; missing groups receive
  generated colors. Unnamed colors are recycled in group order. If
  `NULL`, valid color labels are used directly, otherwise a palette is
  generated automatically.

- group_ncol:

  Positive integer retained for backward compatibility. Grouped score
  plots are currently overlaid by group color rather than arranged as
  separate group-specific panels, so this option is not used by
  [`plot.dp_score()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score.md).
  The default is `3`.

- show_group_legend:

  Logical value indicating whether to display the group legend in
  grouped score plots. The default is `TRUE`. In combined grouped
  displays, each histogram method uses one legend; standalone panels
  retain their own legends. Groups with no sampled points are omitted
  from sampled-point legends.

## Value

A list of plotting options for
[`plot.dp_score()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score.md).

## Examples

``` r
score_plot_control()
#> $color
#> [1] "#6A5ACD"
#> 
#> $scatter_alpha
#> [1] 0.6
#> 
#> $scatter_size
#> [1] 1.8
#> 
#> $hist_alpha_range
#> [1] 0 1
#> 
#> $base_size
#> [1] 12
#> 
#> $title_size
#> [1] 14
#> 
#> $scatter_title
#> [1] "Non-Private Scatter"
#> 
#> $nonprivate_title
#> [1] "Non-Private Histogram"
#> 
#> $add_title
#> [1] "Add DP Histogram"
#> 
#> $sparse_title
#> [1] "Sparse DP Histogram"
#> 
#> $xlab
#> NULL
#> 
#> $ylab
#> NULL
#> 
#> $group_colors
#> NULL
#> 
#> $group_ncol
#> [1] 3
#> 
#> $show_group_legend
#> [1] TRUE
#> 

score_plot_control(
  scatter_size = 2.2,
  title_size = 13,
  nonprivate_title = "Non-Private score histogram"
)
#> $color
#> [1] "#6A5ACD"
#> 
#> $scatter_alpha
#> [1] 0.6
#> 
#> $scatter_size
#> [1] 2.2
#> 
#> $hist_alpha_range
#> [1] 0 1
#> 
#> $base_size
#> [1] 12
#> 
#> $title_size
#> [1] 13
#> 
#> $scatter_title
#> [1] "Non-Private Scatter"
#> 
#> $nonprivate_title
#> [1] "Non-Private score histogram"
#> 
#> $add_title
#> [1] "Add DP Histogram"
#> 
#> $sparse_title
#> [1] "Sparse DP Histogram"
#> 
#> $xlab
#> NULL
#> 
#> $ylab
#> NULL
#> 
#> $group_colors
#> NULL
#> 
#> $group_ncol
#> [1] 3
#> 
#> $show_group_legend
#> [1] TRUE
#> 
```
