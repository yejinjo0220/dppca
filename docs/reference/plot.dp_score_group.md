# Plot stored group-wise private score histograms

Uses the same display options as
[`plot.dp_score()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score.md).
Group sample sizes follow observed proportions if non-private counts
were retained; otherwise they are allocated equally across public group
labels.

## Usage

``` r
# S3 method for class 'dp_score_group'
plot(
  x,
  type = c("histogram", "sample"),
  sampling_control = NULL,
  plot_control = NULL,
  center = TRUE,
  ...
)

# S3 method for class 'dp_score_plot_group'
plot(x, ...)

# S3 method for class 'dp_score_plot_group'
summary(object, ...)

# S3 method for class 'dp_score_plot_group'
print(x, ...)
```

## Arguments

- x:

  A grouped result from
  [`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
  with `group` supplied.

- type:

  `"histogram"`, `"sample"`, or both (`c("histogram", "sample")` by
  default). Samples are synthetic for DP histograms; ordinary-reference
  sample plots show centered observed scores unless an explicit
  `sample_size` is requested.

- sampling_control:

  Optional
  [`sampling_control()`](https://yejinjo0220.github.io/dppca/reference/sampling_control.md)
  object used for `"sample"`; `NULL` uses defaults.

- plot_control:

  Optional
  [`score_plot_control()`](https://yejinjo0220.github.io/dppca/reference/score_plot_control.md)
  appearance settings; `NULL` uses defaults.

- center:

  Whether to shift displayed private coordinates by the midpoint of the
  private frame (`TRUE` by default). Ordinary reference coordinates
  remain centered at zero. This changes only the display.

- ...:

  Passed to the summary print method.

- object:

  A plotting result returned by
  [`plot.dp_score()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score.md)
  or `plot.dp_score_group()`.

## Value

Invisibly returns a `dp_score_plot_group` object while displaying
`plot$all`.
