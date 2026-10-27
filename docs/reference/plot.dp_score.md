# Plot a stored private score histogram

Reuses the stored frame and histograms. Synthetic sampling is
post-processing; no direction, frame or histogram mechanism is run
again.

## Usage

``` r
# S3 method for class 'dp_score'
plot(
  x,
  type = c("histogram", "sample"),
  sampling_control = NULL,
  plot_control = NULL,
  center = TRUE,
  ...
)

# S3 method for class 'dp_score_plot'
plot(x, ...)

# S3 method for class 'dp_score_plot'
summary(object, ...)

# S3 method for class 'dp_score_plot'
print(x, ...)
```

## Arguments

- x:

  A plotting result returned by `plot.dp_score()` or
  [`plot.dp_score_group()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score_group.md).

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

  A plotting result returned by `plot.dp_score()` or
  [`plot.dp_score_group()`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score_group.md).

## Value

Invisibly returns a `dp_score_plot` object with `score`, optional
`sample`, and `plot` components, while displaying `plot$all`.

Invisibly returns the original or rebuilt comparison result.

A score-histogram summary with stored synthetic point counts.

Invisibly returns `x`.

## Details

Non-private reference panels appear automatically when retained at
fitting. Combined grouped displays show one legend per histogram method;
individual plots retain their legends. Groups with no sampled points are
omitted from sampled-point legends.

## Examples

``` r
# After estimating once:
# summary(score_fit)
# plots <- plot(score_fit, type = "histogram")
# plot(score_fit, type = "sample",
#      sampling_control = sampling_control(sample_size = 500))
```
