# Plotting frame for 2D scores

Internal helper that constructs a plotting frame for two-dimensional
score data. If `frame` is supplied, it is used directly. Otherwise a
differentially private frame is constructed from privatized lower and
upper quantiles computed from the pooled score coordinates.

## Usage

``` r
dp_frame(
  X,
  eps_q = NULL,
  delta_q = NULL,
  inflate = 0.1,
  q = NULL,
  frame = NULL
)
```

## Arguments

- X:

  Numeric matrix with exactly two columns.

- eps_q:

  Privacy budget for DP quantile estimation. Required only when
  `frame = NULL`.

- delta_q:

  Privacy parameter for DP quantile estimation. Required only when
  `frame = NULL`.

- inflate:

  Non-negative numeric value controlling frame expansion.

- q:

  Optional symmetric quantile level in \\(0, 1)\\. If `NULL`, extreme
  order statistics are used.

- frame:

  Optional user-specified frame. This can be either a numeric vector of
  length 2, interpreted as a common range for both axes, or a list with
  components `xlim` and `ylim`.

## Value

A list with components `xlim` and `ylim`.
