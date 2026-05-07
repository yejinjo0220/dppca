# Control parameters for the Huber DP scree estimator

Creates a control list for `method = "huber"` in
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
and
[`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md).

## Usage

``` r
huber_control(
  mu0 = 0,
  eta0 = 1,
  T = NULL,
  M = NULL,
  k_min_m2 = -40,
  k_max_m2 = 40,
  m2_frac = 1/4
)
```

## Arguments

- mu0:

  Initial value used by the Huber noisy gradient descent.

- eta0:

  Fixed step size used by the Huber noisy gradient descent.

- T:

  Optional integer number of Huber gradient descent iterations.

- M:

  Optional integer number of blocks used in
  [`dp_m2()`](https://yejinjo0220.github.io/dppca/reference/dp_m2.md).

- k_min_m2:

  Integer lower bound for dyadic histogram bins used in the Huber
  scale-proxy step.

- k_max_m2:

  Integer upper bound for dyadic histogram bins used in the Huber
  scale-proxy step.

- m2_frac:

  Fraction of the Huber scree privacy budget allocated to the private
  scale-proxy step.

## Value

A list of Huber-estimator tuning parameters.
