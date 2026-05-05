# DP Huber noisy gradient descent for a scalar mean

Internal helper that implements a scalar Huber-type differentially
private mean estimator via noisy gradient descent. The optional argument
`scale_factor` is useful when the released quantity should be a scaled
version of the mean, such as the \\n/(n-1)\\-adjusted scree convention.

## Usage

``` r
dp_huber_noisy_gd(w, eps_gd, delta_gd, tau, T, mu0 = 0, eta0 = 1)
```

## Arguments

- w:

  Numeric vector representing the one-dimensional sample.

- eps_gd:

  Positive privacy epsilon allocated to the noisy-gradient-descent
  stage.

- delta_gd:

  Positive privacy delta allocated to the noisy-gradient-descent stage.

- tau:

  Positive Huber threshold.

- T:

  Positive integer number of gradient-descent iterations.

- mu0:

  Initial value for the iterative procedure.

- eta0:

  Base step size.

- step_schedule:

  Character string specifying the step-size schedule. Allowed values are
  `"fixed"` and `"1/sqrt(t)"`.

- scale_factor:

  Positive multiplicative factor applied to the released estimate and to
  the noise calibration.

## Value

Numeric scalar representing the final private estimate.
