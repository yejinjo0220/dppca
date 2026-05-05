# Convert a private scale proxy into a Huber threshold

Internal helper used by the Huber scree routine. Given a private scale
proxy `m2_hat`, this function produces a scalar Huber robustification
threshold.

## Usage

``` r
tau_from_m2(m2_hat, eps_tau, n_tau)
```

## Arguments

- m2_hat:

  Nonnegative private scale proxy.

- eps_tau:

  Positive privacy epsilon allocated to the Huber noisy-gradient-descent
  step for one component.

- n_tau:

  Effective sample size used in the threshold calculation.

## Value

Positive numeric scalar Huber threshold.
