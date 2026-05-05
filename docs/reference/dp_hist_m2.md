# Private scale proxy from a noisy dyadic histogram

Internal helper used in the Huber scree routine. Starting from a
nonnegative scalar sample `u`, this function places values into dyadic
bins indexed by powers of two, perturbs the resulting histogram counts
with Laplace noise, and returns the dyadic scale corresponding to the
largest noisy bin.

## Usage

``` r
dp_hist_m2(u, eps_m2, k_min_m2 = -20, k_max_m2 = 40)
```

## Arguments

- u:

  Numeric nonnegative vector whose scale is to be summarized.

- eps_m2:

  Positive privacy epsilon used for the noisy histogram step.

- k_min_m2:

  Integer lower bound for the dyadic bin index.

- k_max_m2:

  Integer upper bound for the dyadic bin index.

## Value

A nonnegative numeric scalar corresponding to the selected dyadic scale
level.
