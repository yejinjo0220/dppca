# Differentially private scree estimation via Huber-type robust mean

Internal helper that estimates the leading scree values using a
Huber-type differentially private mean estimator applied to the centered
squared scores of each principal component.

The procedure is as follows:

1.  preprocess the data matrix for PCA,

2.  compute a non-private or private direction matrix,

3.  form the projected score matrix,

4.  for each component \\\ell\\, compute centered squared scores
    \\w\_{i\ell} = (y\_{i\ell} - \bar y\_\ell)^2\\,

5.  privately estimate a scale proxy from \\w\_{1\ell}, \dots,
    w\_{n\ell}\\ using
    [`dp_m2()`](https://yejinjo0220.github.io/dppca/reference/dp_m2.md),

6.  convert that scale proxy into a Huber threshold using
    [`tau_from_m2()`](https://yejinjo0220.github.io/dppca/reference/tau_from_m2.md),

7.  run noisy gradient descent via
    [`dp_huber_noisy_gd()`](https://yejinjo0220.github.io/dppca/reference/dp_huber_noisy_gd.md)
    to obtain a private mean estimate of \\w\_{i\ell}\\,

8.  apply the \\n/(n-1)\\ correction to obtain the final scree estimate.

If `mono = TRUE`, the final scree vector is post-processed so that it is
nonnegative and nonincreasing.

## Usage

``` r
dp_scree_huber(
  X,
  k,
  eps_total,
  delta_total,
  g_dppca = FALSE,
  cpp.option = FALSE,
  center = TRUE,
  standardize = FALSE,
  mu0 = 0,
  eta0 = 1,
  T = NULL,
  M = NULL,
  k_min_m2 = -20,
  k_max_m2 = 40,
  m2_frac = 1/4,
  mono = TRUE
)
```

## Arguments

- X:

  Numeric data matrix with observations in rows.

- k:

  Integer number of leading principal components.

- eps_total:

  Total privacy epsilon allocated to the scree routine.

- delta_total:

  Total privacy delta allocated to the scree routine.

- g_dppca:

  Logical; whether to privatize the PCA direction matrix.

- cpp.option:

  Logical passed to
  [`mech_tau_sph()`](https://yejinjo0220.github.io/dppca/reference/mech_tau_sph.md)
  when private directions are computed.

- center:

  Logical; whether to center columns before PCA.

- standardize:

  Logical; whether to standardize columns before PCA.

- mu0:

  Initial value for noisy gradient descent.

- eta0:

  Fixed step size used in noisy gradient descent.

- T:

  Optional integer number of gradient descent iterations. If `NULL`, a
  default based on `ceiling(log(n))` is used.

- M:

  Optional integer number of blocks used in
  [`dp_m2()`](https://yejinjo0220.github.io/dppca/reference/dp_m2.md).
  If `NULL`, a default based on `floor(sqrt(n / 2))` is used.

- k_min_m2:

  Integer lower bound for dyadic histogram bins used in
  [`dp_hist_m2()`](https://yejinjo0220.github.io/dppca/reference/dp_hist_m2.md).

- k_max_m2:

  Integer upper bound for dyadic histogram bins used in
  [`dp_hist_m2()`](https://yejinjo0220.github.io/dppca/reference/dp_hist_m2.md).

- m2_frac:

  Fraction of the scree privacy budget allocated to the private
  scale-proxy step
  [`dp_m2()`](https://yejinjo0220.github.io/dppca/reference/dp_m2.md).

- mono:

  Logical; whether to post-process the final scree vector so that it is
  nonnegative and nonincreasing.

## Value

A list with components:

- `scree`:

  Private scree estimates.

- `scree_np`:

  Non-private scree estimates computed as \\(n/(n-1)) \times
  mean(w\_{i\ell})\\.

- `evr`:

  Explained variance ratios computed from the private scree estimates.
